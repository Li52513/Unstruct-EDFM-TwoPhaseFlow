#pragma once
#include <iostream>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "Solver_AssemblerCOO.h"

namespace GeneralTools
{
    //  复制标量场：dst <- src
    inline bool copyField
    (
        FieldRegistry& reg, // 需要确保 src 和 dst 都存在
        const std::string& src_name, //被复制的场
        const std::string& dst_name  //目标场
    )
    {
        auto src = reg.get<volScalarField>(src_name);
        auto dst = reg.get<volScalarField>(dst_name);
        if (!src || !dst) return false;
        if (dst->data.size() != src->data.size()) dst->data.resize(src->data.size(), 0.0);
        std::copy(src->data.begin(), src->data.end(), dst->data.begin()); // 复制
        return true;
    }
    //  外迭代步开始
    inline bool startOuterIteration_scatter
    (
        FieldRegistry& reg,
        const std::string& x_name,
        const std::string& x_prev_name
    )
    {
        bool ok = true;
        ok = ok && copyField(reg, x_prev_name, x_name);
        return ok;

    }

    //  建立 “cellId -> 方程自由度编号”
    inline std::vector<int> buildUnknownMap(Mesh& mesh, int& nUnknowns)
    {
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        std::vector<int> lid_of_cell(cells.size(), -1); // -1 表示该 cellId 不在方程中（如 ghost）
        nUnknowns = 0;
        for (const auto& c : cells)
        {
            if (c.id < 0) continue; 
            const size_t i = id2idx.at(c.id);
            lid_of_cell[i] = nUnknowns++;
        }
        return lid_of_cell;
    }

    // 将网格场变量转化为向量
    inline std::vector<double> gatherFieldToVec(const FieldRegistry& reg, Mesh& mesh, const std::string& fld, const std::vector<int>& lid_of_cell, int N)
    {
        std::vector<double> x(N, 0.0);
        auto f = reg.get<volScalarField>(fld);
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        if (!f) return x;
        for (const auto& c : cells) {
            if (c.id < 0)continue;
            size_t i = id2idx.at(c.id);
            int r = lid_of_cell[i];
            if (r >= 0 && r < N) x[r] = (*f)[i]; // 只赋值自由度
        }
        return x;
    }
    
    // 将向量映射到网格场变量中
    inline void scatterVecToField(FieldRegistry& reg, Mesh& mesh, const std::string& fld, const std::vector<int>& lid_of_cell, const std::vector<double>& x)
    {
        auto f = reg.get<volScalarField>(fld);
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        if (!f) return;
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            size_t i = id2idx.at(c.id);
            int r = lid_of_cell[i];
            if (r >= 0) (*f)[i] = x[r];
        }
    }

    // 将线性系统维度从 Nc 扩展到 Ntot（保留已有 A 与 b 条目）
    inline void extend_linear_system_size(SparseSystemCOO& sys, int Ntot)
    {
        if (Ntot <= sys.n) return;
        // 扩展 RHS 到 Ntot（保留已有值，新增填 0）
        if ((int)sys.b.size() < Ntot) sys.b.resize(Ntot, 0.0);
        sys.n = Ntot; // 注意：在调用 addA/addb 前先设 n，使断言通过
    }
    
    //  欠松弛：x <- x_prev + alpha * (x - x_prev), 在 x 场内就地更新
    inline bool underRelaxInPlace
    (
        FieldRegistry& reg,
        const std::string& x_name,
        const std::string& x_prev_name,
        double alpha)  // 0<alpha<=1
    {
        auto x = reg.get<volScalarField>(x_name);
        auto xp = reg.get<volScalarField>(x_prev_name);
        if (!x || !xp) return false;
        if (x->data.size() != xp->data.size()) return false;
        alpha = std::max(0.0, std::min(1.0, alpha));
        for (size_t i = 0; i < x->data.size(); ++i)
        {
            x->data[i] = xp->data[i] + alpha * (x->data[i] - xp->data[i]);
        }
        return true;
    }

    //  更行迭代层
    inline bool updatePrevIterates(FieldRegistry& reg,
        const std::initializer_list<std::pair<std::string, std::string>>& pairs)
    {
        bool ok = true;
        for (const auto& pr : pairs) {
            ok = ok && copyField(reg, pr.first, pr.second);
        }
        return ok;
    }

    //  收敛判据：∞-范数(最大绝对差)
    inline double maxAbsDiff
    (
        const FieldRegistry& reg,
        const std::string& x_name,
        const std::string& x_prev_name)
    {
        auto a = reg.get<volScalarField>(x_name);
        auto b = reg.get<volScalarField>(x_prev_name);
        if (!a || !b || a->data.size() != b->data.size()) return std::numeric_limits<double>::infinity();
        double m = 0.0;
        for (size_t i = 0; i < a->data.size(); ++i) m = std::max(m, std::abs(a->data[i] - b->data[i]));
        return m;
    }

    //  收敛判据：L2-范数(均方根)
    inline double rmsDiff
    (
        const FieldRegistry& reg,
        const std::string& x_name,
        const std::string& x_prev_name)
    {
        auto a = reg.get<volScalarField>(x_name);
        auto b = reg.get<volScalarField>(x_prev_name);
        if (!a || !b || a->data.size() != b->data.size() || a->data.empty()) return std::numeric_limits<double>::infinity();
        double sum2 = 0.0;
        for (size_t i = 0; i < a->data.size(); ++i) {
            double d = a->data[i] - b->data[i];
            sum2 += d * d;
        }
        return std::sqrt(sum2 / static_cast<double>(a->data.size()));
    }

    //  收敛判据：最大相对变化（∞-范数）
    inline double maxRelChange(const FieldRegistry& reg,
        const std::string& x_name,
        const std::string& x_prev_name,
        double eps = 1e-20)
    {
        auto x = reg.get<volScalarField>(x_name);
        auto xp = reg.get<volScalarField>(x_prev_name);
        if (!x || !xp || x->data.size() != xp->data.size()) return std::numeric_limits<double>::infinity();
        double m = 0.0;
        for (size_t i = 0; i < x->data.size(); ++i) {
            double denom = std::max(std::abs(xp->data[i]), eps);
            m = std::max(m, std::abs(x->data[i] - xp->data[i]) / denom);
        }
        return m;
    }

}