#pragma once
#include <iostream>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "Solver_AssemblerCOO.h"

namespace GeneralTools
{
    //  复制标量场：dst <- src
    inline bool copyField1
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

    inline bool ensureTransientFields_scalar1
    (
        Mesh& mesh,
        FieldRegistry& reg,
        const std::string& var_name,
        const std::string& var_old_name,
        const std::string& var_prev_name
    )
    {
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
        const size_t n = cells.size();

        // —— 基础场：若不存在就创建，若尺寸不符就调整 —— //
        auto var = reg.get<volScalarField>(var_name);
        bool base_created_or_resized = false;
        if (!var) {
            var = reg.getOrCreate<volScalarField>(var_name, n, 0.0);
            base_created_or_resized = true;
        }
        else if (var->data.size() != n) {
            var->data.resize(n, 0.0);
            base_created_or_resized = true;
        }

        // 小工具：确保某个场存在且尺寸为 n
        auto ensureSized = [&](const std::string& name, std::shared_ptr<volScalarField>& out)->bool {
            out = reg.get<volScalarField>(name);
            bool changed = false;
            if (!out) { out = reg.getOrCreate<volScalarField>(name, n, 0.0); changed = true; }
            else if (out->data.size() != n) { out->data.resize(n, 0.0); changed = true; }
            return changed;
            };

        std::shared_ptr<volScalarField> var_old, var_prev;
        const bool need_init_old = ensureSized(var_old_name, var_old);
        const bool need_init_prev = ensureSized(var_prev_name, var_prev);

        // 若新建/尺寸变化（包括基础场）→ 用当前 var 值初始化 *_old/*_prev
        if (need_init_old || need_init_prev || base_created_or_resized) {
            for (const auto& c : cells) {
                if (c.id < 0) continue; // 跳过 ghost
                const size_t i = id2idx.at(c.id);
                if (need_init_old || base_created_or_resized) (*var_old)[i] = (*var)[i];
                if (need_init_prev || base_created_or_resized) (*var_prev)[i] = (*var)[i];
            }
        }
        return true;
    }

    /// 时间步迭代开始
    inline bool startTimeStep_scalar(Mesh& mesh, FieldRegistry& reg, const std::string& x_name, const std::string& x_old_name, const std::string& x_prev_name)
    {
        if (!ensureTransientFields_scalar1(mesh, reg, x_name, x_old_name, x_prev_name)) return false;
        // x^n
        if (!copyField1(reg, x_name, x_old_name)) return false; //将x 复制到_old
        // prev ← old (给 k=0 的外迭代作为初值)
        if (!copyField1(reg, x_old_name, x_prev_name)) return false; //将x_old 复制到 x_prev
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
        ok = ok && copyField1(reg, x_prev_name, x_name);
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
            ok = ok && copyField1(reg, pr.first, pr.second);
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