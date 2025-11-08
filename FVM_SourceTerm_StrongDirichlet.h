#pragma once
#include <vector>
#include <cassert>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "VolField.h"
#include "Solver_AssemblerCOO.h"  

namespace FVM
{
	namespace SourceTerm
	{
        // 将 cell 掩码场转换为 0/1 数组（>th 记 1）
        inline std::vector<char> buildMask01(const ::Mesh& mesh,
            const ::volScalarField& mask, double th = 0.5)
        {
            const size_t n = mesh.getCells().size();
            std::vector<char> mk(n, 0);
            for (size_t i = 0; i < n; ++i) mk[i] = (mask.data[i] > th) ? 1 : 0;
            return mk;
        }

        /**
         * @brief 在 COO 线性系统中对若干未知施加强 Dirichlet: x[i]=target
         * 做法：一次遍历 triplets
         *  (1) 对每个被固定列 c，将 A(r,c) 贡献搬到 RHS：b[r] -= A(r,c)*target[c]，并把 A(r,c)=0（r!=c）
         *  (2) 对每个被固定行 r，将该行全部清零（含对角，稍后另行置对角=1）
         *  (3) 最后对每个固定未知 i，追加/设置 (i,i)=1, b[i]=target[i]
         *
         * @param mesh  网格（仅用于尺寸校验）
         * @param lid_cell2unk  cell->unknown 的局部索引映射（负数表示无未知，例如 ghost）
         * @param cellMask01    |cells|，1 表示对该 cell 的未知施加强制
         * @param targetVal     |cells|，给被强制 cell 的目标值（未被强制的元素忽略）
         * @param sys           稀疏系统（就地修改 triplets 与 rhs）
        */

        inline void applyStrongDirichletCOO(const ::Mesh& mesh,
            const std::vector<int>& lid_cell2unk,   // |cells| -> unknown id (>=0)
            const std::vector<char>& cellMask01,    // |cells|, 1=强制
            const std::vector<double>& targetVal,   // |cells|
            ::SparseSystemCOO& sys)
        {
            // 0) 基本校验
            const int nUnknowns = sys.n;
            assert(nUnknowns >= 0);
            assert((size_t)nUnknowns == sys.b.size());

            // 1) 标记被固定的未知及其目标值（按 unknown 索引存储，便于 O(1) 查询）
            std::vector<char>   fixed(nUnknowns, 0);
            std::vector<double> tval(nUnknowns, 0.0);

            const size_t nC = mesh.getCells().size();
            for (size_t cid = 0; cid < nC; ++cid) {
                if (!cellMask01[cid]) continue;
                const int i = lid_cell2unk[cid];
                if (i < 0 || i >= nUnknowns) continue; // 该 cell 无对应未知（可能是 ghost/裁剪）
                fixed[i] = 1;
                tval[i] = targetVal[cid];
            }

            // 2) 单次扫描三元组：清列并搬 RHS；清行（含对角）
            for (auto& t : sys.A) {
                const int r = t.r, c = t.c;
                // 列被固定：把列贡献搬到 RHS（除对角位置外），然后把该条目清为 0
                if (c >= 0 && c < nUnknowns && fixed[c] && r != c) {
                    sys.b[r] -= t.v * tval[c];
                    t.v = 0.0;
                }
                // 行被固定：整行清零（包括对角，后面再单独补 (i,i)=1）
                if (r >= 0 && r < nUnknowns && fixed[r]) {
                    t.v = 0.0;
                }
            }

            // 3) 为每个固定未知 i：补 (i,i)=1，并设置 RHS
            for (int i = 0; i < nUnknowns; ++i) {
                if (!fixed[i]) continue;
                sys.addA(i, i, 1.0);
                sys.b[i] = tval[i];
            }
        }
       
        /** 求解后（便于通量重建/可视化）把掩码上的 cell 标量场置常数值 */
        inline void assignCellValueOnMask(::FieldRegistry& reg,
            const std::string& cellFieldName,
            const std::string& maskName,
            double value)
        {
            auto* fld = reg.get< ::volScalarField >(cellFieldName).get();
            auto* msk = reg.get< ::volScalarField >(maskName).get();
            if (!fld || !msk) return;
            const size_t n = fld->data.size();
            for (size_t i = 0; i < n; ++i) if (msk->data[i] > 0.5) fld->data[i] = value;
        }


	}
}