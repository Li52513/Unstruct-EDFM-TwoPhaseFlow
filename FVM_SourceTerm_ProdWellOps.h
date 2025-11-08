// FVM_SourceTerm_ProdWellOps.h
#pragma once
#include <vector>
#include <cassert>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "VolField.h"      // volScalarField

namespace FVM {
    namespace SourceTerm {

        // —— 小工具：计算“流入指定单元”的面质量通量（约定 mf>0: owner→neighbour）——
        static inline double influxIntoCell(const ::Cell& c, const ::Face& f, double mf_face)
        {
            if (f.ownerCell == c.id) return -mf_face; // 流入 owner 记正
            if (f.neighborCell == c.id) return +mf_face; // 流入 neighbour 记正
            return 0.0;
        }

        /**
         * @brief 统计井区净质量流入（kg/s；2D 默认单位厚度）
         * 规则：对 mask=1 的所有单元，把“指向该单元”的面质量通量累加；掩码内部共享面会抵消，仅剩边界净流入。
         * @param maskName   井区掩码名（cell 标量场，>0.5 视为 1）
         * @param mf_name    面质量通量场名（如 "mf_g"）
         * @param thickness  几何厚度（2D=1.0；3D 传实际厚度）
         * @return           掩码域的净流入（采出井应为正）
         */
        inline double accumulateMassIntoMask(::MeshManager& mgr,
            ::FieldRegistry& reg,
            ::FaceFieldRegistry& freg,
            const std::string& maskName,   // "prod_mask"
            const std::string& mf_name,    // "mf_g" (kg/s; owner->neigh 为正)
            double thickness = 1.0)
        {
            ::Mesh& mesh = mgr.mesh();

            auto* mask = reg.get< ::volScalarField >(maskName).get();
            auto* mf = freg.get< ::faceScalarField >(mf_name).get();
            if (!mask || !mf) return 0.0;

            const auto& faces = mesh.getFaces();
            const auto& id2idx = mesh.getCellId2Index();

            double Qm = 0.0;
            for (size_t fid = 0; fid < faces.size(); ++fid) {
                const ::Face& f = faces[fid];
                if (f.neighborCell < 0) continue;                 // 真边界面跳过
                const int io = id2idx.at(f.ownerCell);
                const int in = id2idx.at(f.neighborCell);
                const bool mo = mask->data[io] > 0.5;
                const bool mn = mask->data[in] > 0.5;
                if (mo == mn) continue;                           // 不是“掩码边界面”

                const double mdot = mf->data[fid];                // ★ 用 fid 访问
                // 流入 mask 的条件：mask 在 owner 侧 => mdot<0；mask 在 neigh 侧 => mdot>0
                if (mo && !mn && mdot < 0.0) Qm += (-mdot) * thickness;
                if (!mo && mn && mdot > 0.0) Qm += (mdot)*thickness;
            }
            return Qm;                                            // kg/s
        }

        /**
         * @brief 为井区添加“温度抽能源”（半隐式）：a_src_T += (∑入井面 ṁ) * cp，b_src_T 不变
         * 物理：把进入井区的质量流率按 cp·T 抽走，避免井内“囤热”。等价于在 T 方程增加 a*T 的耗散项。
         * @param maskName   井区掩码（cell 标量场，>0.5 视为 1）
         * @param mf_name    面质量通量名（如 "mf_g"）
         * @param cp_name    单元比热名（如 "cp_g"）
         * @param tagT       温度场标识（默认 "T"；写入 a_src_T / b_src_T）
         * @param accumulate 是否在原有源项上累加（false 则先清零）
         * @param thickness  几何厚度（2D=1.0；3D 传实际厚度）
         * @return           true 表示成功写入/累加
         */
        inline bool addTemperatureSinkFromFlux(::MeshManager& mgr,
            ::FieldRegistry& reg,
            ::FaceFieldRegistry& freg,
            const std::string& maskName,
            const std::string& mf_name,
            const std::string& cp_name,
            const std::string& tagT = "T",
            bool accumulate = true,
            double thickness = 1.0)
        {
            ::Mesh& mesh = mgr.mesh();

            ::volScalarField* mask = reg.get< ::volScalarField >(maskName).get();
            ::faceScalarField* mf = freg.get< ::faceScalarField >(mf_name).get();
            ::volScalarField* cp = reg.get< ::volScalarField >(cp_name).get();
            if (!mask || !mf || !cp) return false;

            const auto& faces = mesh.getFaces();
            const auto& id2idx = mesh.getCellId2Index();
            const size_t nCell = mesh.getCells().size();

            // 源项场：a_src_T / b_src_T
            std::shared_ptr< ::volScalarField > a_sh = reg.getOrCreate< ::volScalarField >("a_src_" + tagT, nCell, 0.0);
            std::shared_ptr< ::volScalarField > b_sh = reg.getOrCreate< ::volScalarField >("b_src_" + tagT, nCell, 0.0);
            ::volScalarField* a = a_sh.get();
            ::volScalarField* b = b_sh.get();
            if (!accumulate) { a->data.assign(nCell, 0.0); b->data.assign(nCell, 0.0); }

            for (const ::Cell& c : mesh.getCells()) {
                const int ci = id2idx.at(c.id);
                if (mask->data[ci] <= 0.5) continue;

                double qm_in = 0.0; // 流入井单元的质量流率
                for (int fid : c.CellFaceIDs) {
                    const ::Face& f = faces[fid];
                    const double mdot = mf->data[fid];
                    qm_in += influxIntoCell(c, f, mdot);
                }
                if (qm_in <= 0.0) continue;
                const double V = std::max(0.0, c.volume);
                const double cpv = cp->data[ci];
                // 半隐式抽能：等价于在方程左端加 K*T（把能量抽向“0 K”的虚拟库），K = qm*cp
                 // 与你的 PI 注入源保持同一量纲：把“边界通量(W)”作为“体积分布源”写入 => 乘 V
                a->data[ci] += (qm_in * thickness) * cpv; // 半隐式：a += qm*cp；b 不变
            }
            (void)b; // 当前策略不写 b，避免未使用警告
            return true;
        }


        // 严格抽走“流入井区的焓流”：b_src_T -= Σ_f( m_in * cp_up * T_up ) * V
       // 仅跨“掩码边界面”抽能：b_src_T -= Σ_f (m_in * cp_up * T_up) * V_mask
        inline bool addTemperatureExtractionFromInflowFaces(::MeshManager& mgr,
            ::FieldRegistry& reg,
            ::FaceFieldRegistry& freg,
            const std::string& maskName,   // "prod_mask"
            const std::string& mf_name,    // "mf_g"
            const std::string& T_name,     // "T"
            const std::string& cp_name,    // "cp_g"
            const std::string& tagT = "T",
            bool accumulate = true,
            double thickness = 1.0)
        {
            ::Mesh& mesh = mgr.mesh();

            auto* mask = reg.get< ::volScalarField >(maskName).get();
            auto* mf = freg.get< ::faceScalarField >(mf_name).get();
            auto* T = reg.get< ::volScalarField >(T_name).get();
            auto* cp = reg.get< ::volScalarField >(cp_name).get();
            if (!mask || !mf || !T || !cp) return false;

            const auto& faces = mesh.getFaces();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            auto  b_sh = reg.getOrCreate< ::volScalarField >("b_src_" + tagT, cells.size(), 0.0);
            auto* b = b_sh.get();
            if (!accumulate) b->data.assign(cells.size(), 0.0);

            for (size_t fid = 0; fid < faces.size(); ++fid) {
                const ::Face& f = faces[fid];
                if (f.neighborCell < 0) continue;
                const int io = id2idx.at(f.ownerCell);
                const int in = id2idx.at(f.neighborCell);
                const bool mo = mask->data[io] > 0.5;
                const bool mn = mask->data[in] > 0.5;
                if (mo == mn) continue; // 只要掩码边界面

                const double mdot = mf->data[fid];                // kg/s
                int cell_in_mask = -1, upCell = -1; double m_in = 0.0;
                if (mo && !mn && mdot < 0.0) { cell_in_mask = f.ownerCell;    upCell = f.neighborCell; m_in = -mdot; }
                if (!mo && mn && mdot > 0.0) { cell_in_mask = f.neighborCell; upCell = f.ownerCell;    m_in = mdot; }
                if (m_in <= 0.0) continue;

                const int ui = id2idx.at(upCell);
                const int ci = id2idx.at(cell_in_mask);
                const double Tup = T->data[ui];
                const double cpup = cp->data[ui];
                const double Qin = m_in * cpup * Tup * thickness;         // W
                b->data[ci] -= Qin;       // 负号=抽能（体积分布）
            }
            return true;
        }



    }
} // namespace FVM::SourceTerm
