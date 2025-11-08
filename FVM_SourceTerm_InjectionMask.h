// FVM_SourceTerm_InjectionMask.h
#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "VolField.h"

// 所有实现都放在 FVM::SourceTerm 下
namespace FVM {
    namespace SourceTerm {

        // 圆形注入几何
        struct wellCircle {
            double cx{ 0.0 }, cy{ 0.0 }, r{ 0.0 };
        };

        /**
         * @brief 基于“单元质心落在圆内”的稳健掩码构建（2D）
         * 生成体标量场 out_name（通常叫 "inj_mask"），长度=Cells.size()。
         * - mask[i]=1：单元质心在圆内；
         * - mask[i]=0：否则。
         * 说明：
         *   1) 该版本用质心判定，极其稳定；2D 下你可以把它当作面积权重的“第一近似”（后续可升级为求积占比）。
         *   2) 3D 时此函数同样可用（忽略 z，仅用 x,y 距离），例如圆柱形竖井在水平截面上的投影选择。
         */
        inline bool build_Well_mask_circle_centroid(
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const wellCircle& well,
            const std::string& out_name = "inj_mask",
			bool verbose = true
        ) {
            ::Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const size_t nC = cells.size();

            // 输出字段：体标量（每个单元一个值）
            auto maskPtr = reg.getOrCreate<::volScalarField>(out_name, nC, 0.0);
            auto& mask = maskPtr->data;
            mask.assign(nC, 0.0);

            const double r2 = std::max(well.r, 0.0) * std::max(well.r, 0.0);
            size_t n_inside = 0;
            double covered_area = 0.0; // 2D: 总覆盖面积近似（仅用于粗略检查）

            for (size_t i = 0; i < nC; ++i) {
                const ::Cell& c = cells[i];
                const double dx = c.center.m_x - well.cx;
                const double dy = c.center.m_y - well.cy;
                const double d2 = dx * dx + dy * dy;

                if (d2 <= r2) {
                    mask[i] = 1.0;
                    ++n_inside;
                    // 2D 时 c.volume == 单元面积；3D 你也可用它统计井柱投影的“带权个数”
                    covered_area += c.volume;
                }
                else {
                    mask[i] = 0.0;
                }
            }

            if (verbose) {
                std::cout << "[InjectionMask] '" << out_name << "' built by centroid rule: "
                    << n_inside << "/" << nC << " cells selected. "
                    << "Approx. covered area = " << covered_area << "\n";
            }
            return true;
        }


        inline bool build_Injec_mask_circle_centroid
        (
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const wellCircle& inj,            // 圆心(cx,cy)、半径r
            const std::string& out_name = "inj_mask",
            bool verbose = true
        )
        {
            return build_Well_mask_circle_centroid(mgr, reg, inj, "inj_mask", true);
        }

        inline bool build_Prod_mask_circle_centroid
        (
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const wellCircle& prod,            // 圆心(cx,cy)、半径r
            const std::string& out_name = "prod_mask",
            bool verbose = true

        )
        {
            return build_Well_mask_circle_centroid(mgr, reg, prod, "prod_mask", true);

        }

    } // namespace SourceTerm
} // namespace FVM
