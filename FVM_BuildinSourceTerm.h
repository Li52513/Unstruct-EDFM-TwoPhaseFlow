#pragma once
#include <string>
#include <memory>
#include <cmath>
#include <iostream>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "VolField.h"
#ifdef _WIN32
#  ifndef NOMINMAX
#    define NOMINMAX   // 禁掉 Windows 的 min/max 宏
#  endif
#endif
#include <algorithm>

namespace FVM {
    namespace SourceTerm {

        // 工具：按装配器约定生成字段名
        inline std::string name_a_src(const std::string& tag) { return "a_src_" + tag; }
        inline std::string name_b_src(const std::string& tag) { return "b_src_" + tag; }

        /**
         * @brief 压力方程 Su–Sp 井源（Peaceman PI → 质量源）
         *        q_m = PI*(p_inj - p)
         *   =>   a_src = (+PI) * V,   b_src = (PI*p_inj) * V
         * @param tag 与装配时的 OperatorFieldNames 一致，比如压力用 "p_g"
         */
        inline bool build_pressure_sources_from_PI(
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const std::string& PI_name,          // "PI"
            const std::string& p_field_name,     // "p_g"（此处仅读取，为与接口统一）
            double p_inj,
            const std::string& tag,              // ==> 生成 a_src_<tag> / b_src_<tag>
            bool verbose,
            bool accumulate = false
        ) {
            (void)p_field_name; // 当前 a,b 不直接用 p；保留接口方便扩展
            ::Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const size_t nC = cells.size();

            // 读 PI
            auto PIF = reg.get<::volScalarField>(PI_name);
            if (!PIF) {
                std::cerr << "[Source:p] field '" << PI_name << "' not found; sources=0.\n";
                auto aZ = reg.getOrCreate<::volScalarField>(name_a_src(tag), nC, 0.0);
                auto bZ = reg.getOrCreate<::volScalarField>(name_b_src(tag), nC, 0.0);
                aZ->data.assign(nC, 0.0); bZ->data.assign(nC, 0.0);
                return true;
            }
            const auto& PI = PIF->data;

            // 输出 a,b
            auto aPtr = reg.getOrCreate<::volScalarField>(name_a_src(tag), nC, 0.0);
            auto bPtr = reg.getOrCreate<::volScalarField>(name_b_src(tag), nC, 0.0);
            auto& a = aPtr->data;
            auto& b = bPtr->data;
            
            if (!accumulate) { a.assign(nC, 0.0); b.assign(nC, 0.0); }

            size_t nAct = 0;
            for (size_t i = 0; i < nC; ++i) {
                const double V = std::max(0.0, cells[i].volume);  // 2D: 面积；3D: 体积
                if (V <= 0.0) continue;
                const double PIi = PI[i];
                if (!(PIi > 0.0)) continue;

                a[i] += PIi * V;          // 加到对角
                b[i] += PIi * p_inj * V;  // 加到 RHS
                ++nAct;
            }

            if (verbose) {
                std::cout << "[Source:p] built a='" << name_a_src(tag)
                    << "', b='" << name_b_src(tag) << "'  active=" << nAct << "\n";
            }
            return true;
        }

        /**
         * @brief 温度/能量方程 Su–Sp 井源（半隐式，冻结 q_m*）
         *        q_m* = PI*(p_inj - p) 来自最新压力场
         *        q_E  = q_m* c_p* (T_inj - T)
         *   =>   a_src = ( q_m* c_p* ) * V,   b_src = ( q_m* c_p* T_inj ) * V
         * @param cp_const 若传 cp_field=nullptr，则使用该常量
         * @param tag 与装配时一致，温度用 "T"
         */
        inline bool build_temperature_sources_from_PI
(
            ::MeshManager& mgr,
            ::FieldRegistry& reg,
            const std::string& PI_name,          // "PI"
            const std::string& p_field_name,     // 用于 q_m* 的 p，如 "p_g"
            double p_inj,
            double cp_const,
            const ::volScalarField* cp_field,    // 若为 nullptr 则用 cp_const
            double T_inj,
            const std::string& tag,              // ==> 生成 a_src_<tag> / b_src_<tag>
            bool verbose,
            bool accumulate = false
        ) {
            ::Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const size_t nC = cells.size();

            // 读字段：PI、p
            auto PIF = reg.get<::volScalarField>(PI_name);
            auto pF = reg.get<::volScalarField>(p_field_name);
            if (!PIF || !pF) {
                std::cerr << "[Source:T] missing '" << PI_name << "' or '" << p_field_name
                    << "'; sources=0.\n";
                auto aZ = reg.getOrCreate<::volScalarField>(name_a_src(tag), nC, 0.0);
                auto bZ = reg.getOrCreate<::volScalarField>(name_b_src(tag), nC, 0.0);
                aZ->data.assign(nC, 0.0); bZ->data.assign(nC, 0.0);
                return true;
            }
            const auto& PI = PIF->data;
            const auto& p = pF->data;

            // 输出 a,b
            auto aPtr = reg.getOrCreate<::volScalarField>(name_a_src(tag), nC, 0.0);
            auto bPtr = reg.getOrCreate<::volScalarField>(name_b_src(tag), nC, 0.0);
            auto& a = aPtr->data;
            auto& b = bPtr->data;
            if (!accumulate) { a.assign(nC, 0.0); b.assign(nC, 0.0); }

            size_t nAct = 0;
            for (size_t i = 0; i < nC; ++i) {
                const double V = std::max(0.0, cells[i].volume);
                if (V <= 0.0) continue;

                const double qm = PI[i] * (p_inj - p[i]);     // q_m*
                if (qm == 0.0) continue;

                const double cp = (cp_field ? cp_field->data[i] : cp_const);
                const double K = qm * cp;                    // 系数
                a[i] += K * V;                                 // 对角
                b[i] += K * T_inj * V;                         // RHS
                ++nAct;
            }

            if (verbose) {
                std::cout << "[Source:T] built a='" << name_a_src(tag)
                    << "', b='" << name_b_src(tag) << "'  active=" << nAct
                    << "  (semi-implicit)\n";
            }
            return true;
        }


    } // namespace SourceTerm
} // namespace FVM
