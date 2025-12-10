#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FVM_WellDOF_TwoPhase.h"
#include "CapRelPerm.h"            
#include "TwoPhaseWells_StrictRate.h" 

namespace IMPES_Iteration {

    struct WellMassReport
    {
        std::string name;
        WellDOF_TwoPhase::Role role;
        WellDOF_TwoPhase::Mode mode;
        int hostCellId;   ///< 目前只考虑每井一个完井单元的简单情形
        double p_bh;      ///< 井底压力 [Pa]
        double Mw_tot;    ///< 该井总水相质量流量 [kg/s]（注入为+，采出为-）
        double Mg_tot;    ///< 该井总气相质量流量 [kg/s]
        double fw_mass;   ///< Mw_tot / (Mw_tot+Mg_tot)，若总流量≈0用0处理
        double Tin;       ///< 对注入井：WellConfig_TwoPhase::Tin；对采出井可填占位值
    };
    /// \brief 边界质量通量的简单汇总信息（按 outflow 汇总）.
    struct BoundaryMassReport
    {
        double Mw_out = 0.0;  ///< 边界总流出水相质量流量 [kg/s] (owner->外侧 & outflow)
        double Mg_out = 0.0;  ///< 边界总流出气相质量流量 [kg/s]
        double fw_out = 0.0;  ///< 边界流出混相的质量分数 f_w = Mw_out / (Mw_out + Mg_out)

        double Mw_in = 0.0;   ///< (可选) 边界总流入水相质量流量 [kg/s] (外侧->owner & inflow)
        double Mg_in = 0.0;   ///< (可选) 边界总流入气相质量流量 [kg/s]
        double fw_in = 0.0;   ///< (可选) 边界流入混相 f_w
    };

    /// \brief 单个时间步的井 + 边界质量流量记录（一行CSV）
    struct WellBoundaryCSVRow
    {
        double time = 0.0;

        // 注入井
        double pbh_inj = 0.0;
        double pcell_inj = 0.0;
        double Mw_inj = 0.0;
        double Mg_inj = 0.0;
        double fw_inj = 0.0;
        double Tin_inj = 0.0;

        // 采出井
        double pbh_prod = 0.0;
        double pcell_prod = 0.0;
        double Mw_prod = 0.0;
        double Mg_prod = 0.0;
        double fw_prod = 0.0;
        double Tout_prod = 0.0; ///< 先占位，可填 NaN 或 -1

        // 边界流出
        double Mw_bdry_out = 0.0;
        double Mg_bdry_out = 0.0;
        double fw_bdry_out = 0.0;
    };

    /// \brief 如果文件不存在，则写入一行表头；如果已存在则不动。
    inline void ensureWellBoundaryCSVHeader(const std::string& filename)
    {
        std::ifstream fin(filename);
        if (fin.good()) return; // 已存在，不重复写

        std::ofstream fout(filename);
        if (!fout.good()) {
            std::cerr << "[WellCSV] Cannot open file for header: " << filename << "\n";
            return;
        }

        fout << "time,"
            // inj
            << "pbh_inj,pcell_inj,Mw_inj,Mg_inj,fw_inj,Tin_inj,"
            // prod
            << "pbh_prod,pcell_prod,Mw_prod,Mg_prod,fw_prod,Tout_prod,"
            // boundary
            << "Mw_bdry_out,Mg_bdry_out,fw_bdry_out"
            << "\n";
    }

    /// \brief 追加一行记录到 CSV 文件（假定表头已经存在或稍后会手动调用 ensureHeader）
    inline void appendWellBoundaryCSVRow(
        const std::string& filename,
        const WellBoundaryCSVRow& row,
        int time_precision = 8,
        int value_precision = 8
    )
    {
        std::ofstream fout(filename, std::ios::app);
        if (!fout.good()) {
            std::cerr << "[WellCSV] Cannot open file for append: " << filename << "\n";
            return;
        }

        fout << std::setprecision(time_precision) << std::scientific
            << row.time << ",";

        auto printVal = [&](double v, bool last = false)
            {
                fout << std::setprecision(value_precision) << std::scientific << v;
                if (!last) fout << ",";
            };

        // inj
        printVal(row.pbh_inj);
        printVal(row.pcell_inj);
        printVal(row.Mw_inj);
        printVal(row.Mg_inj);
        printVal(row.fw_inj);
        printVal(row.Tin_inj);

        // prod
        printVal(row.pbh_prod);
        printVal(row.pcell_prod);
        printVal(row.Mw_prod);
        printVal(row.Mg_prod);
        printVal(row.fw_prod);
        printVal(row.Tout_prod);

        // boundary
        printVal(row.Mw_bdry_out);
        printVal(row.Mg_bdry_out, /*last=*/false);
        fout << std::setprecision(value_precision) << std::scientific << row.fw_bdry_out;

        fout << "\n";
    }

    /// \brief 从 wellReports 向量中挑出单个注入井 / 采出井，并根据当前场变量补齐 p_cell。
    inline void buildWellBoundaryRowFromReports(
        MeshManager& mgr,
        FieldRegistry& reg,
        const std::vector<WellMassReport>& wellReports,
        const std::string& p_field,   ///< 压力场名，一般是 pressureCtrl.assembly.pressure_field
        const BoundaryMassReport& bdry,
        double time,
        WellBoundaryCSVRow& row
    )
    {
        row.time = time;

        auto& mesh = mgr.mesh();
        const auto& id2idx = mesh.getCellId2Index();

        auto getCellPressure = [&](int cellId)->double
            {
                if (cellId < 0) return 0.0;
                const auto pF = reg.get<volScalarField>(p_field.c_str());
                if (!pF) return 0.0;
                const size_t i = id2idx.at(cellId);
                return (*pF)[i];
            };

        // 注入井 / 采出井：这里假定每类最多 1 口（当前测试是1注1采）
        for (const auto& wrep : wellReports)
        {
            if (wrep.role == WellDOF_TwoPhase::Role::Injector)
            {
                row.pbh_inj = wrep.p_bh;
                row.pcell_inj = getCellPressure(wrep.hostCellId);
                row.Mw_inj = wrep.Mw_tot;
                row.Mg_inj = wrep.Mg_tot;
                row.fw_inj = wrep.fw_mass;
                row.Tin_inj = wrep.Tin;
            }
            else if (wrep.role == WellDOF_TwoPhase::Role::Producer)
            {
                row.pbh_prod = wrep.p_bh;
                row.pcell_prod = getCellPressure(wrep.hostCellId);
                row.Mw_prod = wrep.Mw_tot;
                row.Mg_prod = wrep.Mg_tot;
                row.fw_prod = wrep.fw_mass;
                row.Tout_prod = -1.0; // 当前未耦合温度，可后续替换
            }
        }

        // 边界
        row.Mw_bdry_out = bdry.Mw_out;
        row.Mg_bdry_out = bdry.Mg_out;
        const double denom = std::max(std::abs(bdry.Mw_out) + std::abs(bdry.Mg_out), 1e-30);
        row.fw_bdry_out = (denom > 0.0 ? bdry.Mw_out / denom : 0.0);
    }

    /// \brief 诊断边界上的流入/流出总质量通量，并打印对应的 f_w（若存在），
    ///        同时汇总边界总流出/流入质量（用于守恒检查/CSV输出）.
    ///
    /// 约定：
    ///  - 正通量 mf_total > 0: owner -> 外侧（outflow）
    ///  - 负通量 mf_total < 0: 外侧 -> owner（inflow）
    ///  - 只在 |mf_total| > flux_sign_epsilon 时认为是真正的流入/流出
    ///
    /// 参数：
    ///  - total_mass_flux_name : total Darcy mass flux 面场名（通常来自 FaceMassRateConfig::total_mass_flux）
    ///  - fw_face_name         : f_w 面场名（通常来自 FluxSplitConfig::fractional_flow_face，例 "fw_face"）
    ///  - flux_sign_epsilon    : 判定流向的阈值，应与 splitTwoPhaseMassFlux 中一致
    ///  - report               : (可选 out) 边界质量通量汇总结果
    inline void diagnoseBoundaryFluxWithFw(
        MeshManager& mgr,
        FaceFieldRegistry& freg,
        const std::string& total_mass_flux_name,
        const std::string& fw_face_name,
        double flux_sign_epsilon,
        BoundaryMassReport* report = nullptr
    )
    {
        auto mf_total = freg.get<faceScalarField>(total_mass_flux_name.c_str());
        if (!mf_total) {
            std::cerr << "[FluxCheck][boundary] total mass flux field '"
                << total_mass_flux_name << "' not found.\n";
            return;
        }

        std::shared_ptr<faceScalarField> fw_face = nullptr;
        if (!fw_face_name.empty()) {
            fw_face = freg.get<faceScalarField>(fw_face_name.c_str());
            if (!fw_face) {
                std::cerr << "[FluxCheck][boundary] fw face field '"
                    << fw_face_name << "' not found.\n";
            }
        }

        const Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();

        BoundaryMassReport dummy;       // 若用户不关心汇总，用本地变量吞掉
        BoundaryMassReport& rep = report ? *report : dummy;
        rep = BoundaryMassReport{};     // 清零

        std::cout << "[FluxCheck][boundary] listing inflow/outflow faces with mf_total & fw\n";

        for (const auto& F : faces)
        {
            if (!F.isBoundary()) continue;

            const int iF = F.id - 1;
            if (iF < 0 || iF >= static_cast<int>(mf_total->data.size())) continue;

            const double flux = (*mf_total)[iF];   // [kg/s]
            if (std::abs(flux) <= flux_sign_epsilon) continue; // 视为近似无流

            const bool isInflow = (flux < -flux_sign_epsilon); // 边界 -> owner
            const bool isOutflow = (flux > flux_sign_epsilon); // owner -> 边界

            if (!isInflow && !isOutflow) continue;

            double fw_val = -1.0;
            if (fw_face && iF < static_cast<int>(fw_face->data.size())) {
                fw_val = (*fw_face)[iF];
            }

            // ---- 1) 打印单个面的信息（原有功能） ----
            std::cout << "  Face " << F.id
                << "  type = " << (isInflow ? "INFLOW " : "OUTFLOW")
                << ", mf_total = " << flux;

            if (fw_face) {
                std::cout << ", fw = " << fw_val;
            }
            else {
                std::cout << ", fw = (not available)";
            }

            std::cout << "\n";

            // ---- 2) 进行边界质量流量汇总 ----
            if (fw_face) {
                // 有两相信息：按质量分流函数拆 Mw, Mg
                const double fw_mass = fw_val;
                const double Mw = flux * fw_mass;
                const double Mg = flux * (1.0 - fw_mass);

                if (isOutflow) {
                    rep.Mw_out += Mw;
                    rep.Mg_out += Mg;
                }
                else { // inflow
                    rep.Mw_in += Mw;
                    rep.Mg_in += Mg;
                }
            }
            else {
                // 没有 f_w 面场：可以选择全部算到 Mw_out/Mw_in，或只统计总质量不拆相.
                // 这里示例：全部算作水相（可按需要修改逻辑）
                if (isOutflow) {
                    rep.Mw_out += flux;
                }
                else {
                    rep.Mw_in += flux;
                }
            }
        }

        // ---- 3) 汇总出 f_w_out / f_w_in（可选）----
        {
            const double denom_out = std::abs(rep.Mw_out) + std::abs(rep.Mg_out);
            rep.fw_out = (denom_out > 0.0) ? (rep.Mw_out / denom_out) : 0.0;

            const double denom_in = std::abs(rep.Mw_in) + std::abs(rep.Mg_in);
            rep.fw_in = (denom_in > 0.0) ? (rep.Mw_in / denom_in) : 0.0;
        }

        // 你也可以在这里追加一行整体 summary 的打印，便于快速看：
        std::cout << "[FluxCheck][boundary-summary] "
            << "Mw_out = " << rep.Mw_out << " kg/s, "
            << "Mg_out = " << rep.Mg_out << " kg/s, "
            << "fw_out = " << rep.fw_out << " | "
            << "Mw_in = " << rep.Mw_in << " kg/s, "
            << "Mg_in = " << rep.Mg_in << " kg/s, "
            << "fw_in = " << rep.fw_in
            << "\n";
    }



}