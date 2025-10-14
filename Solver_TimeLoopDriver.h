#pragma once
#include <functional>
#include <iostream>
#include <string>
#include <sstream>   // ★ 新增
#include <iomanip>   // ★ 新增
#include "Solver_TimeLoopSkeleton.h"  // outerIter_OneStep_singlePhase(...)
#include "Solver_PostChecks.h"        // PostChecks::export_* / dumpCOO_to_matrix_market
#include "TecplotExporter.h"          // TecplotExport::export_pT_tecplot_after_step
#include "Solver_AssemblerCOO.h"      // assemblePressure_* / assembleTemperature_*

// 可选导出回调：用户自定义额外输出
using WriteCallback = std::function<void(int step, double time)>;

// 小工具：根据相态返回压力场名
inline const char* pNameForPhase(const std::string& phase) {
    return (phase == "CO2" || phase == "co2") ? "p_g" : "p_w";
}

/**
 * @brief 单相（CO2 / water）渗流-传热：多步时间推进驱动
 */
inline bool runTransient_singlePhase(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const PressureBCAdapter& Pbc,
    const TemperatureBCAdapter& Tbc,
    const GravUpwind& gu,
    const RockDefaults& rock,
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    const std::string& phase = "CO2",
    int writeEvery = 0,
    WriteCallback onWrite = nullptr,
    bool exportCSV = false,
    bool exportTXT = true,
    bool exportMM = false,
    bool exportTecplotP = false,
    bool exportTecplotT = false
) {
    if (nSteps <= 0 || dt <= 0.0) {
        std::cerr << "[runTransient] invalid nSteps/dt.\n";
        return false;
    }

    const char* pField = pNameForPhase(phase);
    std::string outFolderCSV = std::string("out_pT_") + phase;
    std::string outFolderTXT = std::string("out_txt_") + phase;
    std::string mmFolder = "mm";

    // 初始时刻 t=0, step=0 的输出（可选）
    if (exportTXT) {
        PostChecks::export_pT_txt_after_step(mgr, reg, pField, "T",
            /*step=*/0, /*t=*/0.0, outFolderTXT);
    }
    if (exportCSV) {
        PostChecks::export_pT_csv_after_step(mgr, reg, pField, "T",
            /*step=*/0, /*t=*/0.0, outFolderCSV);
    }
    if (exportTecplotP || exportTecplotT) {
        TecplotExport::TecplotOptions opt;
        opt.folder = std::string("out_tec_") + phase;
        opt.filePrefix = phase;
        if (exportTecplotP) {
            TecplotExport::export_cellField_to_tecplot_after_step(
                mgr, reg, freg, opt,
                /*cellField*/ pField,
                /*tmpFace*/   std::string(pField) + "_face_tmp",   // ★ 修正拼接
                /*varName*/   "P",
                /*step=*/0, /*t=*/0.0,
                /*PBC*/ &Pbc, /*TBC*/ nullptr);
        }
        if (exportTecplotT) {
            TecplotExport::export_cellField_to_tecplot_after_step(
                mgr, reg, freg, opt,
                /*cellField*/ "T",
                /*tmpFace*/   "T_face_tmp",
                /*varName*/   "T",
                /*step=*/0, /*t=*/0.0,
                /*PBC*/ nullptr, /*TBC*/ &Tbc);
        }
    }

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step) 
    {
        const int step1 = step + 1;
        t += dt;

        // ——单步推进（含：startTimeStep/外迭代/提交 n+1）——
        bool ok = outerIter_OneStep_singlePhase(
            mgr, reg, freg, ppm, Pbc, Tbc, gu, rock, dt, ctrl, phase
        );
        if (!ok) {
            std::cerr << "[runTransient] step " << step1 << " failed.\n";
            return false;
        }

        // ——每步导出：TXT/CSV——
        if (exportTXT) {
            PostChecks::export_pT_txt_after_step(mgr, reg, pField, "T", step1, t, outFolderTXT);
        }
        if (exportCSV) {
            PostChecks::export_pT_csv_after_step(mgr, reg, pField, "T", step1, t, outFolderCSV);
        }


        const bool doThisStep = (writeEvery <= 0) || (step1 % writeEvery == 0);
        if ((exportTecplotP || exportTecplotT) && doThisStep) {
            TecplotExport::TecplotOptions opt;
            opt.folder = std::string("out_tec_") + phase;  // 目录可自定义
            opt.filePrefix = phase;
            opt.precision = 16;

            if (exportTecplotP) {
                const bool okP = TecplotExport::export_cellField_to_tecplot_after_step_cell2face2node(
                    mgr, reg, freg, opt,
                    /*cellField*/ pField,
                    /*tmpFace */ std::string(pField) + "_face_tmp",
                    /*varName */ "P",
                    step1, t,
                    /*PBC*/ &Pbc, /*TBC*/ nullptr);
                if (!okP) std::cerr << "[runTransient] Tecplot(P) export failed at step " << step1 << "\n";
            }
            if (exportTecplotT) {
                const bool okT = TecplotExport::export_cellField_to_tecplot_after_step_cell2face2node(
                    mgr, reg, freg, opt,
                    /*cellField*/ "T",
                    /*tmpFace */ "T_face_tmp",
                    /*varName */ "T",
                    step1, t,
                    /*PBC*/ nullptr, /*TBC*/ &Tbc);
                if (!okT) std::cerr << "[runTransient] Tecplot(T) export failed at step " << step1 << "\n";
            }
        }
        //if ((exportTecplotP || exportTecplotT) && doThisStep) {
        //    TecplotExport::TecplotOptions opt;
        //    opt.folder = std::string("out_tec_") + phase;
        //    opt.filePrefix = phase;

        //    if (exportTecplotP) {
        //        const bool okP = TecplotExport::export_cellField_to_tecplot_after_step(
        //            mgr, reg, freg, opt,
        //            /*cellField*/ pField,
        //            /*tmpFace*/   std::string(pField) + "_face_tmp", // ★ 修正拼接
        //            /*varName*/   "P",
        //            step1, t,
        //            /*PBC*/ &Pbc, /*TBC*/ nullptr
        //        );
        //        if (!okP) std::cerr << "[runTransient] Tecplot(P) export failed at step " << step1 << "\n";
        //    }
        //    if (exportTecplotT) {
        //        const bool okT = TecplotExport::export_cellField_to_tecplot_after_step(
        //            mgr, reg, freg, opt,
        //            /*cellField*/ "T",
        //            /*tmpFace*/   "T_face_tmp",
        //            /*varName*/   "T",
        //            step1, t,
        //            /*PBC*/ nullptr, /*TBC*/ &Tbc
        //        );
        //        if (!okT) std::cerr << "[runTransient] Tecplot(T) export failed at step " << step1 << "\n";
        //    }
        //}



        // ——可选：导出 MatrixMarket（重新装配一次用于外部核查）——
        if (exportMM && ( (writeEvery <= 0) || (step1 % writeEvery == 0) )) {  // ★ 修正优先级
            SparseSystemCOO sysP, sysT;

            if (phase == "CO2" || phase == "co2") {
                assemblePressure_CO2_singlePhase_COO(mgr, reg, freg, &sysP);
            } else {
                assemblePressure_water_singlePhase_COO(mgr, reg, freg, &sysP);
            }
            sysP.compressInPlace(0.0);

            assembleTemperature_singlePhase_COO(
                mgr, reg, freg,
                /*a_f_Diff_T*/ "a_f_Diff_T",
                /*s_f_Diff_T*/ "s_f_Diff_T",
                /*aPP_conv*/   "aPP_conv",
                /*aPN_conv*/   "aPN_conv",
                /*bP_conv */   "bP_conv",
                /*a_time*/     "aC_time_T",
                /*b_time*/     "bC_time_T",
                &sysT
            );
            sysT.compressInPlace(0.0);

            std::ostringstream fnA, fnB, fnAT, fnBT;
            fnA  << mmFolder << "/A_P_" << phase << "_step_" << std::setw(6) << std::setfill('0') << step1 << ".mtx";
            fnB  << mmFolder << "/b_P_" << phase << "_step_" << std::setw(6) << std::setfill('0') << step1 << ".txt";
            fnAT << mmFolder << "/A_T_" << phase << "_step_" << std::setw(6) << std::setfill('0') << step1 << ".mtx";
            fnBT << mmFolder << "/b_T_" << phase << "_step_" << std::setw(6) << std::setfill('0') << step1 << ".txt";

            PostChecks::dumpCOO_to_matrix_market(sysP,  fnA.str(),  fnB.str(),  /*compress*/false);
            PostChecks::dumpCOO_to_matrix_market(sysT, fnAT.str(), fnBT.str(), /*compress*/false);

            // 装配体检（可视化到控制台）
            auto RP = PostChecks::reportAssembly(sysP, /*force_compress*/false);
            auto RT = PostChecks::reportAssembly(sysT, /*force_compress*/false);
            PostChecks::printAssemblyReport(RP, "P");
            PostChecks::printAssemblyReport(RT, "T");
        }

        // ——可选自定义回调——
        if (writeEvery > 0 && (step1 % writeEvery == 0) && onWrite) {
            onWrite(step1, t);
            std::cout << "[runTransient] wrote step " << step1 << " at t=" << t << "\n";
        }





    }

    if (writeEvery > 0 && onWrite && (nSteps % writeEvery != 0)) {
        onWrite(nSteps, nSteps * dt);
        std::cout << "[runTransient] wrote final step " << nSteps
                  << " at t=" << (nSteps * dt) << "\n";
    }
    return true;
}
