#pragma once
#include <functional>
#include <iostream>
#include <string>
#include <sstream>   // ★ 新增
#include <iomanip>   // ★ 新增
#include <filesystem>   
#include "Solver_TimeLoopSkeleton.h"  // outerIter_OneStep_singlePhase(...)
#include "Solver_PostChecks.h"        // PostChecks::export_* / dumpCOO_to_matrix_market
#include "TecplotExporter.h"          // TecplotExport::export_pT_tecplot_after_step
#include "Solver_AssemblerCOO.h"      // assemblePressure_* / assembleTemperature_*
#include "PostProcess_.h"


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
    if (exportTecplotP || exportTecplotT)
    {
        TecplotExport::TecplotOptions opt;
        opt.folder = std::string("out_tec_") + phase;
        opt.filePrefix = phase;
        if (exportTecplotP) 
        {
            TecplotExport::export_cellField_to_tecplot_after_step
            (
                mgr, reg, freg, opt,
                /*cellField*/ pField,
                /*tmpFace*/   std::string(pField) + "_face_tmp",   // ★ 修正拼接
                /*varName*/   "P",
                /*step=*/0, /*t=*/0.0,
                /*PBC*/ &Pbc, /*TBC*/ nullptr
            );
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
        //if ((exportTecplotP || exportTecplotT) && doThisStep) {
        //    TecplotExport::TecplotOptions opt;
        //    opt.folder = std::string("out_tec_") + phase;  // 目录可自定义
        //    opt.filePrefix = phase;
        //    opt.precision = 16;

        //    if (exportTecplotP) {
        //        const bool okP = TecplotExport::export_cellField_to_tecplot_after_step_cell2face2node(
        //            mgr, reg, freg, opt,
        //            /*cellField*/ pField,
        //            /*tmpFace */ std::string(pField) + "_face_tmp",
        //            /*varName */ "P",
        //            step1, t,
        //            /*PBC*/ &Pbc, /*TBC*/ nullptr);
        //        if (!okP) std::cerr << "[runTransient] Tecplot(P) export failed at step " << step1 << "\n";
        //    }
        //    if (exportTecplotT) {
        //        const bool okT = TecplotExport::export_cellField_to_tecplot_after_step_cell2face2node(
        //            mgr, reg, freg, opt,
        //            /*cellField*/ "T",
        //            /*tmpFace */ "T_face_tmp",
        //            /*varName */ "T",
        //            step1, t,
        //            /*PBC*/ nullptr, /*TBC*/ &Tbc);
        //        if (!okT) std::cerr << "[runTransient] Tecplot(T) export failed at step " << step1 << "\n";
        //    }
        //}
        if ((exportTecplotP || exportTecplotT) && doThisStep) {
            TecplotExport::TecplotOptions opt;
            opt.folder = std::string("out_tec_") + phase;
            opt.filePrefix = phase;

            if (exportTecplotP) {
                const bool okP = TecplotExport::export_cellField_to_tecplot_after_step(
                    mgr, reg, freg, opt,
                    /*cellField*/ pField,
                    /*tmpFace*/   std::string(pField) + "_face_tmp", // ★ 修正拼接
                    /*varName*/   "P",
                    step1, t,
                    /*PBC*/ &Pbc, /*TBC*/ nullptr
                );
                if (!okP) std::cerr << "[runTransient] Tecplot(P) export failed at step " << step1 << "\n";
            }
            if (exportTecplotT) {
                const bool okT = TecplotExport::export_cellField_to_tecplot_after_step(
                    mgr, reg, freg, opt,
                    /*cellField*/ "T",
                    /*tmpFace*/   "T_face_tmp",
                    /*varName*/   "T",
                    step1, t,
                    /*PBC*/ nullptr, /*TBC*/ &Tbc
                );
                if (!okT) std::cerr << "[runTransient] Tecplot(T) export failed at step " << step1 << "\n";
            }
        }



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


inline bool runTransient_test_constProperties_singlePhase_CO2_T_diffusion
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const GravUpwind& gu,
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    int writeEvery = 1,                         // 每步/每隔几步输出
    const std::string& outPrefix = "T_CO2_diff" // 输出前缀
)
{

    //创建并打开文件夹
#if __cplusplus >= 201703L
        // 保证父目录存在（只用一次，没必要每步都建）
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefix).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient] cannot create directory for: " << outPrefix << "\n";
        return false;
    }
#endif

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step)
    {
        
        //时间步推进
        const int step1 = step + 1;
        t += dt;
        
        if (!startTimeStep_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, "T", "T_old", "T_prev")) return false;

		// ——单步推进（含：startTimeStep/外迭代/提交 n+1）——
        bool ok = outerIter_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg, freg, ppm, Tbc, gu, dt, ctrl);
		
        if (!ok) {
			std::cerr << "[runTransient] step " << step1 << " failed.\n";
			return false;
		}

        //导出该时间步结果
        const bool doThisStep = (writeEvery <= 0) || (step1 % writeEvery == 0);
		// ——每步导出：TXT——
        if (doThisStep)
        {
            const std::vector<Vector> gradT  = computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "T", /*gradSmoothIters=*/0);


            std::ostringstream fn;
            fn << outPrefix << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPlt = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ &Tbc, /*Pbc*/ nullptr,
                /*cell*/ "T",
                /*face*/ "T_face_tmp",
                /*gradBuf*/ &gradT,       // ★ 零拷贝传梯度缓冲
                /*out*/  fn.str()
            );
            if (!okPlt) {
                std::cerr << "[Transient(T-diff)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }
        }
    }
	return true;
}


// 测试案例：2D-常物性-单相-CO2-热扩散问题
inline bool runTransient_test_varyingProperties_singlePhase_CO2_T_diffusion
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const GravUpwind& gu,
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    int writeEvery = 1,                         // 每步/每隔几步输出
    const std::string& outPrefix = "T_CO2_diff" // 输出前缀
)
{
    //创建并打开文件夹
#if __cplusplus >= 201703L
        // 保证父目录存在（只用一次，没必要每步都建）
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefix).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient] cannot create directory for: " << outPrefix << "\n";
        return false;
    }
#endif

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step)
    {

        //时间步推进
        const int step1 = step + 1;
        t += dt;

        if (!startTimeStep_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, "T", "T_old", "T_prev")) return false;

        // ——单步推进（含：startTimeStep/外迭代/提交 n+1）——
		bool ok = outerIter_test_varyProperties_singlePhase_CO2_T_diffusion(mgr, reg, freg, ppm, Tbc, gu, dt, ctrl);

        if (!ok) {
            std::cerr << "[runTransient] step " << step1 << " failed.\n";
            return false;
        }

        //导出该时间步结果
        const bool doThisStep = (writeEvery <= 0) || (step1 % writeEvery == 0);
        // ——每步导出：TXT——
        if (doThisStep)
        {
            const std::vector<Vector> gradT = computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "T", /*gradSmoothIters=*/0);


            std::ostringstream fn;
            fn << outPrefix << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPlt = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ &Tbc, /*Pbc*/ nullptr,
                /*cell*/ "T",
                /*face*/ "T_face_tmp",
                /*gradBuf*/ &gradT,       // ★ 零拷贝传梯度缓冲
                /*out*/  fn.str()
            );
            if (!okPlt) {
                std::cerr << "[Transient(T-diff)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }
        }
    }
    return true;

}


// 测试案例：2D-常物性-单相-CO2-达西流
// ================== Transient driver: 单相·常物性·CO2·达西渗流（压力） ==================
inline bool runTransient_test_constProperties_singlePhase_CO2_p_flow
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const  PressureBCAdapter& Pbc,
    const Vector& g,                       // 例如 {0,-9.81,0}
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    int writeEvery = 1,                    // 每步/每隔几步输出
    const std::string& outPrefix = "P_CO2_flow" // 输出前缀（可含路径）
)
{
#if __cplusplus >= 201703L
    // 保证父目录存在
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefix).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient(P)] cannot create directory for: " << outPrefix << "\n";
        return false;
    }
#endif

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step) {
        const int step1 = step + 1;
        t += dt;

        if (!startTimeStep_scalar (mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev")) return false;

        // —— 单步外迭代—— //
        bool ok = outerIter_constProperties_singlePhase_CO2_Pressure(
            mgr, reg, freg, ppm, Pbc, g, dt, ctrl
        );
        if (!ok) {
            std::cerr << "[runTransient(P)] step " << step1 << " failed.\n";
            return false;
        }

        // —— 导出 —— //
        const bool doThisStep = (writeEvery <= 0) || (step1 % writeEvery == 0);
        if (doThisStep) {
            const std::vector<Vector> gradP =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "p_g", /*gradSmoothIters=*/0);

            std::ostringstream fn;
            fn << outPrefix << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPlt = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ nullptr, /*Pbc*/ &Pbc,
                /*cell*/ "p_g",
                /*face*/ "p_g_face_tmp",
                /*gradBuf*/ &gradP,
                /*out*/ fn.str()
            );
            if (!okPlt) {
                std::cerr << "[Transient(P)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }
        }
    }
    return true;
}



//测试案例：2D-常物性-单相-CO2-达西渗流-传热耦合问题
inline bool runTransient_constProperties_singlePhase_CO2_T_H
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const PressureBCAdapter& Pbc,
    const Vector& g,                   // 例如 {0,-9.81,0}
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    // —— 输出控制（可分别设置频率与路径）——
    int writeEveryP = 1,
    int writeEveryT = 1,
    const std::string& outPrefixP = "P_CO2_pTH/p",
    const std::string& outPrefixT = "T_CO2_pTH/t"
)
{
    if (dt <= 0.0) { std::cerr << "[runTransient(p+T)] invalid dt.\n"; return false; }

#if __cplusplus >= 201703L
    // 保证父目录存在
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefixP).parent_path());
        std::filesystem::create_directories(std::filesystem::path(outPrefixT).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient(p+T)] cannot create directory for: "
            << outPrefixP << " or " << outPrefixT << "\n";
        return false;
    }
#endif

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step) {
        const int step1 = step + 1;
        t += dt;

        // —— 时间层切换：压力与温度 —— //
        if (!startTimeStep_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev")) return false;
        if (!startTimeStep_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev")) return false;

        // —— 单步外迭代：p→mf→T（常物性）—— //
        bool ok = outerIter_constProperties_singlePhase_CO2_T_H(
            mgr, reg, freg, ppm, Tbc, Pbc, g, dt, ctrl
        );
        if (!ok) {
            std::cerr << "[runTransient(p+T)] step " << step1 << " failed.\n";
            return false;
        }

        // —— 导出（可分别控制频率与路径）—— //
        const bool dumpP = (writeEveryP <= 0) || (step1 % writeEveryP == 0);
        const bool dumpT = (writeEveryT <= 0) || (step1 % writeEveryT == 0);

        if (dumpP) {
            const std::vector<Vector> gradP =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "p_g", /*smoothIters=*/0);

            std::ostringstream fnP;
            fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ nullptr, /*Pbc*/ &Pbc,
                /*cell*/ "p_g",
                /*face*/ "p_g_face_tmp",
                /*gradBuf*/ &gradP,
                /*out*/ fnP.str()
            );
            if (!okPltP) 
            {
                std::cerr << "[Transient(P)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }

        }

        if (dumpT) {
            const std::vector<Vector> gradT =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "T", /*smoothIters=*/0);

            std::ostringstream fnT;
            fnT << outPrefixT << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPltT = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ &Tbc, /*Pbc*/ nullptr,
                /*cell*/ "T",
                /*face*/ "T_face_tmp",
                /*gradBuf*/ &gradT,
                /*out*/ fnT.str()
            );
            if (!okPltT) {
                std::cerr << "[Transient(T)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }

            // CO2 相：压力场名是 p_g；如切到水相就改为 "p_w"
            const char* pField = "p_g";
            const char* TField = "T";
            // 自定义一个你喜欢的文件夹路径
            const std::string outFolderTXT = "./Postprocess_Data/TXT";

            // 导出：# columns: cell_id  cx  cy  cz  volume  p  T
            PostChecks::export_pT_txt_after_step(
                mgr, reg,
                pField, TField,
                step1, t,
                outFolderTXT
            );
        }
    }
    return true;
}


// 测试案例：2D-变物性-单相-CO2-达西渗流-传热耦合-带井源问题
//测试案例：2D-常物性-单相-CO2-达西渗流-传热耦合问题
inline bool runTransient_constProperties_singlePhase_CO2_T_H_withWell
(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const PressureBCAdapter& Pbc,
    const Vector& g,                   // 例如 {0,-9.81,0}
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    const std::vector<WellConfig>& wellsCfg_in,
    // —— 输出控制（可分别设置频率与路径）——
    int writeEveryP = 1,
    int writeEveryT = 1,
    const std::string& outPrefixP = "P_CO2_pTH/p",
    const std::string& outPrefixT = "T_CO2_pTH/t"
)
{
    if (dt <= 0.0) { std::cerr << "[runTransient(p+T)] invalid dt.\n"; return false; }

#if __cplusplus >= 201703L
    // 保证父目录存在
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefixP).parent_path());
        std::filesystem::create_directories(std::filesystem::path(outPrefixT).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient(p+T)] cannot create directory for: "
            << outPrefixP << " or " << outPrefixT << "\n";
        return false;
    }
#endif

    double t = 0.0;
    for (int step = 0; step < nSteps; ++step) {
        const int step1 = step + 1;
        t += dt;

        // —— 时间层切换：压力与温度 —— //
        if (!startTimeStep_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev")) return false;
        if (!startTimeStep_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev")) return false;

        // —— 单步外迭代：p→mf→T（常物性）—— //
        bool ok = outerIter_constProperties_singlePhase_CO2_T_H_withWell( mgr, reg, freg, ppm, Tbc, Pbc, g, wellsCfg_in,dt, ctrl);
        if (!ok) 
        {
            std::cerr << "[runTransient(p+T)] step " << step1 << " failed.\n";
            return false;
        }

        // —— 导出（可分别控制频率与路径）—— //
        const bool dumpP = (writeEveryP <= 0) || (step1 % writeEveryP == 0);
        const bool dumpT = (writeEveryT <= 0) || (step1 % writeEveryT == 0);

        if (dumpP) {
            const std::vector<Vector> gradP =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "p_g", /*smoothIters=*/0);

            std::ostringstream fnP;
            fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPltP = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ nullptr, /*Pbc*/ &Pbc,
                /*cell*/ "p_g",
                /*face*/ "p_g_face_tmp",
                /*gradBuf*/ &gradP,
                /*out*/ fnP.str()
            );
            if (!okPltP)
            {
                std::cerr << "[Transient(P)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }

        }

        if (dumpT) {
            const std::vector<Vector> gradT =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "T", /*smoothIters=*/0);

            std::ostringstream fnT;
            fnT << outPrefixT << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";

            const bool okPltT = outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg,
                /*Tbc*/ &Tbc, /*Pbc*/ nullptr,
                /*cell*/ "T",
                /*face*/ "T_face_tmp",
                /*gradBuf*/ &gradT,
                /*out*/ fnT.str()
            );
            if (!okPltT) {
                std::cerr << "[Transient(T)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }

            // CO2 相：压力场名是 p_g；如切到水相就改为 "p_w"
            const char* pField = "p_g";
            const char* TField = "T";
            // 自定义一个你喜欢的文件夹路径
            const std::string outFolderTXT = "./Postprocess_Data/TXT";

            // 导出：# columns: cell_id  cx  cy  cz  volume  p  T
            PostChecks::export_pT_txt_after_step(
                mgr, reg,
                pField, TField,
                step1, t,
                outFolderTXT
            );
        }
    }
    return true;
}
