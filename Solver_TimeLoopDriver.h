#pragma once
#include <functional>
#include <iostream>
#include <string>
#include <sstream>   // ★ 新增
#include <iomanip>   // ★ 新增
#include <filesystem>   
#include "Solver_TimeLoopSkeleton.h"  // outerIter_OneStep_singlePhase(...)
#include "Solver_Accel.h"
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


inline bool runTransient_constProperties_singlePhase_CO2_T_H_withWell_accel(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const PressureBCAdapter& Pbc,
    const Vector& g,
    int nSteps,
    double dt,
    const SolverControls& ctrl,
    const std::vector<WellConfig>& wellsCfg_in,
    int writeEveryP = 1,
    int writeEveryT = 1,
    const std::string& outPrefixP = "P_CO2_pTH/p",
    const std::string& outPrefixT = "T_CO2_pTH/t",
    Solver::Accel::OuterIterRuntime* runtimeIn = nullptr)
{
    if (dt <= 0.0) {
        std::cerr << "[runTransient_accel] invalid dt.\n";
        return false;
    }

#if __cplusplus >= 201703L
    try {
        std::filesystem::create_directories(std::filesystem::path(outPrefixP).parent_path());
        std::filesystem::create_directories(std::filesystem::path(outPrefixT).parent_path());
    }
    catch (...) {
        std::cerr << "[runTransient_accel] cannot create directory for: "
            << outPrefixP << " or " << outPrefixT << "\n";
        return false;
    }
#endif

    Solver::Accel::OuterIterRuntime ownedRuntime;
    Solver::Accel::OuterIterRuntime& runtime = runtimeIn ? *runtimeIn : ownedRuntime;
    runtime.lastDt = dt;
    Solver::Accel::ensureOuterRuntimeCapacity(mgr.mesh(), runtime);

    double t = 0.0;
    double dt_step = dt;
    for (int step = 0; step < nSteps; ++step) {
        const int step1 = step + 1;
        const double dt_this = dt_step;
        t += dt_this;

        runtime.lastDt = dt_this;
        Solver::Accel::ensureOuterRuntimeCapacity(mgr.mesh(), runtime);

        if (!startTimeStep_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev")) return false;
        if (!startTimeStep_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev"))   return false;

        bool ok = Solver::Accel::outerIter_constProperties_singlePhase_CO2_T_H_withWell_accel(
            mgr, reg, freg, ppm, Tbc, Pbc, g, wellsCfg_in, dt_this, ctrl, runtime);
        if (!ok) {
            std::cerr << "[runTransient_accel] step " << step1 << " failed.\n";
            return false;
        }

        const bool dumpP = (writeEveryP <= 0) || (step1 % writeEveryP == 0);
        const bool dumpT = (writeEveryT <= 0) || (step1 % writeEveryT == 0);

        if (dumpP) {
            const std::vector<Vector> gradP =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "p_g", /*smoothIters=*/0);
            std::ostringstream fnP;
            fnP << outPrefixP << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";
            if (!outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg, /*Tbc*/ nullptr, /*Pbc*/ &Pbc,
                "p_g", "p_g_face_tmp", &gradP, fnP.str())) {
                std::cerr << "[Transient_accel(P)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }
        }

        if (dumpT) {
            const std::vector<Vector> gradT =
                computeCellGradients_LSQ_with_GG(mgr.mesh(), reg, "T", /*smoothIters=*/0);
            std::ostringstream fnT;
            fnT << outPrefixT << "_step_" << std::setw(5) << std::setfill('0') << step1 << ".plt";
            if (!outputTecplot_cellToFaceToNode_BC(
                mgr, reg, freg, /*Tbc*/ &Tbc, /*Pbc*/ nullptr,
                "T", "T_face_tmp", &gradT, fnT.str())) {
                std::cerr << "[Transient_accel(T)] Tecplot export failed at step " << step1 << ".\n";
                return false;
            }
        }

        if (ctrl.enable_dt_adapt) {
            double cfl = std::max(0.0, runtime.lastCFL_T);
            double target = (ctrl.dt_adapt_CFL_target > 0.0) ? ctrl.dt_adapt_CFL_target : ctrl.CFL_T_threshold;
            double hyst = std::max(0.0, ctrl.dt_adapt_CFL_hysteresis);
            double lower = target * std::max(0.0, 1.0 - hyst);
            double upper = target * (1.0 + hyst);
            double newDt = dt_this;
            bool changed = false;

            if (cfl > upper && ctrl.dt_adapt_shrink > 0.0 && ctrl.dt_adapt_shrink < 1.0) {
                newDt = dt_this * ctrl.dt_adapt_shrink;
                changed = true;
            }
            else if (cfl > 0.0 && cfl < lower && ctrl.dt_adapt_grow > 1.0) {
                newDt = dt_this * ctrl.dt_adapt_grow;
                changed = true;
            }

            newDt = std::min(std::max(newDt, ctrl.dt_min), ctrl.dt_max);
            dt_step = newDt;

            if (ctrl.reportPerIter && changed) {
                std::cout << "[dt-adapt] step " << step1
                    << "  maxCFL_T=" << cfl
                    << "  dt: " << dt_this << " -> " << dt_step << "\n";
            }
        }
    }
    return true;
}
