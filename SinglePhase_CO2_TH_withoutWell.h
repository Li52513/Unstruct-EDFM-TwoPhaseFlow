#pragma once
#include <chrono>   //用于记录时间
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include "PhysicalPropertiesManager.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "Solver_TimeLoopUtils.h" 
#include "Solver_TimeLoopDriver.h"

int SinglePhase_CO2_TH_withoutWell()
{

    //定义几何区域参数
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    //确定网格划分策略及参数
    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    // 1) 构造并预处理网格
    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";

    // 2）生成变量场
    FieldRegistry reg;  //基岩变量场生成器
    FaceFieldRegistry freg; //基岩面变量场生成器
    InitFields ic;

    Initializer::createPrimaryFields_singlePhase_CO2_P(mgr.mesh(), reg);
    Initializer::createPrimaryFields_singlePhase_CO2_T(mgr.mesh(), reg);
    Initializer::fillBaseDistributions_singlePhase_CO2_P(mgr.mesh(), reg, ic);
    Initializer::fillBaseDistributions_singlePhase_CO2_T(mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    //3）物性参数设置,当前给定常物性参数
    PhysicalPropertiesManager ppm;
    ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
    ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
    ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

    //4) 边界条件设置
    const auto& bfaces = mgr.boundaryFaces();
    /// 压力边界条件设置
    //给定参数：2D情况，给定基岩四个边界条件系数
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0,0,6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0,0,5.5e6 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0,1,0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0,1,0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter PbcA{ pbc_pw };

    // /  温度边界条件设置
    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 1.0,0,673.15 };
    TemperatureBC::BoundaryCoefficient T_Right{ 0.0,1,0.0 };
    TemperatureBC::BoundaryCoefficient T_Down{ 0.0,1,0.0 };
    TemperatureBC::BoundaryCoefficient T_Up{ 0.0,1,0.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    tbc.printInformationofBoundarySetting(mgr);
    TemperatureBCAdapter TbcA{ tbc };

    // 5) 求解器与时间推进设置
    Vector g = { 0.0, 0.0, 0.0 };

    SolverControls sc;
    sc.maxOuter = 300;

    // 压力
    sc.tol_p_abs = 1e-6;     // 绝对容差（Pa）
    sc.tol_p_rel = 1e-6;     // 相对容差
    sc.urf_p = 1;     // 欠松弛
    sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_p.maxIters = 5000;
    sc.lin_p.tol = sc.tol_p_abs;
    sc.lin_p.iluFill = 10;
    sc.lin_p.iluDrop = 1e-4;

    // 温度
    sc.tol_T_abs = 1e-6;
    sc.tol_T_rel = 1e-6;
    sc.urf_T = 1;
    sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_T.maxIters = 5000;
    sc.lin_T.tol = sc.tol_T_abs;
    sc.lin_T.iluFill = 10;
    sc.lin_T.iluDrop = 1e-4;

    sc.useJacobi = false;

    // 时间步
    int    nSteps = 100;
    double dt = 10000 / nSteps;

    // 6) 运行：p→mf→T（分别输出到两个路径）
    bool ok = runTransient_constProperties_singlePhase_CO2_T_H
    (
        mgr, reg, freg, ppm,
        TbcA, PbcA, g,
        nSteps, dt, sc,
        /*writeEveryP*/ 10, /*writeEveryT*/ 10,
        /*outPrefixP*/ "./Postprocess_Data/SinglePhase/P_CO2/p",
        /*outPrefixT*/ "./Postprocess_Data/SinglePhase/T_CO2/T"
    );

    if (!ok) {
        std::cerr << "[MAIN] Transient run failed.\n";
        return 1;
    }

    std::cout << "[MAIN] Done.\n";
    return 0;

}