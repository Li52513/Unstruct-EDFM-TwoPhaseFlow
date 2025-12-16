#pragma once
#include <chrono>   //用于记录时间
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "Solver_Loop.h"

int SinglePhase_CO2_TH_withWell_reviese()
{
    //定义几何区域参数
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    //确定网格划分策略及参数
    int sectionNumX = 50, sectionNumY = 50, sectionNumZ = 0;
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

    //4) 边界条件设置
    const auto& bfaces = mgr.boundaryFaces();
    /// 压力边界条件设置
    //给定参数：2D情况，给定基岩四个边界条件系数
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0,1.0,0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值，3D是按照左右 下上 前后的顺序
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter PbcA{ pbc_pw };

    // /  温度边界条件设置
    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Right{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Down{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Up{ 0.0,1.0,0.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    tbc.printInformationofBoundarySetting(mgr);
    TemperatureBCAdapter TbcA{ tbc };

    // 全局/统一 Peaceman 参数（如需每井不同，可逐井修改 cfg.pm）
    PeacemanParams prm;
    prm.mu = 1.48e-5; prm.rho = 800; prm.reFactor = 1.12; prm.fallbackKh = 1e-14; prm.PI_is_mass = true;

    // —— 构造两个注入井，一个采出井（示例） —— //
    std::vector<WellConfig> wellsCfg;
    wellsCfg.reserve(3);

    // INJ1
    WellConfig inj1;
    inj1.name = "inj_1";
    inj1.role = WellDOF::Role::Injector;
    inj1.geom.name = inj1.name;               // 便于日志
    inj1.geom.pos = Vector(0.25, 0.8, 0.0);
    inj1.geom.rw = 0.001;
    inj1.geom.skin = 0.0;
    inj1.geom.H = 1.0;
    inj1.geom.perfRadius =0;
    inj1.geom.maxHitCells = 8;
    inj1.pm = prm;                           // 拷贝全局 Peaceman 参数
    inj1.mode = WellDOF::Mode::Rate;       // 固定井底压
    inj1.target = 2;                         // BHP
    inj1.Tin = 373.15;                      // 注入温度（只对注井用）
    inj1.derive_names_if_empty();              // mask_inj_1 / PI_inj_1 / p_w_inj_1
    wellsCfg.push_back(inj1);

    // INJ2
    WellConfig inj2;
    inj2.name = "inj_2";
    inj2.role = WellDOF::Role::Injector;
    inj2.geom.name = inj2.name;
    inj2.geom.pos = Vector(0.25, 0.2, 0.0);
    inj2.geom.rw = 0.001;
    inj2.geom.skin = 0.0;
    inj2.geom.H = 1.0;
    inj2.geom.perfRadius = 0;
    inj2.geom.maxHitCells = 8;
    inj2.pm = prm;
    inj2.mode = WellDOF::Mode::Rate; //Rate Pressure
    inj2.target = 2;
    inj2.Tin = 373.15;
    inj2.derive_names_if_empty();              // mask_inj_2 / PI_inj_2 / p_w_inj_2
    wellsCfg.push_back(inj2);

    // PROD1
    WellConfig prod1;
    prod1.name = "prod_1";
    prod1.role = WellDOF::Role::Producer;
    prod1.geom.name = prod1.name;
    prod1.geom.pos = Vector(0.75, 0.5, 0.0);
    prod1.geom.rw = 0.001;
    prod1.geom.skin = 0.0;
    prod1.geom.H = 1.0;
    prod1.geom.perfRadius = 0.0;
    prod1.geom.maxHitCells = 1;
    prod1.pm = prm;
    prod1.mode = WellDOF::Mode::Pressure;
    prod1.target = 5e6;
    prod1.Tin = 0.0;                       
    prod1.derive_names_if_empty();             // mask_prod_1 / PI_prod_1 / p_w_prod_1
    wellsCfg.push_back(prod1);

    // 6) 求解器与时间推进设置
    Vector g = { 0.0, 0.0, 0.0 };

    SinglePhase::FaceMassRateConfig mf_ctrl;
    SinglePhase::PressureSolveControls p_sol_ctrl;
    SinglePhase::TemperatureSolveControls T_sol_ctrl;
    const int writeEveryP = 1;
    const int writeEveryT = 1;

    const std::string outPrefixP = "./Postprocess_Data/SinglePhase_revised/Case1/P";
    const std::string outPrefixT = "./Postprocess_Data/SinglePhase_revised/Case1/T";
    const int snapshotEveryCsv = 1;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/Case13/pT_state";

    // —— 时间步 —— //
    const int    nSteps = 100;
    const double dt = 1000.0 / nSteps;

    bool ok = SinglePhase::runTransient_SinglePhase_HT_Iteration(mgr, reg, freg, ppm, TbcA, PbcA, wellsCfg, p_sol_ctrl, mf_ctrl, T_sol_ctrl, dt, nSteps, g, writeEveryP, writeEveryT, outPrefixP, outPrefixT, snapshotEveryCsv, snapshotPrefix);

    std::cout << "[MAIN] Accelerated run finished.\n";
    return 0;
}