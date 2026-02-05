#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>

// =========================================================
// 引入核心模块头文件
// =========================================================
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "3D_PhysicalPropertiesManager.h"
#include "3D_VariableInitializer.h"
#include "3D_PostProcess.h"

using namespace PhysicalProperties_string_op;

// =========================================================
// 辅助结构：测试工况
// =========================================================
struct PropertyTestCase {
    std::string name;
    double P_target; // [Pa]
    double T_target; // [K]
    std::string desc;
};

// =========================================================
// 辅助函数：构建标准测试几何
// =========================================================
void SetupStandardGeometry1(MeshManager_3D& meshMgr)
{
    std::cout << " -> [Geometry] Building 100x100x100 Matrix with 10x10x10 Grid..." << std::endl;
    // 1. 基岩网格 (100x100x100 m, 10x10x10 cells)
    meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Bench_Geom");

    std::cout << " -> [Geometry] Embedding Fracture at x=52..." << std::endl;
    // 2. 裂缝 (x=52, YZ平面)
    std::vector<Vector> fracPts = {
        Vector(52, 0, 0), Vector(52, 100, 0), Vector(52, 100, 100), Vector(52, 0, 100)
    };
    Fracture_2D frac(0, fracPts);
    frac.aperture = 1e-3;       // 1 mm
    frac.conductivity = 1000.0; // 1000 mD*m

    meshMgr.addFracturetoFractureNetwork(frac);

    // 3. 裂缝网格化 (10x10)
    meshMgr.meshAllFracturesinNetwork(10, 10);

    // 4. 重建索引 (关键)
    meshMgr.rebuildFractureGlobalIndex();
}

// =========================================================
// 测试 B: 物性计算准确性扫描 (One-by-One Verify)
// =========================================================
void RunTest_Property_Accuracy_Sweep()
{
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "   Test B: Fluid Property Accuracy Sweep (CSV Output)" << std::endl;
    std::cout << "========================================================================" << std::endl;

    // 1. 环境准备 (同上)
    MeshManager_3D meshMgr(100.0, 100.0, 100.0, 10, 10, 10, true, false);
    SetupStandardGeometry1(meshMgr);

    FieldManager_3D fieldMgr;
    size_t nMat = meshMgr.mesh().getCells().size();
    size_t nMatFaces = meshMgr.mesh().getFaces().size();
    size_t nFrac = meshMgr.fracture_network().fracElemIndex.total;
    size_t nFracEdges = meshMgr.fracture_network().getAllFractureEdges().size();
    fieldMgr.InitSizes(nMat, nFrac, 0, 0, nMatFaces, nFracEdges);

    VariableInitializer_3D initMgr(meshMgr, fieldMgr);
    auto pConfig = PressureEquation_String::SinglePhase();
    auto tConfig = TemperatureEquation_String::SinglePhase();

    PhysicalPropertiesManager_3D propsMgr(meshMgr, fieldMgr, pConfig, tConfig);

    // 2. 定义 6 个测试工况
    std::vector<PropertyTestCase> cases = {
        {"Case1_LowLow",      5.0e6,  285.0, "Boundary (Liq CO2)"},
        {"Case2_HighHigh",    70.0e6, 690.0, "Deep Supercritical"},
        {"Case3_NearCrit",    8.0e6,  310.0, "Near CO2 Critical Point"},
        {"Case4_SubCritW",    20.0e6, 600.0, "Sub-critical Water"},
        {"Case5_Reservoir",   20.0e6, 350.0, "Typical Reservoir"},
        {"Case6_MidRange",    35.0e6, 450.0, "Mid-Range Condition"}
    };

    // 3. CSV 输出准备
    std::ofstream csv("Test/PropertyTest/Property_Accuracy_Sweep.csv");
    csv << "CaseID,Name,P_Set[MPa],T_Set[K],"
        << "Mat_Rho_W,Mat_Mu_W,Mat_Cp_W,Mat_Lam_W,"
        << "Mat_Rho_C,Mat_Mu_C,Mat_Cp_C,Mat_Lam_C,"
        << "Frac_Rho_C,Frac_Mu_C" // 抽样裂缝数据验证一致性
        << "\n";

    // 获取场引用 (用于读取结果)
    PhysicalProperties_string_op::CO2 co2Tags;
    PhysicalProperties_string_op::Water watTags;

    // 4. 循环测试
    for (size_t i = 0; i < cases.size(); ++i) {
        const auto& tc = cases[i];
        std::cout << " -> Running " << tc.name << " (P=" << tc.P_target / 1e6 << "MPa, T=" << tc.T_target << "K)...";

        // 4.1 重置全场状态 (Uniform Initialization)
        // 使用 VariableInitializer 确保全场(Matrix & Fracture)都被正确赋值
        LinearInitParams pInit(tc.P_target);
        LinearInitParams tInit(tc.T_target);
        initMgr.InitSinglePhaseState(pConfig, tConfig, pInit, tInit);

        // 4.2 执行物性更新
        propsMgr.UpdateAllFluidEOS();

        // 4.3 抽样读取 (Sampling)
        // 选取 Matrix Cell [0] 和 Fracture Cell [0]
        // 注意：FieldManager 的 get 接口返回 shared_ptr，需解引用访问 data
        double mat_rho_w = fieldMgr.getMatrixScalar(watTags.rho_tag)->data[0];
        double mat_mu_w = fieldMgr.getMatrixScalar(watTags.mu_tag)->data[0];
        double mat_cp_w = fieldMgr.getMatrixScalar(watTags.cp_tag)->data[0];
        double mat_lam_w = fieldMgr.getMatrixScalar(watTags.k_tag)->data[0];

        double mat_rho_c = fieldMgr.getMatrixScalar(co2Tags.rho_tag)->data[0];
        double mat_mu_c = fieldMgr.getMatrixScalar(co2Tags.mu_tag)->data[0];
        double mat_cp_c = fieldMgr.getMatrixScalar(co2Tags.cp_tag)->data[0];
        double mat_lam_c = fieldMgr.getMatrixScalar(co2Tags.k_tag)->data[0];

        // 裂缝数据 (验证初始化和更新是否覆盖到了裂缝)
        double frac_rho_c = 0.0, frac_mu_c = 0.0;
        if (nFrac > 0) {
            frac_rho_c = fieldMgr.getFractureScalar(co2Tags.rho_tag)->data[0];
            frac_mu_c = fieldMgr.getFractureScalar(co2Tags.mu_tag)->data[0];
        }

        // 4.4 写入 CSV
        csv << i << "," << tc.name << "," << tc.P_target / 1e6 << "," << tc.T_target << ","
            << mat_rho_w << "," << mat_mu_w << "," << mat_cp_w << "," << mat_lam_w << ","
            << mat_rho_c << "," << mat_mu_c << "," << mat_cp_c << "," << mat_lam_c << ","
            << frac_rho_c << "," << frac_mu_c
            << "\n";

        std::cout << " Done." << std::endl;
    }

    csv.close();
    std::cout << " -> [Result] File 'Property_Accuracy_Sweep.csv' generated." << std::endl;
    std::cout << "    Action: Verify values against CoolProp." << std::endl;
}