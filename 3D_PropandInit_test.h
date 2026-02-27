/**
 * @file 3D_PropandInit_test.h
 * @brief 3D EDFM 热流耦合基准测试模块 (合并单头文件版)
 * @details 负责串联网格管理、场管理、物性初始化、流体相更新及后处理。
 * 验证阶段一（The Physics Kernel）落地后的完整正确性。
 * 提供独立的内联函数供 main() 直接调用。
 */

#pragma once

#include <string>
#include <iostream>
#include <chrono>

 // 引入相关的头文件
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "3D_VariableInitializer.h"
#include "3D_PhysicalPropertiesManager.h"
#include "3D_PostProcess.h"
#include "CapRelPerm_HD.h"
#include "SolverContrlStrName_op.h"


/**
 * @brief 运行全流程热流耦合基准测试
 * @param exportPrefix 导出 Tecplot 时的文件前缀，默认为 "3D_PropTest_Output"
 */
void RunBenchmark_3D_PropTest(const std::string& exportPrefix = "3D_PropTest_Output")
{
    std::cout << "========== [3D EDFM Benchmark Started] ==========" << std::endl;

    // =========================================================
    // Step 1: 几何重构与拓扑建立 (3D_MeshManager)
    // =========================================================

    double lengthX = 100.0, lengthY = 100.0, lengthZ = 100.0;
    int nx = 50, ny = 50, nz = 50;
    int nU = 10, nV = 10;

    std::cout << "\n[Step 1] Constructing Mesh and Topology..." << std::endl;
    std::cout << " 1. Generating 3D Matrix Mesh (" << nx << "x" << ny << "x" << nz << ")..." << std::endl;
    MeshManager_3D meshMgr(lengthX, lengthY, lengthZ, nx, ny, nz, true, false);
    meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "3D_Matrix_Mesh");

    std::vector<Fracture_2D> inputFractures;
    std::vector<Vector> frac1_pts = { Vector(20, 10, 10), Vector(80, 10, 10), Vector(80, 90, 90), Vector(20, 90, 90) };
    inputFractures.push_back(Fracture_2D(1, frac1_pts, 100.0, 0.005));

    std::vector<Vector> frac2_pts = { Vector(10, 50, 20), Vector(90, 50, 20), Vector(90, 50, 80), Vector(10, 50, 80) };
    inputFractures.push_back(Fracture_2D(2, frac2_pts, 50.0, 0.002));

    std::cout << " 2. Adding " << inputFractures.size() << " Fractures and Meshing (" << nU << "x" << nV << ")..." << std::endl;
    for (const auto& f : inputFractures) {
        meshMgr.addFracturetoFractureNetwork(f);
    }

    meshMgr.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

    std::cout << " 3. Setting up Global Solver Indices..." << std::endl;
    meshMgr.setupGlobalIndices();

    MeshManager_3D::IntersectionStrategy strategy = MeshManager_3D::IntersectionStrategy::Octree_Optimized;
    std::cout << " 4. Running M-F Intersection (Strategy: " << (int)strategy << ")..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    meshMgr.SolveIntersection3D_improved_twist_accleration(strategy);

    std::cout << " 5. Running F-F Intersection..." << std::endl;
    meshMgr.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);

    std::cout << " 6. Rebuilding Edge Properties for FVM..." << std::endl;
    meshMgr.fracture_network().rebuildEdgeProperties();

    meshMgr.removeDuplicateInteractions();
    meshMgr.resolveCoplanarInteractions();

    std::cout << " 7. Building Topology Maps..." << std::endl;
    meshMgr.buildTopologyMaps();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    size_t totalFracCells = 0;
    for (const auto& f : meshMgr.fracture_network().getFractures()) {
        totalFracCells += f.fracCells.size();
    }

    std::cout << "    -> Geometry & Topology Processing Time: " << elapsed.count() << " s" << std::endl;
    std::cout << "    -> Matrix Cells: " << meshMgr.mesh().getCells().size() << std::endl;
    std::cout << "    -> Fractures Elements: " << totalFracCells << std::endl;
    std::cout << "    -> NNC Pairs: " << meshMgr.getInteractionPairs().size() << std::endl;

    // =========================================================
    // Step 2: 场内存预分配 (3D_FieldManager)
    // =========================================================
    std::cout << "\n[Step 2] Initializing Field Manager Sizes..." << std::endl;
    FieldManager_3D fieldMgr;

    fieldMgr.InitSizes(
        meshMgr.mesh().getCells().size(),
        totalFracCells,
        meshMgr.getInteractionPairs().size(),
        meshMgr.fracture_network().ffIntersections.size(),
        meshMgr.countBoundaryFaces(),
        meshMgr.fracture_network().getGlobalEdges().size()
    );

    // =========================================================
    // Step 3: 物性计算与静态更新 (3D_PhysicalPropertiesManager)
    // =========================================================
    std::cout << "\n[Step 3] Initializing Static Properties (Rock & Fracture)..." << std::endl;

    PhysicalProperties_string_op::PressureEquation_String pConfig = PhysicalProperties_string_op::PressureEquation_String::FIM();
    PhysicalProperties_string_op::TemperatureEquation_String tConfig = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
    PhysicalProperties_string_op::SaturationEquation_String sConfig = PhysicalProperties_string_op::SaturationEquation_String::FIM();

    PhysicalPropertiesManager_3D propMgr(meshMgr, fieldMgr, pConfig, tConfig);

    RockPropertyParams bgProps(10.0, 10.0, 1.0, 0.2, 2600.0, 1000.0, 2.5, 1e-9);
    std::vector<std::pair<std::string, std::pair<BoundingBox3D, RockPropertyParams>>> regions;
    BoundingBox3D highPermBox = { 30.0, 70.0, 30.0, 70.0, 30.0, 70.0 };
    RockPropertyParams highPermProps(100.0, 100.0, 50.0, 0.3, 2600.0, 1000.0, 2.5, 1e-9);
    regions.push_back({ "HighPermLayer", {highPermBox, highPermProps} });
    propMgr.InitRockProperties(bgProps, regions);

    FractureGlobalParams fProps(1.0, 2500.0, 1200.0, 1.5, 0.5);
    propMgr.InitFractureProperties(fProps);

    // =========================================================
    // Step 4: 变量场初始化 (3D_VariableInitializer)
    // =========================================================
    std::cout << "\n[Step 4] Seeding FIM Variables (Pressure, Temp, Saturation)..." << std::endl;
    VariableInitializer_3D varInit(meshMgr, fieldMgr);

    LinearInitParams pInit(10.0e6, 0.0, 0.1e6);
    LinearInitParams tInit(300.0);
    tInit.grad_x = 1.0; tInit.x_ref = 0.0;
    LinearInitParams sInit(0.2);
    sInit.grad_y = 0.006; sInit.y_ref = 0.0;

    varInit.InitFIMState<3>(pConfig, tConfig, sConfig, pInit, tInit, sInit, 0, 1, 2);

    // =========================================================
    // Step 5: 基于饱和度计算相渗与毛管力 (VG & RelPerm)
    // =========================================================
    std::cout << "\n[Step 5] Deriving Relative Permeability & Capillary Pressure fields..." << std::endl;
    CapRelPerm::VGParams vgParams;
    CapRelPerm::RelPermParams rpParams;

    std::string krw_tag = PhysicalProperties_string_op::Water().k_rw_tag;
    std::string krg_tag = PhysicalProperties_string_op::CO2().k_rg_tag;
    std::string pc_tag = PhysicalProperties_string_op::TwoPhaseState_String().Pc_field;

    auto krw_m = fieldMgr.createMatrixADScalar<3>(krw_tag);
    auto krg_m = fieldMgr.createMatrixADScalar<3>(krg_tag);
    auto pc_m = fieldMgr.createMatrixADScalar<3>(pc_tag);
    auto sw_m = fieldMgr.getMatrixADScalar<3>(sConfig.saturation);

    auto krw_f = fieldMgr.createFractureADScalar<3>(krw_tag);
    auto krg_f = fieldMgr.createFractureADScalar<3>(krg_tag);
    auto pc_f = fieldMgr.createFractureADScalar<3>(pc_tag);
    auto sw_f = fieldMgr.getFractureADScalar<3>(sConfig.saturation);

    if (sw_m && krw_m && krg_m && pc_m) {
        for (size_t i = 0; i < sw_m->data.size(); ++i) {
            double sw_val = sw_m->data[i].val;
            double kw, kg;
            CapRelPerm::kr_Mualem_vG(sw_val, vgParams, rpParams, kw, kg);
            double pc = CapRelPerm::pc_vG(sw_val, vgParams);

            krw_m->data[i] = ADVar<3>(kw);
            krg_m->data[i] = ADVar<3>(kg);
            pc_m->data[i] = ADVar<3>(pc);
        }
    }

    if (sw_f && krw_f && krg_f && pc_f) {
        for (size_t i = 0; i < sw_f->data.size(); ++i) {
            double sw_val = sw_f->data[i].val;
            double kw, kg;
            CapRelPerm::kr_Mualem_vG(sw_val, vgParams, rpParams, kw, kg);
            double pc = CapRelPerm::pc_vG(sw_val, vgParams);

            krw_f->data[i] = ADVar<3>(kw);
            krg_f->data[i] = ADVar<3>(kg);
            pc_f->data[i] = ADVar<3>(pc);
        }
    }

    // =========================================================
    // Step 6: 动态流体物性更新 (Fluid EOS AD Update)
    // =========================================================
    std::cout << "\n[Step 6] Updating Dynamic Fluid EOS (AD Formulation)..." << std::endl;
    propMgr.UpdateFluidEOS_Water_AD<3>(pConfig.pressure_field, tConfig.temperatue_field, false);
    propMgr.UpdateFluidEOS_CO2_AD<3>(pConfig.pressure_field, tConfig.temperatue_field, false);

    // =========================================================
    // Step 7: 结果降维导出与可视化 (3D_PostProcess)
    // =========================================================
    std::cout << "\n[Step 7] Decoupling AD fields and Exporting to Tecplot..." << std::endl;

    // 【核心修复】：为导出的纯标量场统一添加 "_out" 后缀，绝对避免与已存在的同名 ADVar 场产生强转错误与越界。

    // 主变量降维同步
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, pConfig.pressure_field, pConfig.pressure_field + "_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, tConfig.temperatue_field, tConfig.temperatue_field + "_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, sConfig.saturation, sConfig.saturation + "_out");

    // 注意：岩石骨架属性(phi, K_xx 等)在 InitRockProperties 中本身就是作为纯标量场(volScalarField)创建的！
    // 它们根本不是 AD 场，因此已将它们从 SyncADFieldToScalar 中剔除。ExportTecplot 会自动识别这些纯标量场并输出。

    // 水相同步
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, PhysicalProperties_string_op::Water().rho_tag, "Rho_Water_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, PhysicalProperties_string_op::Water().mu_tag, "Mu_Water_out");

    // CO2相同步
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, PhysicalProperties_string_op::CO2().rho_tag, "Rho_CO2_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, PhysicalProperties_string_op::CO2().mu_tag, "Mu_CO2_out");

    // 派生相渗与毛管力同步
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, krw_tag, "Krw_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, krg_tag, "Krg_out");
    PostProcess_3D::SyncADFieldToScalar<3>(fieldMgr, pc_tag, "Pc_out");

    PostProcess_3D postProcess(meshMgr, fieldMgr);
    std::string filename = exportPrefix + ".dat";
    postProcess.ExportTecplot(filename, 0.0);
	postProcess.ExportVTK(exportPrefix+".vtk", 0.0);

    std::cout << "\n========== [3D EDFM Benchmark Completed Successfully] ==========" << std::endl;
}