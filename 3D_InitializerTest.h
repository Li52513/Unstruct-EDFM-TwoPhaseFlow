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
// 辅助函数：构建标准测试几何
// =========================================================
void SetupStandardGeometry(MeshManager_3D& meshMgr)
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
// 测试 A: 初始化场的可视化验证 (Tecplot)
// =========================================================
void RunTest_Initialization_And_Viz()
{
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "   Test A: Initialization & Visualization (Gradient Verification)" << std::endl;
    std::cout << "========================================================================" << std::endl;

    // 1. 环境准备
    MeshManager_3D meshMgr(100.0, 100.0, 100.0, 10, 10, 10, true, false);
    SetupStandardGeometry(meshMgr);

    FieldManager_3D fieldMgr;
    size_t nMat = meshMgr.mesh().getCells().size();
    size_t nMatFaces = meshMgr.mesh().getFaces().size();
    size_t nFrac = meshMgr.fracture_network().fracElemIndex.total;
    size_t nFracEdges = meshMgr.fracture_network().getAllFractureEdges().size();
    fieldMgr.InitSizes(nMat, nFrac, 0, 0, nMatFaces, nFracEdges);

    // 2. 管理器实例化
    VariableInitializer_3D initMgr(meshMgr, fieldMgr);
    // 配置
    auto pConfig = PressureEquation_String::SinglePhase();
    auto tConfig = TemperatureEquation_String::SinglePhase();
    auto sConfig = SaturationEquation_String();

    // 3. 设定具有梯度的物理场 (模拟真实油藏)
    std::cout << " -> [Init] Applying Gradients..." << std::endl;

    // P: 底部(Z=0) 30MPa, 顶部(Z=100) 20MPa (梯度 -100,000 Pa/m)
    LinearInitParams pInit(30.0e6, 0.0, -100000.0);

    // T: 底部 360K, 顶部 340K (梯度 -0.2 K/m)
    LinearInitParams tInit(360.0, 0.0, -0.2);

    // Sw: 底部 1.0 (水), 顶部 0.2 (油/气) (梯度 -0.008 /m)
    LinearInitParams sInit(1.0, 0.0, -0.008);

    // VG & RelPerm 参数 (用于生成 Pc 和 Kr)
    CapRelPerm::VGParams vg; vg.alpha = 1e-5; vg.n = 2.0;
    CapRelPerm::RelPermParams rp; rp.L = 0.5;

    // 4. 执行初始化 (IMPES)
    initMgr.InitIMPESState(pConfig, tConfig, sConfig, pInit, tInit, sInit, vg, rp);

    // 5. 后处理可视化导出 (Tecplot)
    std::cout << " -> [Viz] Exporting Tecplot file..." << std::endl;
    PostProcess_3D post(meshMgr, fieldMgr);
    std::string filename = "Test/FieldOperator/Viz_Init_Gradients.dat";
    post.ExportTecplot(filename, 0.0);

    std::cout << " -> [Result] File '" << filename << "' generated." << std::endl;
    std::cout << "    Action: Open in Tecplot. Verify Z-gradients for P, T, Sw, Pc." << std::endl;
}