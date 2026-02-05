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

/**
 * @brief 构建复杂裂缝网络 (十字锁 + 浮动板)
 * @details 包含 3 条互相正交或倾斜的裂缝，用于测试拓扑稳定性和坐标映射精度。
 */
void SetupComplexFractureNetwork(MeshManager_3D& meshMgr)
{
    std::cout << " -> [Geometry] Building 20x20x20 Matrix (High Resolution)..." << std::endl;
    // 使用 20x20x20 的基岩网格 (总计 8000 cells)，保证足够的分辨率观察梯度
    meshMgr.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "Network_Geom");

    std::cout << " -> [Geometry] Generating Fracture Network (3 Fracs)..." << std::endl;

    // 1. Fracture A: 垂直裂缝 (YZ平面, x=50) - "The Wall"
    // 特征：常数 X，变量 Y, Z
    std::vector<Vector> ptsA = {
        Vector(50, 20, 0), Vector(50, 80, 0), Vector(50, 80, 100), Vector(50, 20, 100)
    };
    Fracture_2D fracA(0, ptsA);
    fracA.conductivity = 1000.0;

    // 2. Fracture B: 水平裂缝 (XY平面, z=50) - "The Floor"
    // 特征：常数 Z，变量 X, Y (与 A 正交)
    std::vector<Vector> ptsB = {
        Vector(20, 20, 50), Vector(80, 20, 50), Vector(80, 80, 50), Vector(20, 80, 50)
    };
    Fracture_2D fracB(1, ptsB);
    fracB.conductivity = 500.0;

    // 3. Fracture C: 倾斜裂缝 (Slanted) - "The Slope"
    // 特征：X, Y, Z 全变 (非轴对齐测试)
    std::vector<Vector> ptsC = {
        Vector(20, 20, 10), Vector(40, 20, 10), Vector(40, 40, 40), Vector(20, 40, 40)
    };
    Fracture_2D fracC(2, ptsC);
    fracC.conductivity = 200.0;

    // 加入管理器
    meshMgr.addFracturetoFractureNetwork(fracA);
    meshMgr.addFracturetoFractureNetwork(fracB);
    meshMgr.addFracturetoFractureNetwork(fracC);

    // 统一网格化
    std::cout << " -> [Mesh] Meshing all fractures (10x10)..." << std::endl;
    meshMgr.meshAllFracturesinNetwork(10, 10);

    // 重建全局索引
    meshMgr.rebuildFractureGlobalIndex();
}

void RunTest_DFN_Initialization_And_Viz()
{
    std::cout << "\n========================================================================" << std::endl;
    std::cout << "   Test C: Fracture Network Stability & Full-Field Visualization" << std::endl;
    std::cout << "========================================================================" << std::endl;

    // 1. 环境准备 (20x20x20 网格)
    MeshManager_3D meshMgr(100.0, 100.0, 100.0, 20, 20, 20, true, false);
    SetupComplexFractureNetwork(meshMgr);

    FieldManager_3D fieldMgr;
    size_t nMat = meshMgr.mesh().getCells().size();
    size_t nMatFaces = meshMgr.mesh().getFaces().size();
    size_t nFrac = meshMgr.fracture_network().fracElemIndex.total;
    size_t nFracEdges = meshMgr.fracture_network().getAllFractureEdges().size();
    fieldMgr.InitSizes(nMat, nFrac, 0, 0, nMatFaces, nFracEdges);

    std::cout << " -> [Info] Matrix Cells: " << nMat << std::endl;
    std::cout << " -> [Info] Matrix Faces: " << nMatFaces << std::endl;
    std::cout << " -> [Info] Fracture Cells: " << nFrac << std::endl;
    std::cout << " -> [Info] Fracture Faces: " << nFracEdges << std::endl;

    VariableInitializer_3D initMgr(meshMgr, fieldMgr);

    // 配置字符串
    auto pConfig = PressureEquation_String::SinglePhase();
    auto tConfig = TemperatureEquation_String::SinglePhase();
    auto sConfig = SaturationEquation_String();

    // 2. 设定正交梯度场 (Orthogonal Gradients)
    std::cout << " -> [Init] Applying Orthogonal Gradients (P->Z, T->X, Sw->Y)..." << std::endl;

    // [P] Pressure: 随 Z 增加而增加 (模拟静水压力)
    // P = 10 MPa + 0.1 MPa/m * Z
    LinearInitParams pInit;
    pInit.refVal = 10.0e6;
    pInit.grad_z = 1.0e5;

    // [T] Temperature: 随 X 增加而增加 (横向温差)
    // T = 300 K + 1.0 K/m * X
    LinearInitParams tInit;
    tInit.refVal = 300.0;
    tInit.grad_x = 1.0;

    // [Sw] Saturation: 随 Y 增加而增加 (纵向水驱)
    // Sw = 0.2 + 0.006 * Y (从 Y=0 的 0.2 到 Y=100 的 0.8)
    LinearInitParams sInit;
    sInit.refVal = 0.3;
    sInit.grad_y = 0.006;

    // VG & RelPerm 参数 (用于计算 Pc 和 Kr)
    // 这些参数将基于 Sw 场自动生成 Pc 和 Kr 场，验证派生变量计算的稳定性
    VGParams vg; vg.alpha = 1e-4; vg.n = 2.0;
    RelPermParams rp; rp.L = 0.5;

    // 3. 执行完整初始化 (IMPES Pipeline)
    // 这将自动计算 Matrix 和 Fracture 的 P, T, Sw, Pc, Krw, Krg
    initMgr.InitIMPESState(pConfig, tConfig, sConfig, pInit, tInit, sInit, vg, rp);

    // 4. 导出 Tecplot
    std::string filename = "Test/FieldOperator/Network_Viz_Full.dat";
    std::cout << " -> [Viz] Exporting " << filename << "..." << std::endl;

    PostProcess_3D post(meshMgr, fieldMgr);
    post.ExportTecplot(filename, 0.0);

    std::cout << " -> [Result] Export Done." << std::endl;
    std::cout << "    [Action Checklist for Tecplot]:" << std::endl;
    std::cout << "    1. Check Matrix Contour:" << std::endl;
    std::cout << "       - Pressure should vary along Z axis (Vertical Layers)." << std::endl;
    std::cout << "       - Temperature should vary along X axis (Left-Right)." << std::endl;
    std::cout << "       - Sw should vary along Y axis (Front-Back)." << std::endl;
    std::cout << "    2. Check Fracture B (The Floor, Z=50):" << std::endl;
    std::cout << "       - Should have CONSTANT Pressure (approx 15 MPa)." << std::endl;
    std::cout << "       - Should have VARYING Sw (along Y)." << std::endl;
    std::cout << "    3. Check Fracture A (The Wall, X=50):" << std::endl;
    std::cout << "       - Should have CONSTANT Temperature (approx 350 K)." << std::endl;
    std::cout << "       - Should have VARYING Sw (along Y)." << std::endl;
}