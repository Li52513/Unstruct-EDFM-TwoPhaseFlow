#pragma once
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include "3D_MeshManager.h"

#include "Mesh.h"
#include "2D_Fracture.h"
#include "2D_FractureNetwork.h"
#include "FractureNetwork.h"
#include "MeshManager.h"
#include "PhysicalPropertiesManager.h"
#include "Initializer.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "FaceFieldRegistry.h"
#include "Solver_TimeLoopDriver.h"
#include "CouplingAssembler.h" 
#include "PostProcessor.h"

int EDFM_withFracture_Geomtry()
{

    /***************************3D基岩网格全局参数定义*******************************/
    double lengthX = 1, lengthY = 1, lengthZ = 1;
    int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 5;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
    
    /**************************3D基岩网格绘制并导出******************************/
    std::cout << "===========================================================" << std::endl;
    std::cout << "          3D-EDFM Matrix mesh is generating          " << std::endl;
    std::cout << "===========================================================" << std::endl;

    //【基岩网格划分】
    MeshManager_3D mgr_3D(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase); //构建基岩-裂缝网格管理器
    mgr_3D.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection,"3D_Mesh_test");
    mgr_3D.mesh().printMeshInfo();
    mgr_3D.exportMeshInfortoTxt("check3D");

    auto groups = BoundaryFaceClassify_byTag::ClassifyBoundaryFacesByTag_3D(mgr_3D.mesh());
    std::cout << "Boundary Faces Classification Result:" << std::endl;
    std::cout << "  X=0 (Left)   : " << groups.x0.size() << std::endl;
    std::cout << "  X=L (Right)  : " << groups.xL.size() << std::endl;
    std::cout << "  Y=0 (Bottom) : " << groups.y0.size() << std::endl;
    std::cout << "  Y=L (Top)    : " << groups.yL.size() << std::endl;
    std::cout << "  Z=0 (Front)  : " << groups.z0.size() << std::endl;
    std::cout << "  Z=L (Back)   : " << groups.zL.size() << std::endl;
   

    /**************************2D裂缝添加并独立进行网格划分******************************/
    std::cout << "===========================================================" << std::endl;
    std::cout << "          3D-EDFM Fracture Network is generating            " << std::endl;
    std::cout << "===========================================================" << std::endl;

    // =========================================================
    // 场景 A: 标准相交的一对平面裂缝 (基准对照组)
    // =========================================================

    //【裂缝四个角点坐标设置】 逆时针，左下 右下 右上 左上
    //// 裂缝 0
    std::vector<Vector> corners0 = 
    {
    Vector(0.0, 1.0, 0.0), Vector(4.0, 1.0, 0.0),
    Vector(4.0, 1.2, 2.0), Vector(0.0, 0.8, 2.0)
    };
    
    Fracture_2D frac0(0, corners0, 100.0, 0.001); //分别对应 编号  裂缝边界顶点  导流能力 开度

    //// 裂缝 1
    std::vector<Vector> corners1 = {
        Vector(1.0, 0.0, 0.5), Vector(1.0, 2.0, 0.5),
        Vector(1.0, 2.0, 1.5), Vector(1.0, 0.0, 1.5)
    };
    Fracture_2D frac1(1, corners1, 50.0, 0.002);

    //// 裂缝 2:
    std::vector<Vector> corners2 = {
        Vector(2.0, 0.5, 0.2), // P0 (BL) - Z=0.0
        Vector(2.5, 1.5, 0.2), // P1 (BR) - Z=0.5
        Vector(2.5, 1.5, 1.8), // P2 (TR) - Z=0.0
        Vector(2.0, 0.5, 1.8)  // P3 (TL) - Z=0.5
    };
    //// 加扭曲
    corners2[2].m_x += 0.2;
    Fracture_2D frac2(2, corners2, 10.0, 0.001);

    //// Frac 3
    std::vector<Vector> corners3 = {
        Vector(3.0, 0.0, 1.0), // P0 (Left-Bottom)  - Z=0.0
        Vector(3.0, 1.1, 1.0), // P1 (Right-Bottom) - Z=0.6 (比Frac2高一点)
        Vector(3.5, 1.1, 1.0), // P2 (Right-Top)    - Z=0.0
        Vector(3.5, 0.0, 1.0)  // P3 (Left-Top)     - Z=0.6
    };
    Fracture_2D frac3(3, corners3, 10.0, 0.001);

    //【添加裂缝至裂缝网络】
    mgr_3D.addFracturetoFractureNetwork(frac0);
    mgr_3D.addFracturetoFractureNetwork(frac1);
    mgr_3D.addFracturetoFractureNetwork(frac2);
    mgr_3D.addFracturetoFractureNetwork(frac3);
    std::cout << "[Main] Added 4 fractures to the network." << std::endl;

    //【裂缝网格独立划分】 生成网格后会自动重建索引，内部自动调用rebuildGlobalIndex()
    // #Note:增加网格密度以更好地逼近空间曲线
    int nU = 5;
    int nV = 5;
    mgr_3D.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

    //【处理裂缝-裂缝相交】
    //#Note 支持非共面裂缝，生成微观线段集合 
    // strategy 策略开关 Octree_Optimized 八叉树优化  BruteForce 暴力遍历
    mgr_3D.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::BruteForce);

    //【后处理】导出数据
    ///#导出基岩网格信息
    mgr_3D.exportMeshInfortoTxt("EDFM_test_");

    ///#导出裂缝网格信息至Txt文件
    mgr_3D.exportFracturesNetworkInfortoTxt("EDFM_test_");

    ///#导出特定宏观裂缝裂缝微观裂缝的网格面非正交信息至.csv，以#3号裂缝为例
    mgr_3D.inspectFractureEdges_non_orthogonalInfor(3, "EDFM_test_");

    ///#导出 F-F 交线及微观网格对应关系至 CSV
    mgr_3D.inspectIntersections_FracturetoFracture("EDFM_test_");

    std::cout << "===========================================================" << std::endl;
    std::cout << "          Test Completed. Check output files.              " << std::endl;
    std::cout << "===========================================================" << std::endl;








    /*----------------------------------------------------------------------*/
    // 1) 构造并预处理网格
    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";

    // 2) 添加裂缝 & 裂缝网格划分 & 计算几何耦合系数
    mgr.addFracture({ 0.1,0.2,0 }, { 0.3,0.9,0 });
    mgr.addFracture({ 0.7,0.1,0 }, { 0.1,0.8,0 });
    mgr.addFracture({ 0.0725,0.1825,0 }, { 0.4025,0.3925,0 });
    auto t2 = std::chrono::high_resolution_clock::now(); // 计时开始

    //mgr.setDFNRandomSeed(12345); // 设置随机种子
    //mgr.generateDFN (/*N=*/2,/*minPoint=*/{ 0.0,0.0,0.0 },/*maxPoint=*/{ 1.0,1.0,0.0 },/*Lmin=*/0.5,/*Lmax=*/1.4,/*alpha=*/0,/*kappa=*/0,/*avoidOverlap=*/true);
    auto t8 = std::chrono::high_resolution_clock::now();

    mgr.setDistanceMetric(DistanceMetric::CrossAwareGauss); // 设置距离度量方式AreaWeight CellCenter NodeAverage CrossAwareGauss
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GlobalAABB);

    auto t9 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms10 = std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count();
    std::cout << "FractureDetect in " << ms10 << " ms.\n";
    mgr.ComputeFractureGeometryCouplingCoefficient();   //明确一点： 因为这里还没有给出物性参数，这里计算的只是几何耦合系数 CI_geo和 geomAlpha 各自的表达式分别为：CI_geo = L*1/d_avg geomAlpha = 2 / L;
    // 导出网格、裂缝信息
    mgr.exportMesh("mesh");
    mgr.exportFractures("fractures");
    mgr.printFractureInfo();
    /**************************主变量初始化以及物性参数设置模块******************************/
    /*----------------------------------------------------------------------*/
    auto t4 = std::chrono::high_resolution_clock::now(); // 计时开始

    // 3) 基岩和裂缝介质区域划分
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry
    (
        mgr,
        {
            // 高渗区
          { Cell::RegionType::High,

             {
                {0.0, 0.5, 0.0},
                {0.0, 1.0, 0.0},
                {0.5, 1.0, 0.0}
             }
          },
        // 低渗区
          { Cell::RegionType::Low,
            {
            {0.5, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {1.0, 0.5, 0.0}
            }
          }
        },
        /* defaultRegion = */
        Cell::RegionType::Medium
    );
    ppm.classifyFractureElementsByGeometry(mgr, 0, { 0.1, 0.2, 0 }, { 0.3, 0.9, 0 }, FractureElementType::Blocking, FractureElementType::Conductive);


    // 4) 生成场变量以及设置初始化参数和VG模型（温度T,水相饱和度Sw,水相压力Pw)
    FieldRegistry reg;      // 基岩场
    FieldRegistry reg_fr;   // 裂缝场（与你当前实现一致，用同一类型的注册表管理）
    InitFields ic;          // p0/T0/Sw0 及其梯度（默认为均匀）
    VGParams vg;            // vG 参数 =vg模型
    RelPermParams rp;       // 相对渗透率参数（默认 L=0.5）
    RockDefaults rock;      // 基岩默认热物性
    InitDiagnostics diag;   // 诊断统计

    // 5) 基岩“主变量场”创建与填充：p_w, S_w, T (+ p_c, p_g, kr_w, kr_g)

    //创建基岩的主变量场
    Initializer::createPrimaryFields(mgr.mesh(), reg,"T");

    // 测试案例-2D-常物性-温度扩散

    //填充基岩主变量场
	Initializer::fillBaseDistributions1(mgr.mesh(), reg, ic, "T");

    //7）裂缝主变量初始化 （当前还未涉及裂缝）**
    Initializer::initFracturePrimaries(mgr.mesh(), mgr.fracture_network(), reg, reg_fr);

    //8）基岩物性参数设置
    //   a.固相参数
    ppm.UpdateRockProperties(mgr, reg, "p_g", "T");

    //9）裂缝物性参数设置
    //   a.固相参数


    //   b.流体参数
  


    /**************************边界条件设置模块******************************/

    //边界面识别
    const auto& bfaces = mgr.boundaryFaces();

    // 你可以先自检打印（便于验证识别是否准确）
    std::cout << "[BC] x0=" << bfaces.x0.size()
        << " xL=" << bfaces.xL.size()
        << " y0=" << bfaces.y0.size()
        << " yL=" << bfaces.yL.size()
        << " z0=" << bfaces.z0.size()
        << " zL=" << bfaces.zL.size() << "\n";


    /// 压力边界条件设置

    //给定参数：2D情况，给定基岩四个边界条件系数
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0,8e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0,6e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Down{ 1.0, 0.0,6e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Up{ 1.0, 0.0,8e6 }; // p = 2e5 Pa
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter bcAdapter{ pbc_pw };

    // /  温度边界条件设置
    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 450.0 };
    TemperatureBC::BoundaryCoefficient T_Right{ 1.0, 0.0, 373.15 };
    TemperatureBC::BoundaryCoefficient T_Down{ 1.0, 0.0, 373.15 };
    TemperatureBC::BoundaryCoefficient T_Up{ 1.0, 0.0, 450.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    tbc.printInformationofBoundarySetting(mgr);
    TemperatureBCAdapter TbcA{ tbc };

    /**************************基岩裂缝面系数计算 当前系数不考虑裂缝的离散**************************/



    //生成储存储存网格面上离散系数和源项的面场
    FaceFieldRegistry freg;

    // 打印初始化诊断
    Initializer::printDiag(diag);

    auto t5 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count();
    std::cout << "Fields & properties initialized in " << ms2 << " ms.\n";
    /***********************************************************************/

    /**************************耦合系数（CI/TI）初始化@t=0  #include "CouplingAssembler.h" *******************/
    // 1) 先给裂缝注册 CI 场
    {
        // 计算全局段数，用于创建场
        size_t Nseg = 0;
        for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
        ensureFracCouplingFields(reg_fr, Nseg); // 创建 CIw / CIg
    }

    // 2) 计算 CI（矩阵?裂缝）
    {
        bool upwind = true;            // 采用迎风（按相位势比较）
        bool include_gravity = false;  // 2D 平面暂不考虑重力项；需要的话置 true
        double g = 9.80665;
        updateMatrixFractureCI(mgr, reg, reg_fr, vg, upwind, include_gravity, g);
    }

    // 3) 计算 TI（裂缝?裂缝，多臂交点 Star–Delta）
    {
        bool include_gravity = false;
        double g = 9.80665;
        updateFractureFractureTI(mgr, reg_fr, vg, rp, include_gravity, g);
    }
    /***********************************************************************/


    /****************** 调试 & 输出 ******************/
    // 打印所有 Cell 和 FractureElement 的当前物性
    ppm.debugPrintProperties(mgr, reg, reg_fr, /*maxPrint=*/10);
    // 导出网格、裂缝信息
    mgr.exportMesh("mesh");
    mgr.exportFractures("fractures");
    mgr.printFractureInfo();
    printCI_Diagnostics(mgr, reg, reg_fr,  /*maxFracs*/(size_t)-1, /*maxSegs*/(size_t)-1);
    printTI_Diagnostics(mgr, reg_fr,       /*maxFracs*/(size_t)-1);
    //mgr.printCISourceTerms();


    cout << "Finished initial setup and property assignment.\n";

    // —— 导出 t=0、step=0 的场，用于初始时刻可视化 ——
    const double time0 = 0.0;
    const int    step0 = 0;

    PostProcessor outM("out\\matrix");
    PostProcessor outF("out\\fracture");

    // 可选自检：确保裂缝三大主变量长度与全局段数一致
    {
        size_t Nseg = 0; for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
        auto pf = reg_fr.get<volScalarField>("pf_w");
        auto sf = reg_fr.get<volScalarField>("Sf_w");
        auto Tf = reg_fr.get<volScalarField>("Tf");
        std::cout << "[fracture] Nseg=" << Nseg
            << " pf_w=" << (pf ? pf->data.size() : 0)
            << " Sf_w=" << (sf ? sf->data.size() : 0)
            << " Tf=" << (Tf ? Tf->data.size() : 0) << "\n";
    }

    // 分开导出：矩阵 / 裂缝
    outM.exportMatrixValue(mgr.mesh(), reg, time0, step0);
    outF.exportFractureValue(mgr.mesh(), mgr.fracture_network(), reg_fr, time0, step0);

    std::cout << "Finished initial setup and wrote t=0 state for MATLAB.\n";


    return 0;

}


