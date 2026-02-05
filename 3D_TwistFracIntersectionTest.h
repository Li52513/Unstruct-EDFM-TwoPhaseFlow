#pragma once
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include "3D_MeshManager.h"

int TwistFracIntersectionTest_3D()
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
    mgr_3D.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, "3D_Mesh_test");
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
    mgr_3D.exportMeshInfortoTxt("Test/MeshTest/3D-EDFM/TwistFractureIntersectionTest/EDFM_test_");

    ///#导出裂缝网格信息至Txt文件
    mgr_3D.exportFracturesNetworkInfortoTxt("Test/MeshTest/3D-EDFM/TwistFractureIntersectionTest/EDFM_test_");

    ///#导出特定宏观裂缝裂缝微观裂缝的网格面非正交信息至.csv，以#3号裂缝为例
    mgr_3D.inspectFractureEdges_non_orthogonalInfor(3, "Test/MeshTest/3D-EDFM/TwistFractureIntersectionTest/EDFM_test_");

    ///#导出 F-F 交线及微观网格对应关系至 CSV
    mgr_3D.inspectIntersections_FracturetoFracture("Test/MeshTest/3D-EDFM/TwistFractureIntersectionTest/EDFM_test_");

    std::cout << "===========================================================" << std::endl;
    std::cout << "          Test Completed. Check output files.              " << std::endl;
    std::cout << "===========================================================" << std::endl;
    return 0;

}


