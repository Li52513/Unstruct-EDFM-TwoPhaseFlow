/**
 * @file main.cpp
 * @brief Improved 3D-EDFM 系统集成测试入口
 * @details 包含四个标准测试用例：
 * 1. 单一平面裂缝
 * 2. 单一扭曲(Twisted)裂缝
 * 3. 多根非相交裂缝
 * 4. 多根相交裂缝
 * @author Professor (AI Assistant)
 * @date 2026-01-13
 */

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <chrono>

 // 引入核心头文件
#include "UserDefineVarType.h" // 包含 Vector 定义
#include "3D_MeshManager.h"
#include "2D_Fracture.h"


/**
 * @brief 快速创建裂缝对象的辅助函数
 * @param id 裂缝 ID
 * @param p1, p2, p3, p4 裂缝的四个角点 (逆时针顺序需满足围成环)
 * @return Fracture_2D 对象
 */
Fracture_2D CreateFracture(int id, const Vector& p1, const Vector& p2, const Vector& p3, const Vector& p4)
{
    std::vector<Vector> boundaryVertices = { p1, p2, p3, p4 };
    // [Assumption]: Fracture_2D 构造函数接受 ID 和顶点列表
    return Fracture_2D(id, boundaryVertices);
}
/**
 * @brief 运行单个测试用例的通用流程
 * @param caseName 测试用例名称 (用于日志和文件名)
 * @param fractures 该用例包含的裂缝列表
 */

void RunTestCase(const std::string& caseName, const std::vector<Fracture_2D>& fractures)
{
    std::cout << "\n#########################################################" << std::endl;
    std::cout << "Running Test Case: " << caseName << std::endl;
    std::cout << "#########################################################" << std::endl;

    auto startTotal = std::chrono::high_resolution_clock::now();

    // 1. 初始化网格管理器 (3D Mesh Manager)
    // -----------------------------------------------------
    // 域范围: 100m x 100m x 100m
    // 网格数: 10 x 10 x 10 = 1000 Cells
    // 网格类型
    double lengthX = 100, lengthY = 100, lengthZ = 100;
    int sectionNumX = 10, sectionNumY = 10, sectionNumZ = 10;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    std::cout << "-> Initializing MeshManager_3D (" << lengthX << "x" << lengthY << "x" << lengthZ << ")..." << std::endl;
    MeshManager_3D manager(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ,usePrism,useQuadBase);

    // 2. 生成基岩网格 (Matrix Mesh)
    // -----------------------------------------------------
    // 使用正交修正 (Orthogonal)
    std::string meshName = caseName + "_Matrix";
    manager.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, meshName);

    // 3. 添加裂缝 (Add Fractures)
    // -----------------------------------------------------
    std::cout << "-> Adding " << fractures.size() << " fractures to network..." << std::endl;
    for (const auto& frac : fractures)
    {
        manager.addFracturetoFractureNetwork(frac);
    }

    // 4. 裂缝离散化 (Mesh Fractures)
    // -----------------------------------------------------
    // 将每条裂缝划分为 5x5 的微元 (Fracture Elements)
    int nU = 5, nV = 5;
    std::cout << "-> Meshing fractures (" << nU << "x" << nV << " elements per fracture)..." << std::endl;
    manager.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

    // === 初始化全局索引 ===
    std::cout << "-> Setting up global solver indices..." << std::endl;
    manager.setupGlobalIndices();

    // 5. 裂缝-裂缝求交 (F-F Intersection)
    // -----------------------------------------------------
    // 使用八叉树加速
    // [Assumption]: FFIntersectionStrategy::Octree_Optimized 是有效的枚举值
    std::cout << "-> Detecting Fracture-Fracture Intersections..." << std::endl;
    manager.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);

    // 6. 裂缝-基岩求交 (F-M Intersection - The Core)
    // -----------------------------------------------------
    // 使用我们刚落地的 SolveIntersection3D
    std::cout << "-> Solving Fracture-Matrix Intersections (Improved 3D-EDFM)..." << std::endl;
    manager.SolveIntersection3D_improved_twist_accleration(MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    std::cout << "-> Post-processing intersection pairs..." << std::endl;
    manager.removeDuplicateInteractions();     // 清除幽灵交互和重复项
    manager.resolveCoplanarInteractions();     // 解决共面重叠
    manager.buildTopologyMaps();               // 构建双向映射

    // 7. 导出结果 (Export Results)
    // -----------------------------------------------------
    std::string prefix = "Test/MeshTest/3D-EDFM/Result_" + caseName + "_";

    std::cout << "-> Exporting data to files with prefix '" << prefix << "'..." << std::endl;

    // 导出网格几何信息 (.txt)
    manager.exportMeshInfortoTxt(prefix);

    // 导出裂缝网络信息 (.txt / .csv)
    manager.exportFracturesNetworkInfortoTxt(prefix);

    // 导出最重要的交互对信息 (.csv)
    manager.exportInteractionPairsToCSV(prefix + "_InteractionPairs.csv");

    // 导出交互多边形.txt,用于可视化
    manager.exportInteractionPolygonsToTxt(prefix + "_InteractionPolygons.txt");
    
    // 0. 导出 BruteForce 视角 (全红，填满整个基岩)
    manager.exportSearchSpaceToTxt(prefix + "_Debug_BruteForce.txt", 1,
        MeshManager_3D::IntersectionStrategy::BruteForce);

    // 1. 导出光栅化视角
    manager.exportSearchSpaceToTxt(prefix + "_Debug_Raster.txt", 1,
        MeshManager_3D::IntersectionStrategy::Rasterization_14DOP);

    // 2. 导出 Octree+14DOP 视角
    manager.exportSearchSpaceToTxt(prefix + "_Debug_Octree14DOP.txt", 1,
        MeshManager_3D::IntersectionStrategy::Octree_With_14DOP);

    // 3. 导出纯 Octree 视角 (基准)
    manager.exportSearchSpaceToTxt(prefix + "_Debug_OctreePure.txt", 1,
        MeshManager_3D::IntersectionStrategy::Octree_Optimized);

    auto endTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTotal - startTotal;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "Test Case [" << caseName << "] Completed in " << elapsed.count() << " seconds." << std::endl;
    std::cout << "#########################################################\n" << std::endl;
}

void RunTestCase_for_notwist(const std::string& caseName, const std::vector<Fracture_2D>& fractures)
{
    std::cout << "\n#########################################################" << std::endl;
    std::cout << "Running Test Case: " << caseName << std::endl;
    std::cout << "#########################################################" << std::endl;

    auto startTotal = std::chrono::high_resolution_clock::now();

    // 1. 初始化网格管理器 (3D Mesh Manager)
    // -----------------------------------------------------
    // 域范围: 100m x 100m x 100m
    // 网格数: 10 x 10 x 10 = 1000 Cells
    // 网格类型
    double lengthX = 100, lengthY = 100, lengthZ = 100;
    int sectionNumX = 3, sectionNumY = 3, sectionNumZ = 3;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    std::cout << "-> Initializing MeshManager_3D (" << lengthX << "x" << lengthY << "x" << lengthZ << ")..." << std::endl;
    MeshManager_3D manager(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);

    // 2. 生成基岩网格 (Matrix Mesh)
    // -----------------------------------------------------
    // 使用正交修正 (Orthogonal)
    std::string meshName = caseName + "_Matrix";
    manager.BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod::OrthogonalCorrection, meshName);

    // 3. 添加裂缝 (Add Fractures)
    // -----------------------------------------------------
    std::cout << "-> Adding " << fractures.size() << " fractures to network..." << std::endl;
    for (const auto& frac : fractures)
    {
        manager.addFracturetoFractureNetwork(frac);
    }

    // 4. 裂缝离散化 (Mesh Fractures)
    // -----------------------------------------------------
    // 将每条裂缝划分为 5x5 的微元 (Fracture Elements)
    int nU = 5, nV = 5;
    std::cout << "-> Meshing fractures (" << nU << "x" << nV << " elements per fracture)..." << std::endl;
    manager.meshAllFracturesinNetwork(nU, nV, NormalVectorCorrectionMethod::OrthogonalCorrection);

    // 5. 裂缝-裂缝求交 (F-F Intersection)
    // -----------------------------------------------------
    // 使用八叉树加速
    // [Assumption]: FFIntersectionStrategy::Octree_Optimized 是有效的枚举值
    std::cout << "-> Detecting Fracture-Fracture Intersections..." << std::endl;
    manager.DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy::Octree_Optimized);

    // 6. 裂缝-基岩求交 (F-M Intersection - The Core)
    // -----------------------------------------------------
    // 使用我们刚落地的 SolveIntersection3D
    std::cout << "-> Solving Fracture-Matrix Intersections (Improved 3D-EDFM)..." << std::endl;
    manager.SolveIntersection3D_improved(MeshManager_3D::IntersectionStrategy::Octree_With_14DOP);

    // 7. 导出结果 (Export Results)
    // -----------------------------------------------------
    std::string prefix = "Result_" + caseName + "_";

    std::cout << "-> Exporting data to files with prefix '" << prefix << "'..." << std::endl;

    // 导出网格几何信息 (.txt)
    manager.exportMeshInfortoTxt(prefix);

    // 导出裂缝网络信息 (.txt / .csv)
    manager.exportFracturesNetworkInfortoTxt(prefix);

    // 导出最重要的交互对信息 (.csv)
    manager.exportInteractionPairsToCSV(prefix + "InteractionPairs.csv");

    // 导出交互多边形.txt,用于可视化
    manager.exportInteractionPolygonsToTxt(prefix + "InteractionPolygons.txt");


    auto endTotal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTotal - startTotal;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "Test Case [" << caseName << "] Completed in " << elapsed.count() << " seconds." << std::endl;
    std::cout << "#########################################################\n" << std::endl;
}

//用于测试 不同策略以及不同裂缝网络下Improved 3D-EDFM 核心求交性能变化
int Improved_EDFM_test_3D()
{
    try {
        std::cout << "=========================================================" << std::endl;
        std::cout << "       Improved 3D-EDFM System Verification Suite        " << std::endl;
        std::cout << "=========================================================" << std::endl;

        // -----------------------------------------------------
        // Case 1: 单一非扭曲裂缝 (Single Planar)
        // -----------------------------------------------------
        // 垂直平面 x=50, y=[20,80], z=[10,90]
        std::vector<Fracture_2D> case1_fracs;
        case1_fracs.push_back(CreateFracture(1,
            Vector(50, 20, 10), Vector(50, 80, 10),
            Vector(50, 80, 90), Vector(50, 20, 90)));

        RunTestCase("Case1_SinglePlanar", case1_fracs);

        // -----------------------------------------------------
        // Case 2: 单一扭曲裂缝 (Single Twisted)
        // -----------------------------------------------------
        // 基础为平面 z=50，但将一个角点拉高到 z=90，形成双曲抛物面
        // P1(20,20,50), P2(80,20,50), P3(80,80,90) <--- Twisted, P4(20,80,50)
        std::vector<Fracture_2D> case2_fracs;
        case2_fracs.push_back(CreateFracture(1,
            Vector(20, 20, 22), Vector(80, 20, 22),
            Vector(80, 80, 90), Vector(20, 80, 22)));

        RunTestCase("Case2_SingleTwisted", case2_fracs);

        // -----------------------------------------------------
        // Case 3: 多根非相交裂缝 (Multi Non-Intersecting)
        // -----------------------------------------------------
        // 两条平行的垂直裂缝: x=30 和 x=70
        std::vector<Fracture_2D> case3_fracs;
        case3_fracs.push_back(CreateFracture(1,
            Vector(30, 20, 10), Vector(30, 80, 10),
            Vector(30, 80, 90), Vector(30, 20, 90)));
        case3_fracs.push_back(CreateFracture(2,
            Vector(70, 20, 10), Vector(70, 80, 10),
            Vector(70, 80, 90), Vector(70, 20, 90)));

        RunTestCase("Case3_MultiParallel", case3_fracs);

        // -----------------------------------------------------
        // Case 4: 多根相交裂缝 (Multi Intersecting)
        // -----------------------------------------------------
        // 十字交叉: Fracture A (x=50), Fracture B (y=50)
        std::vector<Fracture_2D> case4_fracs;
        case4_fracs.push_back(CreateFracture(1,
            Vector(50, 10, 10), Vector(50, 90, 10),
            Vector(50, 90, 90), Vector(50, 10, 90))); // Plane X=50
        case4_fracs.push_back(CreateFracture(2,
            Vector(10, 50, 10), Vector(90, 50, 10),
            Vector(90, 50, 90), Vector(10, 50, 90))); // Plane Y=50

        RunTestCase("Case4_MultiIntersecting", case4_fracs);

        

        // -----------------------------------------------------
        // Case 5: 两根扭转相交裂缝 (Multi Intersecting)
        // -----------------------------------------------------
        std::vector<Fracture_2D> case5_fracs;
        case5_fracs.push_back(CreateFracture(1,
            Vector(20, 10, 10), Vector(80,10,75),
            Vector(80, 70, 10), Vector(20,80,45)));
        case5_fracs.push_back(CreateFracture(2,
            Vector(55, 10, 10), Vector(55, 70, 10),
            Vector(30, 70, 88), Vector(70, 10, 68)));
        RunTestCase("Case5_TwoVetrix", case5_fracs);



        // -----------------------------------------------------
        // Case 6: 多根扭转相交裂缝 (Multi Intersecting)
        // -----------------------------------------------------
        std::vector<Fracture_2D> case6_fracs;
		case6_fracs.push_back(CreateFracture(1,
			Vector(0.0, 0.0, 50.0), Vector(100.0, 0.0, 20.0),
			Vector(100.0, 100.0, 20.0), Vector(0.0, 100.0, 0.0)));
        case6_fracs.push_back(CreateFracture(2,
            Vector(25.0, 0.0, 50.0), Vector(75.0, 0.0, 30.0),
            Vector(75.0, 100.0, 70.0), Vector(25.0, 100.0, 50.0)));
        case6_fracs.push_back(CreateFracture(3,
            Vector(50.0, 0.0, 20.0), Vector(50.0, 100.0, 20.0),
            Vector(50.0, 100.0, 80.0), Vector(50.0, 0.0, 80.0)));
        case6_fracs.push_back(CreateFracture(4,
            Vector(50.0, 0.0, 20.0), Vector(50.0, 100.0, 20.0),
            Vector(70.0, 100.0, 80.0), Vector(35.0, 0.0, 80.0)));
        case6_fracs.push_back(CreateFracture(5,
            Vector(0.0, 0.0, 80.0), Vector(100.0, 0.0, 80.0),
            Vector(100.0, 100.0, 80.0), Vector(0.0, 100.0, 80.0)));
        RunTestCase("Case6_MultiVetrix", case6_fracs);
    }
    catch (const std::exception& e) {
        std::cerr << "[Fatal Error]: " << e.what() << std::endl;
        return -1;
    }
    catch (...) {
        std::cerr << "[Fatal Error]: Unknown exception occurred." << std::endl;
        return -2;
    }

    return 0;
}