#pragma once
#include <chrono>
#include <string>
#include <vector>
#include <iostream>

#include "MeshManager.h"

//* @brief 运行单个2D-EDFM测试用例的通用流程
//* @param caseName 测试用例名称 (用于日志和文件名)
//* @param fractures 该用例包含的裂缝列表
//* @param faceVector_method 裂缝面法向修正方法
//* @param dis_method 距离度量方法
//* @param search_strategy 裂缝求交搜索策略

void run2D_EDFM_test(
    const std::string& caseName, 
    const std::vector<Fracture>& fractures,
    const NormalVectorCorrectionMethod& faceVector_method,
	const DistanceMetric& dis_method,
	const IntersectionSearchStrategy_2D& search_strategy
)
{
    std::cout << "\n\n#########################################################" << std::endl;
    std::cout << "Running Test Case: " << caseName << std::endl;

	
    /**************************基岩网格几何及网格划分******************************/
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
    /**************************3D基岩网格绘制并导出******************************/
    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "          3D-EDFM Matrix mesh is generating               " << std::endl;
    std::cout << " --------------------------------------------------------" << std::endl;
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
	mgr.BuildSolidMatrixGrid_2D(faceVector_method);
	mgr.exportMesh(caseName + "_MatrixMesh");

    /**************************添加裂缝******************************/
    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "                   1D-Fractions is adding                 " << std::endl;
    std::cout << " --------------------------------------------------------" << std::endl;
	int frac_num = fractures.size();
	for (int i = 0; i < frac_num; ++i)
	{
        mgr.addFracture(fractures[i].start, fractures[i].end);
	}

    /**************************基于基岩网格实现裂缝网格划分******************************/
    std::cout << "\n--------------------------------------------------------" << std::endl;
    std::cout << "                   1D-Fractions is meshing                 " << std::endl;
    std::cout << " --------------------------------------------------------" << std::endl;
    mgr.setDistanceMetric(dis_method);
	auto startTotal = std::chrono::high_resolution_clock::now();
    mgr.DetectAndSubdivideFractures(search_strategy);
	auto endTotal = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = endTotal - startTotal;

    /**************************EDFM-几何耦合系数计算******************************/
	std::cout << "\n--------------------------------------------------------" << std::endl;
	std::cout << "              2D-EDFM Geometric Coupling is calculating            " << std::endl;
	std::cout << " --------------------------------------------------------" << std::endl;
    mgr.ComputeFractureGeometryCouplingCoefficient();

	/**************************裂缝网格结果导出******************************/
    std::cout << "\n--------------------------------------------------------" << std::endl;
	std::cout << "              1D-Fractions mesh information is exporting            " << std::endl;
	std::cout << " --------------------------------------------------------" << std::endl;
	mgr.exportFractures(caseName + "_FractureMesh");

	
	
	std::cout << "--------------------------------------------------------" << std::endl;
	std::cout << "Test Case [" << caseName << "] Completed in " << elapsed.count() << " seconds." << std::endl;
	std::cout << "#########################################################\n" << std::endl;
}

int EDFM_test_2D()
{	
	// Case 1: 单一裂缝-方案A-暴力遍历
    std::vector<Fracture> fractures1;
	fractures1.push_back(Fracture(Vector(0.25, 0.24, 0), Vector(0.75, 0.76, 0)));
	run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case1_SingleFracture_BruteForce", fractures1, NormalVectorCorrectionMethod::OrthogonalCorrection, 
		DistanceMetric::CrossAwareGauss, IntersectionSearchStrategy_2D::BruteForce);

	//Case 2: 单一裂缝-方案B-GlobalAABB
    std::vector<Fracture> fractures2;
    fractures2.push_back(Fracture(Vector(0.25, 0.24, 0), Vector(0.75, 0.76, 0)));
	run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case2_SingleFracture_GlobalAABB", fractures2, NormalVectorCorrectionMethod::OrthogonalCorrection, 
		DistanceMetric::CrossAwareGauss, IntersectionSearchStrategy_2D::GlobalAABB);

	//Case 3:单一裂缝-方案C-GridIndexing
	std::vector<Fracture> fractures3;
	fractures3.push_back(Fracture(Vector(0.25, 0.24, 0), Vector(0.75, 0.76, 0)));
	run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case3_SingleFracture_GridIndexing", fractures3, NormalVectorCorrectionMethod::OrthogonalCorrection,
		DistanceMetric::CrossAwareGauss, IntersectionSearchStrategy_2D::GridIndexing);

	//Case 4:单一裂缝-方案D-GridIndexing_BasedOn8DOP
	std::vector<Fracture> fractures4;
	fractures4.push_back(Fracture(Vector(0.25, 0.24, 0), Vector(0.75, 0.76, 0)));
	run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case4_SingleFracture_GridIndexing_BasedOn8DOP", fractures4, NormalVectorCorrectionMethod::OrthogonalCorrection,
		DistanceMetric::CrossAwareGauss, IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP);

	//Case 5:单一裂缝-方案B-GridIndexing_BasedOn8DOP_DDA
	std::vector<Fracture> fractures5;
	fractures5.push_back(Fracture(Vector(0.25, 0.24, 0), Vector(0.75, 0.76, 0)));
	run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case5_SingleFracture_GridIndexing_BasedOn8DOP_DDA", fractures5, NormalVectorCorrectionMethod::OrthogonalCorrection,
		DistanceMetric::CrossAwareGauss, IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);

	return 0;
}

int EDFM_DFN_test_2D()
{
    // =========================================================
    // Step 1: 预生成统一的随机裂缝网络 (DFN)
    // =========================================================
    std::cout << ">>> Generating Random DFN Network..." << std::endl;

    // 创建一个临时的 Manager 仅用于生成裂缝数据
    // 尺寸必须与 run2D_EDFM_test 中的 lengthX, lengthY 一致 (1.0 x 1.0)
    double Lx = 1.0, Ly = 1.0;
    MeshManager genMgr(Lx, Ly, 0.0, 1, 1, 0, false, false);

    // DFN 参数设置
    int N_fractures = 50;           // 裂缝数量 (建议 50-200 以测试性能)
    unsigned int seed = 2024;       // 固定种子，保证每次运行结果一致
    double Lmin = 0.1, Lmax = 0.4;  // 裂缝长度范围
    double alpha = 2.5;             // 幂律指数
    double kappa = 0.0;             // 0.0 表示完全随机方向
    bool avoidOverlap = true;       // 【关键】开启避免重合检测

    genMgr.setDFNRandomSeed(seed);
    genMgr.generateDFN(
        N_fractures,
        Vector(0, 0, 0),   // 最小坐标 (基岩左下角)
        Vector(Lx, Ly, 0), // 最大坐标 (基岩右上角)
        Lmin, Lmax, alpha, kappa, avoidOverlap
    );

    // 提取生成的裂缝列表
    // 注意：这里我们拷贝一份 fracture_network 中的裂缝数据
    std::vector<Fracture> dfnFractures = genMgr.fracture_network().fractures;

    if (dfnFractures.empty()) {
        std::cerr << "[Error] DFN generation failed or generated 0 fractures." << std::endl;
        return -1;
    }
    std::cout << ">>> DFN Generated. Total fractures: " << dfnFractures.size() << "\n" << std::endl;


    // =========================================================
    // Step 2: 使用同一套 DFN 数据运行对比测试
    // =========================================================

    // Case 1: 暴力遍历 (BruteForce)
    // 这是基准，用于验证正确性，但在裂缝多时会非常慢
    run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case1_DFN_BruteForce", dfnFractures,
        NormalVectorCorrectionMethod::OrthogonalCorrection,
        DistanceMetric::CrossAwareGauss,
        IntersectionSearchStrategy_2D::BruteForce);

    // Case 2: GlobalAABB (全局包围盒预筛)
    run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case2_DFN_GlobalAABB", dfnFractures,
        NormalVectorCorrectionMethod::OrthogonalCorrection,
        DistanceMetric::CrossAwareGauss,
        IntersectionSearchStrategy_2D::GlobalAABB);

    // Case 3: GridIndexing (背景网格索引 + AABB)
    run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case3_DFN_GridIndexing", dfnFractures,
        NormalVectorCorrectionMethod::OrthogonalCorrection,
        DistanceMetric::CrossAwareGauss,
        IntersectionSearchStrategy_2D::GridIndexing);

    // Case 4: GridIndexing + 8DOP (背景网格 + 8-DOP 精筛)
    run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case4_DFN_GridIndexing_8DOP", dfnFractures,
        NormalVectorCorrectionMethod::OrthogonalCorrection,
        DistanceMetric::CrossAwareGauss,
        IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP);

    // Case 5: GridIndexing + DDA + 8DOP (终极方案：DDA 游走 + 8-DOP 精筛)
    run2D_EDFM_test("Test/MeshTest/2D-EDFM/Case5_DFN_GridIndexing_DDA_8DOP", dfnFractures,
        NormalVectorCorrectionMethod::OrthogonalCorrection,
        DistanceMetric::CrossAwareGauss,
        IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);

    return 0;
}
