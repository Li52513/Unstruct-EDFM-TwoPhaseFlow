/*
 
 * @file    main.cpp
 * @brief   当前实现基岩与裂缝网格划分与相关参数的计算
 * @date    2025-04-29
 
 * --------------------------------------------------------------------------------------
 * @note 历史版本  修改日期    修改内容                                   待修改内容         
 * @note   v1.0   2025-04-28   1.添加了边界虚拟单元 GhostCells            1.整理现在代码将mesh场捋清楚，并解决特殊工况
 *                             2.计算了非结构化网格非正交参数（A  T  E）
 * 
 * @note   v2.0   2025-04-29   1.对基岩和裂缝网格划分功能封装进了MeshManager中  1.子函数解决特殊工况
 *                             2.分离了裂缝CI_geo和CI_phys                      2. 构建物性参数类，实现基岩和裂缝物性参数的赋值
 * 
 * 
 */



#include <chrono>
#include <iostream>
#include "Mesh.h"
#include "FractureNetwork.h"
#include "Fluid.h"
#include "Matrix.h"
#include "pressureinit.h"
#include "MeshManager.h"
#include "PhysicalPropertiesManager.h"



int main()
{

    /***************************全局参数定义区*******************************/
    /*----------------------------------------------------------------------*/
    double lengthX = 1, lengthY = 1, lengthZ = 0;  
    int sectionNumX =1, sectionNumY =1, sectionNumZ =0;
	bool usePrism = true;       /// true-使用棱柱单元 false-使用四面体单元
	bool useQuadBase = false;    /// true-是用四边形底面 false -使用三角形底面
    /*----------------------------------------------------------------------*/

    /**************************网格设置模块******************************/
    /*----------------------------------------------------------------------*/
    // 1) 构造并预处理网格
	auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法
	auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
	auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	std::cout << "MatrixMesh built in " << ms0 << " ms.\n";

	// 2) 添加裂缝 & 裂缝网格划分 & 计算几何耦合系数
    mgr.addFracture({ 0.1,0.2,0 }, { 0.3,0.9,0 });
    mgr.addFracture({ 0.7,0.1,0 }, { 0.1,0.8,0 });
    mgr.addFracture({ 0.0725,0.1825,0 }, { 0.4025,0.3925,0 });

    //@brief 生成随机 DFN 裂缝网络
 /**
  * @brief 基于 DFN 方法随机生成裂缝
  * @param N            要生成的裂缝数量
  * @param minPoint     裂缝中心坐标下限 (x,y,z)
  * @param maxPoint     裂缝中心坐标上限 (x,y,z)
  * @param Lmin         裂缝最小长度
  * @param Lmax         裂缝最大长度
  * @param alpha        长度幂律指数 (p(L) ∝ L^{-α})
  * @param kappa        von Mises 浓度 (κ≈0→均匀取向)
  * @param avoidOverlap 是否简单避让已有裂缝重叠
  */
	auto t2 = std::chrono::high_resolution_clock::now(); // 计时开始

	mgr.setDFNRandomSeed(12345); // 设置随机种子
	mgr.generateDFN
	(
		/*N=*/100,
		/*minPoint=*/{ 0.0,0.0,0.0 },
		/*maxPoint=*/{ 1.0,1.0,0.0 },
		/*Lmin=*/0.5,
		/*Lmax=*/1.4,
		/*alpha=*/0,
		/*kappa=*/0,
		/*avoidOverlap=*/true
	);
    auto t8 = std::chrono::high_resolution_clock::now();
    mgr.DetectAndSubdivideFractures();
    auto t9 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms10 = std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count();
    std::cout << "FractureDetect in " << ms10 << " ms.\n";
    mgr.ComputeFractureGeometryCouplingCoefficient();
  
	//auto t3 = std::chrono::high_resolution_clock::now(); // 计时结束
	//auto ms1 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	//std::cout << "FractureNetwork built in " << ms1 << " ms.\n";

    /**************************物性参数设置模块******************************/
    /*----------------------------------------------------------------------*/
	auto t4 = std::chrono::high_resolution_clock::now(); // 计时开始
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
                {0.5, 0.5, 0.0}
            }
            },
        // 低渗区
        { Cell::RegionType::Low,
        {
            {0.5, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.5, 0.0}
        }
        }
        },
        /* defaultRegion = */ Cell::RegionType::Medium
    );
	ppm.classifyFractureElementsByGeometry(mgr, 0, { 0.1, 0.2, 0 }, { 0.3, 0.9, 0 }, FractureElementType::Blocking, FractureElementType::Conductive);
    //初始化固相、液相和气相参数（基岩和裂缝内：初始温度为303.15K，初始压力为1e6，初始水相饱和度为0.8，初始CO2饱和度为0.2） 目前的初始化温度是通过构造函数初始化的
    ppm.InitializeRockMatrixProperties(mgr);
    ppm.InitializeFractureElementsProperties(mgr);
	auto t5 = std::chrono::high_resolution_clock::now(); // 计时结束
	auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count();
	std::cout << "Physical properties initialized in " << ms2 << " ms.\n";
    /****************** 调试 & 输出 ******************/
    // 打印所有 Cell 和 FractureElement 的当前物性
    ppm.debugPrintProperties(mgr);
    // 导出网格、裂缝信息
    mgr.exportMesh("mesh");
    mgr.exportFractures("fractures");
    mgr.printFractureInfo();
	mgr.printCISourceTerms();
    cout << "Finished initial setup and property assignment.\n";
    return 0; 

    //
//  用 PhysicalPropertiesManager 中的物性来做物性耦合
//    （暂时接口还没改，这里先把它拆成两步调用原来函数）
//
////     // 5) 物性参数赋值
//    mgr.assignRockProperties(matrix);
//    mgr.assignFractureProperties(0, porosity_f_1, permeability_f_1, compressibility_f_1, aperture_1);
//    mgr.assignFractureProperties(1, porosity_f_2, permeability_f_2, compressibility_f_2, aperture_2);
//
//    // 6) 计算 CI（裂缝–基岩） & TI（裂缝–裂缝）
//   // ) 再做物性耦合（CI_phys、aW/aE/aP、alpha_phys）
//    mgr.computeFracturePhysicalCoupling(fluid, matrix);
//    mgr.computeTI(fluid);

}

//   // 构造 BC 列表
//   vector<BoundaryCondition> bcs;
   //for (int fid : groups.bottom) bcs.push_back({ BoundaryType::Dirichlet, fid, 2e5 }); //Directlet边界条件的VALUE  就是给定值
//   for (int fid : groups.left) bcs.push_back({ BoundaryType::Dirichlet, fid, 1e5 });
//   for (int fid : groups.right) bcs.push_back({ BoundaryType::Dirichlet, fid, 1e5 });

//   // Neumann（第二类）在下边界，指定通量值 5e3
//   for (int fid : groups.bottom)
//       bcs.push_back({ BoundaryType::Neumann, fid, 5e3 });
//   

//   // **调试输出—检查边界条件到底写到哪个 Cell 上了**
//   cout << "\n=== Applied boundary conditions ===\n";
//   for (const auto& bc : bcs) 
//   {
//       const auto& face = mesh.faces[bc.faceId - 1];
//       int owner = face.ownerCell;
//       const auto& cell = mesh.cells[mesh.cellId2index.at(owner)];
//           cout
//           << "Face " << bc.faceId
//           << " (" << (bc.type == BoundaryType::Dirichlet ? "Dirichlet" :
//               bc.type == BoundaryType::Neumann ? "Neumann" : "Robin")
//           << "), value=" << bc.value
//           << "  to  Cell " << owner
//           << ": pressure=" << cell.pressure
//           << ", sourceTerm=" << cell.sourceTerm
//           << ", faceDiscreCoef=" << cell.faceDiscreCoef
//           << "\n";
//   }

/*未封装前的网格划分代码*/
 //   // ========== 1. 构造网格 ==========
 //   // 构造 Mesh 对象，并构建网格
 //   Mesh mesh;
 //   mesh.BuildMesh(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    //// 根据每个单元所包含的面判断内部单元与边界单元
 //   mesh.ClassifySolidMatrixCells();
 //   // 拿到所有边界面并分配各自位置*******  这行代码在处理边界条件是要用
 //   auto groups = BoundaryClassify::ClassifySolidMatrixMeshBoundaryFaces(mesh, lengthX, lengthY);
 //   //auto bfaceIDs = mesh.getBoundaryFaceIDs(); 
 //   mesh.CreateSolidMatrixGhostCells();
 //   // 4) 这时候再进行 A=E+T 分解
 //   mesh.ComputeSolidMatrixMeshFaceGeometricInfor(CorrectionMethod::OrthogonalCorrection); //MinimumCorrection,OrthogonalCorrection,OverRelaxed
 //   // **调试输出—检查面分解结果**
 //   cout << "\n=== Face decomposition (A, |E|, |T|) ===\n";
 //   for (const auto& face : mesh.faces) 
 //   {
 //  
 //       Vector A_vec = face.normal * face.length;
 //       double A_mag = A_vec.Mag();
 //       double E_mag = face.vectorE.Mag();
 //       double T_mag = face.vectorT.Mag();
 //       std::cout
 //           << "Face " << face.id
 //           << ": |A|=" << A_mag
 //           << ", |E|=" << E_mag
 //           << ", |T|=" << T_mag
 //           << "  (E+A·e dot/|d|?)\n";
 //   }
    ////给基岩网格分配物性参数
 //   mesh.assignRockProperties(matrix);
 //   // 输出分配后的基岩物性参数
 //   cout << "\n―― 基岩物性参数分配 ――" << endl;
 //   for (const auto& cell : mesh.cells)
 //   {
 //       cout << "Cell " << cell.id << " 孔隙度: " << cell.materialProps.matrix_Porosity
 //           << ", 渗透率: " << cell.materialProps.matrix_Permeability
 //           << ", 压缩系数: " << cell.materialProps.matrix_Beta << endl;
 //   }
    //// 输出及可视化网格基本信息
 //   mesh.printMeshInfo();
 //   mesh.exportToTxt("mesh");
    // 输出部分网格信息（节点、单元、面）
 //========== 2. 构造裂缝网络 ==========
 //   FractureNetwork fractureNetwork;
 //   fractureNetwork.addFracture(Vector(0.1, 0.2, 0.0), Vector(0.3, 0.8, 0.0));  // 裂缝 1
 //   fractureNetwork.addFracture(Vector(0.7, 0.1, 0.0), Vector(0.1, 0.8, 0.0));  // 裂缝 2
 //   fractureNetwork.addFracture(Vector(0.1, 0.5, 0.0), Vector(0.8, 0.5, 0.0));  // 裂缝 3
 //    2) 构造裂缝几何：交点 + 分段（但不计算离散系数）
 //   fractureNetwork.detectFractureIntersections();
 //   for (auto& F : fractureNetwork.fracture_network)
 //       F.DetectFracturetoMeshFaceIntersections(mesh.faces);
 //   fractureNetwork.DeduplicateAndRenumberFractureToFractureIntersections();
 //   for (auto& F : fractureNetwork.fracture_network)
 //       F.sortAndRenumberIntersections();
 //   for (auto& F : fractureNetwork.fracture_network)
 //       F.subdivide(mesh.cells, mesh.nodesMap);
 //    ==========7. 基岩与裂缝初始化 ==========
 //   initializePressure(mesh.cells, mesh.cellId2index, fractureNetwork.fracture_network, 1e5);
 //    ==========5. 裂缝物性参数分配 ==========
 //   fractureNetwork.setFractureProperties(0, porosity_f_1, permeability_f_1, compressibility_f_1, aperture_1);
 //   fractureNetwork.setFractureProperties(1, porosity_f_2, permeability_f_2, compressibility_f_2, aperture_2);
 //    4) 现在做离散系数计算
 //   for (auto& F : fractureNetwork.fracture_network)
 //   {
 //       F.computeFractureDiscreteCoefficients(fluid, matrix, mesh);
 //   }
 //    ========== 4. 输出裂缝信息 ==========
 //   fractureNetwork.printFractureInfo();
 //    ==========6. 可视化裂缝信息 ==========
 //   fractureNetwork.exportToTxt("fracture_network");