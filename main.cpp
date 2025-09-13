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
#include "Initializer.h"
#include "PostProcessor.h"
#include "CouplingAssembler.h" 
#include "BoundaryFaceClassify.h"


int main()
{

    /***************************全局参数定义区*******************************/
    /*----------------------------------------------------------------------*/
    double lengthX = 1, lengthY = 1, lengthZ = 0;  
    int sectionNumX =3, sectionNumY =3, sectionNumZ =0;
	bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
	bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
    /*----------------------------------------------------------------------*/

    /**************************网格设置模块******************************/
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

	//mgr.setDFNRandomSeed(12345); // 设置随机种子
	//mgr.generateDFN
	//(
	//	/*N=*/2,
	//	/*minPoint=*/{ 0.0,0.0,0.0 },
	//	/*maxPoint=*/{ 1.0,1.0,0.0 },
	//	/*Lmin=*/0.5,
	//	/*Lmax=*/1.4,
	//	/*alpha=*/0,
	//	/*kappa=*/0,
	//	/*avoidOverlap=*/true
	//);
    auto t8 = std::chrono::high_resolution_clock::now();

	mgr.setDistanceMetric(DistanceMetric::CrossAwareGauss); // 设置距离度量方式AreaWeight CellCenter NodeAverage CrossAwareGauss
    mgr.DetectAndSubdivideFractures();

    auto t9 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms10 = std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count();
    std::cout << "FractureDetect in " << ms10 << " ms.\n";
	mgr.ComputeFractureGeometryCouplingCoefficient();   //明确一点： 因为这里还没有给出物性参数，这里计算的只是几何耦合系数 CI_geo和 geomAlpha 各自的表达式分别为：CI_geo = L*1/d_avg geomAlpha = 2 / L;

    /**************************物性参数设置模块******************************/
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


	// 4) 主变量场初始化（温度T,水相饱和度Sw,水相压力Pw)
    FieldRegistry reg;      // 基岩场
    FieldRegistry reg_fr;   // 裂缝场（与你当前实现一致，用同一类型的注册表管理）
    InitFields ic;          // p0/T0/Sw0 及其梯度（默认为均匀）
    VGParams vg;            // vG 参数
    RelPermParams rp;       // 相对渗透率参数（默认 L=0.5）
    RockDefaults rock;      // 基岩默认热物性
    InitDiagnostics diag;   // 诊断统计
    // 5) 基岩“主变量场”创建与填充：p_w, S_w, T (+ p_c, p_g, kr_w, kr_g)

    //创建基岩的主变量场
    Initializer::createPrimaryFields(mgr.mesh(), reg); 

    //填充基岩主变量场
    Initializer::fillBaseDistributions(mgr.mesh(), reg, ic); 

    //对水相饱和度进行限幅并记录限制修改的次数
    Initializer::enforceSaturationBounds(reg, vg, diag); 

    //计算闭合关系包括毛细压力和相对渗透率
    Initializer::computeClosure(mgr.mesh(), reg, vg, rp, diag);

    // 6) 基岩有效热物性（C_eff, lambda_eff）
    Initializer::computerEffectiveThermals(mgr.mesh(), reg, rock, diag);

    // 7) 裂缝段主变量场：pf_w, Sf_w, Tf（从宿主基岩单元拷贝）
    Initializer::initFracturePrimaries(mgr.mesh(), mgr.fracture_network(), reg, reg_fr);

    // 8) 基岩/裂缝固相场：由 region/type + (P,T) 通过你的 computeSolidProperties() 得到
    ppm.InitializeRockMatrixProperties(mgr, reg);
    ppm.InitializeFractureElementsProperties(mgr, reg_fr);

    // 9) 基岩/裂缝流体物性场：查 WaterPropertyTable / CO2PropertyTable
    //    - 基岩用 p_w/p_g/T
    //    - 裂缝用 pf_w/Tf 与 pg_f = pf_w + pc_vG(Sf_w)

	ppm.MatrixFluidPropertiesTest(315.54930556, 5027664.1668); // 测试水和 CO2 物性表

    ppm.InitializeMatrixFluidProperties(mgr, reg, vg);
    ppm.InitializeFractureFluidProperties(mgr, reg, reg_fr, vg);

    //边界面识别
    const auto& bfaces = mgr.boundaryFaces();

    // 你可以先自检打印（便于验证识别是否准确）
    std::cout << "[BC] x0=" << bfaces.x0.size()
        << " xL=" << bfaces.xL.size()
        << " y0=" << bfaces.y0.size()
        << " yL=" << bfaces.yL.size()
        << " z0=" << bfaces.z0.size()
        << " zL=" << bfaces.zL.size() << "\n";


    // 打印初始化诊断
    Initializer::printDiag(diag);

    auto t5 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count();
    std::cout << "Fields & properties initialized in " << ms2 << " ms.\n";
    /***********************************************************************/

    /**************************耦合系数（CI/TI）初始化@t=0*******************/
    // 1) 先给裂缝注册 CI 场
    {
        // 计算全局段数，用于创建场
        size_t Nseg = 0;
        for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
        ensureFracCouplingFields(reg_fr, Nseg); // 创建 CIw / CIg
    }
    
    // 2) 计算 CI（矩阵↔裂缝）
    {
        bool upwind = true;            // 采用迎风（按相位势比较）
        bool include_gravity = false;  // 2D 平面暂不考虑重力项；需要的话置 true
        double g = 9.80665;
        updateMatrixFractureCI(mgr, reg, reg_fr, vg, upwind, include_gravity, g);
    }

    // 3) 计算 TI（裂缝↔裂缝，多臂交点 Star–Delta）
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



