#include <chrono>
#include <iostream>
#include "Mesh.h"
#include "FractureNetwork.h"
#include "pressureinit.h"
#include "MeshManager.h"
#include "PhysicalPropertiesManager.h"
#include "Initializer.h"
#include "PostProcessor.h"
#include "CouplingAssembler.h" 
#include "BoundaryFaceClassify.h"
#include "PressureBC.h"
#include "Diff_TPFA_Operators.h"
#include "TemperatureBC.h"
#include "TemperatureBCAdapter.h"
#include "Conv_FirstOrder_Operators.h"
#include "TimeIterm_Euler_SinglePhase_PressureEq.h"
#include "TimeIterm_Euler_SinglePhase_TemperatureEq.h"
#include "Solver_TimeLoopDriver.h"
#include "Solver_TimeLoopUtils.h" 



int main()
{

    /***************************全局参数定义区*******************************/
    /*----------------------------------------------------------------------*/
    double lengthX = 1, lengthY = 1, lengthZ = 0;  
    int sectionNumX =50, sectionNumY =50, sectionNumZ =0;
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
        /* defaultRegion = */ Cell::RegionType::Medium
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
    Initializer::createPrimaryFields(mgr.mesh(), reg); 

    //填充基岩主变量场
    Initializer::fillBaseDistributions(mgr.mesh(), reg, ic); 

    //对水相饱和度进行限幅并记录限制修改的次数
    Initializer::enforceSaturationBounds(reg, vg, diag); 

    //计算闭合关系包括毛细压力和相对渗透率 //两相流时启用
    Initializer::computeClosure(mgr.mesh(), reg, vg, rp, diag);

    //6）创建 transient 辅助场，后面时间推进会用
    //单相渗流-传热问题
	ensureTransientFields(mgr.mesh(), reg, /*p_name=*/"p_w", /*T_name=*/"T", /*p_old_name=*/"p_w_old", /*T_old_name=*/"T_old", /*p_prev_name=*/"p_w_prev", /*T_prev_name=*/"T_prev"); //p_w T代表当前场变量，*_old代表上一时步变量，*_prev代表上一迭代步的变量,水相
	ensureTransientFields(mgr.mesh(), reg, /*p_name=*/"p_g", /*T_name=*/"T", /*p_old_name=*/"p_g_old", /*T_old_name=*/"T_old", /*p_prev_name=*/"p_g_prev", /*T_prev_name=*/"T_prev");//p_g T代表当前场变量，*_old代表上一时步变量，*_prev代表上一迭代步的变量,气相

	//7）裂缝主变量初始化 （当前还未涉及裂缝）**
    Initializer::initFracturePrimaries(mgr.mesh(), mgr.fracture_network(), reg, reg_fr);

    //8）基岩物性参数设置
    //   a.固相参数
    ppm.UpdateMatrixRockAt(mgr, reg, "p_g", "T");

	//   b.流体参数
    ppm.UpdateMatrixFluidAt(mgr, reg, "p_g", "T", "CO2");

    //   c.有效热物性参数
    ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_g", "T", "CO2", 1e-12);

	//9）裂缝物性参数设置
    //   a.固相参数
	ppm.UpdateFractureRockAt(mgr, reg, reg_fr, "pf_w", "Tf");

    //   b.流体参数
	ppm.UpdateFractureFluidAt(mgr, reg, reg_fr, "pf_w", "Tf", "CO2");

	//   c.有效热物性参数 (待补充)
	
    ppm.MatrixFluidPropertiesTest(315.54930556, 5027664.1668); // 测试水和 CO2 物性表

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
    PressureBC::BoundaryCoefficient P_Left { 1.0, 0.0,8e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0,6e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Down   { 1.0, 0.0,6e6 }; // p = 2e5 Pa
    PressureBC::BoundaryCoefficient P_Up { 1.0, 0.0,8e6 }; // p = 2e5 Pa
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up );   //按照左 右 下 上 的顺序赋值
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

	// 下面进入时间推进模块

    //********************求解单相渗流-传热问题***********************//
    //控制方程包括：质量守恒（单相达西）+能量守恒（基于单相达西速度的对流扩散方程）
    //质量守恒（单相达西，代求变量为: P_w）：时间项+扩散项=源项
    //能量守恒（基于单相达西速度的对流扩散方程为: T）：时间项+对流项+扩散项=源项

    //1 设置重力大小及重力方向
//2 在对密度进行迎风性判别时，是否需要考虑重力势能
    GravUpwind gu;
    gu.g = Vector(0.0, -9.80665, 0.0);
    //gu.g = Vector(0.0,0.0, 0.0);
    gu.use_potential = true;


    // 求解器参数
    SolverControls sc;
    sc.maxOuter = 300;
    sc.tol_p_abs = 10.0;     // 举例：绝对容差 10 Pa
    sc.tol_T_abs = 1e-6;
    sc.tol_p_rel = 1e-6;     // 新增：压力相对容差
    sc.tol_T_rel = 1e-6;     // 新增：温度相对容差
	sc.urf_p = 0.4; //欠松弛因子
	sc.urf_T = 0.25; //欠松弛因子
    sc.c_phi_const = 1e-9;
    sc.jac_p = { 500, 0.8, 1e-8 };
    sc.jac_T = { 500, 0.8, 1e-8 };
    
    //sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB; // 压力：也可选 CG
    //sc.lin_p.maxIters = 5000;
    //sc.lin_p.tol = sc.tol_p_abs;
    //sc.lin_p.iluFill = 10;
    //sc.lin_p.iluDrop = 1e-4;
    //sc.lin_p.equil = true;   // 建议保持 on
    //sc.lin_p.reusePreconditioner = 3;


    //sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
    //sc.lin_T.maxIters = 5000;
    //sc.lin_T.tol = sc.tol_T_abs;
    //sc.lin_T.iluFill = 10;
    //sc.lin_T.iluDrop = 1e-4;

    sc.lin_p.type = LinearSolverOptions::Type::SparseLU;  // 压力用直接法
    sc.lin_p.equil = true;   // 默认就是 true，可不写
    // 直接法不关心 maxIters/tol/ILU 设置

    sc.lin_T.type = LinearSolverOptions::Type::SparseLU;  // 先试直接法
    sc.lin_T.equil = true;

    //sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB;
    //sc.lin_T.maxIters = 10000;
    //sc.lin_T.tol = sc.tol_T_abs;
    //sc.lin_T.iluFill = 20;
    //sc.lin_T.iluDrop = 1e-5;
    //sc.lin_T.equil = true;   // 建议保持 on
    //sc.lin_T.reusePreconditioner = 2;


    sc.useJacobi = false;
    int    nSteps = 20000;
    double dt = 10.0 / nSteps;

    // 只导出 TXT，每步一份；不导出 CSV / MM：
    runTransient_singlePhase(
        mgr, reg, freg, ppm,
        bcAdapter, TbcA, gu, rock,
        nSteps, dt, sc,
        /*phase=*/"CO2",
        /*writeEvery=*/20,
        /*onWrite=*/nullptr,
        /*exportCSV=*/false,
        /*exportTXT=*/true,
        /*exportMM=*/true,
        /*exportTecplotP=*/false,
        /*exportTecplotT=*/true
    );
    

    return 0;

}




