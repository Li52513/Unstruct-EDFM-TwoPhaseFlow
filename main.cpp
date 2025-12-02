#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "Mesh.h"
#include "FractureNetwork.h"
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
#include "BCAdapter.h"
#include "Conv_FirstOrder_Operators.h"
#include "TimeIterm_Euler_SinglePhase_PressureEq.h"
#include "TimeIterm_Euler_SinglePhase_TemperatureEq.h"
#include "Solver_TimeLoopDriver.h"
#include "Solver_TimeLoopDriver_IMPES.h"
#include "Solver_TimeLoopUtils.h" 


//#include "FVM_SourceTerm_InjectionMask.h"
//#include "FVM_SourceTerm_PeacemanModel.h"

#include "FVM_WellDOF.h"
#include "FVM_Peaceman.h"
#include "WellConfig.h"

#include "WellConfig_TwoPhase.h"


#include    "MultiPhaseProperties.h"
#include    "0_PhysicalParametesCalculateandUpdata.h"
#include    "TimeLoopDriver.h"   //分析解
#include    "IMPES_Iteration_Loop.h"


//int main1()
//{
//
//    /***************************全局参数定义区*******************************/
//    /*----------------------------------------------------------------------*/
//    double lengthX = 1, lengthY = 1, lengthZ = 0;  
//    int sectionNumX =5, sectionNumY =5, sectionNumZ =0;
//	bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//	bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
//    /*----------------------------------------------------------------------*/
//
//    /**************************网格设置模块******************************/
//    /*----------------------------------------------------------------------*/
//    // 1) 构造并预处理网格
//	auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//	auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
//	auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//	std::cout << "MatrixMesh built in " << ms0 << " ms.\n";
//
//	// 2) 添加裂缝 & 裂缝网格划分 & 计算几何耦合系数
//    mgr.addFracture({ 0.1,0.2,0 }, { 0.3,0.9,0 });
//    mgr.addFracture({ 0.7,0.1,0 }, { 0.1,0.8,0 });
//    mgr.addFracture({ 0.0725,0.1825,0 }, { 0.4025,0.3925,0 });
//
//    /*mgr.addFracture({ 3.9925,0.0066356,0 }, { 3.9925,5.0066,0 });
//    mgr.addFracture({ 3.9925,1.3439,0 }, { 1.4925,2.985,0 });
//    mgr.addFracture({ 3.9925,1.7748,0 }, { 6.4925,3.6727,0 });
//    mgr.addFracture({ 2.371,0.97711,0 }, { 2.7235,1.6617,0 });
//    mgr.addFracture({ 2.5472,1.3194,0 }, { 2.7563,1.3194,0 });
//    mgr.addFracture({ 1.3467,1.3503,0 }, { 1.4925,1.5089,0 });
//    mgr.addFracture({ 2.1825,1.5945,0 }, { 1.4925,2.004,0 });
//    mgr.addFracture({ 1.4925,1.5089,0 }, { 3.33836,1.7436,0 });
//    mgr.addFracture({ 1.4925,1.5089,0 }, { 3.33836,1.7436,0 });
//    mgr.addFracture({ 1.4925,1.7381,0 }, { 1.718,1.8702,0 });
//    mgr.addFracture({ 1.102,3.182,0 }, { 1.5179,3.1829,0 });
//    mgr.addFracture({ 1.4925,2.985,0 }, { 1.5433,3.3789,0 });
//    mgr.addFracture({ 2.4291,2.3702,0 }, { 2.6646,3.3334,0 });
//    mgr.addFracture({ 2.3067,3.3701,0 }, { 2.5469,2.8518,0 });
//    mgr.addFracture({ 2.3067,3.3701,0 }, { 2.5469,2.8518,0 });
//    mgr.addFracture({ 2.6139,3.126,0 }, { 3.0591,3.2417,0 });
//    mgr.addFracture({ 4.6805,2.2971,0 }, { 6.4925,1.8756,0 });
//    mgr.addFracture({ 5.318,2.1488,0 }, { 5.8117,1.2797,0 });
//    mgr.addFracture({ 5.8028,2.0361,0 }, {6.4925,2.2699,0 });
//    mgr.addFracture({ 6.4925,1.8756,0 }, { 6.6653,1.632,0 });
//    mgr.addFracture({ 6.4925,1.8756,0 }, { 6.7387,2.0502,0 });
//    mgr.addFracture({ 6.2888,2.25266,0 }, { 6.3178,2.2106,0 });
//    mgr.addFracture({ 5.2425,2.7237,0 }, { 5.564,4.6445,0 });
//    mgr.addFracture({ 4.6556,4.2411,0 }, { 5.3717,3.4953,0 });
//    mgr.addFracture({ 5.2887,5.0066,0 }, { 5.5234,4.4022,0 });
//	mgr.addFracture({ 5.2887,5.0066,0 }, { 5.5234,4.4022,0 });
//    mgr.addFracture({ 5.452,3.9755,0 }, { 6.3347,4.452,0 });
//    mgr.addFracture({ 6.1004,4.3255,0 }, { 6.1178,4.7011,0 });
//    mgr.addFracture({ 5.8934,4.2138,0 }, { 6.4925,4.0853,0 });
//    mgr.addFracture({ 6.3347,4.452,0 }, { 6.5304,4.8415,0 });
//    mgr.addFracture({ 6.4326,4.6467,0 }, {6.7749,4.5888,0 });
//    mgr.addFracture({2.9949,5.0066,0 }, { 3.9982,3.6268,0 });
//    mgr.addFracture({ 2.3893,4.3328,0 }, { 3.4937,4.3167,0 });
//    mgr.addFracture({ 2.7472,3.856,0 }, { 3.1882,4.3212,0 });
//    mgr.addFracture({ 3.2771,4.6164,0 }, { 3.3986,5.0066,0 });
//    mgr.addFracture({ 3.9925,4.507,0 }, { 4.4446,5.0066,0 });
//    mgr.addFracture({ 4.2186,4.7568,0 }, { 4.6672,4.7825,0 });
//    mgr.addFracture({3.6395,5.2578,0 }, { 3.9925,5.0066,0 });
//    mgr.addFracture({ 3.7324,5.1917,0 }, { 3.9312,5.3935,0 });*/
//	auto t2 = std::chrono::high_resolution_clock::now(); // 计时开始
//
//	//mgr.setDFNRandomSeed(12345); // 设置随机种子
//	//mgr.generateDFN (/*N=*/2,/*minPoint=*/{ 0.0,0.0,0.0 },/*maxPoint=*/{ 1.0,1.0,0.0 },/*Lmin=*/0.5,/*Lmax=*/1.4,/*alpha=*/0,/*kappa=*/0,/*avoidOverlap=*/true);
//    auto t8 = std::chrono::high_resolution_clock::now();
//
//	mgr.setDistanceMetric(DistanceMetric::CrossAwareGauss); // 设置距离度量方式AreaWeight CellCenter NodeAverage CrossAwareGauss
//    mgr.DetectAndSubdivideFractures();
//
//    auto t9 = std::chrono::high_resolution_clock::now(); // 计时结束
//    auto ms10 = std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count();
//    std::cout << "FractureDetect in " << ms10 << " ms.\n";
//	mgr.ComputeFractureGeometryCouplingCoefficient();   //明确一点： 因为这里还没有给出物性参数，这里计算的只是几何耦合系数 CI_geo和 geomAlpha 各自的表达式分别为：CI_geo = L*1/d_avg geomAlpha = 2 / L;
//
//    /**************************主变量初始化以及物性参数设置模块******************************/
//    /*----------------------------------------------------------------------*/
//	auto t4 = std::chrono::high_resolution_clock::now(); // 计时开始
//    
//	// 3) 基岩和裂缝介质区域划分
//    PhysicalPropertiesManager ppm;
//    ppm.classifyRockRegionsByGeometry
//    (
//        mgr,
//        {
//            // 高渗区
//          { Cell::RegionType::High,
//        
//             {
//                {0.0, 0.5, 0.0},
//                {0.0, 1.0, 0.0},
//                {0.5, 1.0, 0.0}
//             }
//          },
//        // 低渗区
//          { Cell::RegionType::Low,
//            {
//            {0.5, 0.0, 0.0},
//            {1.0, 1.0, 0.0},
//            {1.0, 0.5, 0.0}
//            }
//          }
//        },
//        /* defaultRegion = */ 
//        Cell::RegionType::Medium
//    );
//	ppm.classifyFractureElementsByGeometry(mgr, 0, { 0.1, 0.2, 0 }, { 0.3, 0.9, 0 }, FractureElementType::Blocking, FractureElementType::Conductive);
//
//
//	// 4) 生成场变量以及设置初始化参数和VG模型（温度T,水相饱和度Sw,水相压力Pw)
//    FieldRegistry reg;      // 基岩场
//    FieldRegistry reg_fr;   // 裂缝场（与你当前实现一致，用同一类型的注册表管理）
//    InitFields ic;          // p0/T0/Sw0 及其梯度（默认为均匀）
//    VGParams vg;            // vG 参数 =vg模型
//    RelPermParams rp;       // 相对渗透率参数（默认 L=0.5）
//    RockDefaults rock;      // 基岩默认热物性
//    InitDiagnostics diag;   // 诊断统计
//    
//    // 5) 基岩“主变量场”创建与填充：p_w, S_w, T (+ p_c, p_g, kr_w, kr_g)
//
//    //创建基岩的主变量场
//    Initializer::createPrimaryFields(mgr.mesh(), reg); 
//    // 测试案例-2D-常物性-温度扩散
//
//    //填充基岩主变量场
//    Initializer::fillBaseDistributions(mgr.mesh(), reg, ic); 
//
//    //对水相饱和度进行限幅并记录限制修改的次数
//    Initializer::enforceSaturationBounds(reg, vg, diag); 
//
//    //计算闭合关系包括毛细压力和相对渗透率 //两相流时启用
//    Initializer::computeClosure(mgr.mesh(), reg, vg, rp, diag);
//
//    //6）创建 transient 辅助场，后面时间推进会用
//    //单相渗流-传热问题
//	ensureTransientFields(mgr.mesh(), reg, /*p_name=*/"p_w", /*T_name=*/"T", /*p_old_name=*/"p_w_old", /*T_old_name=*/"T_old", /*p_prev_name=*/"p_w_prev", /*T_prev_name=*/"T_prev"); //p_w T代表当前场变量，*_old代表上一时步变量，*_prev代表上一迭代步的变量,水相
//	ensureTransientFields(mgr.mesh(), reg, /*p_name=*/"p_g", /*T_name=*/"T", /*p_old_name=*/"p_g_old", /*T_old_name=*/"T_old", /*p_prev_name=*/"p_g_prev", /*T_prev_name=*/"T_prev");//p_g T代表当前场变量，*_old代表上一时步变量，*_prev代表上一迭代步的变量,气相
//
//	//7）裂缝主变量初始化 （当前还未涉及裂缝）**
//    Initializer::initFracturePrimaries(mgr.mesh(), mgr.fracture_network(), reg, reg_fr);
//
//    //8）基岩物性参数设置
//    //   a.固相参数
//    ppm.UpdateMatrixRockAt(mgr, reg, "p_g", "T");
//
//	//   b.流体参数
//    ppm.UpdateMatrixFluidAt(mgr, reg, "p_g", "T", "CO2");
//
//    //   c.有效热物性参数
//    ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_g", "T", "CO2", 1e-12);
//
//	//9）裂缝物性参数设置
//    //   a.固相参数
//	ppm.UpdateFractureRockAt(mgr, reg, reg_fr, "pf_w", "Tf");
//
//    //   b.流体参数
//	ppm.UpdateFractureFluidAt(mgr, reg, reg_fr, "pf_w", "Tf", "CO2");
//
//	//   c.有效热物性参数 (待补充)
//	
//    ppm.MatrixFluidPropertiesTest(315.54930556, 5027664.1668); // 测试水和 CO2 物性表
//
//    /**************************边界条件设置模块******************************/
//
//    //边界面识别
//    const auto& bfaces = mgr.boundaryFaces();
//
//    // 你可以先自检打印（便于验证识别是否准确）
//    std::cout << "[BC] x0=" << bfaces.x0.size()
//        << " xL=" << bfaces.xL.size()
//        << " y0=" << bfaces.y0.size()
//        << " yL=" << bfaces.yL.size()
//        << " z0=" << bfaces.z0.size()
//        << " zL=" << bfaces.zL.size() << "\n";
//
//
//    /// 压力边界条件设置
//
//    //给定参数：2D情况，给定基岩四个边界条件系数
//    PressureBC::Registry pbc_pw;
//    PressureBC::BoundaryCoefficient P_Left { 1.0, 0.0,8e6 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0,6e6 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Down   { 1.0, 0.0,6e6 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Up { 1.0, 0.0,8e6 }; // p = 2e5 Pa
//    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up );   //按照左 右 下 上 的顺序赋值
//    pbc_pw.printInformationofBoundarySetting(mgr);
//    PressureBCAdapter bcAdapter{ pbc_pw };
//
//	// /  温度边界条件设置
//    TemperatureBC::Registry tbc;
//    TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 450.0 };
//    TemperatureBC::BoundaryCoefficient T_Right{ 1.0, 0.0, 373.15 };
//    TemperatureBC::BoundaryCoefficient T_Down{ 1.0, 0.0, 373.15 };
//    TemperatureBC::BoundaryCoefficient T_Up{ 1.0, 0.0, 450.0 };
//    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
//    tbc.printInformationofBoundarySetting(mgr);
//    TemperatureBCAdapter TbcA{ tbc };
//
//    /**************************基岩裂缝面系数计算 当前系数不考虑裂缝的离散**************************/
//
//
//
//	//生成储存储存网格面上离散系数和源项的面场
//	FaceFieldRegistry freg;
//
//    // 打印初始化诊断
//    Initializer::printDiag(diag);
//
//    auto t5 = std::chrono::high_resolution_clock::now(); // 计时结束
//    auto ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count();
//    std::cout << "Fields & properties initialized in " << ms2 << " ms.\n";
//    /***********************************************************************/
//
//    /**************************耦合系数（CI/TI）初始化@t=0*******************/
//    // 1) 先给裂缝注册 CI 场
//    {
//        // 计算全局段数，用于创建场
//        size_t Nseg = 0;
//        for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
//        ensureFracCouplingFields(reg_fr, Nseg); // 创建 CIw / CIg
//    }
//    
//    // 2) 计算 CI（矩阵?裂缝）
//    {
//        bool upwind = true;            // 采用迎风（按相位势比较）
//        bool include_gravity = false;  // 2D 平面暂不考虑重力项；需要的话置 true
//        double g = 9.80665;
//        updateMatrixFractureCI(mgr, reg, reg_fr, vg, upwind, include_gravity, g);
//    }
//
//    // 3) 计算 TI（裂缝?裂缝，多臂交点 Star–Delta）
//    {
//        bool include_gravity = false;
//        double g = 9.80665;
//        updateFractureFractureTI(mgr, reg_fr, vg, rp, include_gravity, g);
//    }
//    /***********************************************************************/
// 
//
//    /****************** 调试 & 输出 ******************/
//    // 打印所有 Cell 和 FractureElement 的当前物性
//    ppm.debugPrintProperties(mgr, reg, reg_fr, /*maxPrint=*/10);
//    // 导出网格、裂缝信息
//    mgr.exportMesh("mesh");
//    mgr.exportFractures("fractures");
//    mgr.printFractureInfo();
//    printCI_Diagnostics(mgr, reg, reg_fr,  /*maxFracs*/(size_t)-1, /*maxSegs*/(size_t)-1);
//    printTI_Diagnostics(mgr, reg_fr,       /*maxFracs*/(size_t)-1);
//	//mgr.printCISourceTerms();
//
//
//    cout << "Finished initial setup and property assignment.\n";
//
//    // —— 导出 t=0、step=0 的场，用于初始时刻可视化 ——
//    const double time0 = 0.0;
//    const int    step0 = 0;
//
//    PostProcessor outM("out\\matrix");
//    PostProcessor outF("out\\fracture");
//
//    // 可选自检：确保裂缝三大主变量长度与全局段数一致
//    {
//        size_t Nseg = 0; for (auto& F : mgr.fracture_network().fractures) Nseg += F.elements.size();
//        auto pf = reg_fr.get<volScalarField>("pf_w");
//        auto sf = reg_fr.get<volScalarField>("Sf_w");
//        auto Tf = reg_fr.get<volScalarField>("Tf");
//        std::cout << "[fracture] Nseg=" << Nseg
//            << " pf_w=" << (pf ? pf->data.size() : 0)
//            << " Sf_w=" << (sf ? sf->data.size() : 0)
//            << " Tf=" << (Tf ? Tf->data.size() : 0) << "\n";
//    }
//
//    // 分开导出：矩阵 / 裂缝
//    outM.exportMatrixValue(mgr.mesh(), reg, time0, step0);
//    outF.exportFractureValue(mgr.mesh(), mgr.fracture_network(), reg_fr, time0, step0);
//
//    std::cout << "Finished initial setup and wrote t=0 state for MATLAB.\n";
//
//	// 下面进入时间推进模块
//
//    //********************求解单相渗流-传热问题***********************//
//    //控制方程包括：质量守恒（单相达西）+能量守恒（基于单相达西速度的对流扩散方程）
//    //质量守恒（单相达西，代求变量为: P_w）：时间项+扩散项=源项
//    //能量守恒（基于单相达西速度的对流扩散方程为: T）：时间项+对流项+扩散项=源项
//
//    //1 设置重力大小及重力方向
////2 在对密度进行迎风性判别时，是否需要考虑重力势能
//    GravUpwind gu;
//    gu.g = Vector(0.0, -9.80665, 0.0);
//    //gu.g = Vector(0.0,0.0, 0.0);
//    gu.use_potential = true;
//
//
//    // 求解器参数
//    SolverControls sc;
//    sc.maxOuter = 300;
//    sc.tol_p_abs = 10.0;     // 举例：绝对容差 10 Pa
//    sc.tol_T_abs = 1e-6;
//    sc.tol_p_rel = 1e-6;     // 新增：压力相对容差
//    sc.tol_T_rel = 1e-6;     // 新增：温度相对容差
//	sc.urf_p = 0.2; //欠松弛因子
//	sc.urf_T = 0.25; //欠松弛因子
//    sc.c_phi_const = 1e-9;
//    sc.jac_p = { 500, 0.8, 1e-8 };
//    sc.jac_T = { 500, 0.8, 1e-8 };
//    sc.theta_p = 0.5;
//	sc.theta_T = 0.5;
//    
//    //sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB; // 压力：也可选 CG
//    //sc.lin_p.maxIters = 5000;
//    //sc.lin_p.tol = sc.tol_p_abs;
//    //sc.lin_p.iluFill = 10;
//    //sc.lin_p.iluDrop = 1e-4;
//    //sc.lin_p.equil = true;   // 建议保持 on
//    //sc.lin_p.reusePreconditioner = 3;
//
//
//    //sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
//    //sc.lin_T.maxIters = 5000;
//    //sc.lin_T.tol = sc.tol_T_abs;
//    //sc.lin_T.iluFill = 10;
//    //sc.lin_T.iluDrop = 1e-4;
//
//    sc.lin_p.type = LinearSolverOptions::Type::SparseLU;  // 压力用直接法
//    sc.lin_p.equil = true;   // 默认就是 true，可不写
//    // 直接法不关心 maxIters/tol/ILU 设置
//
//    sc.lin_T.type = LinearSolverOptions::Type::SparseLU;  // 先试直接法
//    sc.lin_T.equil = true;
//
//    //sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB;
//    //sc.lin_T.maxIters = 10000;
//    //sc.lin_T.tol = sc.tol_T_abs;
//    //sc.lin_T.iluFill = 20;
//    //sc.lin_T.iluDrop = 1e-5;
//    //sc.lin_T.equil = true;   // 建议保持 on
//    //sc.lin_T.reusePreconditioner = 2;
//
//
//    sc.useJacobi = false;
//    int    nSteps = 1000;
//    double dt = 1000.0 / nSteps;
//
//    // 只导出 TXT，每步一份；不导出 CSV / MM：
//    runTransient_singlePhase(
//        mgr, reg, freg, ppm,
//        bcAdapter, TbcA, gu, rock,
//        nSteps, dt, sc,
//        /*phase=*/"CO2",
//        /*writeEvery=*/10,
//        /*onWrite=*/nullptr,
//        /*exportCSV=*/false,
//        /*exportTXT=*/true,
//        /*exportMM=*/true,
//        /*exportTecplotP=*/false,
//        /*exportTecplotT=*/true
//    );
//    
//
//    return 0;
//
//}
//
////测试案例：2D-常物性-温度扩散且不考虑裂隙及重力作用
//int main333()
//{
//	//定义几何区域参数
//    double lengthX = 1, lengthY = 1, lengthZ = 0;
//    //确定网格划分策略及参数
//    int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 0;
//    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
//
//    // 1) 构造并预处理网格
//    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
//    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";
//
//	// 2）生成变量场
//    FieldRegistry reg;  //基岩变量场生成器
//	FaceFieldRegistry freg; //基岩面变量场生成器
//    InitFields ic;
//    Initializer::createPrimaryFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg);
//	Initializer::fillBaseDistributions_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg,ic);
//	ensureTransientFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, "T", "T_old", "T_prev");
//
//	//3）物性参数设置,当前给定常物性参数 
//    PhysicalPropertiesManager ppm;
//	ppm.RockProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
//    ppm.CO2Properties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
//	ppm.ComputeEffectiveThermalProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
//
//	//4）边界条件设置
//	const auto& bfaces = mgr.boundaryFaces();
//	TemperatureBC::Registry tbc;
//    TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 573.15 };
//    TemperatureBC::BoundaryCoefficient T_Right{ 0.0, 1.0, 0 };
//    TemperatureBC::BoundaryCoefficient T_Down{ 1.0, 0.0, 373.15 };
//    TemperatureBC::BoundaryCoefficient T_Up{ 1.0, 0.0, 573.15 };
//	TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
//	tbc.printInformationofBoundarySetting(mgr);
//	TemperatureBCAdapter TbcA{ tbc };
//
//    //5 求解器设置
//    GravUpwind gu;
//    //gu.g = Vector(0.0, -9.80665, 0.0);
//    gu.g = Vector(0.0,0.0, 0.0);
//    gu.use_potential = false;
//
//    SolverControls sc;
//    sc.maxOuter = 300;
//    sc.tol_T_abs = 1e-6;
//    sc.tol_T_rel = 1e-6;     // 新增：温度相对容差
//    sc.urf_T = 0.25; //欠松弛因子
//    sc.c_phi_const = 1e-9;
//    sc.jac_T = { 500, 0.8, 1e-8 };
//    sc.theta_T = 0.5;
//
//    sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
//    sc.lin_T.maxIters = 5000;
//    sc.lin_T.tol = sc.tol_T_abs;
//    sc.lin_T.iluFill = 10;
//    sc.lin_T.iluDrop = 1e-4;
//
//	sc.useJacobi = false;
//
//    int    nSteps = 1000;
//    double dt = 100000.0 / nSteps;
//
//	// 6）只导出 TXT，每步一份；不导出 CSV / MM：
//    runTransient_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg, freg, ppm, TbcA, gu, nSteps, dt, sc, 1, "./T_CO2_diff/T/T");
//	return 0;
//}
//
////测试案例：2D-变物性-温度扩散且不考虑裂隙-利用COMSOL中内嵌模型对CO2物性进行计算与更新
//
//int main5555()
//{
//    //定义几何区域参数
//    double lengthX = 1, lengthY = 1, lengthZ = 0;
//    //确定网格划分策略及参数
//    int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 0;
//    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
//
//    // 1) 构造并预处理网格
//    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
//    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";
//
//    // 2）生成变量场
//    FieldRegistry reg;  //基岩变量场生成器
//    FaceFieldRegistry freg; //基岩面变量场生成器
//    InitFields ic;
//    Initializer::createPrimaryFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg);
//    Initializer::fillBaseDistributions_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, ic);
//    ensureTransientFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, "T", "T_old", "T_prev");
//
//    //3）物性参数设置,当前给定常物性参数 
//    PhysicalPropertiesManager ppm;
//    ppm.RockProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
//	ppm.CO2Properties_test_varProperties_singlePhase_CO2_T_diffusion(mgr, reg,"T");
//	ppm.ComputeEffectiveThermalProperties_test_varProperties_singlePhase_CO2_T_diffusion(mgr, reg, "T");   
//
//    // 4）边界条件设置
//    const auto & bfaces = mgr.boundaryFaces();
//    TemperatureBC::Registry tbc;
//    TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 573.15 };
//    TemperatureBC::BoundaryCoefficient T_Right{ 0.0, 1.0, 0 };
//    TemperatureBC::BoundaryCoefficient T_Down{ 1.0, 0.0, 473.15 };
//    TemperatureBC::BoundaryCoefficient T_Up{ 1.0, 0.0, 573.15 };
//    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
//    tbc.printInformationofBoundarySetting(mgr);
//    TemperatureBCAdapter TbcA{ tbc };
//
//    //5 求解器设置
//    GravUpwind gu;
//    //gu.g = Vector(0.0, -9.80665, 0.0);
//    gu.g = Vector(0.0, 0.0, 0.0);
//    gu.use_potential = false;
//
//    SolverControls sc;
//    sc.maxOuter = 300;
//    sc.tol_T_abs = 1e-6;
//    sc.tol_T_rel = 1e-6;     // 新增：温度相对容差
//    sc.urf_T = 0.25; //欠松弛因子
//    sc.c_phi_const = 1e-9;
//    sc.jac_T = { 500, 0.8, 1e-8 };
//    sc.theta_T = 0.5;
//
//    sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
//    sc.lin_T.maxIters = 5000;
//    sc.lin_T.tol = sc.tol_T_abs;
//    sc.lin_T.iluFill = 10;
//    sc.lin_T.iluDrop = 1e-4;
//
//    sc.useJacobi = false;
//
//    int    nSteps = 1000;
//    double dt = 100000.0 / nSteps;
//    runTransient_test_varyingProperties_singlePhase_CO2_T_diffusion(mgr, reg, freg, ppm, TbcA, gu, nSteps, dt, sc, 1, "./T_CO2_diff/T2/T");
//    return 0;
//}
//
//
////int main1111()
////{
//// //   //定义几何区域参数
//// //   double lengthX = 1, lengthY = 1, lengthZ = 0;
//// //   //确定网格划分策略及参数
//// //   int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 0;
//// //   bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//// //   bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
////
//// //   // 1) 构造并预处理网格
//// //   auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
//// //   MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//// //   mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//// //   auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
//// //   auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//// //   std::cout << "MatrixMesh built in " << ms0 << " ms.\n";
////
//// //   // 2）生成变量场
//// //   FieldRegistry reg;  //基岩变量场生成器
//// //   FaceFieldRegistry freg; //基岩面变量场生成器
//// //   InitFields ic;
//// //   Initializer::createPrimaryFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg);
//// //   Initializer::fillBaseDistributions_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, ic);
//// //   ensureTransientFields_scalar (mgr.mesh(), reg, "T", "T_old", "T_prev");
////
////	////3）物性参数设置,当前给定常物性参数
//// //   PhysicalPropertiesManager ppm;
//// //   ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
////	//ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);  
//// //   ppm.ComputeEffectiveThermalProperties_test_constProperties_singlePhase_CO2_T_diffusion(mgr, reg);
////
////	////4) 边界条件设置
//// //       //边界面识别
//// //   const auto& bfaces = mgr.boundaryFaces();
////
//// //   //给定参数：2D情况，给定基岩四个边界条件系数
//// // 
//// //   TemperatureBC::Registry tbc;
//// //   TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 573.15 };
//// //   TemperatureBC::BoundaryCoefficient T_Right{ 0.0, 1.0, 0 };
//// //   TemperatureBC::BoundaryCoefficient T_Down{ 0.0, 1.0, 0 };
//// //   TemperatureBC::BoundaryCoefficient T_Up{ 1.0, 0.0, 573.15 };
//// //   TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
//// //   tbc.printInformationofBoundarySetting(mgr);
//// //   TemperatureBCAdapter TbcA{ tbc };
////
////	////5 求解器设置
////	//Vector g = { 0.0, 0.0, 0.0 };
////
//// //   SolverControls sc;
//// //   sc.maxOuter = 300;
//// //   sc.tol_T_abs = 1e-9;
//// //   sc.tol_T_rel = 1e-9;     // 新增：温度相对容差
//// //   sc.urf_T = 0.25; //欠松弛因子
//// //   sc.c_phi_const = 1e-9;
//// //   sc.jac_T = { 500, 0.8, 1e-8 };
//// //   sc.theta_T = 0.5;
////
//// //   sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
//// //   sc.lin_T.maxIters = 5000;
//// //   sc.lin_T.tol = sc.tol_T_abs;
//// //   sc.lin_T.iluFill = 10;
//// //   sc.lin_T.iluDrop = 1e-4;
////
//// //   sc.useJacobi = false;
//// //   int    nSteps = 1000;
//// //   double dt = 100000.0 / nSteps;
////	//runTransient_test_constProperties_singlePhase_CO2_p_flow(mgr, reg, freg, ppm, TbcA, g, nSteps, dt, sc, 100, "./p_CO2/p/p");
////	//return 0;
////
////}
//
///// 实现单相-达西方程的求解
//int main11111()
//{
//    //定义几何区域参数
//    double lengthX = 1, lengthY = 1, lengthZ = 0;
//    //确定网格划分策略及参数
//    int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 0;
//    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
//
//    // 1) 构造并预处理网格
//    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
//    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
//    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";
//
//    // 2）生成变量场
//    FieldRegistry reg;  //基岩变量场生成器
//    FaceFieldRegistry freg; //基岩面变量场生成器
//    InitFields ic;
//    Initializer::createPrimaryFields_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg);
//    Initializer::fillBaseDistributions_test_singlePhase_CO2_T_diffusion(mgr.mesh(), reg, ic);
//    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");
//
//    //3）物性参数设置,当前给定常物性参数
//    PhysicalPropertiesManager ppm;
//    ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
//    ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
// 
//
//    //4) 边界条件设置
//        //边界面识别
//    const auto& bfaces = mgr.boundaryFaces();
//
//    //给定参数：2D情况，给定基岩四个边界条件系数
//
//    PressureBC::Registry pbc_pw;
//    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0,8e6 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0,0.0 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Down{ 1.0, 0.0,6e6 }; // p = 2e5 Pa
//    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0,0.0 }; // p = 2e5 Pa
//    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值
//    pbc_pw.printInformationofBoundarySetting(mgr);
//    PressureBCAdapter pbcA{ pbc_pw };
//
//    //5 求解器设置
//    Vector g = { 0.0, 0.0, 0.0 };
//
//    SolverControls sc;
//    sc.maxOuter = 300;
//    sc.tol_p_abs = 10.0;     // 举例：绝对容差 10 Pa
//    sc.tol_p_rel = 1e-6;     // 新增：压力相对容差
//    sc.urf_p = 0.2; //欠松弛因子
//
//    sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB; // 温度：非对称，建议 BiCGSTAB
//    sc.lin_p.maxIters = 5000;
//    sc.lin_p.tol = sc.tol_p_abs;
//    sc.lin_p.iluFill = 10;
//    sc.lin_p.iluDrop = 1e-4;
//
//    sc.useJacobi = false;
//    int    nSteps = 100;
//    double dt = 1.0 / nSteps;
//    runTransient_test_constProperties_singlePhase_CO2_p_flow(mgr, reg, freg, ppm, pbcA, g, nSteps, dt, sc, 1, "./p_CO2/p/p");
//    return 0;
//
//}

//实现CO2单相-常物性-渗流传热控制方程的离散与求解
int SinglePhase_CO2_TH_withoutWell()
{
    
    //定义几何区域参数
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    //确定网格划分策略及参数
    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    // 1) 构造并预处理网格
    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";

    // 2）生成变量场
    FieldRegistry reg;  //基岩变量场生成器
    FaceFieldRegistry freg; //基岩面变量场生成器
    InitFields ic;

	Initializer::createPrimaryFields_singlePhase_CO2_P(mgr.mesh(), reg);
	Initializer::createPrimaryFields_singlePhase_CO2_T(mgr.mesh(), reg);
	Initializer::fillBaseDistributions_singlePhase_CO2_P(mgr.mesh(), reg, ic);
	Initializer::fillBaseDistributions_singlePhase_CO2_T(mgr.mesh(), reg, ic);
  
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");
	ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    //3）物性参数设置,当前给定常物性参数
    PhysicalPropertiesManager ppm;
    ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
    ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
	ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

    //4) 边界条件设置
    const auto& bfaces = mgr.boundaryFaces();
    /// 压力边界条件设置
    //给定参数：2D情况，给定基岩四个边界条件系数
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0,0,6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0,0,5.5e6 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0,1,0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0,1,0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter PbcA{ pbc_pw };

    // /  温度边界条件设置
    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 1.0,0,673.15 };
    TemperatureBC::BoundaryCoefficient T_Right{ 0.0,1,0.0 };
    TemperatureBC::BoundaryCoefficient T_Down{ 0.0,1,0.0 };
    TemperatureBC::BoundaryCoefficient T_Up{ 0.0,1,0.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    tbc.printInformationofBoundarySetting(mgr);
    TemperatureBCAdapter TbcA{ tbc };

    // 5) 求解器与时间推进设置
    Vector g = { 0.0, 0.0, 0.0 };

    SolverControls sc;
    sc.maxOuter = 300;

    // 压力
    sc.tol_p_abs = 1e-6;     // 绝对容差（Pa）
    sc.tol_p_rel = 1e-6;     // 相对容差
    sc.urf_p = 0.20;     // 欠松弛
    sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_p.maxIters = 5000;
    sc.lin_p.tol = sc.tol_p_abs;
    sc.lin_p.iluFill = 10;
    sc.lin_p.iluDrop = 1e-4;

    // 温度
    sc.tol_T_abs = 1e-6;
    sc.tol_T_rel = 1e-6;
    sc.urf_T = 0.20;
    sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_T.maxIters = 5000;
    sc.lin_T.tol = sc.tol_T_abs;
    sc.lin_T.iluFill = 10;
    sc.lin_T.iluDrop = 1e-4;

    sc.useJacobi = false;

    // 时间步
    int    nSteps = 100;
    double dt = 10000 / nSteps;

    // 6) 运行：p→mf→T（分别输出到两个路径）
    bool ok = runTransient_constProperties_singlePhase_CO2_T_H
    (
        mgr, reg, freg, ppm,
        TbcA, PbcA, g,
        nSteps, dt, sc,
        /*writeEveryP*/ 10, /*writeEveryT*/ 10,
        /*outPrefixP*/ "./Postprocess_Data/SinglePhase/P_CO2/p",
        /*outPrefixT*/ "./Postprocess_Data/SinglePhase/T_CO2/T"
    );

    if (!ok) {
        std::cerr << "[MAIN] Transient run failed.\n";
        return 1;
    }

    std::cout << "[MAIN] Done.\n";
    return 0;

}

//实现CO2单相-常物性-渗流传热带井源控制方程的离散与求解 不带裂缝
int SinglePhase_CO2_TH_withWell()
{
    //定义几何区域参数
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    //确定网格划分策略及参数
    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    // 1) 构造并预处理网格
    auto t0 = std::chrono::high_resolution_clock::now(); // 计时开始
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
    auto t1 = std::chrono::high_resolution_clock::now(); // 计时结束
    auto ms0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "MatrixMesh built in " << ms0 << " ms.\n";

    // 2）生成变量场
    FieldRegistry reg;  //基岩变量场生成器
    FaceFieldRegistry freg; //基岩面变量场生成器
    InitFields ic;

    Initializer::createPrimaryFields_singlePhase_CO2_P(mgr.mesh(), reg);
    Initializer::createPrimaryFields_singlePhase_CO2_T(mgr.mesh(), reg);
    Initializer::fillBaseDistributions_singlePhase_CO2_P(mgr.mesh(), reg, ic);
    Initializer::fillBaseDistributions_singlePhase_CO2_T(mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    //3）物性参数设置,当前给定常物性参数
    PhysicalPropertiesManager ppm;
    ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
    ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
    ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

    //4) 边界条件设置
    const auto& bfaces = mgr.boundaryFaces();
    /// 压力边界条件设置
    //给定参数：2D情况，给定基岩四个边界条件系数
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0,1.0,0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0,1.0,0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值，3D是按照左右 下上 前后的顺序
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter PbcA{ pbc_pw };

    // /  温度边界条件设置
    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Right{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Down{ 0.0,1.0,0.0 };
    TemperatureBC::BoundaryCoefficient T_Up{ 0.0,1.0,0.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    tbc.printInformationofBoundarySetting(mgr);
    TemperatureBCAdapter TbcA{ tbc };

    // 全局/统一 Peaceman 参数（如需每井不同，可逐井修改 cfg.pm）
    PeacemanParams prm;
    prm.mu = 1.48e-5; prm.rho = 800; prm.reFactor = 1.12; prm.fallbackKh = 1e-14; prm.PI_is_mass = true;

    // —— 构造两个注入井，一个采出井（示例） —— //
    std::vector<WellConfig> wellsCfg;
    wellsCfg.reserve(3);

    // INJ1
    WellConfig inj1;
    inj1.name = "inj_1";
    inj1.role = WellDOF::Role::Injector;
    inj1.geom.name = inj1.name;               // 便于日志
    inj1.geom.pos = Vector(0.25, 0.25, 0.0);
    inj1.geom.rw = 0.001;
    inj1.geom.skin = 0.0;
    inj1.geom.H = 1.0;
    inj1.geom.perfRadius = 0.0;
    inj1.geom.maxHitCells = 3;
    inj1.pm = prm;                           // 拷贝全局 Peaceman 参数
    inj1.mode = WellDOF::Mode::Rate;       // 固定井底压
    inj1.target = 10;                         // BHP
    inj1.Tin = 373.15;                      // 注入温度（只对注井用）
    inj1.derive_names_if_empty();              // mask_inj_1 / PI_inj_1 / p_w_inj_1
    wellsCfg.push_back(inj1);

    // INJ2
    WellConfig inj2;
    inj2.name = "inj_2";
    inj2.role = WellDOF::Role::Injector;
    inj2.geom.name = inj2.name;
    inj2.geom.pos = Vector(0.25, 0.85, 0.0);
    inj2.geom.rw = 0.001;
    inj2.geom.skin = 0.0;
    inj2.geom.H = 1.0;
    inj2.geom.perfRadius = 0.0;
    inj2.geom.maxHitCells = 3;
    inj2.pm = prm;
    inj2.mode = WellDOF::Mode::Rate; //Rate Pressure
    inj2.target = 8;
    inj2.Tin = 373.15;
    inj2.derive_names_if_empty();              // mask_inj_2 / PI_inj_2 / p_w_inj_2
    wellsCfg.push_back(inj2);

    // PROD1
    WellConfig prod1;
    prod1.name = "prod_1";
    prod1.role = WellDOF::Role::Producer;
    prod1.geom.name = prod1.name;
    prod1.geom.pos = Vector(0.75, 0.55, 0.0);
    prod1.geom.rw = 0.001;
    prod1.geom.skin = 0.0;
    prod1.geom.H = 1.0;
    prod1.geom.perfRadius = 0.0;
    prod1.geom.maxHitCells = 1;
    prod1.pm = prm;
    prod1.mode = WellDOF::Mode::Pressure;
    prod1.target = 5e6;
    prod1.Tin = 0.0;                        // 采井不用
    prod1.derive_names_if_empty();             // mask_prod_1 / PI_prod_1 / p_w_prod_1
    wellsCfg.push_back(prod1);

    // —— 生成掩码与 PI（每井一组字段） —— //
    build_masks_and_PI_for_all(mgr, reg, wellsCfg);

    // —— 注册 DoF —— //
    const int Nc = (int)mgr.mesh().getCells().size();
    std::vector<WellDOF> wells;
    const int Ntot = register_well_dofs_for_all(Nc, wellsCfg, wells);

    // 6) 求解器与时间推进设置
    Vector g = { 0.0, 0.0, 0.0 };

    SolverControls sc;
    sc.maxOuter = 200;              // 允许的最大外迭代
    sc.useJacobi = false;

    // —— 压力收敛与线性求解 —— //
    sc.tol_p_abs = 1e-6;
    sc.tol_p_rel = 1e-6;
    sc.urf_p = 0.35;

    sc.lin_p.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_p.maxIters = 4000;
    sc.lin_p.tol = sc.tol_p_abs;
    sc.lin_p.iluFill = 15;
    sc.lin_p.iluDrop = 5e-5;

    // —— 温度收敛与线性求解 —— //
    sc.tol_T_abs = 1e-8;
    sc.tol_T_rel = 1e-6;
    sc.urf_T = 1;

    sc.lin_T.type = LinearSolverOptions::Type::BiCGSTAB;
    sc.lin_T.maxIters = 4000;
    sc.lin_T.tol = sc.tol_T_abs;
    sc.lin_T.iluFill = 25;
    sc.lin_T.iluDrop = 1e-5;

    // —— Step 1/3: 压力子循环与增量通量控制 —— //
    sc.NsweepP_max = 3;       // 每次外迭代内的压力子扫上限
    sc.tol_p_inner = 2e-7;    // 子扫提前退出阈值
    sc.enable_incremental_convection = false;
    sc.kappa_p_for_flux_update = 0.5;   // 压力进度阈值
    sc.incremental_flip_ratio = 0.02;    // 面符号翻转占比阈值
    sc.flux_sign_epsilon = 1e-14;

    // —— Step 4: Aitken 加速 —— //
    sc.enable_Aitken_p = true;
    sc.enable_Aitken_T = true;
    sc.aitken_omega_min = 0.05;
    sc.aitken_omega_max = 1.10;

    // —— Step 5: 线性求解器复用 —— //
    sc.reuse_linear_pattern = true;
    sc.rebuild_precond_every = 3;    // 每 3 次子扫刷新一次 ILUT 预条件器

    // —— Step 6: CFL Guard —— //
    sc.enable_CFL_guard = true;
    sc.CFL_T_threshold = 0.8;      // 若 maxCFL_T > 0.8 即触发子步
    sc.CFL_dt_scale_min = 0.25;     // 子步时间不得小于 0.25*dt

    // —— 时间步 —— //
    const int    nSteps = 100;
    const double dt = 1000.0 / nSteps;

    bool ok = runTransient_constProperties_singlePhase_CO2_T_H_withWell(mgr, reg, freg, ppm, TbcA, PbcA, g, nSteps, dt, sc, wellsCfg, 1, 1, "./Postprocess_Data/P_CO2_withWell_rate/p", "./Postprocess_Data/T_CO2_withWell_rate/T");
    // —— 调用加速版时间推进 —— //
    //bool ok = runTransient_constProperties_singlePhase_CO2_T_H_withWell_accel(
    //    mgr, reg, freg, ppm,
    //    TbcA, PbcA, g,
    //    nSteps, dt,
    //    sc, wellsCfg,
    //    /*writeEveryP*/ 1,
    //    /*writeEveryT*/ 1,
    //    "./Postprocess_Data/P_CO2_accel/p",
    //    "./Postprocess_Data/T_CO2_accel/T",
    //    /*runtimeIn*/ nullptr   // 若想跨步复用，可传入持久化的 OuterIterRuntime
    //);

    if (!ok) {
        std::cerr << "[MAIN] Accelerated transient run failed.\n";
        return 1;
    }

    std::cout << "[MAIN] Accelerated run finished.\n";
    return 0;
}





/// <summary>
/// /2D_TwoPhase_TH_constProperties_FVM test
/// </summary>
/// <returns></returns>

int legacy_singlePhase_main()
{
    //定义几何区域参数
    double lengthX = 100, lengthY = 100, lengthZ = 0;
    //确定网格划分策略及参数
    int sectionNumX = 5, sectionNumY = 5, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型

    // 1) 构造并预处理网格
    cout << "--- Step 1: Building and Pre-processing Mesh ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性

    // 2) 网格单元场和网格面场注册
    cout << "--- Step 2: Registering Field Registries ---\n";
    FieldRegistry reg;
    FaceFieldRegistry freg;

    //3) 给定地层初始参数
    cout << "--- Step 3: Setting Initial Conditions ---\n";
    InitFields ic;
    {
        //==压力==//
        ic.p_w0 = 1e7;
        ic.dp_wdx = 0;
        ic.dp_wdy = 0;
        ic.dp_wdz = 0;

        //==温度==//
        ic.T0 = 593.15;
        ic.dTdx = 0;
        ic.dTdy = 0;
        ic.dTdz = 0;

        //==初始水相饱和度==//
        ic.s_w = 0.8;   //水相初始饱和度
    }

    //4.1) 基岩主变量场注册名称
    {
        cout << "--- Step 4.1: Creating and Initializing Primary Variable Fields ---\n";
        Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg,"p_w","s_w","T");
        Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);
    }

    //4.2) 时层拷贝
    cout << "--- Step 4.2: Ensuring Transient Fields (old, prev) ---\n";
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
  

    // 5) 物性参数设置
    PhysicalPropertiesManager ppm;
    {
        // 5.1) 基岩区域划分及基岩物性参数初始化，定物性参数工况
        cout << "  - Updating Matrix Rock Properties...\n";
        {
            ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
            ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");
        }

        // 5.2)裂缝区域划分及裂缝物性参数初始化,函数已给出，这里先不调用，因为还没有裂缝 的主变量场
        {

        }

        //5.3) 基岩内流体物性参数初始化
        cout << "  - Updating Fluid Properties in Rock...\n";
        {
            // 在计算流体物性前，需要有 p_g 场。p_g 是由 updateRockTwoPhaseProperties_IMPES 计算的
            // 因此，我们先进行一次初步的两相属性计算，以获得初始的 p_g 场
            VGParams temp_vg_params; // 临时参数，仅为获得p_g
            temp_vg_params.Swr = 0.25; temp_vg_params.Sgr = 0.05;
            temp_vg_params.alpha = 2.0e-6; temp_vg_params.n = 2.5;
            RelPermParams temp_rp_params; temp_rp_params.L = 0.5;
            multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, temp_vg_params, temp_rp_params);
           
            //CO2
            ppm.UpdateCO2inRockAt(mgr, reg, "p_w", "T");

            //water
            ppm.UpdateWaterinRockAt(mgr, reg, "p_w", "T");
        }


        //5.4) 裂缝内流体物性参数初始化
        {
            //CO2 函数已内置，由于还没有裂缝 的主变量场这里不计算

            //water 函数已内置，由于还没有裂缝 的主变量场这里不计算
        }
    }
        

    // 6) 两相参数计算及有效热物性参数计算
    cout << "\n--- Step 6: Two-Phase & Effective Property Calculation ---\n";
    {
        // 6.1) 定义两相流模型参数 (VG & Mualem)
        // 这些参数通常与岩石类型相关，后续可以扩展为按区域赋值
        VGParams vg_rock_params;
        vg_rock_params.Swr = 0.25;      // 水相残余饱和度
        vg_rock_params.Sgr = 0.05;      // 气相残余饱和度
        vg_rock_params.alpha = 2.0e-6;  // Van Genuchten alpha (进气压力相关), [1/Pa]
        vg_rock_params.n = 2.5;         // Van Genuchten n (孔隙尺寸分布指数)
        
        RelPermParams rp_rock_params;
        rp_rock_params.L = 0.5;         // Mualem 模型孔隙连通性参数

        // 6.2) 计算基岩的派生两相属性
        // 这个函数会根据主变量 s_w 和 p_w，以及基础物性，计算出基岩内的 Pc, p_g, k_r, lambda 等场

        bool ok = multiPhase::updateRockTwoPhaseProperties_IMPES( mgr,reg,vg_rock_params,rp_rock_params);

        if (ok)
        {
            auto pc_field = reg.get<volScalarField>("Pc");
            auto pg_field = reg.get<volScalarField>("p_g");
            auto krw_field = reg.get<volScalarField>("k_rw");

            if (pc_field && pg_field && krw_field && mgr.mesh().getCells().size() > 5)
            {
                cout << "==========================================================================" << endl;
                cout << "  - Check cell 5: Pc = " << (*pc_field)[5] << " Pa, p_g = " << (*pg_field)[5] << " Pa, k_rw = " << (*krw_field)[5] << std::endl;
                cout << "  - Check cell 10: Pc = " << (*pc_field)[10] << " Pa, p_g = " << (*pg_field)[10] << " Pa, k_rw = " << (*krw_field)[10] << std::endl;
                cout << "  - Check cell 15: Pc = " << (*pc_field)[15] << " Pa, p_g = " << (*pg_field)[15] << " Pa, k_rw = " << (*krw_field)[15] << std::endl;
                cout << "==========================================================================" << endl;
            }
        }
        else
        {
            cerr << "ERROR: Failed to update two-phase properties for the rock matrix. Aborting.\n";
            return -1; // 返回错误码
        }
        // 6.3) 裂缝的两相属性计算 (当前未激活)
        {
            // VGParams vg_frac_params = ... ;
            // RelPermParams rp_frac_params = ... ;
            // multiPhase::updateFractureTwoPhaseProperties_IMPES(mgr, reg_fr, vg_frac_params, rp_frac_params);
        }
    }

    // 7) 边界条件设置
    cout << "\n--- Step 7: Setting Boundary Conditions and Source Terms ---\n";
    {
        // 7.1) 获取边界面的分组信息
        const auto& bfaces = mgr.boundaryFaces();

        // 7.2) 压力边界条件 (p_w)
        PressureBC::Registry pbc_pw;
        PressureBC::BoundaryCoefficient P_Left{ 0.0,1.0,0.0 };
        PressureBC::BoundaryCoefficient P_Right{ 0.0,1.0,0.0 };
        PressureBC::BoundaryCoefficient P_Down{ 0.0,1.0,0.0 };
        PressureBC::BoundaryCoefficient P_Up{ 0.0,1.0,0.0 };
        PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);   //按照左 右 下 上 的顺序赋值，3D是按照左右 下上 前后的顺序
        pbc_pw.printInformationofBoundarySetting(mgr);
        
        pbc_pw.printInformationofBoundarySetting(mgr);

        // 创建压力边界条件适配器，供求解器使用
        PressureBCAdapter PbcA{ pbc_pw };

        // /  温度边界条件设置
        TemperatureBC::Registry tbc;
        TemperatureBC::BoundaryCoefficient T_Left{ 0.0,1.0,0.0 };
        TemperatureBC::BoundaryCoefficient T_Right{ 0.0,1.0,0.0 };
        TemperatureBC::BoundaryCoefficient T_Down{ 0.0,1.0,0.0 };
        TemperatureBC::BoundaryCoefficient T_Up{ 0.0,1.0,0.0 };
        TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
        tbc.printInformationofBoundarySetting(mgr);
        
        tbc.printInformationofBoundarySetting(mgr);
        
        // 创建温度边界条件适配器，供求解器使用
        TemperatureBCAdapter TbcA{ tbc };

    }

    // 8）注入井和采出井设置
  
    cout << "\n--- Step 8: Setting up Two-Phase Well Configurations ---\n";

    // Create a vector to hold all two-phase well configurations
    std::vector<WellConfig_TwoPhase> wells_cfg_2p;
    wells_cfg_2p.reserve(2);

    // --- Configure Injection Well INJ1 ---
    {
        WellConfig_TwoPhase inj1;
        inj1.name = "INJ1_2P";
        inj1.role = WellDOF_TwoPhase::Role::Injector;

        // Geometry
        inj1.geom.pos = Vector(25.0, 25.0, 0.0);
        inj1.geom.rw = 0.1;
        inj1.geom.H = (lengthZ > 0) ? lengthZ : 1.0;
        inj1.geom.maxHitCells = 1;

        // Peaceman parameters (can be default)
        // inj1.pm_2p.reFactor = 0.28;

         // Operational control
        inj1.mode = WellDOF_TwoPhase::Mode::Pressure;
        inj1.target = 1.2e7; // 12 MPa

        // Two-phase injection parameters
        inj1.Tin = 373.15;   // 100 C
        inj1.s_w_bh = 1.0;     // Initial injection is pure water
        inj1.mu_w_inj = 2.82e-4;  // Viscosity of water
        inj1.rho_w_inj = 958.4;    // Density of water
        inj1.cp_w_inj = 4216.0;   // Heat capacity of water
        // Set CO2 properties to zero or default, as they are not used for pure water injection
        inj1.mu_g_inj = 0.0;
        inj1.rho_g_inj = 0.0;
        inj1.cp_g_inj = 0.0;

        wells_cfg_2p.push_back(inj1);
    }

    // --- Configure Production Well PROD1 ---
    {
        WellConfig_TwoPhase prod1;
        prod1.name = "PROD1_2P";
        prod1.role = WellDOF_TwoPhase::Role::Producer;

        // Geometry
        prod1.geom.pos = Vector(75.0, 75.0, 0.0);
        prod1.geom.rw = 0.1;
        prod1.geom.H = (lengthZ > 0) ? lengthZ : 1.0;

        // Operational control
        prod1.mode = WellDOF_TwoPhase::Mode::Pressure;
        prod1.target = 0.9e7; // 9 MPa

        wells_cfg_2p.push_back(prod1);
    }
    // --- Pre-process all wells with a single function call ---
    build_masks_and_WI_for_all(mgr, reg, wells_cfg_2p);

    // --- Register all well DOFs with a single function call ---
    const int Nc = (int)mgr.mesh().getCells().size();
    cout << "Mesh Size: " << Nc << endl;

    std::vector<WellDOF_TwoPhase> wells_dof_2p; // This vector will be used by the solver
    const int N_total = register_well_dofs_for_all_TwoPhase(Nc, wells_cfg_2p, wells_dof_2p);

    cout << "Total unknowns registered: " << N_total << endl;

    cout << "\n============================================\n";
    cout << "=== Step 0 (Initialization) Fully Complete ===\n";
    cout << "============================================\n";

    // --- Time Loop Example ---
    double current_time = 0.0;
    double end_time = 100.0;
    double dt = 10.0;

    while (current_time < end_time) {
        // ...

        // Find the injection well config and update its s_w_bh for phase switching
        for (auto& cfg : wells_cfg_2p) {
            if (cfg.name == "INJ1_2P") {
                if (current_time + dt <= 60.0) {
                    cfg.s_w_bh = 1.0; // Water
                }
                else {
                    cfg.s_w_bh = 0.0; // CO2
                }
                cout << "  - Well " << cfg.name << " injecting with s_w_bh = " << cfg.s_w_bh << endl;

                // CRITICAL: Update the corresponding DOF object for the solver
                for (auto& dof : wells_dof_2p) {
                    if (dof.name == cfg.pw_name) { // Match using the unique pw_name
                        dof.s_w_bh = cfg.s_w_bh;
                        break;
                    }
                }
                break;
            }
        }

        // Now, call your solvers. The solver will use the `wells_dof_2p` vector
        // which contains the up-to-date s_w_bh for the current time step.
        // ...

        current_time += dt;
    }







    return 0;

}

int runFullIMPESCase()
{
    using std::cout;
    using std::endl;

    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 30;
    const int sectionNumY = 30;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    cout << "--- IMPES: building mesh ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry reg;
    FaceFieldRegistry freg;

    cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;
    ic.p_w0 = 6.5e6;
    ic.T0 = 560.0;
    ic.s_w = 1;
    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");



    cout << "--- IMPES: updating rock/fluid properties ---\n";
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    VGParams vg_params;
    vg_params.Swr = 0.00;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;
    RelPermParams rp_params;
    rp_params.L = 0.5;

    ppm.UpdateCO2inRockAt(mgr, reg, "p_w", "T");
    ppm.UpdateWaterinRockAt(mgr, reg, "p_w", "T");
    if (!multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params))
    {
        std::cerr << "[MAIN] failed to build two-phase rock properties (lambda_t, etc.).\n";
        return 1;
    }
   // ppm.ComputeMatrixEffectiveThermalsAt(mgr, reg, "p_w", "T", "both");

    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Down{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Up{ 1.0, 0.0, 6.5e6 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    TemperatureBC::Registry tbc;
    TemperatureBC::BoundaryCoefficient T_Left{ 1.0, 0.0, 680 };
    TemperatureBC::BoundaryCoefficient T_Right{ 0.0, 1.0, 0.0 };
    TemperatureBC::BoundaryCoefficient T_Down{ 0.0, 1.0, 0.0 };
    TemperatureBC::BoundaryCoefficient T_Up{ 0.0, 1.0, 0.0 };
    TemperatureBC::setBoxBCs2D(tbc, bfaces, T_Left, T_Right, T_Down, T_Up);
    TemperatureBCAdapter TbcA{ tbc };

    cout << "--- IMPES: configuring wells ---\n";
    std::vector<WellConfig_TwoPhase> wells_cfg;
    /*{
        WellConfig_TwoPhase inj;
        inj.name = "INJ1_2P";
        inj.role = WellDOF_TwoPhase::Role::Injector;
        inj.geom.pos = Vector(25.0, 50.0, 0.0);
        inj.geom.rw = 0.1;
        inj.geom.H = (lengthZ > 0.0) ? lengthZ : 1.0;
        inj.geom.maxHitCells = 1;
        inj.mode = WellDOF_TwoPhase::Mode::Pressure;
        inj.target = 1.2e7;
        inj.Tin = 380.0;
        inj.s_w_bh = 1.0;
        inj.mu_w_inj = 2.8e-4;
        inj.mu_g_inj = 1.5e-5;
        inj.rho_w_inj = 960.0;
        inj.rho_g_inj = 720.0;
        inj.cp_w_inj = 4200.0;
        inj.cp_g_inj = 1200.0;
        wells_cfg.push_back(inj);
    }
    {
        WellConfig_TwoPhase prod;
        prod.name = "PROD1_2P";
        prod.role = WellDOF_TwoPhase::Role::Producer;
        prod.geom.pos = Vector(75.0, 50.0, 0.0);
        prod.geom.rw = 0.1;
        prod.geom.H = (lengthZ > 0.0) ? lengthZ : 1.0;
        prod.geom.maxHitCells = 1;
        prod.mode = WellDOF_TwoPhase::Mode::Pressure;
        prod.target = 9.0e6;
        wells_cfg.push_back(prod);
    }*/

    build_masks_and_WI_for_all(mgr, reg, wells_cfg);
    const int Nc = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Nc, wells_cfg, wells_dof);

    Solver::IMPES::PressureSolveControls pressureCtrl;
    pressureCtrl.assembly.pressure_field = "p_w";
    pressureCtrl.assembly.pressure_old_field = "p_w_old";
    pressureCtrl.assembly.pressure_prev_field = "p_w_prev";
    pressureCtrl.assembly.phi_field = "phi_r";
    pressureCtrl.assembly.lambda_total_field = "lambda_t";
    pressureCtrl.assembly.kxx_field = "kxx";
    pressureCtrl.assembly.kyy_field = "kyy";
    pressureCtrl.assembly.kzz_field = "kzz";
    pressureCtrl.linear.type = LinearSolverOptions::Type::BiCGSTAB;
    pressureCtrl.linear.maxIters = 2000;
    pressureCtrl.linear.tol = 1e-8;
    pressureCtrl.linear.iluFill = 10;
    pressureCtrl.linear.iluDrop = 1e-4;
    pressureCtrl.under_relax = 0.5;
    pressureCtrl.max_outer_iterations = 80;
    pressureCtrl.assembly.use_constant_phase_density = true;
	pressureCtrl.report_outer_iterations = true;

    Solver::IMPES::FluxSplitConfig fluxCfg;
    fluxCfg.lambda_water = "lambda_w";
    fluxCfg.lambda_gas = "lambda_g";
    fluxCfg.rho_water = "rho_w";
    fluxCfg.rho_gas = "rho_g";

    Solver::IMPES::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.water_flux = fluxCfg.water_mass_flux;
    satCfg.porosity = "phi_r";
    satCfg.rho_water = "rho_w";
    satCfg.pressure = pressureCtrl.assembly.pressure_field;
    satCfg.thickness = (lengthZ > 0.0) ? lengthZ : 1.0;
    satCfg.vg_params = vg_params;
    satCfg.rp_params = rp_params;

    Solver::IMPES::TemperatureSolveControls tempCtrl;
    tempCtrl.assembly.temperature_field = "T";
    tempCtrl.assembly.temperature_old_field = "T_old";
    tempCtrl.assembly.temperature_prev_field = "T_prev";
    tempCtrl.assembly.C_eff_field = "C_eff";
    tempCtrl.assembly.lambda_eff_field = "lambda_eff";
    tempCtrl.assembly.mass_flux_water = fluxCfg.water_mass_flux;
    tempCtrl.assembly.mass_flux_gas = fluxCfg.gas_mass_flux;
    tempCtrl.assembly.cp_water_field = "cp_w";
    tempCtrl.assembly.cp_gas_field = "cp_g";
    tempCtrl.assembly.lambda_water_field = "lambda_w";
    tempCtrl.assembly.lambda_gas_field = "lambda_g";
    tempCtrl.assembly.rho_water_field = "rho_w";
    tempCtrl.assembly.rho_gas_field = "rho_g";
    tempCtrl.assembly.pressure_field = pressureCtrl.assembly.pressure_field;
    tempCtrl.linear.type = LinearSolverOptions::Type::BiCGSTAB;
    tempCtrl.linear.maxIters = 2000;
    tempCtrl.linear.tol = 1e-8;
    tempCtrl.linear.iluFill = 10;
    tempCtrl.linear.iluDrop = 1e-4;
    tempCtrl.under_relax = 0.2;
    tempCtrl.max_outer_iterations = 100;
    tempCtrl.urf_step = 0.2;
	tempCtrl.report_outer_iterations = true;

    /*const int nSteps = 100000000;
    const double totalDays = 6.0;
    const double dt = (totalDays * 24.0 * 3600.0) / static_cast<double>(nSteps);
    const int dumpFrequency = std::max(1, nSteps / 100);*/

    const int    nSteps = 10000000;
    const double dt = 1000.0 / nSteps;
    //const double dt = (totalDays * 24.0 * 3600.0) / static_cast<double>(nSteps);
    //const int dumpFrequency = std::max(1, nSteps / 100);

    cout << "--- IMPES: running transient simulation ---\n";
    //const bool ok = Solver::IMPES::runTransient_IMPES(
    //    mgr, reg, freg,
    //    PbcA, TbcA,
    //    wells_dof,
    //    nSteps,
    //    dt,
    //    pressureCtrl,
    //    fluxCfg,
    //    satCfg,
    //    tempCtrl,
    //    dumpFrequency,
    //    dumpFrequency,
    //    dumpFrequency,
    //    "./Postprocess_Data/p_impes/p",
    //    "./Postprocess_Data/T_impes/t",
    //    "./Postprocess_Data/sw_impes/sw");

    const bool ok = Solver::IMPES::runTransient_IMPES(
        mgr, reg, freg,
        PbcA, TbcA,
        wells_dof,
        nSteps,
        dt,
        pressureCtrl,
        fluxCfg,
        satCfg,
        tempCtrl,
        10,
        10,
        10,
        "./Postprocess_Data/p_impes/p2",
        "./Postprocess_Data/T_impes/t2",
        "./Postprocess_Data/sw_impes/sw2",
        "./Postprocess_Data/time_series/impes_stats.csv",
        10,
        "./Postprocess_Data/csv_snapshots/state",
        std::vector<std::string>{
            pressureCtrl.assembly.pressure_field,
            tempCtrl.assembly.temperature_field,
            satCfg.saturation });

    if (!ok)
    {
        std::cerr << "[MAIN] IMPES transient run failed.\n";
        return 1;
    }

    cout << "[MAIN] IMPES transient run completed successfully.\n";
    return 0;
}

//可运行，且可看到驱替效果
int runCO2Displacement_IMPES_Ps()
{
    using std::cout;
    using std::endl;

    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 30;
    const int sectionNumY = 30;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    cout << "--- IMPES: building mesh (CO2 displacement) ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry reg;
    FaceFieldRegistry freg;

    cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;
    ic.p_w0 = 6.5e6;
    ic.T0 = 373.15;
    ic.s_w = 1.0;
    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);

    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");

    {
        const double co2StripeWidth = std::max(5.0, 0.1 * lengthX);
        auto sw = reg.get<volScalarField>("s_w");
        auto sw_old = reg.get<volScalarField>("s_w_old");
        auto sw_prev = reg.get<volScalarField>("s_w_prev");
        if (sw && sw_old && sw_prev)
        {
            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                if (c.center.m_x <= co2StripeWidth)
                {
                    const size_t idx = id2idx.at(c.id);
                    (*sw)[idx] = 0.0;
                    (*sw_old)[idx] = 0.0;
                    (*sw_prev)[idx] = 0.0;
                }
            }
        }
    }

    cout << "--- IMPES: updating rock/fluid properties ---\n";
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    VGParams vg_params;
    vg_params.Swr = 0.20;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;
    RelPermParams rp_params;
    rp_params.L = 0.5;

    ppm.UpdateCO2inRockAt(mgr, reg, "p_w", "T");
    ppm.UpdateWaterinRockAt(mgr, reg, "p_w", "T");

    if (!multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params))
    {
        std::cerr << "[MAIN][P-s case] failed to build two-phase rock properties.\n";
        return 1;
    }

    const size_t Nc = mgr.mesh().getCells().size();
    auto fillConstField = [&](const std::string& name, double value)
    {
        auto field = reg.getOrCreate<volScalarField>(name, Nc, value);
        std::fill(field->data.begin(), field->data.end(), value);
    };

    fillConstField("rho_w", 1000.0);
    fillConstField("mu_w", 5.0e-4);
    fillConstField("cp_w", 4200.0);
    fillConstField("k_w", 0.6);
    fillConstField("Drho_Dp_w", 0.0);
    fillConstField("rho_g", 650.0);
    fillConstField("mu_g", 1.5e-5);
    fillConstField("cp_g", 1000.0);
    fillConstField("k_g", 0.08);
    fillConstField("Drho_Dp_g", 0.0);

    if (!multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params))
    {
        std::cerr << "[MAIN][P-s case] failed to rebuild two-phase properties after applying constants.\n";
        return 1;
    }

    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 8.0e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    std::vector<WellConfig_TwoPhase> wells_cfg;
    /*{
        WellConfig_TwoPhase inj;
        inj.name = "CO2_INJ";
        inj.role = WellDOF_TwoPhase::Role::Injector;
        inj.geom.pos = Vector(10.0, 50.0, 0.0);
        inj.geom.rw = 0.1;
        inj.geom.H = (lengthZ > 0.0) ? lengthZ : 1.0;
        inj.geom.maxHitCells = 1;
        inj.mode = WellDOF_TwoPhase::Mode::Pressure;
        inj.target = 8.0e6;
        inj.Tin = 373.15;
        inj.s_w_bh = 0.0;
        inj.mu_w_inj = 5.0e-4;
        inj.mu_g_inj = 1.5e-5;
        inj.rho_w_inj = 1000.0;
        inj.rho_g_inj = 650.0;
        inj.cp_w_inj = 4200.0;
        inj.cp_g_inj = 1000.0;
        wells_cfg.push_back(inj);
    }
    {
        WellConfig_TwoPhase prod;
        prod.name = "PROD";
        prod.role = WellDOF_TwoPhase::Role::Producer;
        prod.geom.pos = Vector(90.0, 50.0, 0.0);
        prod.geom.rw = 0.1;
        prod.geom.H = (lengthZ > 0.0) ? lengthZ : 1.0;
        prod.geom.maxHitCells = 1;
        prod.mode = WellDOF_TwoPhase::Mode::Pressure;
        prod.target = 6.3e6;
        wells_cfg.push_back(prod);
    }*/

    build_masks_and_WI_for_all(mgr, reg, wells_cfg);
    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    Solver::IMPES::PressureSolveControls pressureCtrl;
    pressureCtrl.assembly.pressure_field = "p_w";
    pressureCtrl.assembly.pressure_old_field = "p_w_old";
    pressureCtrl.assembly.pressure_prev_field = "p_w_prev";
    pressureCtrl.assembly.phi_field = "phi_r";
    pressureCtrl.assembly.lambda_total_field = "lambda_t";
    pressureCtrl.assembly.kxx_field = "kxx";
    pressureCtrl.assembly.kyy_field = "kyy";
    pressureCtrl.assembly.kzz_field = "kzz";
    pressureCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };
    pressureCtrl.assembly.enable_buoyancy = false;
    pressureCtrl.assembly.use_constant_phase_density = true;
    pressureCtrl.linear.type = LinearSolverOptions::Type::BiCGSTAB;
    pressureCtrl.linear.maxIters = 1500;
    pressureCtrl.linear.tol = 1e-10;
    pressureCtrl.linear.iluFill = 5;
    pressureCtrl.linear.iluDrop = 1e-4;
    pressureCtrl.under_relax = 0.5;
    pressureCtrl.max_outer_iterations = 300;
    pressureCtrl.report_outer_iterations = true;
    pressureCtrl.tol_abs = 1e-7;        // 或更严格
    pressureCtrl.tol_rel = 1e-10;        // 相当于把允许的相对残差降到 0.008 Pa

    Solver::IMPES::FluxSplitConfig fluxCfg;
    fluxCfg.lambda_water = "lambda_w";
    fluxCfg.lambda_gas = "lambda_g";
    fluxCfg.rho_water = "rho_w";
    fluxCfg.rho_gas = "rho_g";

    Solver::IMPES::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.water_flux = fluxCfg.water_mass_flux;
    satCfg.porosity = "phi_r";
    satCfg.rho_water = "rho_w";
    satCfg.pressure = pressureCtrl.assembly.pressure_field;
    satCfg.thickness = (lengthZ > 0.0) ? lengthZ : 1.0;
    satCfg.vg_params = vg_params;
    satCfg.rp_params = rp_params;
    satCfg.include_well_sources = true;

    const int nSteps = 500;
    const double dt = 3600.0;

    cout << "--- IMPES: running P-s transient simulation ---\n";
    const bool ok = Solver::IMPES::runTransient_IMPES_P_s(
        mgr, reg, freg,
        PbcA,
        wells_dof,
        nSteps,
        dt,
        pressureCtrl,
        fluxCfg,
        satCfg,
        [&]() -> bool {
            return multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params);
        },
        25,
        25,
        "./Postprocess_Data/p_impes_ps/p_ps",
        "./Postprocess_Data/sw_impes_ps/sw_ps",
        "./Postprocess_Data/time_series/impes_ps_stats.csv",
        25,
        "./Postprocess_Data/csv_snapshots/ps_state",
        std::vector<std::string>{
            pressureCtrl.assembly.pressure_field,
            satCfg.saturation });

    if (!ok)
    {
        std::cerr << "[MAIN][P-s case] IMPES transient run failed.\n";
        return 1;
    }

    cout << "[MAIN][P-s case] IMPES transient run completed successfully.\n";
    return 0;
}

int run_IMPES_revised_PS()
{
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 30;
    const int sectionNumY = 30;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "--- IMPES: building mesh (CO2 displacement) ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry    reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ---------- //
    std::cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;
    ic.p_w0 = 6.5e6;
    ic.T0 = 373.15;
    ic.s_w = 1.0;

    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    // 时间相关辅助场 p_old/prev, s_old/prev, T_old/prev
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");

    // 左侧给一条 CO2 条带：将 s_w=0
    {
        const double co2StripeWidth = std::max(1.0, 0.1 * lengthX);
        auto sw = reg.get<volScalarField>("s_w");
        auto sw_old = reg.get<volScalarField>("s_w_old");
        auto sw_prev = reg.get<volScalarField>("s_w_prev");
        if (sw && sw_old && sw_prev)
        {
            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                if (c.center.m_x <= co2StripeWidth)
                {
                    const size_t idx = id2idx.at(c.id);
                    (*sw)[idx] = 0.0;
                    (*sw_old)[idx] = 0.0;
                    (*sw_prev)[idx] = 0.0;
                }
            }
        }
    }

    // ---------- 2. 基岩 / 流体物性初始化 ---------- //
    std::cout << "--- IMPES: updating rock/fluid properties ---\n";
    PhysicalPropertiesManager ppm;

    // 区域分类（这里只用一个 RegionType::Medium）
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    VGParams vg_params;
    vg_params.Swr = 0.20;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;

    RelPermParams rp_params;
    rp_params.L = 0.5;

    // 初始 CO2 / 水物性 + 两相派生量（rho_t, c_t, lambda_w/g 等）
    ppm.UpdateCO2inRockAt(mgr, reg, "p_w", "T");
    ppm.UpdateWaterinRockAt(mgr, reg, "p_w", "T");

    if (!multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params))
    {
        std::cerr << "[MAIN][P-s case] failed to build two-phase rock properties.\n";
        return 1;
    }

    // ---------- 3. 压力边界条件 ---------- //
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 8.0e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 4. 井配置（此例可以先为空，将来可 push_back 注采井） ---------- //
    std::vector<WellConfig_TwoPhase> wells_cfg;
    // TODO: 这里填入严格定流量井： wells_cfg.push_back(...);

    build_masks_and_WI_for_all(mgr, reg, wells_cfg);
    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 5. 压力求解控制参数 ---------- //
    IMPES_revised::PressureSolveControls pCtrl;

    // PressureAssemblyConfig 中字段名与当前字段保持一致即可，默认值已对齐
    // 这里只做少量明确设置（重力关掉 -> 纯水平驱动；若要浮力就改成 true）
    pCtrl.assembly.vg_params = vg_params;
    pCtrl.assembly.relperm_params = rp_params;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };
    pCtrl.assembly.enable_buoyancy = false;

    // 外迭代控制
    pCtrl.under_relax = 0.7;
    pCtrl.urf_min = 0.3;
    pCtrl.urf_max = 0.9;
    pCtrl.urf_step = 0.05;
    pCtrl.max_outer_iterations =300;
    pCtrl.tol_abs = 1e-5;
    pCtrl.tol_rel = 1e-5;
    pCtrl.report_outer_iterations = true;

    // 线性求解器参数用缺省初始化即可（LinearSolverOptions 有默认值）

    // ---------- 6. 饱和度输运配置（revised 包装版） ---------- //
    IMPES_revised::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.porosity = "phi_r";
    satCfg.rho_water = "rho_w";
    satCfg.rho_gas = "rho_g";
    satCfg.pressure = "p_w";
    satCfg.thickness = 1.0;      // 2D 模型厚度（和你岩石厚度一致）
    satCfg.Sw_min = 0.0;
    satCfg.Sw_max = 1.0;
    satCfg.CFL_limit = 0.5;
    satCfg.track_cfl = true;
    satCfg.vg_params = vg_params;
    satCfg.rp_params = rp_params;
    satCfg.include_well_sources = true;

    // 相通量拆分相关字段：与压力装配中的名字对齐
    satCfg.recompute_phase_flux = true;
    satCfg.total_mass_flux = pCtrl.assembly.total_mass_flux_name;      // "mf_total"
    satCfg.water_mass_flux = "mf_w";
    satCfg.gas_mass_flux = "mf_g";
    satCfg.fractional_flow_face = "fw_face";
    satCfg.lambda_water = pCtrl.assembly.lambda_w_field;           // "lambda_w"
    satCfg.lambda_gas = pCtrl.assembly.lambda_g_field;           // "lambda_g"
    satCfg.capillary_correction_flux = pCtrl.assembly.capillary_correction_flux_name; // "mf_capillary_corr"
    satCfg.gravity_correction_flux = pCtrl.assembly.gravity_correction_flux_name;   // "mf_gravity_corr"
    satCfg.min_lambda = 1e-30;
    satCfg.flux_sign_epsilon = 1e-15;

    // ---------- 7. 通量拆分配置（老 IMPES::FluxSplitConfig，runTransient 里现在基本不再用） ---------- //
    IMPES_revised::FluxSplitConfig fluxCfg;
    // 与上面 satCfg 中保持一致
    fluxCfg.total_mass_flux = satCfg.total_mass_flux;
    fluxCfg.water_mass_flux = satCfg.water_mass_flux;
    fluxCfg.gas_mass_flux = satCfg.gas_mass_flux;
    fluxCfg.fractional_flow_face = satCfg.fractional_flow_face;
    fluxCfg.lambda_water = satCfg.lambda_water;
    fluxCfg.lambda_gas = satCfg.lambda_gas;
    fluxCfg.rho_water = satCfg.rho_water;
    fluxCfg.rho_gas = satCfg.rho_gas;
    fluxCfg.saturation = satCfg.saturation;
    fluxCfg.capillary_correction_flux = satCfg.capillary_correction_flux;
    fluxCfg.gravity_correction_flux = satCfg.gravity_correction_flux;
    fluxCfg.min_lambda = satCfg.min_lambda;
    fluxCfg.flux_sign_epsilon = satCfg.flux_sign_epsilon;
    fluxCfg.pressure_bc = &PbcA;

    // ---------- 8. 物性冻结更新回调：时间步开始调用一次 ---------- //
    auto updateMobility = [&]() -> bool
        {
            // 1) CO2 / 水物性随着当前 (p_w, T) 更新
            ppm.UpdateCO2inRockAt(mgr, reg, "p_w", "T");
            ppm.UpdateWaterinRockAt(mgr, reg, "p_w", "T");

            // 2) 两相派生量：ct, rho_t, lambda_w/g, Pc(Sw) 等
            if (!multiPhase::updateRockTwoPhaseProperties_IMPES(mgr, reg, vg_params, rp_params))
            {
                std::cerr << "[MAIN][P-s case] updateRockTwoPhaseProperties_IMPES failed in updateMobility.\n";
                return false;
            }
            return true;
        };

    // ---------- 9. 调用总的 IMPES 时间推进 ---------- //
    const int    nSteps = 2000;   // 例子：200 步
    const double dt0 = 1E-15;  // 初始时间步 [s]，后续你可根据 CFL 报告调整

    const int writeEveryP = 5;
    const int writeEverySw = 5;

    const std::string outPrefixP = "./Postprocess_Data/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/s_impes_ps_revised/s_ps";
    const std::string timeSeriesFn = "./Postprocess_Data/time_series/impes_ps_revised_stats.csv";  // 如需写出时间序列可填路径

    const int snapshotEveryCsv = 5;   // 若需 CSV 快照可改成 >0
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/ps_state_reviesed";
    std::vector<std::string> snapshotFields; // 留空 -> runTransient_IMPES 内部默认只导出 p & s_w

    std::cout << "--- IMPES: starting transient run (revised P-S) ---\n";

    bool ok = IMPES_revised::runTransient_IMPES(
        mgr,
        reg,
        freg,
        PbcA,
        wells_dof,
        nSteps,
        dt0,
        pCtrl,
        fluxCfg,
        satCfg,
        updateMobility,
        writeEveryP,
        writeEverySw,
        outPrefixP,
        outPrefixSw,
        timeSeriesFn,
        snapshotEveryCsv,
        snapshotPrefix,
        snapshotFields);

    if (!ok)
    {
        std::cerr << "[MAIN][P-s case] runTransient_IMPES failed.\n";
        return 1;
    }

    std::cout << "--- IMPES: simulation finished successfully ---\n";
    return 0;
}

/// \brief 简单测试：两相物性 + 基本物性更新 + IMPES 时间项装配
///
/// 场景：
/// - 100m x 100m，30x30 剖分，2D 网格；
/// - 初始：全场 p_w = 6.5 MPa，T = 373.15 K；
/// - 左侧一条 CO2 条带：S_w = 0，其余区域 S_w = 1；
/// - 使用 VG-Mualem 计算 Pc, krw, krg；
/// - 常物性：rho_w=1000, rho_g=800, drho/dp=0（在 TwoPhase 里已写死为常数）；
/// - 只组装时间项，不求解线性系统。
int run_IMPES_TwoPhase_TimeTerm_Test()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 20;
    const int sectionNumY = 20;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] IMPES two-phase time-term assembly =====\n";

    std::cout << "--- IMPES: building mesh (CO2 displacement) ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;
    ic.p_w0 = 6.5e6;
    ic.dp_wdx = 1e6;
    ic.T0 = 373.15;
    ic.s_w = 1.0;

    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    // 时间相关辅助场 p_old/prev, s_old/prev, T_old/prev
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");

    // 把当前初场拷贝到 *_old / *_prev，确保时间层一致
    copyField(reg, "p_w", "p_w_old");
    copyField(reg, "p_w", "p_w_prev");
    copyField(reg, "s_w", "s_w_old");
    copyField(reg, "s_w", "s_w_prev");
    copyField(reg, "T", "T_old");
    copyField(reg, "T", "T_prev");

    // 左侧给一条 CO2 条带：将 s_w=0
    {
        const double co2StripeWidth = std::max(1.0, 0.1 * lengthX);
        auto sw = reg.get<volScalarField>("s_w");
        auto sw_old = reg.get<volScalarField>("s_w_old");
        auto sw_prev = reg.get<volScalarField>("s_w_prev");
        if (sw && sw_old && sw_prev)
        {
            const auto& cells = mgr.mesh().getCells();
            const auto& id2idx = mgr.mesh().getCellId2Index();
            for (const auto& c : cells)
            {
                if (c.id < 0) continue;
                if (c.center.m_x <= co2StripeWidth)
                {
                    const size_t idx = id2idx.at(c.id);
                    (*sw)[idx] = 0.0;
                    (*sw_old)[idx] = 0.0;
                    (*sw_prev)[idx] = 0.0;
                }
            }
        }
    }

    // ---------- 2. 基岩 / 流体物性初始化 ----------
    std::cout << "--- IMPES: updating rock/fluid properties (matrix only) ---\n";
    PhysicalPropertiesManager ppm;

    // 区域分类（这里只用一个 RegionType::Medium）
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    // 确保 IMPES 时间项需要的岩石孔隙度/压缩性场存在（phi_r, c_r）
    {
        const std::size_t nCells = mgr.mesh().getCells().size();
        TwoPhase::ensureRockFields(reg, nCells);
    }

    // VG / 相对渗透率参数
    VGParams vg_params;
    vg_params.Swr = 0.20;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;

    RelPermParams rp_params;
    rp_params.L = 0.5;

    // 在时间步开始时：基于 Sw 更新 Pc, krw, krg, dPc_dSw
    TwoPhase::updateTwoPhasePropertiesAtTimeStep(
        mgr, reg, "s_w", vg_params, rp_params);

    // ---------- 3. 压力边界条件（其实对时间项没影响，这里只是保持完整性） ----------
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 4. 井配置（此例为空） ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    build_masks_and_WI_for_all(mgr, reg, wells_cfg);
    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 5. 基本物性更新 & old 层备份 ----------
    std::cout << "--- IMPES: update water/CO2 basic properties at step ---\n";
    // 这里用 p_w 作为“当前迭代的压力评估点”，T 为温度场
    TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, "p_w", "T");
    TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, "p_w", "T");

    // 把当前步的物性复制到 *_old，供时间项使用
    if (!TwoPhase::copyBasicPropertiesToOldLayer(reg))
    {
        std::cerr << "[TEST] copyBasicPropertiesToOldLayer failed.\n";
        return 1;
    }

    // ---------- 6. 时间项装配测试 ----------
    std::cout << "--- IMPES: assemble time term (ITERATIVE formulation) ---\n";

    const double dt = 1.0;   // 测试用时间步长 [s]，只影响 aC, bC 的量级
    const std::string aC_name = "aC_time";
    const std::string bC_name = "bC_time";

    if (!IMPES_Iteration::TimeTerm_IMPES_Pressure(
        mgr, reg, dt,
        /*p_old_name=*/"p_w_old",
        /*p_eval_name=*/"p_w",      // 当前迭代猜测压力，这里用初始 p_w 测一遍
        aC_name,
        bC_name))
    {
        std::cerr << "[TEST] TimeTerm_IMPES_Pressure failed.\n";
        return 1;
    }

    auto aC = reg.get<volScalarField>(aC_name);
    auto bC = reg.get<volScalarField>(bC_name);
    if (!aC || !bC)
    {
        std::cerr << "[TEST] Missing aC or bC field after time-term assembly.\n";
        return 1;
    }

    // ---------- 7. 结果诊断打印 ----------
    double aC_min = 1e300, aC_max = -1e300;
    double bC_min = 1e300, bC_max = -1e300;

    const auto& cells = mgr.mesh().getCells();
    const auto& id2idx = mgr.mesh().getCellId2Index();

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        const double ai = (*aC)[i];
        const double bi = (*bC)[i];

        aC_min = std::min(aC_min, ai);
        aC_max = std::max(aC_max, ai);
        bC_min = std::min(bC_min, bi);
        bC_max = std::max(bC_max, bi);
    }

    std::cout << std::scientific;
    std::cout << "[TEST] aC_time range: [" << aC_min << ", " << aC_max << "]\n";
    std::cout << "[TEST] bC_time range: [" << bC_min << ", " << bC_max << "]\n";

    // 额外打印几格诊断（左侧 CO2 区 vs 右侧水区）
    auto p_w = reg.get<volScalarField>("p_w");
    auto s_w = reg.get<volScalarField>("s_w");
    auto rho_w = reg.get<volScalarField>(TwoPhase::Water().rho_tag);
    auto rho_g = reg.get<volScalarField>(TwoPhase::CO2().rho_tag);

    if (p_w && s_w && rho_w && rho_g)
    {
        std::cout << "\n[TEST] sample cells:\n";
        const double x_mid = 0.5 * lengthX;
        const double x_co2 = 0.05 * lengthX;

        int co2_sample = -1;
        int water_sample = -1;

        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);

            if (co2_sample < 0 && c.center.m_x < x_co2)
            {
                co2_sample = static_cast<int>(i);
            }
            if (water_sample < 0 && c.center.m_x > x_mid)
            {
                water_sample = static_cast<int>(i);
            }
            if (co2_sample >= 0 && water_sample >= 0) break;
        }

        auto print_cell = [&](const char* tag, int idx)
            {
                if (idx < 0) return;
                std::cout << "  [" << tag << "] i=" << idx
                    << "  p_w=" << (*p_w)[idx]
                    << "  s_w=" << (*s_w)[idx]
                    << "  rho_w=" << (*rho_w)[idx]
                    << "  rho_g=" << (*rho_g)[idx]
                    << "  aC=" << (*aC)[idx]
                    << "  bC=" << (*bC)[idx]
                    << "\n";
            };

        print_cell("CO2-zone", co2_sample);
        print_cell("water-zone", water_sample);
    }

    std::cout << "===== [TEST] IMPES two-phase time-term assembly done. =====\n";
    return 0;
}


int run_IMPES_Iteration_TimeTerm_AnalyticalTest()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int sectionNumX = 50;
    const int sectionNumY = 50;
    const int sectionNumZ = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] IMPES two-phase time-term assembly =====\n";

    std::cout << "--- IMPES: building mesh (CO2 displacement) ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;
    ic.p_w0 = 6.5e6;
    ic.dp_wdx =0.0;
    ic.T0 = 373.15;
    ic.s_w = 1.0;

    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    // 时间相关辅助场 p_old/prev, s_old/prev, T_old/prev
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "T", "T_old", "T_prev");

    // 区域分类（这里只用一个 RegionType::Medium）
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    // VG / 相对渗透率参数
    VGParams vg_params;
    vg_params.Swr = 0.20;
    vg_params.Sgr = 0.00;
    vg_params.alpha = 2.0e-6;
    vg_params.n = 2.5;

    RelPermParams rp_params;
    rp_params.L = 0.5;

    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 6.5e6 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };
    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

	// ---------- 6. 时间项装配测试 ----------

    const int    nSteps = 1000;
	const double totalTime = 1000.0;
    const double dt = totalTime / nSteps;

	IMPES_Iteration::PressureSolveControls_analytic pCtrl;
	pCtrl.assembly.pressure_field = "p_w";
	pCtrl.assembly.pressure_old_field = "p_w_old";
	pCtrl.assembly.pressure_prev_field = "p_w_prev";
  

	IMPES_Iteration::SaturationTransportConfig satCfg;
	satCfg.saturation = "s_w";
	satCfg.saturation_old = "s_w_old";
	satCfg.saturation_prev = "s_w_prev";
	satCfg.VG_Parameter.vg_params = vg_params;
	satCfg.VG_Parameter.relperm_params = rp_params;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case3/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case3/s_impes_ps_revised/s_ps";
    const std::string timeSeriesFn = "./Postprocess_Data/IMPES_Iteration_Test/Case3/time_series/impes_ps_revised_stats.csv";  // 如需写出时间序列可填路径

    if (!IMPES_Iteration::runTransient_IMPES_AnalyticTest(
        mgr, reg, freg, PbcA,
        nSteps, dt,
        pCtrl, satCfg,
        10, 10,
        outPrefixP, outPrefixSw,
        10, timeSeriesFn))
    {
        std::cerr << "[TEST] runTransient_IMPES_AnalyticTest failed.\n";
        return 1;
    }

    std::cout << "===== [TEST] IMPES two-phase time-term assembly done. =====\n";
    return 0;
}

/// \brief 测试重力项 + 扩散算子在静水压场下的一致性
/// \details
/// 1) 搭建一个简单 2D 竖直网格；
/// 2) 填充常物性水相/气相 mobility 与密度、渗透率；
/// 3) 调用 IMPES_Iteration::runHydrostaticGravityOperatorTest 生成静水压场 p_w，
///    并用扩散+重力算子 L[p_w] 检查其残差；
/// 4) 返回 0 表示流程执行成功（是否“物理正确”看输出残差大小）。
int TestGravityTerm()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;
    const int    sectionNumX = 50;
    const int    sectionNumY = 50;
    const int    sectionNumZ = 0;
    const bool   usePrism = true;
    const bool   useQuadBase = false;

    std::cout << "\n===== [TEST] Gravity + Diffusion operator (hydrostatic column) =====\n";

    std::cout << "--- building solid matrix mesh ---\n";
    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // 压力瞬态场（IMPES 用的那套）
    ensureTransientFields_scalar(mesh, reg, "p_w", "p_w_old", "p_w_prev");

    const auto& cells = mesh.getCells();
    const size_t nCell = cells.size();

    if (nCell == 0) {
        std::cerr << "[TestGravityTerm] mesh has no cells.\n";
        return 1;
    }

    // ---------------- 1. 常物性两相场 & 渗透率 ----------------
    // 这里直接用 TwoPhase 命名空间里的 tag，避免手写字符串
    const auto water = TwoPhase::Water();
    const auto co2 = TwoPhase::CO2();
    const auto aux = TwoPhase::Auxiliaryparameters();

    // 简单常数：单相水主导，气相 mobility 设 0
    const double rho_w_const = 1000.0;   // kg/m3
    const double rho_g_const = 1.0;      // kg/m3 (几乎不参与)
    const double mu_w_const = 1.0e-3;   // Pa·s
    const double mu_g_const = 1.0e-5;   // Pa·s (无所谓)
    const double kr_w_const = 1.0;
    const double kr_g_const = 0.0;
    const double lambda_w_const = kr_w_const / mu_w_const;
    const double lambda_g_const = kr_g_const / mu_g_const;

    const double k_const = 1.0e-14; // m^2，各向同性

    // 渗透率张量
    auto kxx = reg.getOrCreate<volScalarField>("kxx", nCell, k_const);
    auto kyy = reg.getOrCreate<volScalarField>("kyy", nCell, k_const);
    auto kzz = reg.getOrCreate<volScalarField>("kzz", nCell, k_const);

    // 流体密度 & mobility 场
    auto rho_w = reg.getOrCreate<volScalarField>(water.rho_tag, nCell, rho_w_const);
    auto rho_g = reg.getOrCreate<volScalarField>(co2.rho_tag, nCell, rho_g_const);
    auto lam_w = reg.getOrCreate<volScalarField>(water.lambda_w_tag, nCell, lambda_w_const);
    auto lam_g = reg.getOrCreate<volScalarField>(co2.lambda_g_tag, nCell, lambda_g_const);

    // λ_mass 场：computeDerivedMobilityFields 里用 reg.get() 拿这个场，所以这里先创建
    reg.getOrCreate<volScalarField>(aux.lambda_mass_tag, nCell, 0.0);

    // （如果想更严谨，也可以逐 cell 填值；目前全场常数就够测算子）

    // ---------------- 2. 压力边界条件（这里全零流 Neumann 较简单） ----------------
    const auto& bfaces = mgr.boundaryFaces();

    PressureBC::Registry pbc_pw;

    // a p + b ∂p/∂n = c
    // Neumann: a=0, b=1, c=0 -> zero flux
    PressureBC::BoundaryCoefficient P_Left{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };

    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------------- 3. PressureAssemblyConfig 配置重力 ----------------
    IMPES_Iteration::PressureAssemblyConfig cfg;
    cfg.operator_tag = "hydrostatic_test";
    cfg.pressure_field = "p_w";
    cfg.pressure_old_field = "p_w_old";
    cfg.pressure_prev_field = "p_w_prev";
    cfg.Pc_field = "Pc";          // 这里不会用到，可以留空或确保存在
    cfg.gravity = Vector{ 0.0, -9.81, 0.0 }; // y 轴向上，重力向下
    cfg.enable_buoyancy = true;
    cfg.gradient_smoothing = 0;             // 先关掉平滑，方便看纯离散误差

    // 这些字段名要和 computeDerivedMobilityFields 中一致
    cfg.rho_coeff_field = "rho_coeff_mass";
    cfg.rho_capillary_field = "rho_capillary_mass";
    cfg.rho_gravity_field = "rho_gravity_mass";
    cfg.lambda_gravity_field = "lambda_gravity_mass";
    cfg.gravity_dummy_field = "gravity_dummy_scalar";

    // ---------------- 4. 调用静水算子测试 ----------------
    bool ok = IMPES_Iteration::runHydrostaticGravityOperatorTest(
        mgr, reg, freg, PbcA, cfg);

    if (!ok) {
        std::cerr << "[TestGravityTerm] runHydrostaticGravityOperatorTest failed.\n";
        return 1;
    }

    std::cout << "[TestGravityTerm] finished. Check above residuals for consistency.\n";
    return 0;
}

/**
 * @brief 单独测试毛细扩散算子 L[Pc] = ∇·(k λ_g ρ_g ∇Pc)
 *
 * 测试设置：
 * - kxx=kyy=kzz=1，λ_g=1，ρ_g=1，λ_w=0，ρ_w=0 ⇒ L[Pc] = ∇² Pc；
 * - 构造解析 Pc(x,y) = P0 + α (y - ymin)，满足 ∇²Pc = 0；
 * - 左右边界：Neumann 零通量；上下边界：解析 Dirichlet；
 * - 用 TwoPhaseDiffusionTemplate + assemble_COO + applyCOO 得到 y=A*Pc - b，
 *   检查 max|y|、L2|y| 是否接近 0。
 *
 * @return 0 表示流程执行成功（是否“物理上正确”看终端残差大小），非 0 表示流程错误。
 */
 /// \brief 独立测试“毛细压扩散算子”的零模（常数 Pc）性质.
 ///        期望：在均匀 k、均匀 rho_cap、常数 Pc 条件下，
 ///             离散算子 L[Pc] ≈ 0（只剩机器误差）.
int TestCapillaryTerm_v2()
{
    using namespace IMPES_Iteration;

    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 100.0;
    const double lengthY = 100.0;
    const double lengthZ = 0.0;   // 2D
    const int    sectionNumX = 50;
    const int    sectionNumY = 50;
    const int    sectionNumZ = 0;
    const bool   usePrism = true;
    const bool   useQuadBase = false;

    std::cout << "\n===== [TEST] Capillary diffusion operator (v2, constant Pc) =====\n";

    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    FieldRegistry     reg;
    FaceFieldRegistry freg;

    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    const size_t nCells = cells.size();

    // ---------------- 1. 定义压力与系数字段名 ----------------
    PressureAssemblyConfig cfg;
    cfg.operator_tag = "Pc_capillary_test";  // 用于 OperatorFieldNames 前缀
    cfg.Pc_field = "Pc_test";            // 当前测试用的毛细压场
    cfg.rho_capillary_field = "rho_capillary_mass"; // 可随意命名，和 cfg 一致即可
    cfg.gravity_dummy_field = "gravity_dummy_scalar"; // 这里不用重力，只做占位
    cfg.gravity = Vector{ 0.0, 0.0, 0.0 };
    cfg.enable_buoyancy = false;
    cfg.gradient_smoothing = 0;

    // ---------------- 2. 构建 Pc 场 & k 场 & rho_cap 场 ----------------
    // 2.1: 常数 Pc 场
    auto PcF = reg.getOrCreate<volScalarField>(cfg.Pc_field, nCells, 0.0);

    // 2.2: 渗透率场（各向同性 k = 1.0，实际值只影响尺度）
    auto kxx = reg.getOrCreate<volScalarField>("kxx", nCells, 1.0);
    auto kyy = reg.getOrCreate<volScalarField>("kyy", nCells, 1.0);
    auto kzz = reg.getOrCreate<volScalarField>("kzz", nCells, 1.0);

    // 2.3: “毛细系数” rho_capillary，用来承载 λ_g ρ_g 的等效系数
    auto rho_cap = reg.getOrCreate<volScalarField>(cfg.rho_capillary_field, nCells, 0.0);
    auto rho_dummy = reg.getOrCreate<volScalarField>(cfg.gravity_dummy_field, nCells, 0.0);

    const double Pc0 = 2.0e5;  // 任意常数 Pc [Pa]
    const double rhoCapVal = 1.0;    // 等效 λ_g ρ_g；只要正即可
    const double kVal = 1.0;    // 等效 k；只要正即可

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);

        (*PcF)[i] = Pc0;
        (*kxx)[i] = kVal;
        (*kyy)[i] = kVal;
        (*kzz)[i] = kVal;
        (*rho_cap)[i] = rhoCapVal;
        (*rho_dummy)[i] = 0.0;  // 不使用，只是占位
    }

    // ---------------- 3. 边界条件（毛细项：这里完全不设 BC，相当于自然 Neumann=0） ----------------
    PressureBC::Registry pbc_pc;  // 空 registry，不添加任何边界
    PressureBCAdapter    PbcA{ pbc_pc };

    // ---------------- 4. 构建未知量编号映射 ----------------
    int N = 0;
    auto lid_of_cell = buildUnknownMap(mesh, N);
    if (N <= 0)
    {
        std::cerr << "[CapillaryTest_v2] buildUnknownMap failed: N=" << N << "\n";
        return -1;
    }

    // ---------------- 5. 构造 OperatorFieldNames 和 mobility_tokens ----------------
    // 只需要 kxx/kyy/kzz 这三个渗透率场；rho 部分由 rho_capillary_field 承载
    auto nm = makeNames(cfg.operator_tag); // a_f_diff, s_f_diff 等名字

    std::vector<std::string> mobility_tokens{
        "kxx:kxx",
        "kyy:kyy",
        "kzz:kzz"
    };

    // ---------------- 6. 收集 Pc 场到向量 Pc_vec ----------------
    std::vector<double> Pc_vec, Lvec;
    if (!detailed::gatherFieldToVector(mesh, reg, cfg.Pc_field,
        lid_of_cell, N, Pc_vec))
    {
        std::cerr << "[CapillaryTest_v2] gatherFieldToVector failed for Pc.\n";
        return -2;
    }

    // ---------------- 7. 构建算子并作用：Lvec = A * Pc_vec - b ----------------
    if (!detailed::buildAndApplyOperator(
        mgr, reg, freg,
        PbcA,
        nm,
        mobility_tokens,
        cfg.rho_capillary_field,  // rho_coeff_field → 用毛细系数场
        cfg.gravity_dummy_field,  // rho_buoy_field  → 不用
        cfg.Pc_field,             // x_field → Pc
        cfg,
        /*enable_buoyancy=*/false,
        Pc_vec,
        Lvec))
    {
        std::cerr << "[CapillaryTest_v2] buildAndApplyOperator failed.\n";
        return -3;
    }

    // ---------------- 8. 计算残差范数并输出 ----------------
    double maxAbs = 0.0;
    double l2sum = 0.0;
    int    nUsed = 0;

    int    sampleCellId = -1;
    double samplePc = 0.0;
    double sampleL = 0.0;

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        const int    lid = lid_of_cell[i];
        if (lid < 0) continue;

        const double val = std::fabs(Lvec[lid]);
        maxAbs = std::max(maxAbs, val);
        l2sum += val * val;
        ++nUsed;

        // 记录一个样本单元
        if (sampleCellId < 0)
        {
            sampleCellId = c.id;
            samplePc = (*PcF)[i];
            sampleL = Lvec[lid];
        }
    }

    const double l2norm = (nUsed > 0) ? std::sqrt(l2sum / nUsed) : 0.0;
    const double relMax = (std::fabs(Pc0) > 0.0) ? maxAbs / std::fabs(Pc0) : 0.0;

    std::cout << "[CapillaryTest_v2] N = " << nUsed
        << ", max|L[Pc]| = " << maxAbs
        << ", L2|L[Pc]| = " << l2norm
        << ", relMax = " << relMax << "\n";

    if (sampleCellId >= 0)
    {
        std::cout << "  sample cell id = " << sampleCellId
            << ", Pc = " << samplePc
            << ", L[Pc] = " << sampleL << "\n";
    }

    std::cout << "[CapillaryTest_v2] finished. "
        "Expect residuals near machine precision if capillary operator is assembled correctly.\n";

    return 0;
}


/// \brief 测试 assemblePressureTwoPhase：退化两相(单水相)静水平衡 + 重力
int Test_AssemblePressureTwoPhase_Hydrostatic()
{
    using namespace IMPES_Iteration;

    // ---------------- 0. 网格与场注册 ----------------
    const double Lx = 100.0;
    const double Ly = 100.0;
    const double Lz = 0.0;
    const int nx = 50;
    const int ny = 50;
    const int nz = 0;
    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] assemblePressureTwoPhase: hydrostatic limit =====\n";

    MeshManager mgr(Lx, Ly, Lz,
        nx, ny, nz,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    Mesh& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // 压力时间场：p_w, p_w_old, p_w_prev
    ensureTransientFields_scalar(mesh, reg, "p_w", "p_w_old", "p_w_prev");

    // ---------------- 1. 两相相关体场：退化为单水相 ----------------
    const size_t nCells = cells.size();

    // 饱和度：全水相 Sw = 1
    auto Sw = reg.getOrCreate<volScalarField>("s_w", nCells, 1.0);
    std::fill(Sw->data.begin(), Sw->data.end(), 1.0);

    // 常数物性
    const double rho_w_const = 1000.0;   // [kg/m3]
    const double rho_g_const = 600.0;    // [kg/m3] 这里不会用到
    const double mu_w_const = 1.0e-3;   // [Pa·s]
    const double k_abs = 1.0e-14;  // [m2]
    const double phi_const = 0.1;      // [-]

    // mobility: λ_w = k_rw / μ_w，这里直接当成常数 1 / μ_w
    const double lambda_w_const = 1.0 / mu_w_const;
    const double lambda_g_const = 0.0;   // 退化掉气相

    auto rho_w = reg.getOrCreate<volScalarField>(TwoPhase::Water().rho_tag.c_str(),
        nCells, rho_w_const);
    auto rho_g = reg.getOrCreate<volScalarField>(TwoPhase::CO2().rho_tag.c_str(),
        nCells, rho_g_const);
    auto lambda_w = reg.getOrCreate<volScalarField>(TwoPhase::Water().lambda_w_tag.c_str(),
        nCells, lambda_w_const);
    auto lambda_g = reg.getOrCreate<volScalarField>(TwoPhase::CO2().lambda_g_tag.c_str(),
        nCells, lambda_g_const);

    // lambda_mass = λ_w ρ_w + λ_g ρ_g，将由 computeDerivedMobilityFields 更新
    auto lambda_mass = reg.getOrCreate<volScalarField>(
        TwoPhase::Auxiliaryparameters().lambda_mass_tag.c_str(),
        nCells, 0.0);

    // 岩石孔隙度和压缩性，用于时间项
    auto phi_r = reg.getOrCreate<volScalarField>(TwoPhase::Rock().phi_tag, nCells, phi_const);
    auto c_phi = reg.getOrCreate<volScalarField>(TwoPhase::Rock().c_r_tag, nCells, 0.0); // 不考虑岩石压缩

    // 这里假设 TimeTerm_IMPES_Pressure 需要的 rho_old和 drdp 场：
    // 简化成常数（并且我们用 p_new = p_old → 时间残差为 0）
    auto rho_w_old = reg.getOrCreate<volScalarField>(TwoPhase::Water().rho_old_tag, nCells, rho_w_const);
    auto rho_g_old = reg.getOrCreate<volScalarField>(TwoPhase::CO2().rho_old_tag, nCells, rho_w_const);
    auto drho_dp_w = reg.getOrCreate<volScalarField>(TwoPhase::Water().drho_w_dp_tag, nCells, 0.0);
    auto drho_dp_w_old = reg.getOrCreate<volScalarField>(TwoPhase::Water().drho_w_dp_old_tag, nCells, 0.0);
    auto drho_dp_g = reg.getOrCreate<volScalarField>(TwoPhase::CO2().drho_g_dp_tag, nCells, 0.0);
    auto drho_dp_g_old = reg.getOrCreate<volScalarField>(TwoPhase::CO2().drho_g_dp_old_tag, nCells, 0.0);

    // 渗透率张量：各向同性
    auto kxx = reg.getOrCreate<volScalarField>("kxx", nCells, k_abs);
    auto kyy = reg.getOrCreate<volScalarField>("kyy", nCells, k_abs);
    auto kzz = reg.getOrCreate<volScalarField>("kzz", nCells, k_abs);

    // 毛细压力场：设为常数 0 → 不产生毛细贡献
    auto Pc = reg.getOrCreate<volScalarField>("Pc", nCells, 0.0);
    std::fill(Pc->data.begin(), Pc->data.end(), 0.0);

    // 预留：rho_coeff / rho_capillary / rho_gravity / rho_mix / lambda_gravity / gravity_dummy
    auto rho_coeff = reg.getOrCreate<volScalarField>("rho_coeff_mass", nCells, 0.0);
    auto rho_capillary = reg.getOrCreate<volScalarField>("rho_capillary_mass", nCells, 0.0);
    auto rho_gravity = reg.getOrCreate<volScalarField>("rho_gravity_mass", nCells, 0.0);
    auto rho_mix = reg.getOrCreate<volScalarField>("rho_mix_mass", nCells, 0.0);
    auto lambda_gravity = reg.getOrCreate<volScalarField>("lambda_gravity_mass", nCells, 0.0);
    auto grav_dummy = reg.getOrCreate<volScalarField>("gravity_dummy_scalar", nCells, 0.0);

    // ---------------- 2. 构造静水压力场 p_w, p_w_old ----------------
    auto p_w = reg.get<volScalarField>("p_w");
    auto p_w_old = reg.get<volScalarField>("p_w_old");

    const double g = 9.81;
    const double gy = -g; // y 向上，重力向下
    const double y_top = Ly;
    const double p_top = 9.0e6; // 顶部压力

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        const double y = c.center.m_y;

        const double p_exact = p_top + rho_w_const * gy * (y - y_top);
        (*p_w)[i] = p_exact;
        (*p_w_old)[i] = p_exact;  // 保证时间项残差为 0
  
    }

    // ---------------- 3. 边界条件：左右无流，上下给定静水压力 ----------------
    const auto& bfaces = mgr.boundaryFaces();

    // 这里简单做法：左右 Neumann(无流) → a=0,b=1,c=0；上下 Dirichlet
    // 顶边：p = p_top
    // 底边：p = p_top + rho_w * g_y * (0 - y_top) = p_top - rho_w*g*Ly
    const double p_bottom = p_top + rho_w_const * gy * (0.0 - y_top);

    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 0.0, 1.0, 0.0 };       // no-flow
    PressureBC::BoundaryCoefficient P_Right{ 0.0, 1.0, 0.0 };       // no-flow
    PressureBC::BoundaryCoefficient P_Down{ 1.0, 0.0, p_bottom };  // Dirichlet
    PressureBC::BoundaryCoefficient P_Up{ 1.0, 0.0, p_top };     // Dirichlet

    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    pbc_pw.printInformationofBoundarySetting(mgr);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------------- 4. 组装两相压力方程（退化单相） ----------------
    PressureAssemblyConfig cfg;
    cfg.operator_tag = "p_w_IMPES_test";
    cfg.pressure_field = "p_w";
    cfg.pressure_old_field = "p_w_old";
    cfg.pressure_prev_field = "p_w_prev";
    cfg.Pc_field = "Pc";

    cfg.gravity = Vector{ 0.0, gy, 0.0 };
    cfg.enable_buoyancy = true;
    cfg.gradient_smoothing = 0;

    cfg.rho_coeff_field = "rho_coeff_mass";
    cfg.rho_capillary_field = "rho_capillary_mass";
    cfg.rho_gravity_field = "rho_gravity_mass";
    cfg.lambda_gravity_field = "lambda_gravity_mass";
    cfg.gravity_dummy_field = "gravity_dummy_scalar";
    cfg.rho_mix_field = "rho_mix_mass";

    // 让 assemblePressureTwoPhase 真的去计算总通量
    cfg.total_mass_flux_name = "mf_tot";
    cfg.total_vol_flux_name = "Qf_tot";
    cfg.total_velocity_name = "ufn_tot";

    // wells: 本测试不含井
    std::vector<WellDOF_TwoPhase> wells;

    PressureAssemblyResult result;
    const double dt = 1.0; // 任意正数即可，p_new = p_old → 时间残差为 0

    if (!assemblePressureTwoPhase(mgr, reg, freg, PbcA, wells, dt, cfg, result))
    {
        std::cerr << "[Test_AssemblePressureTwoPhase] assemblePressureTwoPhase failed.\n";
        return -1;
    }

    // ---------------- 5. 用解析解 p_exact 检查离散残差 L[p] ----------------
    const int N = result.system.n;
    std::vector<double> p_vec(N, 0.0), Lvec;

    // cell_lid: 每个 cell 对应的未知量局部编号
    const auto& lid_of_cell = result.cell_lid;

    for (const auto& c : cells)
    {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        const int lid = lid_of_cell[i];
        if (lid < 0 || lid >= N) continue;
        p_vec[lid] = (*p_w)[i];
    }

    IMPES_Iteration::detailed::applyCOO(result.system, p_vec, Lvec);

    double maxAbs = 0.0;
    double sum2 = 0.0;
    int    cnt = 0;
    for (int i = 0; i < static_cast<int>(Lvec.size()); ++i)
    {
        const double v = std::abs(Lvec[i]);
        maxAbs = std::max(maxAbs, v);
        sum2 += v * v;
        ++cnt;
    }
    const double L2 = (cnt > 0) ? std::sqrt(sum2 / cnt) : 0.0;
    const double relMax = (p_top > 0.0 ? maxAbs / p_top : maxAbs);

    std::cout << "[AssemblePressureTwoPhase-Hydrostatic] N = " << N
        << ", max|L[p]| = " << maxAbs
        << ", L2|L[p]| = " << L2
        << ", relMax = " << relMax << "\n";

    // ---------------- 6. 检查总质量通量是否接近 0 ----------------
    auto mf_tot = freg.get<faceScalarField>("mf_tot");
    if (mf_tot)
    {
        double maxMf = 0.0;
        for (double mface : mf_tot->data)
        {
            maxMf = std::max(maxMf, std::abs(mface));
        }
        std::cout << "  max|mf_tot| = " << maxMf << " [kg/s per face]\n";
    }
    else
    {
        std::cout << "  [Warn] mf_tot not found, skip flux check.\n";
    }

    std::cout << "[AssemblePressureTwoPhase-Hydrostatic] finished.\n";
    return 0;
}

/**
 * @brief 1D Buckley–Leverett 风格的两相 IMPES 数值测试（仅数值解，不构造解析解）。
 *
 * 设定：
 *  - 1D 区域：x ∈ [0, Lx]，y 方向只有 1 层单元；
 *  - 左端高压、右端低压，形成近似常总通量的两相驱替；
 *  - 初场为典型 Riemann 型：x≈0 处 S_w ≈ S_w_inj，其余区域 S_w ≈ S_w_i；
 *  - VG 模型设置为 Swr = Sgr = 0，且 alpha 较大 → Pc 数值很小，对应经典 Buckley–Leverett 假设。
 *
 * 流程：
 *  1) 构建 1D 网格 + 主变量场 (p_w, s_w, T)；
 *  2) 基于 p_w, T 更新基岩物性；
 *  3) 施加 BL 风格饱和度初场；
 *  4) 设置 VG/Mualem 参数、压力边界条件与 IMPES 控制参数；
 *  5) 调用 IMPES_Iteration::runTransient_IMPES_Iteration 做主时间步推进；
 *  6) 输出 Tecplot 结果用于后处理。
 */
int run_IMPES_Iteration_TwoPhase_BL_Numerical()
{
    // ---------------- 0. 网格与场注册 ----------------
    const double lengthX = 50.0;   // [m]
    const double lengthY = 1.0;    // 几乎 1D
    const double lengthZ = 0.0;

    const int sectionNumX = 50;    // 细一点，方便看前沿
    const int sectionNumY = 5;
    const int sectionNumZ = 0;

    const bool usePrism = true;
    const bool useQuadBase = false;

    std::cout << "\n===== [TEST] IMPES two-phase (1D Buckley–Leverett, numerical only) =====\n";
    std::cout << "--- IMPES: building mesh (1D BL test) ---\n";

    MeshManager mgr(lengthX, lengthY, lengthZ,
        sectionNumX, sectionNumY, sectionNumZ,
        usePrism, useQuadBase);
    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);

    FieldRegistry     reg;
    FaceFieldRegistry freg;

    // ---------- 1. 初始主变量场 ----------
    std::cout << "--- IMPES: initializing primary fields ---\n";
    InitFields ic;  // 默认 p0 / T0 / Sw0，如有需要可以修改 ic.p0 / ic.T0 等
    // 压力：线性从左 9 MPa 过渡到右 8 MPa
    ic.p_w0 = 9.0e6;                        // x=0 处的压力
    ic.dp_wdx = (8.0e6 - 9.0e6) / lengthX;  // 与 lengthX 匹配的梯度
    ic.dp_wdy = 0.0;
    ic.dp_wdz = 0.0;
    // BL: 背景初始饱和度 = 右侧初始值 Sw_i
    ic.s_w = 0.8;   // Sw_i

    // 主变量场：水相压力 p_w、水相饱和度 s_w、温度 T
    Initializer::createPrimaryFields_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, "p_w", "s_w", "T");
    Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(
        mgr.mesh(), reg, ic);

    // 时间层：old / prev（p_w、s_w、p_g 都准备好）
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_w", "p_w_old", "p_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "s_w", "s_w_old", "s_w_prev");
    ensureTransientFields_scalar(mgr.mesh(), reg, "p_g", "p_g_old", "p_g_prev");

    // ---------- 2. 基岩区域分类与固相物性 ----------
    PhysicalPropertiesManager ppm;
    ppm.classifyRockRegionsByGeometry(mgr, {}, Cell::RegionType::Medium);
    ppm.UpdateMatrixRockAt(mgr, reg, "p_w", "T");

    // 3. Buckley–Leverett 型饱和度初场：避免 0/1 极端，留出 0.1 缓冲
    {
        auto sw = reg.get<volScalarField>("s_w");
        auto sw_old = reg.get<volScalarField>("s_w_old");
        auto sw_prev = reg.get<volScalarField>("s_w_prev");
        if (!sw || !sw_old || !sw_prev) { std::cerr << "[BL-TEST] saturation fields not found.\n"; return EXIT_FAILURE; }

        const auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        const double Sw_inj = 0.05; // 注入端，远离 1.0
        const double Sw_i = 0.85; // 背景，远离 0.0
        const double dx = lengthX / static_cast<double>(sectionNumX);
        const double x_threshold = dx;

        for (const auto& c : cells) {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            double Sw_val = (c.center.m_x < x_threshold) ? Sw_inj : Sw_i;
            Sw_val = std::max(0.05 + 1e-6, std::min(1.0 - 0.15 - 1e-6, Sw_val)); // clamp 至 [Swr, 1-Sgr]
            (*sw)[i] = (*sw_old)[i] = (*sw_prev)[i] = Sw_val;
        }
    }

    // 4. VG / 相对渗透率参数：进一步减小 Pc
    VGParams vg_params;
    vg_params.Swr = 0.05;  // 少量束缚水
    vg_params.Sgr = 0.15;  // 适当残余气，避免气相端点过头
    vg_params.alpha = 1e4;   // 保持 Pc 小
    vg_params.n = 1.3;   // 曲线更平缓，前沿不至于过陡

    RelPermParams rp_params;
    rp_params.L = 0.5;       // Mualem 默认，减缓相对渗透率陡峭度


    // ---------- 5. 压力边界条件（左高右低，其他 no-flow） ----------
    const auto& bfaces = mgr.boundaryFaces();
    PressureBC::Registry pbc_pw;
    PressureBC::BoundaryCoefficient P_Left{ 1.0, 0.0, 9.0e6 };
    PressureBC::BoundaryCoefficient P_Right{ 1.0, 0.0, 8.0e6 };
    PressureBC::BoundaryCoefficient P_Down{ 0.0, 1.0, 0.0 };
    PressureBC::BoundaryCoefficient P_Up{ 0.0, 1.0, 0.0 };

    PressureBC::setBoxBCs2D(pbc_pw, bfaces, P_Left, P_Right, P_Down, P_Up);
    PressureBCAdapter PbcA{ pbc_pw };

    // ---------- 6. 井配置（BL 测试：无井源） ----------
    std::vector<WellConfig_TwoPhase> wells_cfg;
    build_masks_and_WI_for_all(mgr, reg, wells_cfg);

    const int Ncells = static_cast<int>(mgr.mesh().getCells().size());
    std::vector<WellDOF_TwoPhase> wells_dof;
    register_well_dofs_for_all_TwoPhase(Ncells, wells_cfg, wells_dof);

    // ---------- 7. 饱和度输运配置 ----------
    IMPES_Iteration::SaturationTransportConfig satCfg;
    satCfg.saturation = "s_w";
    satCfg.saturation_old = "s_w_old";
    satCfg.saturation_prev = "s_w_prev";
    satCfg.VG_Parameter.vg_params = vg_params;
    satCfg.VG_Parameter.relperm_params = rp_params;

    // ---------- 8. 压力求解控制参数 ----------
    IMPES_Iteration::PressureSolveControls pCtrl;
    pCtrl.max_outer = 1;        // 常物性可保持很小
    pCtrl.tol_abs = 1e15;
    pCtrl.tol_rel = 1e-4;
    pCtrl.under_relax = 1;
    pCtrl.verbose = true;

    pCtrl.assembly.enable_buoyancy = false;
    pCtrl.assembly.gradient_smoothing = 1;
    pCtrl.assembly.gravity = Vector{ 0.0, 0.0, 0.0 };

    // ---------- 9. 两相通量拆分配置 ----------
    IMPES_Iteration::FluxSplitConfig fluxCfg;
    fluxCfg.rho_water = TwoPhase::Water().rho_tag; // "rho_w"
    fluxCfg.rho_gas = TwoPhase::CO2().rho_tag;   // "rho_g"
    fluxCfg.pressure_bc = &PbcA;

    // ---------- 10. 初场后处理：同步 Pc / p_g / 物性到 old 层 ----------
    {
        // 更新 Pc, krw, krg, dPc/dSw
        TwoPhase::updateTwoPhasePropertiesAtTimeStep(mgr, reg, "s_w", vg_params, rp_params);

        // 设置 p_g = p_w + Pc，并同步到 old/prev
        auto p_w = reg.get<volScalarField>("p_w");
        auto p_g = reg.get<volScalarField>("p_g");
        auto p_w_old = reg.get<volScalarField>("p_w_old");
        auto p_w_prev = reg.get<volScalarField>("p_w_prev");
        auto p_g_old = reg.get<volScalarField>("p_g_old");
        auto p_g_prev = reg.get<volScalarField>("p_g_prev");
        auto Pc = reg.get<volScalarField>(TwoPhase::Auxiliaryparameters().Pc_tag);
        const auto& cells = mgr.mesh().getCells();
        const auto& id2idx = mgr.mesh().getCellId2Index();
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            const double pg_val = (*p_w)[i] + (*Pc)[i];
            (*p_g)[i] = pg_val;
            (*p_g_old)[i] = pg_val;
            (*p_g_prev)[i] = pg_val;
            // 确保 p_w_old/p_w_prev 也和 p_w 同步
            (*p_w_old)[i] = (*p_w)[i];
            (*p_w_prev)[i] = (*p_w)[i];
        }

        // 更新水/CO2 物性并拷贝到 *_old 层
        TwoPhase::updateWaterBasicPropertiesAtStep(mgr, reg, "p_w", "T");
        TwoPhase::updateCO2BasicPropertiesAtStep(mgr, reg, "p_g", "T");
        TwoPhase::copyBasicPropertiesToOldLayer(reg);
    }

    // ---------- 11. IMPES 主时间推进 ----------
    const int    nSteps = 2000;
    double       dt_initial = 1e-4;   // 更保守的初始时间步

    const int writeEveryP = 1;
    const int writeEverySw = 1;

    const std::string outPrefixP = "./Postprocess_Data/IMPES_Iteration_Test/Case4/p_impes_ps_revised/p_ps";
    const std::string outPrefixSw = "./Postprocess_Data/IMPES_Iteration_Test/Case4/s_impes_ps_revised/s_ps";
    const int snapshotEveryCsv = 1;
    const std::string snapshotPrefix = "./Postprocess_Data/csv_snapshots/Case4/ps_state_reviesed";
    std::vector<std::string> snapshotFields;

    std::cout << "--- IMPES: start transient run (BL numerical test) ---\n";

    const bool ok = IMPES_Iteration::runTransient_IMPES_Iteration
    (
        mgr,
        reg,
        freg,
        PbcA,
        wells_dof,
        nSteps,
        dt_initial,
        pCtrl,
        satCfg,
        fluxCfg,
        writeEveryP,
        writeEverySw,
        outPrefixP,
        outPrefixSw,
        snapshotEveryCsv,
        snapshotPrefix,
        snapshotFields
    );

    if (!ok)
    {
        std::cerr << "[BL-TEST] IMPES transient run failed.\n";
        return EXIT_FAILURE;
    }

    std::cout << "[BL-TEST] IMPES two-phase Buckley–Leverett numerical test finished successfully.\n";
    return EXIT_SUCCESS;
}


int main()
{
    return run_IMPES_Iteration_TwoPhase_BL_Numerical();
}

