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


#include "MultiPhaseProperties.h"


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
int main000000000000()
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
        /*outPrefixP*/ "./Postprocess_Data/P_CO2_BDF2/p",
        /*outPrefixT*/ "./Postprocess_Data/T_CO2_BDF2/T"
    );

    if (!ok) {
        std::cerr << "[MAIN] Transient run failed.\n";
        return 1;
    }

    std::cout << "[MAIN] Done.\n";
    return 0;

}

//实现CO2单相-常物性-渗流传热带井源控制方程的离散与求解 不带裂缝
int main00000000000()
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

int main()
{
    return runCO2Displacement_IMPES_Ps();
}

