/**
 * @file  Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp
 * @brief 独立测试文件：2D 单相 CO2 常物性、无裂缝、无井的压力扩散全流程
 *
 * @details
 * 本文件是一个教学模板 + 可扩展基准，演示完整的 FVM-EDFM 瞬态求解流程：
 *
 *   1. 构建网格
 *   2. 建立全局索引与拓扑
 *   3. 初始化场管理器
 *   4. 设置初始/边界条件
 *   5. 配置求解器模块
 *   6. 运行瞬态求解器
 *   7. 后处理与 VTK 导出
 *
 * 架构：三层结构（名称索引 -> 参数配置 -> 完整管道执行）
 * 独立性：不导入 FullCaseTest 中的任何函数/类型
 * 扩展性：后续可通过修改 TestCaseSpec 添加裂缝、井、变物性等
 */

#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"

// --- 引擎与求解器基础设施 ---
#include "2D_PostProcess.h"              // PostProcess_2D：VTK 导出
#include "BoundaryConditionManager.h"    // BoundarySetting::BoundaryConditionManager
#include "FIM_TransientCaseKit.hpp"      // FIM_CaseKit::InitFieldManager + 传递引入 RunGenericFIMTransient
#include "MeshDefinitions.h"             // MeshTags::LEFT/RIGHT/TOP/BOTTOM
#include "MeshManager.h"                 // MeshManager：2D 网格管理
#include "SolverContrlStrName_op.h"      // PhysicalProperties_string_op：字段名标签
#include "Well_WellControlTypes.h"       // WellScheduleStep（空井列表需要）

// --- 标准库 ---
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// --- 平台相关目录创建 ---
#ifdef _WIN32
#include <direct.h>
#define TEST_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define TEST_MKDIR(path) mkdir(path, 0777)
#endif

// =====================================================================
// 命名空间开始
// =====================================================================
namespace Test_H_CO2_ConstPP {
namespace {

// =====================================================================
// === 本地工具函数（从 FullCaseTest.cpp 独立复制，不引用其命名空间） ===
// =====================================================================

/**
 * @brief 递归创建目录路径
 *
 * 将路径中的反斜杠统一为正斜杠，然后逐级创建。
 * 复制自 FullCaseTest.cpp:153-169，因该函数在匿名命名空间中无法外部访问。
 */
void EnsureDirRecursive(const std::string& rawPath) {
    if (rawPath.empty()) return;
    std::string path = rawPath;
    for (char& ch : path) {
        if (ch == '\\') ch = '/';
    }

    std::stringstream ss(path);
    std::string token;
    std::string current;
    while (std::getline(ss, token, '/')) {
        if (token.empty() || token == ".") continue;
        if (!current.empty()) current += "/";
        current += token;
        TEST_MKDIR(current.c_str());
    }
}

/**
 * @brief 将标量场的所有值设置为均匀常数
 *
 * volScalarField（即 VolField<double>）没有内置的 fill/assign 方法，
 * 因此需要手动遍历 data 向量进行赋值。
 * 复制自 FullCaseTest.cpp:185-188。
 */
void ApplyUniformScalarField(const std::shared_ptr<volScalarField>& field, double value) {
    if (!field) return;
    for (double& v : field->data) v = value;
}

/**
 * @brief 将压力向量和恒温值同步写入 FieldManager 的标量场
 *
 * 功能说明：
 *   - 从求解器输出的 pBlocks 向量中读取每个 DOF 的压力值
 *   - 分别写入基岩（i < nMat）和裂缝（i >= nMat）的压力/温度场
 *   - 同时维护 "P" 可视化场用于 VTK 导出
 *
 * 复制自 FullCaseTest.cpp:224-258（原名 SyncN1SnapshotPressureFields），
 * 仅调用 FieldManager_2D 和 MeshManager 的公共 API。
 */
void SyncPressureFieldsToFM(MeshManager& mgr,
                             FieldManager_2D& fm,
                             const std::vector<double>& pBlocks,
                             double tConst) {
    if (pBlocks.empty()) return;

    // 获取压力/温度方程的字段名标签
    const auto pCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();

    // 基岩场：压力求解场、可视化场、温度场
    auto fPw   = fm.getOrCreateMatrixScalar(pCfg.pressure_field, 0.0);
    auto fPviz = fm.getOrCreateMatrixScalar("P", 0.0);
    auto fT    = fm.getOrCreateMatrixScalar(tCfg.temperatue_field, tConst);

    // 裂缝场（无裂缝时 size=0，循环自动跳过）
    auto fracPw   = fm.getOrCreateFractureScalar(pCfg.pressure_field, 0.0);
    auto fracPviz = fm.getOrCreateFractureScalar("P", 0.0);
    auto fracT    = fm.getOrCreateFractureScalar(tCfg.temperatue_field, tConst);

    const int nMat   = mgr.getMatrixDOFCount();
    const int nTotal = mgr.getTotalDOFCount();
    const int nUse   = std::min(static_cast<int>(pBlocks.size()), nTotal);

    for (int i = 0; i < nUse; ++i) {
        const double p = pBlocks[static_cast<std::size_t>(i)];
        if (i < nMat) {
            // 基岩 DOF
            if (fPw   && i < static_cast<int>(fPw->data.size()))   fPw->data[static_cast<std::size_t>(i)]   = p;
            if (fPviz && i < static_cast<int>(fPviz->data.size())) fPviz->data[static_cast<std::size_t>(i)] = p;
            if (fT    && i < static_cast<int>(fT->data.size()))    fT->data[static_cast<std::size_t>(i)]    = tConst;
        } else {
            // 裂缝 DOF
            const int fi = i - nMat;
            if (fracPw   && fi < static_cast<int>(fracPw->data.size()))   fracPw->data[static_cast<std::size_t>(fi)]   = p;
            if (fracPviz && fi < static_cast<int>(fracPviz->data.size())) fracPviz->data[static_cast<std::size_t>(fi)] = p;
            if (fracT    && fi < static_cast<int>(fracT->data.size()))    fracT->data[static_cast<std::size_t>(fi)]    = tConst;
        }
    }
}

/**
 * @brief 计算网格特征长度 h = sqrt(总体积 / 单元数)
 *
 * 用于后处理指标输出。
 * 复制自 FullCaseTest.cpp:323-330。
 */
double ComputeMeshCharLength(const MeshManager& mgr) {
    const auto& cells = mgr.mesh().getCells();
    if (cells.empty()) return 0.0;
    double totalV = 0.0;
    for (const auto& c : cells) totalV += std::max(c.volume, 0.0);
    if (totalV <= 0.0) return 0.0;
    return std::sqrt(totalV / static_cast<double>(cells.size()));
}

// =====================================================================
// === Layer 2：参数容器与配置 ===
// =====================================================================

/**
 * @brief 测试工况参数容器
 *
 * 包含几何、物性、流体、边界、时间步、求解器等全部可配置参数。
 * 默认值与 FullCaseTest 中 N1CaseSpec 的常物性无裂缝工况一致，
 * 确保结果可对比验证。
 *
 * 扩展方式：后续添加裂缝/井/变物性时，在此结构体中增加相应字段即可。
 */
struct TestCaseSpec {
    // --- 标识 ---
    std::string case_name = "h_co2_constpp_nofrac_nowell";  // 输出子目录名
    std::string output_base_dir = "Test/Transient/FullCaseTest";  // 输出根目录（可外部修改）
    std::string sub_dir = "H_CO2_ConstPP";  // 输出中间目录

    // --- 计算域几何 ---
    double lx = 400.0;   // 域长度 [m]
    double ly = 40.0;    // 域高度 [m]
    int    nx = 48;      // x 方向网格数
    int    ny = 6;       // y 方向网格数

    // --- 基岩物性 ---
    double matrix_phi  = 0.10;     // 孔隙度 [-]
    double matrix_perm = 1.0e-13;  // 渗透率 [m^2]
    double matrix_ct   = 5.0e-9;   // 岩石压缩系数 [1/Pa]

    // --- 流体常物性 ---
    double mu_const  = 6.0e-5;   // 粘度 [Pa*s]
    double rho_const = 700.0;    // 密度 [kg/m^3]

    // --- 初始/边界条件 ---
    double p_init  = 10.0e6;   // 初始压力 [Pa]
    double p_left  = 12.0e6;   // 左边界 Dirichlet 压力 [Pa]
    double p_right =  8.0e6;   // 右边界 Dirichlet 压力 [Pa]
    double t_init  = 360.0;    // 初始温度 [K]（常物性模式下作为 EOS 评估温度）

    // --- 时间步控制 ---
    double dt_init          = 120.0;   // 初始时间步 [s]
    double dt_min           = 1.0;     // 最小时间步 [s]
    double dt_max           = 5000.0;  // 最大时间步 [s]
    double target_end_time_s = 1.0e5;  // 模拟终止时间 [s]
    int    max_steps        = 12000;   // 最大时间步数
    int    max_newton_iter  = 12;      // 每步最大 Newton 迭代数

    // --- 求解器配置（可外部调整） ---

    // 线性求解器
    FIM_Engine::LinearSolverType lin_solver = FIM_Engine::LinearSolverType::AMGCL;  // 线性求解器类型
    bool   amgcl_use_fallback_sparselu     = true;    // AMGCL 失败时回退 SparseLU
    double amgcl_tol                       = 1.0e-6;  // AMGCL 收敛容差
    int    amgcl_maxiter                   = 500;     // AMGCL 最大迭代数

    // 非正交修正与行缩放
    bool enable_non_orthogonal_correction  = true;    // 非正交修正（斜网格时建议开启）
    bool enable_row_scaling                = true;    // Jacobian 行缩放

    // Newton 迭代容差
    double abs_res_tol                     = 1.0e-10; // 绝对残差容差 [Pa]
    double rel_res_tol                     = 1.0e-6;  // 相对残差容差 [-]
    double rel_update_tol                  = 1.0e-8;  // 相对更新容差 [-]

    // 状态限幅
    double max_dP                          = 2.0e7;   // 单次 Newton 最大压力变化 [Pa]

    // 线搜索
    bool enable_armijo_line_search         = false;   // Armijo 线搜索（常物性简单工况可关闭）

    // 自适应时间步控制
    double rollback_shrink_factor          = 0.7;     // 回退时 dt 缩小因子 [-]
    double dt_relres_grow_factor           = 1.08;    // 收敛良好时 dt 增长因子 [-]

    // 重力
    Vector gravity_vector = Vector(0.0, -9.81, 0.0);  // 重力加速度 [m/s^2]

    // 诊断输出
    FIM_Engine::DiagLevel diag_level = FIM_Engine::DiagLevel::Off;  // Off/Summary/Hotspot/Forensic
};

/**
 * @brief 测试运行摘要（求解完成后填充）
 *
 * 暂不含解析解验证字段，后续扩展时添加。
 */
struct TestCaseSummary {
    std::string case_dir;               // 输出目录路径
    std::string convergence_log_path;   // 收敛日志路径
    std::string metrics_csv_path;       // 指标 CSV 路径

    int    nx = 0, ny = 0, n_cells = 0; // 网格信息
    double h_char = 0.0;                // 网格特征长度 [m]

    int    steps = 0;                   // 实际步数
    int    total_rollbacks = 0;         // 时间步回退次数
    double avg_iters = 0.0;            // 平均 Newton 迭代数
    int    max_iters = 0;              // 最大 Newton 迭代数
    double t_end = 0.0;                // 实际终止时间 [s]
};

/**
 * @brief 测试计划（注册表使用）
 */
struct TestCasePlan {
    std::string  plan_key;  // 注册表键
    TestCaseSpec spec;      // 参数（单工况，后续可扩展为 vector）
};

/**
 * @brief 构建求解器参数
 *
 * 对标 FullCaseTest.cpp:388-413 的 BuildN1Params。
 * 物理容差等为固定值，线性求解器类型等从 cfg 读取以支持外部调整。
 */
FIM_Engine::TransientSolverParams BuildSolverParams(const TestCaseSpec& cfg) {
    FIM_Engine::TransientSolverParams p;

    // 时间步控制
    p.max_steps         = cfg.max_steps;
    p.dt_init           = cfg.dt_init;
    p.dt_min            = cfg.dt_min;
    p.dt_max            = cfg.dt_max;
    p.target_end_time_s = cfg.target_end_time_s;

    // Newton 迭代控制
    p.max_newton_iter = cfg.max_newton_iter;
    p.abs_res_tol     = cfg.abs_res_tol;
    p.rel_res_tol     = cfg.rel_res_tol;
    p.rel_update_tol  = cfg.rel_update_tol;

    // 线性求解器
    p.lin_solver                       = cfg.lin_solver;
    p.amgcl_use_fallback_sparselu      = cfg.amgcl_use_fallback_sparselu;
    p.amgcl_tol                        = cfg.amgcl_tol;
    p.amgcl_maxiter                    = cfg.amgcl_maxiter;

    // 非正交修正与行缩放
    p.enable_non_orthogonal_correction = cfg.enable_non_orthogonal_correction;
    p.enable_row_scaling               = cfg.enable_row_scaling;

    // 状态限幅
    p.max_dP    = cfg.max_dP;
    p.max_dT    = 1.0;       // 温度限幅（N=1 压力模式下影响较小）
    p.max_dSw   = 0.1;       // 饱和度限幅（N=1 模式不参与计算）
    p.min_alpha = 1.0e-8;    // 线搜索最小步长

    // 线搜索
    p.enable_armijo_line_search = cfg.enable_armijo_line_search;

    // 自适应时间步
    p.rollback_shrink_factor   = cfg.rollback_shrink_factor;
    p.dt_relres_grow_factor    = cfg.dt_relres_grow_factor;

    // 重力
    p.gravity_vector = cfg.gravity_vector;

    // 诊断输出
    p.diag_level = cfg.diag_level;

    return p;
}

/**
 * @brief 构建默认测试计划（常物性、无裂缝、无井）
 *
 * TestCaseSpec 使用全部默认值，参数与 FullCaseTest 的 N1 模板一致。
 */
TestCasePlan BuildDefaultPlan() {
    TestCasePlan plan;
    plan.plan_key = "h_co2_constpp_nofrac_nowell";
    // plan.spec 使用 TestCaseSpec 的全部默认值
    return plan;
}

/// 函数指针类型：返回 TestCasePlan 的无参构建器
using BuilderFn = TestCasePlan(*)();

/**
 * @brief 获取测试注册表（单例）
 *
 * 当前仅注册一个工况。后续扩展新工况时，
 * 添加新的 BuildPlan 函数并在此注册即可。
 */
const std::unordered_map<std::string, BuilderFn>& GetRegistry() {
    static const std::unordered_map<std::string, BuilderFn> registry = {
        {"h_co2_constpp_nofrac_nowell", &BuildDefaultPlan}
    };
    return registry;
}

// =====================================================================
// === Layer 3：RunCase — 完整管道执行 ===
// =====================================================================

/**
 * @brief 执行一个完整的 2D 常物性无裂缝无井压力扩散测试
 *
 * 本函数是三层架构的核心执行层，包含从几何构建到后处理导出的全部步骤。
 * 每个步骤仅调用当前工程中已有的模块函数，不引入任何新的辅助函数。
 *
 * @param cfg 测试工况参数
 * @return 运行摘要
 */
TestCaseSummary RunCase(const TestCaseSpec& cfg) {
    TestCaseSummary summary;

    // =================================================================
    // === Step 1: 创建输出目录 ===
    // =================================================================
    // 输出路径结构：{output_base_dir}/{sub_dir}/{case_name}/
    // 默认为 Test/Transient/FullCaseTest/H_CO2_ConstPP/h_co2_constpp_nofrac_nowell/
    summary.case_dir = cfg.output_base_dir + "/" + cfg.sub_dir + "/" + cfg.case_name;
    EnsureDirRecursive(summary.case_dir);
    summary.convergence_log_path = summary.case_dir + "/convergence.log";
    summary.metrics_csv_path     = summary.case_dir + "/metrics.csv";
    summary.nx = cfg.nx;
    summary.ny = cfg.ny;

    // 打开收敛日志文件（每步写入 Newton 迭代信息）
    std::ofstream convergenceLog(summary.convergence_log_path, std::ios::out | std::ios::trunc);
    if (!convergenceLog.good()) {
        throw std::runtime_error("[Test_H_CO2] failed to open convergence log: " + summary.convergence_log_path);
    }

    // =================================================================
    // === Step 2: 构建网格 ===
    // =================================================================
    // MeshManager 构造参数：
    //   lx, ly = 域尺寸; lz=0 表示 2D
    //   nx, ny = 网格划分; nz=0
    //   usePrism=true, useQuadBase=false → 生成四边形网格
    // BuildSolidMatrixGrid_2D 完成：
    //   (1) gmsh 生成网格  (2) 面索引映射  (3) 单元分类
    //   (4) 边界面标记 (LEFT/RIGHT/TOP/BOTTOM)  (5) 面法向量计算（OverRelaxed 修正）
    MeshManager mgr(cfg.lx, cfg.ly, 0.0, cfg.nx, cfg.ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod::OverRelaxed);

    // 注意：本工况无裂缝，跳过 addFracture / DetectAndSubdivideFractures

    // =================================================================
    // === Step 3: 建立全局索引与拓扑 ===
    // =================================================================
    // BuildGlobalSystemIndexing：为所有基岩单元（+ 裂缝微元，如有）分配 solverIndex
    // BuildFracturetoFractureTopology：构建裂缝-裂缝邻接图（无裂缝时为空操作）
    // setNumDOFs(1)：N=1 压力单方程模式，每个 block 1 个自由度
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
    mgr.setNumDOFs(1);

    // =================================================================
    // === Step 4: 初始化场管理器 ===
    // =================================================================
    // FieldManager_2D 持有基岩/裂缝/NNC/面的标量场
    // FIM_CaseKit::InitFieldManager 根据 MeshManager 的 DOF 信息分配场存储
    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    const size_t nCells     = mgr.mesh().getCells().size();  // 基岩单元数
    const int    totalBlocks = mgr.getTotalDOFCount();        // 总 DOF 数（基岩 + 裂缝）
    if (nCells == 0) {
        throw std::runtime_error("[Test_H_CO2] matrix cell count is zero");
    }
    summary.n_cells = static_cast<int>(nCells);

    // =================================================================
    // === Step 5: 设置初始条件 ===
    // =================================================================
    // InitialConditions 结构体传入求解器，用于初始化状态向量
    // N=1 模式下只使用 P_init 和 T_init，Sw_init 被忽略但需设置
    FIM_Engine::InitialConditions ic;
    ic.P_init  = cfg.p_init;   // 初始压力场均匀分布 [Pa]
    ic.T_init  = cfg.t_init;   // 初始温度场均匀分布 [K]
    ic.Sw_init = 1.0;          // 单相饱和度（N=1 时不参与计算）

    // =================================================================
    // === Step 6: 设置边界条件 ===
    // =================================================================
    // BoundaryConditionManager 支持 Dirichlet/Neumann/Robin 三种类型
    // 本工况：左右 Dirichlet 压力驱动，上下 Neumann 零通量（不可渗透壁面）
    // 物理意义：左侧注入高压（12 MPa），右侧保持低压（8 MPa），形成 4 MPa 压差驱动流动
    BoundarySetting::BoundaryConditionManager bcP;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   cfg.p_left);    // 左边界：恒定压力 12 MPa
    bcP.SetDirichletBC(MeshTags::RIGHT,  cfg.p_right);   // 右边界：恒定压力  8 MPa
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);              // 下边界：零通量（不可渗透）
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);              // 上边界：零通量（不可渗透）

    // VTK BC-aware 重建上下文（用于导出 *_bc_point / *_bc_mask）
    const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
    VTKBoundaryVisualizationContext bcVizCtx;
    bcVizCtx.bindings.push_back(
        VTKBCVariableBinding{ pEqCfg.pressure_field, &bcP, VTKBCTransportKind::Pressure });

    // =================================================================
    // === Step 7: 初始化压力场并同步到 FieldManager ===
    // =================================================================
    // pBlocksLatest 是求解器输出的压力向量的本地副本，用于在回调中更新 FM
    // SyncPressureFieldsToFM 将压力值写入 FM 的 P_w / P / T 标量场
    std::vector<double> pBlocksLatest(static_cast<std::size_t>(std::max(totalBlocks, 0)), cfg.p_init);
    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);

    // =================================================================
    // === Step 8: 导出初始 VTK ===
    // =================================================================
    // PostProcess_2D 从 MeshManager 和 FieldManager 读取几何与场数据，导出为 VTK 格式
    // initial.vtk 记录 t=0 时刻的均匀压力场分布
    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/initial.vtk", 0.0);

    // =================================================================
    // === Step 9: 配置求解器模块 ===
    // =================================================================
    // TransientOptionalModules 是 RunGenericFIMTransient 的可选配置容器，
    // 通过其 std::function 字段设置物性初始化和步接受回调。
    //
    // 注意：property_initializer 和 on_step_accepted 是模块接口的 std::function 字段，
    // 属于求解器模块配置的一部分，而非管道中的新辅助函数。
    // 这是当前架构中设置这些回调的唯一方式。

    bool midExported = false;         // 标记 mid.vtk 是否已导出
    int  iterSum = 0, iterCount = 0;  // Newton 迭代累计统计
    int  maxIters = 0;                // 单步最大迭代数

    FIM_Engine::TransientOptionalModules<MeshManager, FieldManager_2D> modules;

    // 边界条件指针（求解器内部通过此指针读取 BC 设置）
    modules.pressure_bc = &bcP;

    // 流体模型：CO2 单相（常物性模式下仅用标签，实际物性由 baseline 值决定）
    modules.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::CO2;

    // 物性模式：ConstantBaseline 表示使用用户指定的恒定 rho/mu，不走 EOS 计算
    modules.pressure_only_property_mode = FIM_Engine::PressureOnlyPropertyMode::ConstantBaseline;
    modules.pressure_only_temperature_k = cfg.t_init;      // EOS 评估温度 [K]
    modules.pressure_only_baseline_rho  = cfg.rho_const;   // 恒定密度 [kg/m^3]
    modules.pressure_only_baseline_mu   = cfg.mu_const;    // 恒定粘度 [Pa*s]

    // 禁用引擎默认的 VTK 输出（由本函数在回调中手动控制三帧导出）
    modules.disable_default_vtk_output = true;

    // --- Step 9a: 物性初始化回调 ---
    // 求解器在组装传导率之前调用此回调，将岩石/裂缝/流体的静态物性写入 FieldManager。
    // 各字段名标签由 SolverContrlStrName_op.h 中的命名约定提供。
    modules.property_initializer = [&cfg](MeshManager&, FieldManager_2D& fld) {
        const auto rock  = PhysicalProperties_string_op::Rock();
        const auto frac  = PhysicalProperties_string_op::Fracture_string();
        const auto water = PhysicalProperties_string_op::Water();

        // 基岩渗透率张量（各向同性）: k_xx = k_yy = k_zz = matrix_perm
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_xx_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_yy_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.k_zz_tag, cfg.matrix_perm), cfg.matrix_perm);

        // 基岩孔隙度、压缩系数、骨架密度、比热、热导率
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.phi_tag,    cfg.matrix_phi), cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.c_r_tag,    cfg.matrix_ct),  cfg.matrix_ct);
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.rho_tag,    2600.0), 2600.0);   // 骨架密度 [kg/m^3]
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.cp_tag,     1000.0), 1000.0);   // 比热 [J/(kg*K)]
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(rock.lambda_tag, 2.0),    2.0);      // 热导率 [W/(m*K)]

        // 裂缝物性（无裂缝时 size=0，ApplyUniformScalarField 内部循环不执行）
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_t_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.k_n_tag, cfg.matrix_perm), cfg.matrix_perm);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.phi_tag, cfg.matrix_phi),  cfg.matrix_phi);
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(frac.c_r_tag, cfg.matrix_ct),   cfg.matrix_ct);

        // 水的热导率（基岩 + 裂缝）
        ApplyUniformScalarField(fld.getOrCreateMatrixScalar(water.k_tag,   0.6), 0.6);   // [W/(m*K)]
        ApplyUniformScalarField(fld.getOrCreateFractureScalar(water.k_tag, 0.6), 0.6);
    };

    // --- Step 9b: 时间步接受回调 ---
    // 每个时间步收敛后，求解器调用此回调。用途：
    //   (1) 将最新压力向量同步回 FieldManager（供 VTK 导出使用）
    //   (2) 更新运行摘要统计量
    //   (3) 写入收敛日志
    //   (4) 在模拟时间过半时导出 mid.vtk
    modules.on_step_accepted =
        [&](int step, double timeS, double dtUsedS, int newtonIters, double residualInf,
            int totalRollbacks, const std::string& convergeMode,
            const std::vector<double>& pVec, const std::vector<double>&,
            const std::vector<double>*) {
            if (step <= 0) return;

            // (1) 同步压力场到 FieldManager
            if (!pVec.empty()) {
                pBlocksLatest = pVec;
                SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);
            }

            // (2) 更新运行摘要
            summary.steps           = step;
            summary.t_end           = timeS;
            summary.total_rollbacks = totalRollbacks;
            iterSum  += newtonIters;
            iterCount += 1;
            maxIters  = std::max(maxIters, newtonIters);

            // (3) 写入收敛日志
            convergenceLog << "[Step " << step << "] t=" << std::scientific << std::setprecision(8) << timeS
                           << " dt=" << dtUsedS
                           << " iters=" << newtonIters
                           << " residual_inf=" << residualInf
                           << " rollbacks=" << totalRollbacks
                           << " mode=" << convergeMode << "\n";

            // (4) 在模拟过半时导出 mid.vtk（仅一次）
            if (!midExported && timeS >= 0.5 * cfg.target_end_time_s - 1.0e-12) {
                PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", timeS);
                midExported = true;
            }
        };

    // =================================================================
    // === Step 10: 运行瞬态求解器 ===
    // =================================================================
    // RunGenericFIMTransient<1> 是通用全隐 FIM 求解器模板，N=1 表示压力单方程模式。
    // 内部自动完成：
    //   (1) 静态物性注入（调用 property_initializer 回调）
    //   (2) MM/FI/NNC/FF 传导率计算与连接聚合
    //   (3) 块稀疏 Jacobian 构建与 pattern 冻结
    //   (4) Newton 非线性迭代（AD 自动微分组装 Jacobian）
    //   (5) 线性求解（AMGCL + SparseLU 回退）
    //   (6) 自适应时间步控制（dt 增长/回退）
    //   (7) 每步收敛后调用 on_step_accepted 回调
    const auto params = BuildSolverParams(cfg);
    FIM_Engine::RunGenericFIMTransient<1>(
        cfg.case_name,
        mgr,
        fm,
        ic,
        {},                              // 空井列表（无注采井）
        params,
        FIM_Engine::SolverRoute::FIM,    // 使用 FIM 路线
        modules);

    // =================================================================
    // === Step 11: 后处理导出 ===
    // =================================================================
    // 同步最终压力场到 FieldManager，然后导出 VTK
    SyncPressureFieldsToFM(mgr, fm, pBlocksLatest, cfg.t_init);

    // 如果模拟提前终止（未到达 50% 时间），补导 mid.vtk
    if (!midExported) {
        PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/mid.vtk", summary.t_end);
    }

    // 导出最终时刻 VTK
    PostProcess_2D(mgr, fm, &bcVizCtx).ExportVTK(summary.case_dir + "/final.vtk", summary.t_end);

    // =================================================================
    // === Step 12: 计算指标并写入 metrics.csv ===
    // =================================================================
    summary.max_iters = maxIters;
    summary.avg_iters = (iterCount > 0) ? (static_cast<double>(iterSum) / static_cast<double>(iterCount)) : 0.0;
    summary.h_char    = ComputeMeshCharLength(mgr);

    // 输出指标到 CSV 文件
    std::ofstream metrics(summary.metrics_csv_path, std::ios::out | std::ios::trunc);
    metrics << "case_name,nx,ny,n_cells,h_char,t_end,steps,total_rollbacks,"
               "avg_nonlinear_iters,max_nonlinear_iters,"
               "dt_init,dt_min,dt_max,property_mode,solver_route\n";
    metrics << cfg.case_name << ","
            << cfg.nx << ","
            << cfg.ny << ","
            << nCells << ","
            << std::setprecision(12) << summary.h_char << ","
            << summary.t_end << ","
            << summary.steps << ","
            << summary.total_rollbacks << ","
            << summary.avg_iters << ","
            << summary.max_iters << ","
            << cfg.dt_init << ","
            << cfg.dt_min << ","
            << cfg.dt_max << ","
            << "ConstantBaseline" << ","
            << "RunGenericFIMTransient<1>"
            << "\n";

    return summary;
}

// =====================================================================
// === Layer 1：名称索引与调度 ===
// =====================================================================

/**
 * @brief 通过注册表键名调度执行测试
 *
 * 从注册表中查找构建器函数，构建参数计划，执行 RunCase，输出摘要。
 */
void ExecutePlanByKey(const std::string& key) {
    const auto& registry = GetRegistry();
    auto it = registry.find(key);
    if (it == registry.end()) {
        throw std::runtime_error("[Test_H_CO2] unknown registry key: " + key);
    }

    const TestCasePlan plan = it->second();
    const TestCaseSummary summary = RunCase(plan.spec);

    // 输出运行摘要到控制台
        std::cout << "\n============================================\n";
    std::cout << "[Test_H_CO2] run completed\n";
    std::cout << "  output_dir: " << summary.case_dir << "\n";
    std::cout << "  grid: " << summary.nx << " x " << summary.ny
              << " (" << summary.n_cells << " cells)\n";
    std::cout << "  steps: " << summary.steps
              << "  rollbacks: " << summary.total_rollbacks << "\n";
    std::cout << "  Newton iters: avg=" << std::fixed << std::setprecision(2) << summary.avg_iters
              << "  max=" << summary.max_iters << "\n";
    std::cout << "  final_time: " << std::scientific << std::setprecision(4) << summary.t_end << " s\n";
    std::cout << "============================================\n";
}

} // 匿名命名空间结束

// =====================================================================
// === 公开入口 ===
// =====================================================================

void RunTestCase() {
    ExecutePlanByKey("h_co2_constpp_nofrac_nowell");
}

} // namespace Test_H_CO2_ConstPP

