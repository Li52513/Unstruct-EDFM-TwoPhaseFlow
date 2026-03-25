# PROJECT_CONTEXT

## 0. Metadata

| Item | Value |
|---|---|
| Last Updated | 2026-03-24 (Asia/Shanghai) |
| Source of Truth | 当前工作区代码（`main.cpp` + `FullCaseTest.cpp` + `FIM_TransientEngine/*` + 核心管理器） |
| Git State | dirty working tree（存在用户未提交改动） |
| Project Root | `2D-Unstr-Quadrilateral-EDFM` |
| Entry | `main.cpp`（`--list` / `--case` / `--help`） |
| Default Case | `day2_fvm_ad` |

> 规则：文档与实现冲突时，以 `main.cpp` dispatcher 与函数实现为准。

---

## 1. 工程目标与当前状态

### 1.1 总体目标

1. 构建基于非结构网格 FVM 的 2D/3D EDFM 框架（矩阵-裂缝嵌入、NNC/FF 拓扑、传导率、全隐求解）。
2. 支持单相与两相、常物性与变物性（CO2 EOS-AD）瞬态模拟。
3. 通过可复用测试矩阵（FullCase/Day6/Validation）做稳定性、收敛性、鲁棒性验证。

### 1.2 当前落地程度（关键）

1. 2D/3D EDFM 几何链已可用，含裂缝求交、拓扑清理、映射构建。
2. 通用瞬态入口 `RunGenericFIMTransient<N>` 已工程化，N=1/N=2/N=3 均可走 AD-Jacobian 非线性路径。
3. `FullCaseTest.cpp/.h` 已成为统一测试承载核心（N=1 模板 + 迁移后的 T1..T8 验证入口 + 2D/3D dispatcher）。
4. 新增 N=1 变物性单裂缝入口：`day6l4_2d_sp_co2_varprop_nowell_singlefrac`。

---

## 2. 代码架构总览

```text
2D-Unstr-Quadrilateral-EDFM/
├─ main.cpp                                  # 全局 dispatcher（真实 case 清单）
├─ FullCaseTest.cpp/.h                       # 统一测试承载（N=1模板 + T1..T8迁移入口）
├─ Test_Day6_TransientSolver.cpp/.h          # Day6 campaign 与 3D 单场景 helper
├─ MeshManager.*                             # 2D: matrix + 1D fractures
├─ 3D_MeshManager.*                          # 3D: matrix + 2D fractures
├─ FractureNetwork.* / 2D_FractureNetwork.*  # 裂缝网络与 FF/NNC 拓扑容器
├─ 2D_FieldManager.h / 3D_FieldManager.h     # Matrix/Frac/NNC/Face(Edge) 场
├─ TransmissibilitySolver_2D.* / _3D.*       # MM/FI/NNC/FF 系数计算
├─ FIM_TopologyBuilder2D.h / 3D.h            # 从 T_* 字段装载 Connection
├─ FIM_ConnectionManager.h                   # 连接聚合、去重、一致性校验
├─ FIM_BlockSparseMatrix.h                   # 块稀疏 Jacobian，支持冻结 pattern
├─ FIM_TransientEngine/
│  ├─ Types.hpp
│  ├─ StateSync.hpp
│  ├─ Diagnostics.hpp
│  ├─ StepKernels.hpp/_impl.hpp
│  └─ RunGeneric.hpp/_impl.hpp               # 通用全隐 FIM 主循环
├─ BoundaryAssembler.h                       # 边界/井源 FullJac 装配
├─ WellDOFManager.h                          # 井独立 BHP DOF 管理
└─ 2D_PostProcess.* / 3D_PostProcess.*       # VTK/VTU/Tecplot 导出
```

---

## 3. Dispatcher 与 Case 组织

## 3.1 参数行为

1. `--list`：列出可运行 case。
2. `--case=<name>` 或 `--case <name>`：运行指定 case。
3. 不带 `--case`：默认 `day2_fvm_ad`。

## 3.2 当前关键 case 组（按用途）

1. Day1-Day6 与回归：`day1_*`, `day2_*`, `day3_*`, `day4_*`, `day5_*`, `day6_*`, `issue11_*`, `issue12_*`。
2. N=1 梯度调试链（FullCase）：
   - `day6l1_2d_sp_co2_const_nowell_analytical`
   - `day6l2_2d_sp_co2_const_nowell_grid`
   - `day6l3_2d_sp_co2_const_nowell_solver`
   - `day6l4_2d_sp_co2_varprop_nowell`
   - `day6l4_2d_sp_co2_varprop_nowell_singlefrac`（新增）
   - `day6ladder_2d_sp_co2_const_nowell_all`（L1-L4）
3. N=1 模板：
   - `full_n1_template_const_nowell_nofrac`
   - `full_n1_template_const_nowell_singlefrac`
   - `full_n1_template_const_nowell_crossfrac`
4. Validation 迁移入口（T1..T8 × const/var × topology）：`val_t*_...`。
5. FullCase 批量入口：`2d_t*_nofrac|singlefrac|crossfrac`, `3d_t*_...`, `2d_all`, `3d_all`, `all`。

---

## 4. 关键流程（几何到求解）

## 4.1 2D EDFM 主链

1. `BuildSolidMatrixGrid_2D`
2. `addFracture`
3. `DetectAndSubdivideFractures(...)`
4. `BuildGlobalSystemIndexing()`
5. `BuildFracturetoFractureTopology()`
6. `TransmissibilitySolver_2D::{Matrix,FractureInternal,NNC,FF}`
7. `FIM_TopologyBuilder2D::LoadAllConnections`
8. `FIM_ConnectionManager::FinalizeAndAggregate`

## 4.2 3D EDFM 主链

1. `BuildSolidMatrixGrid_3D`
2. `addFracturetoFractureNetwork`
3. `meshAllFracturesinNetwork`
4. `setupGlobalIndices`
5. `DetectFractureFractureIntersectionsInNetwork`
6. `SolveIntersection3D_improved_twist_accleration`
7. `rebuildEdgeProperties`
8. `removeDuplicateInteractions`
9. `resolveCoplanarInteractions`
10. `buildTopologyMaps`
11. `TransmissibilitySolver_3D::{Matrix,FractureInternal,NNC,FF}`
12. `FIM_TopologyBuilder3D::LoadAllConnections`
13. `FIM_ConnectionManager::FinalizeAndAggregate`

## 4.3 通用瞬态 FIM 主链（`RunGenericFIMTransient<N>`）

1. 注入静态物性与初值。
2. 计算 MM/FI/NNC/FF 连接并聚合。
3. 构建/冻结块稀疏结构。
4. Newton 非线性迭代（AD Jacobian）。
5. 线性求解（AMGCL/CPR/SparseLU fallback 等）。
6. 时间步控制（dt 自适应、rollback、可选线搜索/PTC）。
7. 同步回 FieldManager 并导出 VTK。

---

## 5. FullCaseTest 使用指南（核心）

## 5.1 文件角色

1. `FullCaseTest.cpp/.h`：测试模板、计划、注册表、统一执行入口。
2. `main.cpp`：仅做 case 字符串到 `FullCaseTest` 函数的分发。

## 5.2 N=1 模板架构（推荐二次开发入口）

1. 参数结构：`N1CaseSpec`（网格、裂缝几何、基岩/裂缝物性、边界、时间步、是否 EOS）。
2. 执行器：`RunN1Common(spec)`。
3. case 计划：`BuildPlan...()` 返回 `N1CasePlan`。
4. 注册：`GetN1Registry()` 中 `plan_key -> BuildPlanFn`。
5. 暴露：`RunN1...()` 包装 + `main.cpp` case 分发。

## 5.3 常见“改哪个参数”速查（N=1）

1. 基岩几何与网格：`N1CaseSpec.lx/ly/nx/ny`。
2. 裂缝几何：`frac_x_ratio/frac_y0_ratio/frac_y1_ratio`，具体生成为 `AddTemplateFractures(...)`。
3. 基岩物性：`matrix_phi/matrix_perm/matrix_ct`。
4. 裂缝物性：`fracture_phi/fracture_kt/fracture_kn/fracture_ct`。
5. 变物性开关：`use_co2_eos=true`（走 `PressureOnlyPropertyMode::CO2_EOS`）。
6. 边界压力：`p_left/p_right/p_init`。
7. 求解器与时间步：`dt_init/dt_min/dt_max/max_newton_iter/max_steps`。
8. 观测指标：`enable_barrier_probe=true` 可输出裂缝跨压差指标。

## 5.4 N=2/N=3（T1..T8）修改入口

1. 每个 `RunT*` 内定义本场景的 `Lx, nx, k, phi, t_final` 等。
2. 单相预设：`BuildPreset_SP(...)`。
3. 两相预设：`BuildPreset_TP(...)`。
4. 求解器参数：`BuildValParams(twoPhase, dt_init, dt_max, t_final, max_steps)`。
5. 场景组合包装：`Val_T*_Const/Var_*`。
6. key 到函数映射：`GetLegacyValCaseRegistry()`。

## 5.5 新建“不同工况”推荐步骤

1. 复制最接近的 `BuildPlan...`（例如从 `BuildPlanDay6L4SingleFrac` 复制）。
2. 修改 `plan_key` 与 `case_name`（避免覆盖旧输出目录）。
3. 改 `N1CaseSpec` 参数（几何、物性、求解器）。
4. 在 `GetN1Registry()` 注册新 key。
5. 新增对外 `RunN1...()` 包装。
6. 在 `main.cpp` 添加 `--case` 分发。

> 结论：可以直接复制已有测试 case 进行二次修改，这是当前工程推荐模式。

## 5.6 参数修改速查表（一页版）

| 修改目标 | N=1（FullCase-N1）改动点 | N=2/N=3（Validation T1..T8）改动点 | 影响/注意事项 |
|---|---|---|---|
| 计算域尺寸 | `N1CaseSpec.lx/ly` | 各 `RunT*` 中 `Lx/Ly` | 会改变扩散尺度与解析比对基准 |
| 网格密度 | `N1CaseSpec.nx/ny` | 各 `RunT*` 中 `nx/ny` | 影响收敛阶、线性系统规模、CPU 时间 |
| 裂缝拓扑类型 | `N1CaseSpec.topology` | `RunT*` 调用 `FracConfig::{None,Single,Cross}` | 直接改变 NNC/FF 连接数量 |
| 单裂缝位置 | `frac_x_ratio/frac_y0_ratio/frac_y1_ratio` | 修改 `AddSingleFrac(...)` 几何定义 | 影响压降路径与屏障/导流效应 |
| 交叉裂缝位置 | `AddTemplateFractures(...)` 中 Cross 分支 | 修改 `AddCrossFrac(...)` | 会改变 FF 交汇点位置与强度 |
| 基岩孔隙度 | `matrix_phi` | `phi` + `BuildPreset_SP/TP` | 影响蓄积项、推进速度 |
| 基岩渗透率 | `matrix_perm` | `k` + `BuildPreset_SP/TP` | 影响流动速度与稳定步长 |
| 基岩压缩性 | `matrix_ct` | `c_r` + `BuildPreset_SP/TP` | 影响弱可压扩散系数 |
| 裂缝切向渗透率 | `fracture_kt` | `BuildPreset_SP/TP` 中 `p.frac.permeability` | 高 `k_t` 导流增强，低 `k_t` 接近屏障 |
| 裂缝法向渗透率 | `fracture_kn` | `BuildPreset_SP/TP` 映射到裂缝法向项 | 低 `k_n` 会显著增大跨缝压差 |
| 裂缝开度 | 当前固定在属性注入处（可扩展到 `N1CaseSpec`） | `BuildPreset_SP/TP` 中 `p.frac.aperture` | 影响 NNC 阻力与 FI/FF 强度 |
| 初始压力温度 | `p_init/t_init` | 各 `RunT*` 的 `P0/T0` | 改变初态残差与 Newton 起点 |
| 边界压力 | `p_left/p_right`（`RunN1Common` 中设置） | 各 `RunT*` 的 `bcP` | 改变驱动压差与稳态梯度 |
| 是否变物性 | `use_co2_eos=true/false` | `FluidConfig::VarWater` 或 TP EOS 路径 | N=1 变物性走 `CO2_EOS` 压力单方程模式 |
| 时间步初值/上下限 | `dt_init/dt_min/dt_max` | `BuildValParams(...)` 参数 | 直接影响 rollback 频率与效率 |
| Newton 迭代上限 | `max_newton_iter` | `BuildValParams(...)` 或 `RunT*` 覆写 | 过小可能不收敛，过大可能慢 |
| 线性求解器 | `BuildN1Params` (`AMGCL` + `SparseLU` 回退) | `BuildValParams` (`AMGCL/AMGCL_CPR`) | 影响稳健性与性能 |
| 非正交修正 | `BuildN1Params.enable_non_orthogonal_correction` | `BuildValParams.enable_non_orthogonal_correction` | 建议保持开启，避免斜网格偏差 |
| 输出目录与case名 | `plan_key` + `case_name` + `level_dir` | `cfg.tag`（`RunT*`） | 避免覆盖旧结果，方便批量对比 |
| 解析解验收阈值 | `analytical_check/l2_threshold` | 各 `RunT*` 的 analytical compare 段 | 仅在有解析解场景启用 |
| 屏障定量指标 | `enable_barrier_probe=true` | 需自定义指标逻辑 | 输出 `delta_p_cross_fracture` 到 `metrics.csv` |

> 建议：每次只改 1-2 个参数并保留 `case_name` 差异化命名，便于回归定位。

---

## 6. 核心类与数据结构

## 6.1 网格与裂缝

1. `Mesh`：`nodes_/faces_/cells_`，提供 2D/3D 空间索引。
2. `MeshManager`（2D）：管理基岩单元与 1D 裂缝微元，统一 solverIndex。
3. `MeshManager_3D`：管理 3D 基岩与 2D 裂缝微面，持有 `interactionPairs_` 与双向映射。
4. `FractureNetwork / FractureNetwork_2D`：裂缝微元、交线、边、FF 拓扑容器。

## 6.2 场管理

1. `FieldManager_2D/3D`：matrix/fracture/nnc/face(edge) 场。
2. `FieldRegistry`：字符串标签到场对象映射，支持 create/get/getOrCreate。
3. `SolverContrlStrName_op.h`：字段名约定来源。

## 6.3 代数连接与装配

1. `FIM_TopologyBuilder2D/3D`：从 `T_*` 场加载连接。
2. `FIM_ConnectionManager`：连接聚合、重复校验、顺序化。
3. `FIM_BlockSparseMatrix<N>`：块稀疏 Jacobian 容器，支持 pattern freeze。
4. `BoundaryAssembler`：边界与井源项 FullJac 装配。
5. `WellDOFManager`：井 BHP/Rate 方程耦合。

## 6.4 状态与 AD

1. `ADVar<N>`：自动微分核心类型。
2. `FIM_StateMap<N>`：仅持有 `P/T/(Sw)` 数值状态。
3. `StateSync.hpp`：状态与 FieldManager 同步，含静态物性注入。

---

## 7. 理论公式与代码映射

| 对象 | 公式/思想 | 代码实现 |
|---|---|---|
| 两点串联阻力 | `R = d1/K1 + d2/K2` | `FVM_Ops::Op_Math_SeriesResistance` |
| 传导率 | `T = A / R` | `FVM_Ops::Op_Math_Transmissibility` |
| MM 系数 | 面法向距离 + 投影渗透率 | `TransmissibilitySolver_2D/3D::Calculate_Transmissibility_Matrix` |
| FI 系数 | 裂缝微元串联阻力 | `Calculate_Transmissibility_FractureInternal` |
| NNC 系数 | `A_int / (d_m/K_m + w_f/(2K_f))` | `Calculate_Transmissibility_NNC` |
| FF 系数 | 交汇点 Star-Delta | `Calculate_Transmissibility_FF` |
| 势差 | `ΔΦ = ΔP - ρ g·Δx` | `FVM_Ops::Compute_Potential_Diff` |
| 上游选择 | 按 `ΔΦ` 符号 | `FVM_Ops::Op_Upwind_AD` |
| 非正交修正 | `corr = ∇φ_f · vectorT` | `FVM_Ops::Op_NonOrthogonal_PressureCorr` |
| 两相本构 | VG + Mualem | `CapRelPerm_HD.h` |
| AD-EOS 桥接 | EOS + 鲁棒导数链式传播 | `AD_FluidEvaluator.h` |

---

## 8. 输出规范与目录

1. N=1 FullCase 强制三帧：`initial.vtk`, `mid.vtk`, `final.vtk`。
2. 运行日志：`convergence.log`。
3. 指标：`metrics.csv`（含 `property_mode/solver_route/grad_route/bc_geom_route` 审计字段）。
4. N=1 输出根目录：`Test/Transient/FullCaseTest/N1/<level>/<case_name>/`。
5. Validation(T1..T8) 输出根目录：`Test/Transient/ValidationSuite/<tag>/`。

---

## 9. 当前已知约束与注意事项

1. N=1 路径当前硬约束为无井源调试链（case 计划层面未扩展井控）。
2. `Run3DCase` 通过 FullCase dispatcher 调用 Day6 3D 单场景 helper，保持 CLI 一致。
3. 几何专用 benchmark（如 mesh/intersection）不一定计算 MM/FI/NNC/FF，属于设计目的，不等同瞬态求解链。
4. `test_Transmissibility_2D.h` 已修复索引/拓扑顺序：先 `BuildGlobalSystemIndexing()` 再 `BuildFracturetoFractureTopology()`。
5. `FullCaseTest` 纯度扫描已避免“自命中”误报（禁用 token 动态拼接）。

---

## 10. 构建与运行（示例）

```powershell
# 枚举 case
& $exe --list

# N=1 变物性无井单裂缝（新增）
& $exe --case=day6l4_2d_sp_co2_varprop_nowell_singlefrac

# N=1 模板单裂缝常物性
& $exe --case=full_n1_template_const_nowell_singlefrac

# Validation 单场景
& $exe --case=val_t1_var_singlefrac

# FullCase 批跑
& $exe --case=2d_all
& $exe --case=3d_all
```

---

## 11. 新会话快速接手阅读顺序

1. `main.cpp`
2. `FullCaseTest.cpp`
3. `FIM_TransientEngine/RunGeneric_impl.hpp`
4. `FIM_TransientEngine/Types.hpp`
5. `TransmissibilitySolver_2D.cpp` / `TransmissibilitySolver_3D.cpp`
6. `FIM_TopologyBuilder2D.h` / `FIM_TopologyBuilder3D.h`
7. `FIM_ConnectionManager.h`
8. `MeshManager.h` / `3D_MeshManager.h`

