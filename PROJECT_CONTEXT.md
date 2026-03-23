# PROJECT_CONTEXT

## 0. Metadata

| Item | Value |
|---|---|
| Last Updated | 2026-03-23 (Asia/Shanghai) |
| Source of Truth | 当前工作区代码（`main.cpp` + `FIM_TransientEngine/*` + 核心头文件） |
| Git State | `dirty working tree`（含 Day6/FIM 相关未提交改动） |
| Project Root | `2D-Unstr-Quadrilateral-EDFM` |
| Entry | `main.cpp`（支持 `--list` / `--case` / `--help`） |
| Default Case | `day2_fvm_ad` |

> 说明：若文档与代码不一致，以 `main.cpp` dispatcher 和实现代码为准。

---

## 1. 工程目标与论文映射

### 1.1 章节一：非结构网格 EDFM 快速嵌入与耦合（主链路已形成）

- 2D：2D 基岩 + 1D 裂缝（DFN），支持 AABB/Bin、DDA、8DOP 候选筛选与精确求交。
- 3D：3D 基岩 + 2D 裂缝独立网格，支持 Octree/14-DOP/光栅化预筛，重构基岩-裂缝交互多边形。
- NNC 几何量（交互面积、法向、距离、拓扑）可直接用于传导率和方程装配。

对应代码：

- 2D：`MeshManager.*`, `mesh.*`, `FractureNetwork.*`
- 3D：`3D_MeshManager.*`, `2D_FractureNetwork.*`, `EDFM_Geometry_3D.*`, `FaceIndexedOctree.h`

### 1.2 章节二：单相热-流耦合（框架齐备，持续增强）

- 已联通场管理、物性、边界、静态传导率（MM/FI/NNC/FF）。
- 已具备统一瞬态入口骨架，可承载单相与两相扩展。

对应代码：

- 场管理：`2D_FieldManager.h`, `3D_FieldManager.h`, `FieldRegistry.h`
- 物性：`2D_PhysicalPropertiesManager.*`, `3D_PhysicalPropertiesManager.*`
- 传导率：`TransmissibilitySolver_2D.*`, `TransmissibilitySolver_3D.*`

### 1.3 章节三：CO2-水两相渗流-传热（FIM 主干已工程化）

- 两相本构：VG 毛管压 + Mualem 相对渗透率。
- AD 物性桥接：`AD_FluidEvaluator.h`（表格查值/鲁棒差分，支持 `USE_COOLPROP_EOS`）。
- FIM 链路：连接聚合、块稀疏冻结、Newton+线搜索+自适应 dt、井独立 DOF、矩阵审计。

对应代码：

- 两相本构：`CapRelPerm_HD.h`, `CapRelPerm_HD_AD.h`
- AD/FVM 算子：`ADVar.hpp`, `FVM_Ops.h`, `FVM_Ops_AD.h`
- FIM：`FIM_ConnectionManager.h`, `FIM_TopologyBuilder2D.h`, `FIM_TopologyBuilder3D.h`, `FIM_BlockSparseMatrix.h`, `WellDOFManager.h`, `FIM_TransientEngine/*`

### 1.4 章节四：WAG/闭环策略优化（待系统化）

- 已具备参数化研究基础（几何-物性-求解）。
- 待补齐统一闭环 case 编排、指标统计与自动化报告链路。

---

## 2. 代码架构总览（当前仓库）

```text
2D-Unstr-Quadrilateral-EDFM/
├─ main.cpp                                  # 统一 dispatcher（真实 case 源）
├─ mesh.*, Node.*, Face.*, Cell.*            # 基岩网格、几何拓扑、Bin 索引
├─ MeshManager.*                             # 2D EDFM 管理器（2D matrix + 1D fracture）
├─ 3D_MeshManager.*                          # 3D EDFM 管理器（3D matrix + 2D fracture）
├─ FractureNetwork.*                         # 2D 线裂缝网络（globalFFPts 等）
├─ 2D_FractureNetwork.*                      # 3D 场景下的 2D 裂缝网格网络（ffIntersections/globalEdges）
├─ 2D_FieldManager.h / 3D_FieldManager.h     # Matrix/Fracture/NNC 三域场 + face/edge 场
├─ SolverContrlStrName_op.h                  # 字段名约定
├─ TransmissibilitySolver_2D.* / _3D.*       # MM/FI/NNC/FF 静态传导率
├─ FIM_ConnectionManager.h                   # 连接聚合、去重、统计
├─ FIM_TopologyBuilder2D.h / 3D.h            # 从字段装载代数连接
├─ FIM_BlockSparseMatrix.h                   # 块稀疏矩阵 + Pattern Freeze + CSR cache
├─ FIM_GlobalAssembler.h                     # Acc/Flux/Source 装配接口
├─ BoundaryAssembler.h                       # 边界与井源 FullJac 装配
├─ WellDOFManager.h                          # 井独立 BHP DOF
├─ FIM_StateMap.h                            # 轻量状态容器（P/T/(Sw) 标量）
├─ FIM_TransientEngine/
│  ├─ Types.hpp                              # 求解参数与求解路线
│  ├─ StateSync.hpp                          # state↔field 同步、静态物性注入
│  ├─ StepKernels.hpp/_impl.hpp              # 线性求解、行缩放、trial update
│  ├─ Diagnostics.hpp                        # 预检查、matrix audit
│  └─ RunGeneric.hpp/_impl.hpp               # 通用瞬态 FIM 驱动
├─ 2D_PostProcess.* / 3D_PostProcess.*       # Tecplot/VTK/VTU/PVD 导出
└─ Test_*.h / *_test.h                       # Day1-Day6 与模块测试
```

---

## 3. `main.cpp` dispatcher（真实入口）

### 3.1 参数行为

- `--list`：列出可运行 case。
- `--case=<name>` / `--case <name>`：运行指定 case。
- 无 `--case`：默认 `day2_fvm_ad`。

### 3.2 当前已注册 case（2026-03-23）

- Day1：`day1_arch_conn`, `day1_arch_conn_repro`
- Day2：`day2_fvm_ad`
- Day3：`day3_bc_patch`, `day3_leakoff_switch`, `day3_viz`
- Day4：`day4_well_patch`, `day4_well_viz`
- Day5：`day5_block_matrix_robust`, `day5_global_jac_2d`, `day5_global_jac_3d`
- Issue 回归：`issue11_frozen_matrix`, `issue12_linear_solver`
- Day6 瞬态与 campaign：
  - 常规：`day6_transient_2d_sp_injprod`, `day6_transient_2d_tp_injprod`, `day6_transient_2d_tp_multiwell`
  - 常规：`day6_transient_3d_sp_injprod`, `day6_transient_3d_tp_injprod`, `day6_transient_3d_tp_multiwell`
  - 审计：`day6_matrix_audit_2d_edfm`, `day6_matrix_audit_3d_edfm`
  - Campaign：`day6_campaign_2d_all`, `day6_campaign_3d_all`
  - T01 快速门：`day6_t01_f0`, `day6_t01_f1`, `day6_t01_f2`
  - 解析校核：`day6_t1_2d_sp_nowell_analytical`
- 几何/物性/算子基准：
  - `2d_edfm_single`, `2d_edfm_dfn`, `2d_geom_benchmark_dfn`
  - `3d_twist_intersection`, `3d_edfm_improved`, `3d_distance_accuracy`, `3d_nnc_ff_static`, `3d_mesh_benchmark`
  - `init_single_frac`, `init_dfn`, `property_sweep`, `boundary_export`
  - `grad_2d`, `grad_3d`, `grad_all`
  - `dof_mapper`, `advar`, `prop_init_3d`, `fluid_eval`, `complex_frac_benchmark`
  - `trans_2d`, `trans_3d`, `trans_all`

### 3.3 文档一致性提醒

- 当前 dispatcher **未注册** Day7 `day7_closedloop_*` case。
- 若其他文档出现 Day7 case，应标注为计划项，不可当作现有可执行入口。

---

## 4. 关键流程（几何到求解）

### 4.1 2D EDFM 主链

1. `MeshManager::BuildSolidMatrixGrid_2D`
2. 加裂缝（`addFracture`/DFN）
3. `DetectAndSubdivideFractures(...)`
4. `BuildGlobalSystemIndexing()`
5. `BuildFracturetoFractureTopology()`
6. `TransmissibilitySolver_2D::{Matrix,FI,NNC,FF}`
7. `FIM_TopologyBuilder2D::LoadAllConnections`
8. `FIM_ConnectionManager::FinalizeAndAggregate`

### 4.2 3D EDFM 主链

1. `MeshManager_3D::BuildSolidMatrixGrid_3D`
2. `addFracturetoFractureNetwork`
3. `meshAllFracturesinNetwork`
4. `setupGlobalIndices`
5. `DetectFractureFractureIntersectionsInNetwork`
6. `SolveIntersection3D_improved_twist_accleration`
7. `removeDuplicateInteractions`
8. `resolveCoplanarInteractions`
9. `buildTopologyMaps`
10. `TransmissibilitySolver_3D::{Matrix,FI,NNC,FF}`
11. `FIM_TopologyBuilder3D::LoadAllConnections`
12. `FIM_ConnectionManager::FinalizeAndAggregate`

### 4.3 通用瞬态 FIM 主链（`RunGenericFIMTransient`）

1. `InjectStaticProperties` + 可选 `property_initializer`
2. 构建 MM/FI/NNC/FF 连接并聚合
3. `WellDOFManager::Setup` 扩展井 block
4. `FIM_BlockSparseMatrix` 注册 pattern 后 `FreezePattern`
5. Newton + 时间推进：
   - 累积/通量/边界/井项装配
   - 可选非正交修正（`enable_non_orthogonal_correction`）
   - 可选矩阵审计（`enable_matrix_audit`）
   - 线性求解（`SparseLU`/`BiCGSTAB`/`AMGCL`/`AMGCL_CPR`）
   - 线搜索、PTC、dt 自适应
6. `SyncStateToFieldManager` + VTK/VTU 导出

---

## 5. 核心类与数据结构

### 5.1 网格与裂缝

- `Mesh`
  - 持有 `nodes_ / faces_ / cells_`
  - 2D face-bin（`faceBins_`）与 3D cell-bin（`cellBins_`）空间索引
  - 2D NNC 映射：`CellLocalIndexToFracElemSolverIndexMap_`

- `MeshManager`（2D）
  - 2D matrix + 1D fracture 管理
  - 统一 `solverIndex`（matrix 在前，fracture 在后）
  - `getEquationIndex(solverIndex, dofOffset)`

- `MeshManager_3D`
  - 3D matrix + 2D fracture 管理
  - `interactionPairs_` 为核心几何交互容器
  - `mat2InteractionMap_` / `frac2InteractionMap_` 双向拓扑加速

- `FractureNetwork` / `FractureNetwork_2D`
  - 维护 FF 拓扑、全局边、solver index cache（`getOrderedFractureElements`）

### 5.2 场管理

- `FieldRegistry` / `FaceFieldRegistry`
  - 名称到场对象映射，支持 `create/get/getOrCreate/has`

- `FieldManager_2D` / `FieldManager_3D`
  - `matrixFields`, `fractureFields`, `nncFields`
  - `matrixFaceFields`, `fractureFaceFields`(2D) / `fractureEdgeFields`(3D)
  - `ff_topology`（FF solver index 对）

### 5.3 状态与 AD

- `ADVar<N>`：前向自动微分基础类型。
- `volScalarField`, `volADField<N>`, `faceScalarField`, `faceADField<N>`。
- `FIM_StateMap<N>`：仅保存 `P/T/(Sw)` 标量，不做全局驻留 AD。

### 5.4 连接与全局装配

- `FIM_TopologyBuilder2D/3D`：从 `T_*` 字段构建 MM/FI/NNC/FF 连接。
- `FIM_ConnectionManager`：
  - 连接类型：`Matrix_Matrix`, `Fracture_Internal`, `Matrix_Fracture`, `Fracture_Fracture`
  - NNC/FF 可并联聚合；MM/FI 重复连接直接报错
- `FIM_BlockSparseMatrix<N>`：
  - block Jacobian/Residual 容器
  - `FreezePattern()` 后禁止新增连接
  - `GetFrozenMatrix()` 走 CSR cache（Issue11 关键）
- `FIM_GlobalAssembler`：`AssembleAccumulation`, `AssembleFlux`, `AssembleSource`
- `BoundaryAssembler`：2D/3D 边界与井源 FullJac 装配
- `WellDOFManager<N>`：每口井追加独立 BHP block

### 5.5 瞬态引擎子模块（`FIM_TransientEngine/`）

- `Types.hpp`：`TransientSolverParams`、`LinearSolverType`、诊断开关
- `StateSync.hpp`：`InjectStaticProperties`, `SyncStateToFieldManager`
- `StepKernels.hpp/_impl.hpp`：线性求解与 trial update 内核
- `Diagnostics.hpp`：3D precheck + `RunMatrixAssemblyAudit`
- `RunGeneric.hpp/_impl.hpp`：通用 FIM 瞬态主循环

---

## 6. 理论公式与代码映射

| 对象 | 公式/思想 | 代码实现 |
|---|---|---|
| 两点串联阻力 | `R = d1/K1 + d2/K2` | `FVM_Ops::Op_Math_SeriesResistance` |
| 两点传导率 | `T = A / R` | `FVM_Ops::Op_Math_Transmissibility` |
| MM 传导率 | 面法向距离 + 投影渗透率 | `TransmissibilitySolver_2D/3D::Calculate_Transmissibility_Matrix` |
| FI 传导率 | 邻接裂缝微元串联阻力 | `Calculate_Transmissibility_FractureInternal` |
| NNC 传导率 | `A_int / (d_m/K_m + w_f/(2K_f))` | `Calculate_Transmissibility_NNC` |
| FF 传导率 | 交汇处 Star-Delta 展开 | `Calculate_Transmissibility_FF` + `ff_topology` |
| 势能差 | `ΔΦ = ΔP - ρ g·Δx` | `FVM_Ops::Compute_Potential_Diff` |
| 上游选择 | `ΔΦ` 符号判定迎风侧 | `FVM_Ops::Op_Upwind_AD` |
| 非正交修正 | `corr = ∇φ_f · vectorT` | `FVM_Ops::Op_NonOrthogonal_PressureCorr` |
| 两相本构 | VG `Pc(Sw)` + Mualem `krw/krg` | `CapRelPerm_HD.h` |
| AD 物性桥接 | 表格/解析 EOS + 鲁棒导数 + 链式组装 | `AD_FluidEvaluator.h` |
| 井 BHP 方程 | `R_well = P_wbh - P_target` | `WellDOFManager::AssembleWellEquations` |
| 井 Rate 方程 | `R_well = WI_mob(P_res-P_wbh)-q_spec` | `WellDOFManager::AssembleWellEquations` |

---

## 7. Day6 Campaign 事实快照（当前实现）

### 7.1 场景定义

- 物理场景：`T01..T16`（4 个布尔轴组合）
  - `variable_properties`
  - `with_wells`
  - `two_phase`
  - `heat_enabled`
- 拓扑轴：`F0`（无裂缝）/`F1`（单裂缝）/`F2`（交叉裂缝）

### 7.2 闯关与回归联动

- 每个 `Txx` 必须 `F0 -> F1 -> F2` 全通过（fail-fast）。
- 里程碑 `T04/T08/T12/T16` 自动触发：
  - `day6_matrix_audit_*`
  - `issue11_frozen_matrix`
  - `issue12_linear_solver`

### 7.3 输出结构

- 单 case：`Test/Transient/Day6/<2D|3D>_Txx_Fy/`
  - `convergence.log`
  - `final.vtk`
- Campaign 指标：`Test/Transient/Day6/Campaign/<2D|3D>/Txx_Fy/metrics.csv`
  - 含 `topology` 字段

---

## 8. I/O、命名与依赖

### 8.1 边界 Tag

- `MeshDefinitions.h`：`LEFT=1`, `RIGHT=2`, `BOTTOM=3`, `TOP=4`, `TAG_FRONT=5`, `TAG_BACK=6`

### 8.2 字段命名统一来源

- `SolverContrlStrName_op.h`
- 主变量：`p_w`, `s_w`, `T`
- 典型物性：`rho_w/rho_g`, `mu_w/mu_g`, `h_w/h_g`, `lambda_*`, `K_xx/K_yy/K_zz`, `phi`

### 8.3 后处理能力

- `PostProcess_2D`：`ExportTecplot`, `ExportVTK`, `ExportVTU`, `ExportPVD`, `SyncADFieldToScalar`
- `PostProcess_3D`：`ExportTecplot`, `ExportVTK`, `SyncADFieldToScalar`
- VTK 支持 `DomainID` 分域可视化

### 8.4 构建依赖（VS 工程）

- Visual Studio v143 / MSBuild
- Eigen
- Gmsh (`gmsh.lib`)
- OpenMP

---

## 9. 构建与运行

```powershell
msbuild .\2D-Unstr-Quadrilateral-EDFM.sln /p:Configuration=Debug /p:Platform=x64
.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --list
.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day2_fvm_ad
.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_campaign_2d_all
```

---

## 10. 风险与维护约定

- 风险1：其他文档可能保留 Day7 case 名称，但 dispatcher 目前未接入。
- 风险2：`ff_topology` 与 `T_FF_*` 长度不一致会在 Builder 阶段抛错。
- 风险3：`FreezePattern()` 后新增连接会报错，井耦合 pattern 必须提前注册。
- 风险4：字段字符串变更需同步 `SolverContrlStrName_op.h` 与加载/装配调用端。

维护要求：

1. 新增/删除 case 时，同步更新 `main.cpp --list` 与本文件第 3 节。
2. 关键功能改动后，至少更新 Metadata 与风险/现状描述。
3. 文档冲突时先修正代码事实，再回写文档，不反向“改代码迎合旧文档”。

---

## 11. 新会话快速接手建议（优先阅读）

1. `main.cpp`
2. `FIM_TransientEngine/RunGeneric_impl.hpp`
3. `FIM_TransientEngine/Types.hpp`
4. `WellDOFManager.h`
5. `FIM_BlockSparseMatrix.h`
6. `TransmissibilitySolver_2D.cpp`
7. `TransmissibilitySolver_3D.cpp`
8. `FIM_TopologyBuilder2D.h`
9. `FIM_TopologyBuilder3D.h`
10. `FIM_ConnectionManager.h`
11. `MeshManager.h`
12. `3D_MeshManager.h`
13. `Test_Day6_TransientSolver.cpp`
