## Day4 完整计划（最终工程落地标准版）

### 摘要
本计划将 Day4 从“最小闭环”升级为“最终工程标准预埋”：`2D+3D`、`质量+能量井源耦合`、`配置文件驱动 WAG`、`Matrix+Fracture 完井`。  
实现目标是先把井模型接口、装配链路、调度链路一次性定型，确保 Day5-Day7 只需在同一接口上扩展求解与回归，不再返工井控框架。

### 1) 我对今日任务的理解
1. 井指数统一切换为非结构网格几何等效半径 WI，2D 与 3D 都必须支持。  
2. 井控模式至少支持 `BHP` 与 `Rate`，并可分别作用于压力、两相质量、能量方程。  
3. WAG 调度必须来自外部配置文件，不再硬编码。  
4. 完井对象必须覆盖 `Matrix` 与 `Fracture` 两域，并进入统一装配入口。  
5. Day4 交付应包含新的显式 case 与验收判据，保证可脚本化复现。

### 2) 阶段1执行计划

### A. 问题根因假设与验证路径
1. 根因 H1：当前仓库没有井模型数据入口、井源离散算子、井装配入口。  
2. 根因 H2：当前无配置文件驱动调度器，WAG 仅能靠代码硬写。  
3. 根因 H3：当前无 Matrix+Fracture 统一完井抽象，无法同一流程注入两域。  
4. 验证 V1：静态检查新增接口从 `main.cpp case -> 配置加载 -> 井控求值 -> residual/jacobian` 全链路可达。  
5. 验证 V2：`day4_well_patch` 运行日志出现每类 PASS 关键行（2D/3D、BHP/Rate、WAG 切换、Matrix/Fracture）。  
6. 验证 V3：`trans_2d` 与 `trans_3d` 仍 PASS，证明旧链路未破坏。  

### B. 受影响文件清单（仅必要文件）
1. `FVM_Ops_AD.h`。  
2. `BoundaryAssembler.h`。  
3. `BoundaryAssembler.cpp`。  
4. `Test_FVM_Ops_AD.h`。  
5. `main.cpp`。  
6. `REGRESSION_CASES.md`。  
7. `PROJECT_CONTEXT.md`。  
8. 新增 `Well_WellControlTypes.h`。  
9. 新增 `Well_WellScheduleManager.h`。  
10. 新增 `Well_WellScheduleManager.cpp`。  
11. 新增 `Test/WellSchedule/day4_wag_schedule.csv`。  

### C. 逐步实施计划（函数级）
1. 在 `Well_WellControlTypes.h` 定义最小必需类型：`WellControlMode`、`WellTargetDomain`、`WellComponentMode`、`WellControlRecord`、`WellScheduleStep`，仅包含装配所需字段。  
2. 在 `Well_WellScheduleManager.h/.cpp` 实现配置读取与时间查询接口：`LoadFromCsv(path)`、`GetActiveStep(time)`、`Validate()`。  
3. CSV 固定 schema：`t_start,t_end,well_name,domain,control_mode,target_value,component_mode,rw,skin,well_axis,completion_id,wi_override`。  
4. 在 `FVM_Ops_AD.h` 新增井源 AD 算子：`Op_Well_BHP_Source_AD`、`Op_Well_Rate_Source_AD`、`Op_Well_Energy_Source_AD`，全部遵守 Outflow positive。  
5. 在 `BoundaryAssembler.cpp` 新增 WI helper：`Compute_WI_GeometricEqRadius_2D` 与 `Compute_WI_GeometricEqRadius_3D`。  
6. 2D WI 公式固定：`r_eq=sqrt(A/pi)`、`A=cell.volume`、`k=sqrt(Kxx*Kyy)`、`WI=2*pi*k*h/(ln(r_eq/rw)+S)`、`h=1.0`。  
7. 3D WI 公式固定：`r_eq=sqrt(V/(pi*L))`、`V=cell.volume`、`L=井段在单元内穿透长度`、`WI=2*pi*k*L/(ln(r_eq/rw)+S)`。  
8. 3D 有效渗透率 `k` 采用井轴向主导平面几何均值：`axis=Z -> sqrt(Kxx*Kyy)`、`axis=X -> sqrt(Kyy*Kzz)`、`axis=Y -> sqrt(Kxx*Kzz)`。  
9. 在 `BoundaryAssembler.h/.cpp` 新增统一井装配入口：`Assemble_Wells_2D(...)` 与 `Assemble_Wells_3D(...)`，输入 `FieldManager`、`MeshManager`、`WellScheduleStep`、当前时间、`residual/jacobianDiag`。  
10. 完井域处理策略：`Matrix` 用自动 WI；`Fracture` 优先 `wi_override`，若为空且几何量齐全才自动算，否则报错并跳过该完井。  
11. 质量耦合策略：`component_mode=water/gas/total` 分配到对应方程块；若 `total`，按输入分配因子写入水/气方程。  
12. 能量耦合策略：`q_E = q_mass_w*h_w + q_mass_g*h_g`，焓值来自现有场名 `h_w/h_g`。  
13. 在 `Test_FVM_Ops_AD.h` 新增测试：`Test_Well_WI_2D_3D_Geometric`、`Test_Well_BHP_Rate_MassEnergy`、`Test_WAG_Config_Switching`、`Test_Matrix_Fracture_Completions`、`Run_Day4_Well_Patch`。  
14. 在 `main.cpp` 注册 `day4_well_patch`，并在 `--list` 可见。  
15. 在 `REGRESSION_CASES.md` 增加 Day4 必跑命令与 PASS 关键行。  
16. 在 `PROJECT_CONTEXT.md` 更新元信息与 Day4 能力说明。  

### D. 当日验收标准（对齐 PLAN_7D_CHECKLIST + REGRESSION_CASES）
1. `--case=day4_well_patch` 返回码 0。  
2. 日志出现并通过：`[PASS] WI 2D geometric-eq`、`[PASS] WI 3D geometric-eq`、`[PASS] BHP/Rate mass coupling`、`[PASS] Well energy coupling`、`[PASS] WAG config switching`、`[PASS] Matrix+Fracture completion`。  
3. 日志中无 `Error/Exception/Fatal`。  
4. `--case=trans_2d` PASS。  
5. `--case=trans_3d` PASS。  

### E. 风险点与回退策略
1. 风险：`ln(r_eq/rw)+S` 接近 0。回退：分母阈值保护，记录告警并跳过该条完井。  
2. 风险：3D 穿透长度 `L` 获取困难。回退：若无法稳定求交，配置中允许显式 `L_override`，并在日志标注来源。  
3. 风险：Fracture 自动 WI 几何量不足。回退：Fracture 完井强制使用 `wi_override`，避免伪物理估算。  
4. 风险：质量与能量耦合符号错配。回退：先算子级断言，再装配级守恒断言，任何失败立即阻断 Day4 PASS。  
5. 风险：配置文件格式漂移。回退：`Validate()` 强校验列名与单位，失败即拒绝运行。  

### 3) 重要接口/类型变化（公开面）
1. 新增公开类型：`WellControlMode`、`WellTargetDomain`、`WellComponentMode`、`WellScheduleStep`。  
2. 新增公开类：`WellScheduleManager`（CSV 读取与分段查询）。  
3. 新增公开装配入口：`BoundaryAssembler::Assemble_Wells_2D`、`BoundaryAssembler::Assemble_Wells_3D`。  
4. 新增公开 case：`--case=day4_well_patch`。  
5. 现有接口保持兼容，不改旧 case 名，不改旧字段名。  

### 4) 测试场景清单
1. 2D Matrix BHP 生产井：`P_cell>P_bhp`，`q>0`，`dq/dP>0`。  
2. 2D Matrix BHP 注入井：`P_cell<P_bhp`，`q<0`，`dq/dP>0`。  
3. 2D Matrix Rate 井：`q=q_target`，`dq/dP=0`。  
4. 3D Matrix BHP/Rate 与 2D 同符号逻辑。  
5. Fracture 完井（wi_override）在 2D/3D 均可写入对应方程块。  
6. 质量+能量耦合：同一井在水气质量与能量方程同时产生一致源项。  
7. WAG CSV 三阶段切换：水注入 -> 气注入 -> 采出，切换时刻与目标值一致。  
8. Leakoff OFF/ON 下井控趋势可解释。  
9. 回归保护：`trans_2d`、`trans_3d` 均 PASS。  

### 5) 显式假设与默认值
1. 维度：2D+3D 同时落地。  
2. 完井域：Matrix+Fracture。  
3. 方程耦合：压力+两相质量+能量。  
4. WI 模型：几何等效半径法。  
5. 2D 有效渗透率：`sqrt(Kxx*Kyy)`。  
6. 3D 有效渗透率：按井轴向主导平面几何均值。  
7. 默认参数：`rw=0.1 m`、`S=0`、`h=1.0`。  
8. 符号约定：Outflow positive。  
9. 调度输入：CSV 配置文件驱动。  
10. Fracture 完井默认要求 `wi_override`；自动 WI 仅在几何量完整时启用。  
