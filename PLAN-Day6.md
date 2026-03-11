## Day6 最终落地计划（全 FIM-Newton + 预埋 IMPES 切换接口）

### Summary
1. Day6 六个场景统一走全 `FIM-Newton`，满足 `>=50` 连续时间步稳定、完整非线性日志、物理护栏、回滚与自适应 `dt`。
2. 两相路径预埋 `IMPES` 切换接口，但 Day6 默认与验收均固定为 `FIM`，不启用 IMPES 分支。
3. VTK 输出按你要求执行“每 10 步一帧 + final”。

### 先读取并对齐的文件（实施前必读）
1. [PROJECT_CONTEXT.md](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/PROJECT_CONTEXT.md)
2. [REGRESSION_CASES.md](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/REGRESSION_CASES.md)
3. [PLAN_7D_CHECKLIST.md](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/PLAN_7D_CHECKLIST.md)
4. [main.cpp](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp)
5. [Test_FVM_Ops_AD.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_FVM_Ops_AD.h)
6. [Test_Day5_GlobalAssembly_Jacobian.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_Day5_GlobalAssembly_Jacobian.h)
7. [FIM_BlockSparseMatrix.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_BlockSparseMatrix.h)
8. [FIM_GlobalAssembler.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_GlobalAssembler.h)
9. [FIM_StateMap.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_StateMap.h)
10. [FIM_TopologyBuilder2D.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TopologyBuilder2D.h)
11. [FIM_TopologyBuilder3D.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TopologyBuilder3D.h)
12. [2D_PostProcess.cpp](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/2D_PostProcess.cpp) 与 [3D_PostProcess.cpp](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/3D_PostProcess.cpp)

### 问题根因假设与验证路径
1. 根因 H1：Day6 case 尚未在 dispatcher 中可运行。
2. 根因 H2：已有 Day5 组装基础，但缺少统一的“瞬态推进 + Newton 控制 + 回滚 + 护栏 + Day6 导出”运行器。
3. 根因 H3：两相路径尚无可切换推进器接口，未来扩展 IMPES 风险高。
4. 验证 V1：六个 Day6 case 均可被 `--list` 列出并可执行。
5. 验证 V2：每步打印 `step/iter/residual/dt/rollback/limiter_count`，失败步有 rollback，最后一步收敛。
6. 验证 V3：每 case 至少 50 步稳定、无 `Error/Exception/Fatal`。
7. 验证 V4：`Test/Transient/Day6/...` 下 VTK 每 10 步与 final 文件存在且非空，字段满足 `P/T`（两相含 `S_w`）。
8. 验证 V5：回归护栏 `day2_fvm_ad + trans_2d/trans_3d` 通过。

### 受影响文件（仅必要）
1. [main.cpp](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp)
2. [Test_Day6_TransientSolver.h](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_Day6_TransientSolver.h)（新增）
3. [Test_Day6_TransientSolver.cpp](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_Day6_TransientSolver.cpp)（新增）
4. [2D-Unstr-Quadrilateral-EDFM.vcxproj](d:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/2D-Unstr-Quadrilateral-EDFM.vcxproj)

### 逐步实施计划（函数级，可直接照做）
1. 新增 `Test_Day6_TransientSolver.h/.cpp`，定义 6 个公开入口函数，对应 6 个 Day6 case 名。
2. 在 `main.cpp` 注册 6 个 Day6 case 到 dispatcher 与 `--list`。
3. 实现 2D/3D 场景构建函数，复用 Day4/Transmissibility 的网格、裂缝、DOF、场初始化流程。
4. 实现传导率预处理函数，调用 `TransmissibilitySolver_2D/3D` 生成 `T_Matrix/T_NNC/T_FI/T_FF`。
5. 实现连接聚合函数，使用 `FIM_TopologyBuilder2D/3D + FIM_ConnectionManager` 得到连接列表。
6. 实现状态映射函数，使用 `FIM_StateMap<2/3>` 保存 `P/T/(S_w)` 并与 `FieldManager` 双向同步。
7. 实现系统组装函数，按“积累项 + 通量项 + 井源项”写入 `FIM_BlockSparseMatrix`。
8. 实现 Newton 迭代函数，采用 Eigen `BiCGSTAB` 解线性系统，计算残差范数、更新增量。
9. 实现稳健档控制策略：`tol=1e-6`、`max_iter=8`、失败回滚、`dt` 减半、成功后 `dt` 受限放大。
10. 实现物理护栏与限幅：`S_w` 限制在 `[0,1]`，`P/T` 下限保护，检测 NaN/Inf。
11. 实现推进器分派接口（两相）：`enum SolverRoute { FIM, IMPES }`，Day6 调用固定传 `FIM`；`IMPES` 仅预留函数签名与空实现位。
12. 实现 VTK 导出函数：每 10 步一帧 + final，路径为 `Test/Transient/Day6/<scenario>/`，导出后做非空校验并打印 PASS 行。
13. 每个 case 结束输出汇总：总步数、回滚次数、限幅次数、最终残差、导出文件数。

### 当日验收标准（对齐 PLAN_7D_CHECKLIST + REGRESSION_CASES）
1. Gate A：6 个 Day6 命令全部完成，且各自 `>=50` 步。
2. Gate B：每步有 `iter/residual/dt`；非收敛步有 rollback 记录；最终步收敛。
3. Gate C：无 NaN/Inf；两相 `S_w` 全程在 `[0,1]`；`P/T` 无越下限。
4. Gate D：每个 case 至少 1 个 VTK，实际执行为每 10 步一帧 + final，文件存在且非空，字段正确。
5. Gate E：ParaView 手工验收通过（由你本地完成并回传日志/截图路径）。
6. Gate F：`--case=day2_fvm_ad` PASS。
7. Gate G：受影响维度 `--case=trans_2d` 或 `--case=trans_3d` PASS（建议两者都跑）。

### 风险点与回退策略
1. Newton 发散风险：触发 rollback + `dt` 减半，超过最大回滚次数则该 case 失败并输出最小诊断。
2. 两相震荡风险：对 `ΔS_w` 限幅并裁剪，记录 `limiter_count`，必要时进一步缩小 `dt`。
3. 高频导出 I/O 风险：保留导出开关参数，默认“每 10 步”不变，必要时降为 `initial/mid/final`。
4. 预埋 IMPES 误触发风险：Day6 case 硬编码 `SolverRoute::FIM`，IMPES 分支不参与当日验收。

### 假设与默认
1. Day6 不变更任何既有 case 名与回归命令字符串。
2. 不新增持久化业务数据结构，状态容器沿用 `FIM_StateMap`。
3. IMPES 仅“接口预埋”，不在 Day6 启用。
4. 如后续开放 `--solver=impes`，需同步更新 `REGRESSION_CASES.md` 与 `PROJECT_CONTEXT.md`。