## Day 5 最终方案（稳健通用 + 取热全覆盖）

### Summary
1. Day5 交付生产级全局组装核心，不走最小补丁路线。  
2. 新增两个显式入口：`day5_global_jac_2d`、`day5_global_jac_3d`。  
3. 强制覆盖“超临界 CO2 非线性 + 注采对井取热”场景。  
4. 取热场景按你已确认执行：`单相+两相`、`2D+3D` 全覆盖。  
5. 当日回归门槛：`R1 + R4 + R5` 全通过。

### 一、实现范围（生产代码）
1. 新增 `FIM_BlockSparseMatrix`：统一块稀疏矩阵容器（2x2/3x3）。  
2. 新增 `FIM_GlobalAssembler`：统一组装 `Residual + Jacobian`（accumulation/flux/boundary/well/leakoff）。  
3. 新增 `FIM_JacobianVerifier`：统一 `FD vs AD` 矩阵级比较器。  
4. 复用现有 `BoundaryAssembler`（边界、井、滤失）作为源项/边界贡献输入，不重写旧逻辑。  
5. 在 `main.cpp` 增加 Day5 双 case 路由。

### 二、测试文件与入口
1. 新建测试头文件：`Test_Day5_GlobalAssembly_Jacobian.h`。  
2. 在该文件中提供 Day5 运行入口：  
`Run_Day5_GlobalAssembly_Jacobian_2D<3, ADVar<3>>()`  
`Run_Day5_GlobalAssembly_Jacobian_3D<3, ADVar<3>>()`  
3. `main.cpp` 对应 case：  
`--case=day5_global_jac_2d`  
`--case=day5_global_jac_3d`

### 三、Day5 场景矩阵（强制）
1. 2D 单相：超临界 CO2 工况 + 注采对井取热。  
2. 2D 两相：非混水/CO2 + 注采对井取热。  
3. 3D 单相：超临界 CO2 工况 + 注采对井取热。  
4. 3D 两相：非混水/CO2 + 注采对井取热。  

### 四、验收判据（含你新增需求）
1. `FD vs AD max relative error < 1e-6`（四个场景都必须满足）。  
2. 超临界 CO2 双层门槛：  
稳定超临界点（如 `15MPa/400K`）严格 `<1e-6`。  
近临界点（如 `8MPa/310K`）做鲁棒门禁（无 NaN/Inf、误差可控并输出）。  
3. 注采对井取热门槛：  
生产井能量外流符号正确。  
注井能量注入符号正确。  
系统净取热方向正确（净热采出为正）。  
4. 回归门槛：`R1`、`R4`、`R5` 全 PASS。  
5. 日志必须无 `Error/Exception/Fatal`。

### 五、涉及文件
1. 新增：`FIM_BlockSparseMatrix.h/.cpp`  
2. 新增：`FIM_GlobalAssembler.h/.cpp`  
3. 新增：`FIM_JacobianVerifier.h/.cpp`  
4. 新增：`Test_Day5_GlobalAssembly_Jacobian.h`  
5. 修改：`main.cpp`  
6. 修改：`REGRESSION_CASES.md`  
7. 修改：`PROJECT_CONTEXT.md`

### 六、与 IMPES/Day6 衔接
1. Day5 结果可直接复用到 IMPES（通过 `is_impes` 组装开关退化）。  
2. Day6 仅需接入 Newton/线性求解/时间步控制，不需重写组装核心。  
3. 方程完整求解闭环仍在 Day6 实现，Day5 负责“可求解前正确性基线”。

### 七、风险与回退
1. 近临界 FD 噪声风险：采用自适应扰动与上下限保护。  
2. 2x2/3x3 映射错位风险：先做小系统索引断言再跑全场景。  
3. 复杂度上升风险：旧路径保留，新模块仅由 Day5 case 驱动，不影响已有 Day1-Day4 回归。

### 八、假设与默认
1. 保持现有数据结构优先，不新增持久业务实体。  
2. Day5 允许新增“生产组装核心模块”，不进入完整瞬态求解。  
3. 文档与 case 名同步更新为同一变更集（`REGRESSION_CASES.md` + `PROJECT_CONTEXT.md`）。

