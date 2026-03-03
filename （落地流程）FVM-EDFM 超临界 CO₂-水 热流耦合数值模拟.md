收到！作为这一前沿领域的同行，我非常荣幸能继续协助您推进这一极具挑战性的项目。

您目前的进展非常关键且扎实：**网格拓扑、场数据管理、静态几何交互系数**这“三驾马车”已经就位，这为 FVM 离散奠定了坚实的几何与数据基础。

依据《FVM-EDFM 超临界 CO₂-水 热流耦合数值模拟实施规范（工程落地版）》（以下简称《规范》），为了实现最终的 **超临界 CO₂-水 热流耦合数值模拟**，我们需要从当前的静态几何阶段，跨越到动态物理阶段。

以下是为您制定的**详细可执行工作计划 (Step-to-Step Implementation Plan)**。该计划分为五个核心阶段，旨在确保代码的模块化、鲁棒性和可扩展性。

------

### **阶段一：物理属性计算内核 (The Physics Kernel)**

**目标**：构建独立于网格的高精度物性计算模块，为控制方程提供 $\rho, \mu, h$ 等关键参数。

**依据**：《规范》中的 **Equation of State (EOS)** 章节。

1. **Step 1.1: 基础热力学状态变量管理**
   - **任务**：    定义 `ThermodynamicState` 结构体，封装 $(P, T,  S_w)$ 等状态变量。
   - **关键点**：支持自动微分变量（ADVar）的传递，为后续雅可比矩阵（Jacobian）计算做准备。
2. **Step 1.2: 超临界 CO₂ 物性计算 (EOS_CO2)**
   - **任务**：实现 Span-Wagner 方程（高精度）或其高效查表/拟合版本。
   - **输出**：$\rho_{co2}(P, T)$, $\mu_{co2}(P, T)$, $h_{co2}(P, T)$。
   - **难点**：临界点附近的数值稳定性处理。
3. **Step 1.3: 水/水物性计算 (EOS_Water/Brine)**
   - **任务**：实现 IAPWS-97 标准或相关水物性关联式。
   - **输出**：$\rho_{w}(P, T)$, $\mu_{w}$, $h_{w}$。
4. **Step 1.4: 相对渗透率与毛管力 (RelPerm & Capillary)**
   - **任务**：实现 Brooks-Corey 或 Van Genuchten 模型。
   - **输出**：$k_{rw}(S_w)$, $k_{rg}(S_g)$, $P_c(S_w)$。
   - **验证**：单元测试，确保 $0 \le k_r \le 1$ 且导数连续。

------

### **阶段二：空间离散算子 (Spatial Discretization Operators)**

**目标**：基于当前的 FieldManager，实现从体心到面心的插值与梯度计算，完成 FVM 核心离散。

**依据**：《规范》中的 **Discretization Scheme** 章节。

1. **Step 2.1: 梯度计算算子 (Gradient Operator)**

   - **任务**：实现基于 Green-Gauss 公式（推荐）或最小二乘法的梯度计算模块 `FVM_Grad`。
   - **输入**：体心标量场（如 $P, T$）。
   - **输出**：体心梯度向量场（$\nabla P, \nabla T$）。
   - **用途**：用于非正交网格的扩散项修正（Cross-Diffusion）及高阶重构。

   补充：

   确立了 `FVM_Grad` 为纯几何算子，未来的物理边界处理应遵循以下规范：

   1. **裂缝尖端/边缘**：
      - 在 **组装矩阵 (Assembly)** 遍历 Edge 时，检测到 `neighbor == -1`。
      - **显式** 决定是否添加通量项。
      - 对于 EGS 模拟，通常默认裂缝尖端不与围岩之外交换流体，因此直接 **不累加该 Edge 的通量** 即可。这比在梯度算子离掺杂虚拟点更准确、更可控。
   2. **基岩外边界**：
      - 依靠 `boundaryFaces` 和 `BoundaryConditionManager`。
      - 梯度算子计算出的 $\nabla P$ 仅用于高阶格式（二阶精度）的 **源项修正 (Source Term Correction)** 或 **非正交修正**，而不应直接决定边界通量的主值。

2. **Step 2.2: 面心插值算子 (Interpolation Operator)**

   - **任务**：实现从体心到面心的值插值模块 `FVM_Interp`。
   - **子任务 2.2.1**：**中心差分（Central Difference）**：用于扩散项（导热、压力扩散）。
   - **子任务 2.2.2**：**迎风格式（Upwind Scheme）**：用于对流项（移动能力 $\lambda$, 密度 $\rho$）。需结合流体势能梯度方向判断。

3. **Step 2.3: 动态通量组装 (Dynamic Flux Assembly)**

   - **任务**：结合已有的 `Transmissibility`（静态）和上述插值结果（动态），计算各面的质量/能量通量。
   - **公式**：$F_{ij} = T_{ij}^{static} \cdot \lambda_{up} \cdot \rho_{up} \cdot ( \Delta P - \rho_{avg} g \Delta z )$。

------

### **阶段三：全隐式方程组装 (FIM Assembly Framework)**

**目标**：将物理方程转化为线性代数系统 $Ax=b$。

**依据**：《规范》中的 **Jacobian Assembler** 章节。

1. **Step 3.1: 自动微分变量类 (ADVar Class)**
   - **任务**：实现或引入轻量级 `ADVar` 类（Value + Derivatives）。
   - **功能**：支持加减乘除及常用数学函数的链式法则，自动记录对主变量（$P, T, S$）的偏导数。
2. **Step 3.2: 残差方程构建 (Residual Builder)**
   - **任务**：针对每个网格（基岩 & 裂缝），构建质量守恒和能量守恒残差方程。
   - **Accumulation Term**：$\frac{V}{\Delta t} [(\phi \rho S)^{n+1} - (\phi \rho S)^n]$。
   - **Flux Term**：$\sum F_{face} + \sum F_{NNC} + \sum F_{FF}$。
   - **Source/Sink Term**：井源项处理（Peaceman 模型）。
3. **Step 3.3: 雅可比矩阵填充 (Jacobian Filler)**
   - **任务**：利用 `ADVar` 的导数信息，组装稀疏矩阵（Sparse Matrix）。
   - **难点**：高效处理 CSR/CSC 格式的索引映射（利用 `meshMgr` 的全局索引）。

------

### **阶段四：非线性与线性求解器 (Solver Strategy)**

**目标**：求解耦合的非线性方程组。

**依据**：《规范》中的 **Solver Loop** 章节。

1. **Step 4.1: 线性求解器接口 (Linear Solver Interface)**
   - **任务**：集成 Eigen 或 PETSc。
   - **配置**：BiCGSTAB + ILU/AMG 预处理（针对非对称矩阵）。
2. **Step 4.2: 牛顿-拉夫逊迭代 (Newton-Raphson Loop)**
   - **任务**：实现非线性迭代主循环。
   - **逻辑**：
     1. 更新物性 (EOS)。
     2. 组装 $J$ 和 $R$。
     3. 求解 $J \delta x = -R$。
     4. 更新变量 $x^{k+1} = x^k + \delta x$（含阻尼/限幅）。
     5. 检查收敛性（$|R| < \epsilon$ 且 $|\delta x| < \epsilon$）。

------

### **阶段五：时间步进与整体驱动 (Time Stepping & Driver)**

**目标**：模拟时间演化。

1. **Step 5.1: 时间步长控制 (Time Step Control)**
   - **任务**：实现自适应时间步长策略（基于牛顿迭代次数或变量变化率）。
2. **Step 5.2: 主程序驱动 (Main Driver)**
   - **任务**：整合所有模块，编写 `main.cpp`。
   - **流程**：初始化 -> 时间循环 -> 牛顿循环 -> 输出结果 -> 下一步。

