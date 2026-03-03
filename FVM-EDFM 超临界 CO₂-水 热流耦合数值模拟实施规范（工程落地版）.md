这是一份经过深度整合、逻辑重构并润色后的**《FVM-EDFM 超临界 CO₂-水 热流耦合数值模拟实施规范（工程落地版）》**。

这份文档将物理模型（控制方程）、数值方法（FVM离散与EDFM耦合）以及工程实现方案（基于AD的FIM组装）融为一体，旨在指导开发一套工业级、高稳定性的 EGS 模拟器核心模块。

------

# FVM-EDFM 超临界 CO₂-水 热流耦合数值模拟实施规范（工程落地版）

**文档性质**：工程实施标准 / 核心算法白皮书

**适用框架**：全隐式方法 (Fully Implicit Method, FIM)

**线性求解**：Eigen 库 (SpMat, BiCGSTAB + ILUT/AMG)

**网格类型**：三维非结构化网格 + 嵌入式离散裂缝模型 (EDFM)

**物理场景**：

1. **两相热流耦合**：超临界 CO₂ - 水非等温多相流。
2. **单相热流耦合**：纯超临界 CO₂ 循环（作为两相的退化特例）。

------

## 第一章 求解体系定义

### 1.1 主变量与未知量向量

在 FIM 框架下，每个网格单元 $i$（无论是基岩 Matrix 还是裂缝 Fracture）均包含待求解的独立变量向量 $\mathbf{X}_i$。

- **两相流场景** ($N_{dof}=3$)：

  $$\mathbf{X}_i = \begin{bmatrix} p_{w, i} \\ S_{n, i} \\ T_i \end{bmatrix}$$

  - $p_w$: 水相压力
  - $S_n$: 非润湿相 ($CO_2$) 饱和度 ($S_w = 1 - S_n$)
  - $T$: 温度

- **单相流场景** ($N_{dof}=2$)：

  $$\mathbf{X}_i = \begin{bmatrix} p_{i} \\ T_i \end{bmatrix}$$

  - 退化处理：令 $S_n = 1$ (纯气) 或 $S_n = 0$ (纯水)，移除饱和度方程。

### 1.2 残差方程向量 (Residual Vector)

对应守恒定律，每个单元 $i$ 需满足残差方程 $\mathbf{R}_i \to \mathbf{0}$：

$$\mathbf{R}_i = \begin{bmatrix} R_{water, i} \\ R_{CO2, i} \\ R_{energy, i} \end{bmatrix}$$

*(注：单相时移除 $R_{CO2}$)*

### 1.3 牛顿-拉夫逊迭代 (Newton-Raphson)

将非线性系统线性化为 $A \mathbf{x} = B$ 的形式：

- **$A$ (Jacobian 矩阵)**: 全局稀疏矩阵，由 $3 \times 3$ (或 $2 \times 2$) 的块组成。

  $$J_{ij} = \frac{\partial \mathbf{R}_i}{\partial \mathbf{X}_j}$$

- **$x$ (更新量)**: $\delta \mathbf{X}$

- **$B$ (右端项)**: $-\mathbf{R}(\mathbf{X}^k)$ (当前残差的负值)

迭代更新：$\mathbf{X}^{k+1} = \mathbf{X}^k + \delta \mathbf{X}$，包含限幅处理（如 $S_n \in [0,1]$）。

------

## 第二章 几何基础与通量算子

### 2.1 几何预处理

对于共享面 $f$ 的两个单元 $i$ 和 $j$：

- $\mathbf{d}_{ij}$: 单元中心距离向量 ($\mathbf{x}_j - \mathbf{x}_i$)。

- $\mathbf{A}_f$: 面面积向量 ($Area_f \cdot \mathbf{n}_{ij}$)。

- **正交分解 (Orthogonal Decomposition)**：

  $$\mathbf{A}_f = \mathbf{E}_{ij} + \mathbf{T}_{ij}$$

  - $\mathbf{E}_{ij}$ (正交分量): 平行于 $\mathbf{d}_{ij}$，承载主导通量。
  - $\mathbf{T}_{ij}$ (交叉分量): 垂直于 $\mathbf{d}_{ij}$，用于非正交修正。

### 2.2 离散通量算子 (Flux Operators)

为实现代码复用并支持 FIM 自动微分，定义以下通用算子。

*注意：算子中的 $\Lambda$（流度、密度等）均为含导数的变量。*

1. **正交传导算子 $\mathcal{F}_{Ortho}(\Lambda, \Delta \Psi)$**:

   $$\mathcal{F}_{Ortho} = \mathcal{T}_{ij}^{Geo} \cdot \Lambda_{avg} \cdot (\Psi_i - \Psi_j)$$

   - $\mathcal{T}_{ij}^{Geo} = \frac{|\mathbf{E}_{ij}|}{|\mathbf{d}_{ij}|}$: 纯几何因子。
   - $\Lambda_{avg}$: 界面传递系数（通常采用调和平均或上游加权）。

2. **交叉修正算子 $\mathcal{F}_{Cross}(\Lambda, \nabla \Psi)$**:

   采用 **延迟修正 (Deferred Correction)** 策略，在牛顿迭代内部作为显式源项处理：

   $$\mathcal{F}_{Cross} = \Lambda_{f} \cdot (\overline{\nabla \Psi}_f \cdot \mathbf{T}_{ij})$$

   - $\overline{\nabla \Psi}_f$: 梯度重构值（基于最小二乘法或格林-高斯公式）。

------

## 第三章 基岩系统控制方程离散 (Matrix Discretization)

所有通量项必须包含三部分驱动力：**压力驱动**、**毛细驱动**、**重力驱动**。

### 3.1 质量守恒方程

#### **方程 1: 水相残差 $R_{w,i}$**

$$R_{w,i} = \frac{V_i}{\Delta t} \left[ (\phi \rho_w S_w)^{n+1} - (\phi \rho_w S_w)^n \right]_i + \sum_{j \in Faces} F_{w,ij}^{n+1} - Q_{w,i}^{n+1}$$

**水相通量 $F_{w,ij}$ 表达式** (含重力严格投影)：

$$F_{w,ij} = \lambda_w^{up} \rho_w^{up} \left[ \underbrace{ \mathcal{T}_{ij}^{Geo} K_{harm} (p_{w,i} - p_{w,j}) }_{\text{正交压力通量}} + \underbrace{ J_{w,ij}^{Cross} }_{\text{交叉修正}} - \underbrace{ \mathcal{T}_{ij}^{Geo} K_{harm} \rho_{w}^{avg} g (z_i - z_j) }_{\text{正交重力通量}} - \underbrace{ G_{w,ij}^{Cross} }_{\text{交叉重力修正}} \right]$$

- **交叉重力修正 $G_{w,ij}^{Cross}$**:

  $$G_{w,ij}^{Cross} = \rho_w^{avg} (\mathbf{K}_f \cdot \mathbf{g}) \cdot \mathbf{T}_{ij}$$

  *(注：$\mathbf{g} = (0,0, -9.8)$ 为常向量，直接计算投影，**切勿**通过梯度重构计算)*

#### **方程 2: CO₂相残差 $R_{n,i}$**

$$R_{n,i} = \frac{V_i}{\Delta t} \left[ (\phi \rho_n S_n)^{n+1} - (\phi \rho_n S_n)^n \right]_i + \sum_{j \in Faces} F_{n,ij}^{n+1} - Q_{n,i}^{n+1}$$

**CO₂相通量 $F_{n,ij}$ 表达式** (含毛管力)：

$$\begin{aligned} F\_{n,ij} &= \lambda\_n^{up} \rho\_n^{up} \bigg[ \\ &\quad \underbrace{ \mathcal{T}*{ij}^{Geo} K*{harm} (p\_{w,i} - p\_{w,j}) }*{\text{正交压力项}} + \underbrace{ \mathcal{T}*{ij}^{Geo} K\_{harm} (P\_{c,i} - P\_{c,j}) }*{\text{正交毛细项}} \\ &\quad - \underbrace{ \mathcal{T}*{ij}^{Geo} K\_{harm} \rho\_{n}^{avg} g (z\_i - z\_j) }*{\text{正交重力项}} + \underbrace{ (J*{n,ij}^{Cross\_P} + J\_{n,ij}^{Cross\_Pc} - G\_{n,ij}^{Cross}) }\_{\text{总交叉修正}} \bigg] \end{aligned}$$

- **Jacobian 组装提示**: $\frac{\partial F_n}{\partial S_n}$ 必须包含 $\frac{\partial P_c}{\partial S_n}$ 的链式法则贡献。

### 3.2 能量守恒方程

#### **方程 3: 能量残差 $R_{e,i}$**

$$R_{e,i} = \frac{V_i}{\Delta t} \left[ E_{acc, i}^{n+1} - E_{acc, i}^n \right] + \sum_{j \in Faces} (H_{adv, ij}^{n+1} + H_{cond, ij}^{n+1}) - Q_{E,i}^{n+1}$$

- **累积项 $E_{acc}$** (使用内能 $U$):

  $$E_{acc} = (1-\phi)\rho_r C_{p,r} T + \phi (\rho_w S_w U_w + \rho_n S_n U_n)$$

  *注：$U = h - p/\rho$，需严格区分内能与焓。*

- **热对流 $H_{adv, ij}$** (使用焓 $h$):

  $$H_{adv, ij} = (h_w)^{up} F_{w, ij} + (h_n)^{up} F_{n, ij}$$

  *注：直接复用质量方程计算出的质量通量 $F$，以保证守恒性。*

- **热传导 $H_{cond, ij}$**:

  $$H_{cond, ij} = \mathcal{T}_{ij}^{Geo} \lambda_{eff, harm} (T_i - T_j) + \underbrace{\lambda_{eff, f} (\overline{\nabla T}_f \cdot \mathbf{T}_{ij})}_{\text{热交叉修正}}$$

------

## 第四章 裂缝系统与 EDFM 耦合 (Fracture & Coupling)

裂缝作为独立的一组守恒方程求解，并通过 NNC 与基岩耦合。

### 4.1 裂缝控制方程

形式与基岩完全一致，但有以下简化：

- **无交叉项**：通常假设裂缝面内网格正交 ($\mathbf{T}_{ij} \approx 0$)。
- **重力必留**：裂缝在 3D 空间倾斜，必须包含 $\rho g \Delta z$。
- **孔隙度修正**：基岩需扣除裂缝体积 $\phi_{m}^{corr} = \phi_m (1 - V_{frac}/V_{cell})$。

### 4.2 裂缝-裂缝通量 (FF Flux)

$$F_{\alpha, il}^{FF} = \lambda_{\alpha}^{up} \rho_{\alpha}^{up} \cdot \mathcal{T}_{il}^{FF} \cdot \left[ \Delta p_{\alpha} - \rho_{\alpha}^{avg} g (z_i - z_l) \right]$$

- **几何传导率** (串联电阻模型):

  $$\mathcal{T}_{il}^{FF} = \frac{L_{seg}}{\frac{d_i}{w_i K_{f,i}} + \frac{d_l}{w_l K_{f,l}}}$$

### 4.3 基岩-裂缝耦合 (NNC Flux)

作为源汇项添加到双方的残差方程中。

$$F_{\alpha, ik}^{NNC} = \lambda_{\alpha}^{up} \rho_{\alpha}^{up} \cdot T_{ik}^{NNC} \cdot \left[ (p_{\alpha, i} - p_{\alpha, k}) - \rho_{\alpha}^{avg} g (z_i - z_k) \right]$$

- **NNC 传导率**:

  $$T_{ik}^{NNC} = \frac{A_{ik}}{\frac{d_{i}}{\mathbf{n} \cdot \mathbf{K}_m \cdot \mathbf{n}} + \frac{w_k}{2 K_{f,k}}}$$

- **全局矩阵贡献**：

  该项会产生非对角块 $A_{MF}$ (Matrix方程对Fracture变量的导数) 和 $A_{FM}$。

------

## 第五章 工程落地执行方案 (Implementation Strategy)

为解决从“系数型组装”转向“FIM 组装”的难题，采用 **自定义轻量级自动微分 (AD)** 技术路线。

### 5.1 核心技术：Forward-Mode AD

不依赖庞大的第三方库，自研轻量级模板类 `ADVar`。

- **数据结构**:

  C++

  ```
  template<int N> // N=3 for 2-phase, N=2 for 1-phase
  struct ADVar {
      double val;
      Eigen::Matrix<double, N, 1> grad; // 梯度 [d_dp, d_dS, d_dT]
      // 重载 operator+, -, *, /, exp, ...
  };
  ```

- **优势**: 自动处理链式法则（如 $\rho(p,T)$ 的导数），避免手写导数出错，且对 FIM 的短向量效率极高。

### 5.2 实施阶段规划

#### **阶段一：基础设施建设 (Infrastructure)**

- **任务**: 实现 `ADVar` 类及其数学运算符重载。
- **标准**: 单元测试通过 $f(x) = x^2 + \sin(xy)$ 的导数验证，误差 $< 10^{-14}$。

#### **阶段二：物性库升级 (Property Lib)**

- **任务**: 将 scCO2 (Span & Wagner) 和 Water (IAPWS) 状态方程模板化。

  C++

  ```
  template<typename T> T CalcRho(T p, T Temp);
  ```

- **标准**: 确保在临界点附近的导数连续，无 NaN。

#### **阶段三：组装器重构 (Jacobian Assembler)**

- **任务**: 编写 `BuildJacobian` 函数。
- **逻辑**:
  1. **输入**: 状态向量 $\mathbf{X}^k$。
  2. **累积项循环**: 遍历所有 Cell，使用 `ADVar` 计算 $Mass^{n+1}$，自动得到 $\frac{\partial M}{\partial X}$ 填充对角块。
  3. **通量项循环**: 遍历 Faces/NNC/FF，使用 `ADVar` 计算 Flux。
     - Flux 对 $X_i$ 的导数 $\to$ 填充 $(i,i)$ 和 $(j,i)$ 块。
     - Flux 对 $X_j$ 的导数 $\to$ 填充 $(i,j)$ 和 $(j,j)$ 块。
- **标准**: **静水平衡测试**。在重力场下，静止流体的残差应严格为 0（机器精度）。

#### **阶段四：求解器集成 (Solver Loop)**

- **任务**: 实现牛顿迭代循环。

  C++

  ```
  while (norm(R) > tol) {
      BuildJacobian(X, J, R);
      solve(J * dX = -R);
      X = X + dX; // 含限幅
  }
  ```

- **标准**: 单相热流耦合算例与旧代码结果一致；两相算例前缘推进正常。

### 5.3 资源与目标

- **资源**: C++17, Eigen 3.4+。
- **目标**: 形成一套支持全隐式、EDFM、强非线性热流耦合的通用求解器内核。

------

## 第六章 特别说明与注意事项

1. **IMPES 作为特例**：
   - 虽然代码架构基于 FIM，但可通过在 Jacobian 组装时强制忽略 $S_n$ 的非对角项（即 $\frac{\partial R}{\partial S} = 0$ 且 $\frac{\partial F}{\partial S} = 0$）来退化为 IMPES 模式。
   - **警告**：EDFM 中裂缝体积极小，IMPES 极易导致稳定性问题，仅用于调试或对比。
2. **重力项的陷阱**：
   - 在计算重力通量 $\rho \mathbf{g} \cdot \mathbf{A}$ 时，$\rho$ **必须**是 `ADVar` 类型的变量（即包含 $\frac{\partial \rho}{\partial p}$），否则 Jacobian 矩阵会缺失重力对压力的敏感度，导致迭代不收敛。
3. **交叉项的显隐性**：
   - 由于交叉项包含 $\nabla \Psi$（涉及远处邻居），为保持 Jacobian 稀疏性（仅包含一级邻居），交叉项通常处理为 **RHS Source**（显式或延迟修正），不进入 Jacobian 的左端项。
4. **调试建议**：
   - 实现一个 `NumericalJacobian` 函数（有限差分法），用于在单元测试中验证 `ADVar` 计算出的解析 Jacobian 是否正确。这是排查公式错误的终极手段。