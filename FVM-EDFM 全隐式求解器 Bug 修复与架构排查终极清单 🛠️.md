- - # FVM-EDFM 全隐式求解器 Bug 修复与架构排查终极施工图纸 🛠️

    > **文档说明**：本文档详细剖析了当前 FVM-EDFM 求解器中存在的 3 大致命 Bug 与 28 大深层架构雷区。每个条目均包含**问题根因**、**文件位置**与**修复动作**，以保障 2D/3D 单相与多相流的通用性、极致稳定性以及 HPC 高性能计算的极致效率。

    ---

    ## 一、 控制策略与时间步校准 (解决长时“锯齿失速”)
    *🔍 **根因**：大时间步（跨越数天）下，物理量演化剧烈，若沿用短时严格的收敛标准和激进的步长放大，会频繁触发牛顿法截断和 0.5x 步长砍半，形成“死循环”。*
    *📁 **位置**：`Test_Day6_TransientSolver.cpp` (及其他 Transient 测试文件 Setup 阶段)*

    - [ ] **1. 放宽牛顿迭代上限**
      - **修复动作**：将长时阶段的 `max_newton_iter` 从 12 放宽至 24。给强非线性方程组留足迭代衰减的空间。
    - [ ] **2. 适度松绑相对更新容差**
      - **修复动作**：将长时阶段的 `rel_update_tol` 设定为 `1.0e-3`。防止后期绝对压力极大时，微小的数值震荡引发假性发散。
    - [ ] **3. 抑制激进的步长增长**
      - **修复动作**：将 `dt_relres_grow_factor` 从 1.20 降至 `1.05`。避免步子迈得太大超出牛顿法局部收敛半径。
    - [ ] **4. 引入“软着陆”回退机制**
      - **修复动作**：增加 `dt_relres_soft_shrink_factor = 0.90`。在临界收敛边缘时，仅缩小 10% 步长继续试探，替代底层强硬的 `0.5x` 折半惩罚。
    - [ ] **5. 开启底层防锯齿引擎**
      - **修复动作**：强制设定 `enable_best_iter_guard = true` 与 `enable_stagnation_accept = true`，允许接受停滞但已满足绝对容差的解。

    ---

    ## 二、 核心致命 Bug 修复 (阻断雅可比撕裂)
    *📁 **位置**：全部集中在 `FIM_TransientEngine.hpp`*

    - [ ] **Bug 1：消除岩石“绝对刚体”陷阱（物理属性硬编码）**
      - **问题根因**：残差组装中写死了 `const double phi = 0.2`。因为常数没有压力梯度（$\partial \phi / \partial P = 0$），导致封闭边界注入流体时，岩石孔隙无法膨胀缓冲，压力瞬间无穷大。
      - **修复动作**：删除常数声明，动态从 `FieldManager` 读取 `phi, cp_r, rho_r`。必须引入岩石微压缩系数 $C_r$，构建带 AD 压力的动态孔隙度：`phi_ad = phi_ref * (1.0 + cr * (P - P_ref))`。
    - [ ] **Bug 2：大时间步 Limiter（阻尼器）无限放大失效**
      - **问题根因**：限制牛顿更新步幅的 `damp_scale = dt_eff / dt_ref`。长时阶段 $dt$ 增大上万倍，导致允许的单步压力更新上限（max_dP）暴涨到几十万大气压，完全丧失防爆功能。
      - **修复动作**：强制增加 1.0 的安全封顶：`const double damp_scale = std::min(1.0, std::max(1.0e-12, dt_eff / dt_ref));`。
    - [ ] **Bug 3：多相流 PTC（伪瞬态连续法）正则化符号反转**
      - **问题根因**：气相质量对水饱和度 $S_w$ 的偏导数为**负数**。原代码在增强对角线占优时，无脑加上了绝对正数，导致原有的负数被抵消归零，矩阵瞬间奇异 (`dx_nan_inf`)。
      - **修复动作**：提取原本累积项导数的符号：`double sign_acc = (d_acc < 0.0) ? -1.0 : 1.0;`，然后在对角线加上 `ptc * m_ptc * sign_acc`。

    ---

    ## 三、 深层架构级雷区 (保障多相流与热力学一致性)

    ### ⚠️ A. 物性与相渗评估雷区 (Thermodynamics & PVT)
    - [ ] **雷区 1：相消失与 AD 奇异性陷阱**
      - **问题根因**：盲目的牛顿试探步会让 $S_w$ 变成负数。直接代入相渗幂函数 `pow()` 会产生 `NaN` 并顺着自动微分（AD）链条污染全矩阵。
      - **文件位置**：`BoundaryAssembler.cpp` & `AD_FluidEvaluator.h`
      - **修复动作**：在调用相渗和毛管力前，对 $S_w$ 进行**保留梯度的硬截断**：`if(Sw.val < 1e-8) Sw.val = 1e-8;`。
    - [ ] **雷区 6：毛管力（Capillary Pressure）缺位**
      - **问题根因**：直接用网格主变量 $P_{cell}$ 同时计算水和气的物性。忽略了在致密基质中 $P_c$ 高达数兆帕，导致气相密度与热力学状态计算全错。
      - **文件位置**：任意调用 `evaluateCO2` 的地方
      - **修复动作**：计算超临界 CO₂ 前，必须修正气体压力：`ADVar<N> P_gas = P_cell + Pc(Sw_cell);`。
    - [ ] **雷区 12：异常兜底的“零梯度死亡”**
      - **问题根因**：`AD_FluidEvaluator` 中 `catch` 异常后返回常数物性。因为常数导数为 0，雅可比矩阵瞬间丢失方向，牛顿法无法回调至物理区。
      - **文件位置**：`AD_FluidEvaluator.h` (异常处理块)
      - **修复动作**：赋予微小的虚拟梯度作为“指南针”，例如 `dRho_dP = 1e-6`，`dRho_dT = -1e-3`。

    ### ⚠️ B. 井模型雷区 (Well Model)
    *📁 **位置**：`BoundaryAssembler.cpp` 的 `Assemble_Wells_2D/3D_FullJac`*

    - [ ] **雷区 2：BHP 注气井的“迎风物理闭锁”**
      - **问题根因**：注气井错误使用了当前地层网格的流度。若地层全是水（$S_w=1$），气相流度为 0，导致注气量永远为 0（物理死锁）。
      - **修复动作**：对于注入井，其等效流度必须使用总流度：`mob_g_eff = mob_w + mob_g`。
    - [ ] **雷区 4：定产井（Rate Control）隐式分流梯度丢失**
      - **问题根因**：使用 `MakeConstAD` 抹除了饱和度分流率（`fw`）的 AD 梯度。求解器“看不见” $S_w$ 变化对产水/气比例的影响，引发震荡。
      - **修复动作**：用 AD 变量动态计算：`ADVar<N> fw_ad = mob_w / (mob_w + mob_g);`。
    - [ ] **雷区 5：注入热焓的“虚假偏导数”污染**
      - **问题根因**：使用带有梯度的 `P_cell` 计算外部注入流体的热焓。这在雅可比矩阵中伪造出了“地层压力改变会影响井口注入焓值”的非物理耦合。
      - **修复动作**：剥离梯度：`ADVar<N> P_inj_const = MakeConstAD<N>(P_cell.val);`。

    ### ⚠️ C. 边界与迎风格式雷区 (Boundary & Upwind)
    - [ ] **雷区 3：迎风算子 C0 不连续性**
      - **根因/位置**：`FVM_Ops_AD.h`。势差在 0 附近微小跳变时，迎风直接翻转，导致牛顿极限环不收敛。
      - **修复动作**：热焓迎风**必须严格跟随**对应相质量通量方向。未来引入 Logistic 平滑插值。
    - [ ] **雷区 7：外部边界组装的“标量解耦”**
      - **根因/位置**：`BoundaryAssembler.cpp` (`Assemble_3D`)。按单一自由度轮询组装，彻底丢失了热对流（能量方程）对压力变化的交叉导数（Off-diagonal term）。
      - **修复动作**：重构为一次性处理 $P, Sw, T$ 的 3x3 完整 Jacobian 组装。
    - [ ] **雷区 13：重力势计算的“混合密度”悖论**
      - **根因/位置**：`FVM_Ops_AD.h`。将水气密度平均后算重力，彻底摧毁了超临界 CO₂ 向上运移的纯粹浮力。
      - **修复动作**：水相势差和气相势差必须独立使用各自的相密度计算重力压降。
    - [ ] **雷区 15/20：边界物理通量的“裸奔”**
      - **根因/位置**：`BoundaryAssembler.cpp`。边界通量只用了纯几何 $A/d$。对于达西流动缺失了渗透率 $k$ 和流度 $\lambda$；对于热传导缺失了导热系数 $\lambda_{eff}$。
      - **修复动作**：乘上真实物理量：`T_trans = MakeConstAD(T_geom * k_face) * mob_fluid`。

    ### ⚠️ D. EDFM 专有物理与拓扑雷区 (EDFM Specifics)
    - [ ] **雷区 8：裂缝与基质的相渗模型混用**
      - **根因/位置**：`BoundaryAssembler.cpp` 及内部通量。裂缝被错误套用致密基质的 vG 曲线，CO₂ 流动被严重阻力锁定。
      - **修复动作**：判断 `domain`。裂缝网格强制使用 $X$ 型线性相渗（$k_{rw}=S_w, k_{rg}=1-S_w$）且 $P_c=0$。
    - [ ] **雷区 9：裂缝井的“静默失效”**
      - **根因/位置**：`BoundaryAssembler.cpp`。若未给定裂缝井指数 WI，代码直接 `continue` 跳过，不报错但全盘算错。
      - **修复动作**：将 `continue` 改为抛出异常 `throw std::runtime_error("Fracture WI missing!");`。
    - [ ] **雷区 10：NNC 动态插入带来的灾难级降速**
      - **根因/位置**：`FIM_TopologyBuilder3D.h`。稀疏矩阵 `reserve` 漏掉了基质-裂缝 NNC，导致 Eigen 在迭代中疯狂进行 $O(N^2)$ 的内存重分配，极度拖慢速度。
      - **修复动作**：遍历 `fracture_network()`，将所有 NNC 连接对强制加入拓扑邻接表预分配内存。
    - [ ] **雷区 14：各向异性渗透率在 NNC 上的丢失**
      - **根因/位置**：`TransmissibilitySolver_3D.cpp`。只用标量平均渗透率，导致垂直裂缝算出的基质漏失量完全失真。
      - **修复动作**：使用渗透率张量投影公式：$k_{NNC} = \vec{n}^T \mathbf{K} \vec{n}$。
    - [ ] **雷区 19：裂缝开放边界的“盲区”**
      - **根因/位置**：`BoundaryAssembler.cpp`。末尾对于裂缝元素的边界处理仅写了 `(void)elem;` 空循环。裂缝触及储层边界时等同于“撞墙”。
      - **修复动作**：补全裂缝端面（Faces/Edges）通向外部的 Dirichlet/Neumann 物理组装代码。

    ### ⚠️ E. 框架体系结构雷区 (Framework Architecture)
    - [ ] **雷区 11：全隐式装配接口的“硬编码维度”**
      - **根因/位置**：`BoundaryAssembler.cpp`。函数参数强行写死 `<3>` (即 `std::array<double, 3>`)。
      - **修复动作**：重构为模板泛型 `template<int N>`。保障单相流 ($N=2$) 与多相流 ($N=3$) 无缝切换。
    - [ ] **雷区 16：跨网格迎风导数丢失**
      - **根因/位置**：通量底层架构。`ADVar<3>` 只能存本地梯度，无法将“下游压力导致上游流度改变”的偏导数传递给邻居矩阵块。
      - **修复动作**：引入**双路 AD 计算（Double-Pass）**或升级为 `ADVar<6>`，拯救被削弱的牛顿法二次收敛性。
    - [ ] **雷区 17：状态方程 (EOS) 冗余计算拖垮速度**
      - **根因/位置**：`FIM_TransientEngine.hpp`。同一个网格的 CO₂ 状态方程在一轮牛顿迭代内被井、边界、内部面反复求值。
      - **修复动作**：在迭代最外层建立全局网格的 `Properties Cache`，内部组装直接读取缓存。
    - [ ] **雷区 18：自由度 (DOF) 错位导致单相流段错误崩溃**
      - **根因/位置**：`BoundaryAssembler.cpp`。硬编码了 `T_cell.grad(2) = 1.0`。但在纯水单相热流中只有 2 个变量，访问 index 2 会直接导致 Segfault 内存越界。
      - **修复动作**：利用参数动态赋梯：`if(dofOffset_E >= 0) T_cell.grad(dofOffset_E) = 1.0;`。

    ### ⚠️ F. 硬编码与热力学一致性 (Hardcoding & Thermodynamics)
    *📁 **位置**：`FIM_TransientEngine.hpp` & `FVM_Ops_AD.h`*

    - [ ] **雷区 21：岩石热力学参考基准错位**
      - **问题根因**：计算岩石内能直接用绝对温度 $T$，相当于参考点为 0 K；但流体 EOS 焓值的参考点通常是 273.15 K。基准错位引发虚假巨量热传递。
      - **修复动作**：对齐参考点：$e_{rock} \propto C_r \cdot (T - T_{ref})$。
    - [ ] **雷区 22：初始压力的“标量硬编码”**
      - **问题根因**：计算压缩性孔隙度时硬编码了参考压（如 1 atm），导致深部处于高静水压的网格在第 0 步发生非物理的剧烈膨胀。
      - **修复动作**：从初始场精确读取本地参考压：`P_init[bi]`。
    - [ ] **雷区 23：重力矢量的降维锁死**
      - **问题根因**：硬编码 `g_z = 9.81`。当模型切换为 2D Y-剖面时，Z轴重力导致气体完全丧失浮力。
      - **修复动作**：使用参数配置 `dot_product(params.gravity_vector, dx_vector)`。
    - [ ] **雷区 24：井体积流量未转化标准态密度**
      - **问题根因**：用户输入的 Rate 通常是“地面标准体积”。如果直接乘以地下流体密度或硬编码密度（如 1.8），算出的注入质量完全失真。
      - **修复动作**：调用 `evaluateStdCondition()` 计算精准的标准态密度以换算真实质量流率。

    ### ⚠️ G. 线性求解器与 HPC 高效计算雷区 (Linear Solvers & HPC)
    *📁 **位置**：`FIM_TransientEngine.hpp` & `FIM_BlockSparseMatrix.h`*

    - [ ] **雷区 25：非零元动态剔除导致 SparseLU 符号分解崩溃 (Dynamic Sparsity Crash)**
      - **问题根因**：在 `ExportEigenSparseMatrix()` 中通过 `if (std::abs(val) > 1e-16)` 动态丢弃微小元素。这会导致每次牛顿迭代生成的 CSR 结构随机改变，引发 `SparseLU::factorize()` 指针错乱，导致段错误（Segfault）崩溃，彻底摧毁复用机制。
      - **修复动作**：删除数值大小判断，无条件向 `triplets` 推入最初分配好的所有元素（即使值为 0），确保矩阵拓扑绝对静止。
    - [ ] **雷区 26：重建稀疏矩阵的 $O(N \log N)$ 性能雪崩 (Triplet Reassembly Penalty)**
      - **问题根因**：每次牛顿迭代调用 `ExportEigenSparseMatrix()` 时，都在分配数百万个 `Triplet` 并在 CPU 上进行极耗时的 `std::sort` 排序。
      - **修复动作**：在 `FIM_BlockSparseMatrix` 中引入缓存矩阵 `cached_mat_`，仅在 `FreezePattern()` 时组装一次 CSR。牛顿迭代中直接遍历 `cached_mat_.valuePtr()` 核心数组进行 **$O(N)$ 原位更新（In-place Update）**。
    - [ ] **雷区 27：BiCGSTAB 与 ILUT 的循环内内存搅动 (Memory Thrashing)**
      - **问题根因**：`Eigen::BiCGSTAB` 与预处理器在 `for` 循环内部被反复实例化。不仅引发高频的系统级 GB 级内存申请与释放，而且在执行行缩放时 `A_solve = S * A` 引发了整张大矩阵的深拷贝。
      - **修复动作**：将 `bicgstab_solver` 的声明提升到牛顿循环甚至时间步外部。用 `InnerIterator` 在原矩阵 `A` 上进行**原地行缩放**，彻底实现牛顿循环内零内存分配。
    - [ ] **雷区 28：FIM 两相流代数多重网格的架构盲区 (Missing CPR-AMG)**
      - **问题根因**：直接用标准 AMG 求解 FIM 的全耦合矩阵必然发散，且现有的 SparseLU/BiCGSTAB 无法抗住含 EDFM 的百万级强刚性网格。
      - **修复动作**：在引擎演进规划中引入 AMGCL，并架构 **CPR (Constrained Pressure Residual)** 两阶段预处理器：第一阶用 AMG 求解纯压力解耦块获取全局压力增量，第二阶用 ILU0/ILUT 求解全耦合系统。