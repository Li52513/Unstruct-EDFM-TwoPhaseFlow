# EDFM-FVM-FIM 框架技术全景文档

> **基于代码版本 git HEAD `54e2ec9`，生成日期 2026-03-16**
> 工作目录：`2D-Unstr-Quadrilateral-EDFM/`

---

## 0. 总体定位

本工程是一套基于**非结构化四边形网格**的 EDFM（嵌入式离散裂缝模型）数值模拟器，支持：

| 模式 | 模板参数 | 求解变量 | 应用场景 |
|------|---------|---------|---------|
| **N=2** 单相+热 | `<2, MeshMgr, FieldMgr>` | P（压力）、T（温度） | 超临界CO₂单相渗流传热 |
| **N=3** 两相+热 | `<3, MeshMgr, FieldMgr>` | P（水压）、Sᵥ（水饱和度）、T（温度） | 超临界CO₂-水两相渗流传热 |

几何维度通过 **MeshMgr** 模板参数区分：`MeshManager`（2D）vs `MeshManager_3D`（3D）。

---

## 1. 方程体系与理论基础

### 1.1 N=2：单相流体 + 能量方程

**守恒形式（有限体积积分）**：

$$\frac{V_i}{\Delta t}(\rho\phi - \rho^n\phi^n) + \sum_{j\in \mathcal{N}(i)} F^m_{ij} - Q^m_i = 0 \quad \text{(质量)}$$

$$\frac{V_i}{\Delta t}(e - e^n) + \sum_{j\in \mathcal{N}(i)} F^e_{ij} - Q^e_i = 0 \quad \text{(能量)}$$

**变量定义**：
- 孔隙度（含压缩性）：$\phi = \phi_\text{ref}(1 + c_r(P - P_\text{ref}))$
- 内能密度：$e = \rho(h - P/\rho) + (1-\phi)\rho_r c_{pr} T$（流体焓 + 岩石热容）
- 质量通量：$F^m_{ij} = T_\text{flow} \cdot \frac{\rho}{\mu}\bigg|_\text{upwind} \cdot \Delta\Phi_{ij}$
- 热通量：$F^e_{ij} = T_\text{heat}(T_i - T_j) + F^m_{ij} \cdot h\big|_\text{upwind}$（传导 + 对流）
- 势差（含重力）：$\Delta\Phi_{ij} = P_i - P_j - \bar\rho\, g\,\Delta z$

### 1.2 N=3：水+CO₂两相 + 能量方程

**守恒形式**：

$$\frac{V_i}{\Delta t}(\rho_w\phi S_w - \rho_w^n\phi^n S_w^n) + \sum F^w_{ij} - Q^w_i = 0 \quad \text{(水相质量)}$$

$$\frac{V_i}{\Delta t}(\rho_g\phi S_g - \rho_g^n\phi^n S_g^n) + \sum F^g_{ij} - Q^g_i = 0 \quad \text{(CO₂质量)}$$

$$\frac{V_i}{\Delta t}\bigl[(e_\text{fluid}\phi + e_\text{rock}) - (\cdots)^n\bigr] + \sum F^e_{ij} - Q^e_i = 0 \quad \text{(能量)}$$

**附加关系**：
- $S_w + S_g = 1$
- $P_g = P_w + P_c(S_w)$ — van Genuchten 毛细压力
- $k_{rw}(S_w),\ k_{rg}(S_g)$ — Mualem-van Genuchten 相对渗透率
- 各相通量：$F^w_{ij} = T_\text{flow} \cdot \frac{k_{rw}\rho_w}{\mu_w}\big|_\text{up} \cdot \Delta\Phi_w$，CO₂同理（用 $P_g$ 势差）

### 1.3 流体物性方程（EOS）

| 流体 | EOS | 适用范围 | 代码位置 |
|------|-----|---------|---------|
| 水 | **IAPWS-95** 工业水蒸气方程 | 0°C～800°C，0.1～1000 MPa | `2D/3D_Water_Properties` |
| 超临界CO₂ | **Span-Wagner** 方程 | 超临界区（>7.38 MPa, >31.1°C） | `2D/3D_CO2_Properties` |

**AD 自动微分封装**（`AD_FluidEvaluator.h`）：`evaluateWater<N>(P_ad, T_ad)` 和 `evaluateCO2<N>(P_ad, T_ad)` 返回 `ADFluidProperties<N>` —— 包含 ρ、μ、h、cp、λ 及其关于 P、T 的自动微分导数（4层降级策略：中央差→前向→后向→零梯度保底）。

### 1.4 固体物性

| 物性 | 基岩 | 裂缝 |
|------|------|------|
| 孔隙度 | AABB 分层初始化（`2D_RockSolidProperties`） | 按裂缝或单元赋值（`2D_FractureProperties`） |
| 渗透率 | 各向异性 Kxx/Kyy（/Kzz 3D） | 切向 Kₜ + 法向 Kₙ |
| 压缩系数 cᵣ | 可配置（默认 0） | 可配置 |
| 热导率 λ | 体积混合律 | 体积混合律 |

---

## 2. EDFM 离散框架

### 2.1 连接类型（`FIM_ConnectionManager.h:17`）

```
MM (Matrix-Matrix)       = 0   基岩-基岩相邻面
FI (Fracture-Internal)   = 1   裂缝段内部相邻
NNC (Matrix-Fracture)    = 2   裂缝-基岩非邻连（EDFM核心）
FF (Fracture-Fracture)   = 3   裂缝-裂缝交叉
```

全局求解器索引：`[0, nMat-1]` 基岩，`[nMat, nMat+nFrac-1]` 裂缝段。

### 2.2 导流能力计算（静态预计算，一次性）

核心公式（串联阻力，`FVM_Ops.h:82`）：

$$T = \frac{A}{d_1/K_1 + d_2/K_2}$$

**2D EDFM 各连接类型**：

| 类型 | 面积 A | 距离 d | 渗透率 K |
|------|--------|--------|---------|
| MM | 共享面长×厚度 | 半网格距 | Kₓₓ/Kyy投影 |
| FI | 厚度 h=1.0 | 段半长 L/2 | Kₜ×Wf |
| **NNC** | `L_seg × h` | `avgDistance`（平均垂直距） | `nx²Kxx+ny²Kyy`（法向投影） |
| **FF** | Star-Delta：`G_i×G_j/ΣG_k` | `L/2`（段半长） | Kₜ×Wf |

**3D EDFM 差异**：

| 类型 | 2D | 3D |
|------|----|----|
| NNC面积 | L×h（线×厚） | `A_int`（真实交互多边形面积） |
| NNC距离 | `avgDistance` | `distMatrixToFracPlane`（精确法向） |
| 渗透率投影 | nx²Kxx+ny²Kyy | nx²Kxx+ny²Kyy+nz²Kzz |
| FF聚类 | GlobalFFPoint整数ID | 量化弦KEY（坐标量化，TOL=1e-4） |
| FF距离 | 段半长 | PointToSegmentDistance（质心到交线） |
| FF面积权重 | 交线长=段长 | cluster.length（3D交线长度） |

### 2.3 并联聚合（`FIM_ConnectionManager::FinalizeAndAggregate`）

- **NNC / FF**：同一 `(nodeI, nodeJ, type)` 多路径累加：$T_\text{total}=\sum_k T_k$，`aux_dist` 按面积加权平均
- **MM / FI**：重复则抛出异常（唯一界面约束）
- 排序输出：`MM → FI → NNC → FF`，每类内按 `(nodeI, nodeJ)` 字典序

---

## 3. 方程组装（FIM 全隐式）

### 3.1 ADVar<N> 自动微分机制（`ADVar.hpp`）

```
ADVar<2>: val + grad[2]   → {∂/∂P, ∂/∂T}
ADVar<3>: val + grad[3]   → {∂/∂P, ∂/∂Sw, ∂/∂T}
```

在 Newton 迭代内，`ADVar` 的 `.val` 给出残差，`.grad[j]` 给出 Jacobian 列，**无需有限差分**。

### 3.2 统一通量循环（`RunGeneric_impl.hpp:696`）

所有连接类型（MM/FI/NNC/FF）共用同一段代码，差异仅在预计算的 T_Flow/T_Heat 值：

```cpp
for (const auto& conn : connMgr.GetConnections()) {
    auto f_wrt_i = evalFlux(seed_i);   // ∂F/∂u_i
    auto f_wrt_j = evalFlux(seed_j);   // ∂F/∂u_j
    FIM_GlobalAssembler<N>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, global_mat);
}
```

### 3.3 Jacobian 块写入（`FIM_GlobalAssembler.h`）

每条连接贡献 4 个矩阵块：

| 位置 | 值 | 物理含义 |
|------|----|---------|
| `R[i] += F.val` | 节点i残差 | i流出 |
| `R[j] -= F.val` | 节点j残差 | j流入 |
| `J[i,i] += ∂F/∂u_i` | 对角块 | 自耦 |
| `J[j,j] -= ∂F/∂u_j` | 对角块 | 自耦 |
| `J[i,j] += ∂F/∂u_j` | **非对角块** | j→i耦合（NNC独特贡献） |
| `J[j,i] -= ∂F/∂u_i` | **非对角块** | i→j耦合 |

NNC连接的非对角块直接将基岩-裂缝**完全隐式耦合**，是FIM相对IMPES的核心优势。

### 3.4 块稀疏矩阵（`FIM_BlockSparseMatrix.h`）

```
diag_blocks_[i]             : N×N Eigen对齐矩阵（对角块）
off_diag_blocks_[i][j]      : N×N，unordered_map<int, BlockMat>（稀疏）
residual_[i]                : N×1
```

`FreezePattern()` 在初始化后冻结稀疏结构，后续迭代只填值不改结构。

### 3.5 积累项（含岩石压缩性）

**N=2**（`RunGeneric_impl.hpp:627`）：
- $R_0^\text{acc} = V_i/\Delta t \cdot (\rho_w\phi - \rho_w^n\phi^n)$
- $R_1^\text{acc} = V_i/\Delta t \cdot (e_w - e_w^n)$

**N=3**（`RunGeneric_impl.hpp:647`）：额外引入 $S_w$ 积累项

岩石孔隙度含压缩性：$\phi = \phi_\text{ref}(1+c_r(P-P_\text{ref}))$

### 3.6 边界条件（`BoundaryAssembler.cpp`）

| 类型 | 公式 | 用途 |
|------|------|------|
| Dirichlet | $q = T_\text{geom}(c/a - u_\text{cell})$ | 定压/定温 |
| Neumann | $q = L_\text{face} \cdot c$ | 定通量 |
| Robin | $q = L \cdot C_L(u_\text{cell} - u_\infty)$ | 漏失项 |

边界量组装到 `AssembleSource`（负号：$R_i -= q$）。

### 3.7 状态同步（`FIM_TransientEngine/StateSync.hpp`）

`SyncStateToFieldManager<N>()` 将 `FIM_StateMap(P,T,Sw)` → FieldManager 标量场（ρ、μ、h、krw、krg…），在 Newton 迭代内（井/边界组装前）、VTK输出前调用。

---

## 4. Newton 迭代与线性求解

### 4.1 时间步 + Newton 外层结构

```
for step = 1..max_steps:
    old_state = state
    select_profile (startup / long)
    apply_control_ramp

    for iter = 0..max_newton_iter:
        SetZero(); Accumulation; Flux; Wells; BC
        A = ExportEigenSparseMatrix()
        b = ExportEigenResidual()
        conv_res = ||b||∞

        if converged: break
        if enable_ptc: augment A diagonal
        dx = SolveLinearSystem(A, b)
        alpha = compute_initial_damping(dx)   // max_dP/max_dT/max_dSw 约束
        (state, alpha) = armijo_line_search(state, dx, alpha)

    if converged: t += dt; adapt_dt
    else: state = old_state; dt *= rollback_shrink; step--
```

### 4.2 线性求解器选项

| 求解器 | 适用场景 | 关键参数 |
|--------|---------|---------|
| **SparseLU**（Eigen） | 小/中规模，稳定 | 模式哈希缓存符号分解 |
| **AMGCL**（代数多重网格） | 大规模，迭代 | `amgcl_tol`, `amgcl_maxiter`；失败回退SparseLU |
| **BiCGSTAB** | 备选迭代 | `bicgstab_droptol` |

**行缩放**（`enable_row_scaling`）：缩放分母 = max(|D_acc[r]|, |A[r,r]|, floor)，防止压力/温度量级悬殊导致条件数恶化。

### 4.3 线搜索三级策略

1. **初始阻尼**：按 `max_dP/max_dT/max_dSw` 约束 α，含时间步感知缩放
2. **Armijo 条件**：充分下降 + 指数回退（`armijo_beta` 因子，最多 `armijo_max_backtracks` 次）
3. **非单调窗口**：参考残差取最近 `nonmonotone_window` 步的最大值，避免过保守

**救援机制**：
- **级别1 Fallback接受**：最优试探残差 ≤ 1.1×ref 时强制接受
- **级别2 受控接受**（iter=1）：放宽 `controlled_accept_relax` 倍后接受
- **级别3 PTC救援**：连续 `ls_fail_rescue_threshold` 次失败 → 增大PTC对角增强因子，重试

### 4.4 收敛判断（多模式）

| 模式 | 条件 | 触发 |
|------|------|------|
| `abs_res` | `||b||∞ < abs_res_tol` | 绝对残差足够小 |
| `rel_res_update` | `rel_res ≤ rel_res_tol` AND `rel_update ≤ rel_update_tol`（iter>1） | 双重准则 |
| `stagnation` | 残差已显著下降但未到容差，末次迭代 | 停滞接受（dt×0.85） |
| `best_iter_guard` | 从最优快照反弹 > `best_iter_growth_trigger` | 回退到最优状态 |

### 4.5 自适应时间步

| 收敛模式 | iter_used ≤ N₁ | N₁ < iter ≤ N₂ | iter > N₂ |
|---------|----------------|----------------|-----------|
| `abs_res` | dt × 1.2 | dt × 1.05 | dt × 0.8 |
| `rel_res_update` | dt × grow_factor | dt × neutral_factor | dt × shrink_factor |

回滚时：`dt *= rollback_shrink_factor`，`state = old_state`，重试。

---

## 5. 关键可调参数速查

### 5.1 时间步控制

| 参数 | 默认值 | 作用 |
|------|--------|------|
| `dt_init` | 1.0 s | 初始时间步 |
| `dt_min` | 1e-4 s | 回滚下限 |
| `dt_max` | 86400 s | 上限 |
| `max_steps` | 50 | 总步数上限 |
| `target_end_time_s` | -1（无限） | 目标终止时刻 |
| `rollback_shrink_factor` | 0.7 | 回滚缩减因子 |
| `enable_two_stage_profile` | false | 启动/长期双阶段切换 |

### 5.2 Newton 迭代

| 参数 | 默认值 | 作用 |
|------|--------|------|
| `max_newton_iter` | 8 | 每步最大迭代次数 |
| `abs_res_tol` | 1e-6 | 绝对残差容差 |
| `rel_res_tol` | 1e-3 | 相对残差容差 |
| `rel_update_tol` | 1e-6 | 相对更新量容差 |
| `max_dP` | 1e4 Pa | 单步最大压力变化（阻尼限制） |
| `max_dT` | 2 K | 单步最大温度变化 |
| `max_dSw` | 0.05 | 单步最大饱和度变化（N=3） |

### 5.3 PTC 伪瞬态连续

| 参数 | 默认值 | 作用 |
|------|--------|------|
| `enable_ptc` | false | 开关 |
| `ptc_lambda_init` | 1.0 | 初始对角增强系数 |
| `ptc_lambda_decay` | 0.5 | 每迭代衰减率 |
| `ptc_rescue_boost` | 5.0 | 救援时放大倍率 |

### 5.4 线性求解

| 参数 | 默认值 | 作用 |
|------|--------|------|
| `lin_solver` | AMGCL | SparseLU / AMGCL / BiCGSTAB |
| `enable_row_scaling` | true | 行缩放 |
| `amgcl_tol` | 1e-6 | AMGCL收敛容差 |
| `amgcl_use_fallback_sparselu` | true | AMGCL失败时回退 |

### 5.5 物理模型选项（`TransientOptionalModules`）

| 参数 | 作用 |
|------|------|
| `single_phase_fluid` | N=2 时的流体：Water 或 CO₂ |
| `vg_params` | van Genuchten 毛细压力参数（N=3） |
| `rp_params` | Mualem相对渗透率参数（N=3） |
| `property_initializer` | 用户自定义物性初始化函数 |
| `state_initializer` | 用户自定义初始条件函数 |

### 5.6 井控制（`WellScheduleStep`）

| 参数 | 含义 |
|------|------|
| `control_mode` | BHP（定井底压力） / Rate（定流量） |
| `target_value` | 目标值（Pa 或 kg/s） |
| `injection_temperature` | 注入流温度（>0时生效，K） |
| `wi_override` | 井指数手动覆盖（-1=自动） |
| `frac_w / frac_g` | 两相分流比（N=3） |

---

## 6. 关键模块文件索引

| 模块 | 文件 | 主要内容 |
|------|------|---------|
| **主循环** | `FIM_TransientEngine/RunGeneric_impl.hpp` | Newton迭代、积累/通量/井/边界组装 |
| **状态同步** | `FIM_TransientEngine/StateSync.hpp` | State→FieldManager物性同步 |
| **参数定义** | `FIM_TransientParams.hpp` | TransientSolverParams全字段 |
| **自动微分** | `ADVar.hpp` | N维对偶数（值+梯度） |
| **流体物性AD** | `AD_FluidEvaluator.h` | IAPWS-95 / Span-Wagner + 4层微分降级 |
| **CO₂物性** | `2D/3D_CO2_Properties.h/cpp` | Span-Wagner EOS |
| **水物性** | `2D/3D_Water_Properties.h/cpp` | IAPWS-95 EOS |
| **基岩物性** | `2D_RockSolidProperties.h/cpp` | AABB分层初始化 |
| **裂缝物性** | `2D_FractureProperties.h/cpp` | 裂缝参数（孔隙度、张开度、渗透率） |
| **导流能力2D** | `TransmissibilitySolver_2D.cpp` | MM/FI/NNC/FF 2D导流能力 |
| **导流能力3D** | `TransmissibilitySolver_3D.cpp` | MM/FI/NNC/FF 3D导流能力 |
| **连接管理** | `FIM_ConnectionManager.h` | Connection结构、FinalizeAndAggregate |
| **拓扑构建2D** | `FIM_TopologyBuilder2D.h` | _loadMatrix/FI/NNC/FF → connMgr |
| **拓扑构建3D** | `FIM_TopologyBuilder3D.h` | 3D版本，InteractionPair |
| **组装器** | `FIM_GlobalAssembler.h` | AssembleAccumulation/Flux/Source |
| **块稀疏矩阵** | `FIM_BlockSparseMatrix.h` | N×N块存储（Eigen对齐） |
| **FVM算子** | `FVM_Ops.h` | SeriesResistance、Transmissibility |
| **2D网格管理** | `MeshManager.h/cpp` | getNNCTopologyMap、全局索引分配 |
| **3D网格管理** | `3D_MeshManager.h/cpp` | InteractionPair、双向查询加速 |
| **边界条件** | `BoundaryAssembler.cpp` | Dirichlet/Neumann/Robin施加 |
| **井控制类型** | `Well_WellControlTypes.h` | WellScheduleStep、BHP/Rate控制 |

---

## 7. 2D vs 3D 差异汇总

| 维度 | 2D | 3D |
|------|----|----|
| NNC面积 | L_seg × h（h=1 m） | 真实交互多边形面积 A_int |
| NNC距离 | avgDistance（FractureElement字段） | distMatrixToFracPlane（精确法向） |
| 渗透率投影 | nx²Kxx + ny²Kyy | +nz²Kzz |
| NNC拓扑源 | `getNNCTopologyMap()` | `getInteractionPairs()` |
| NNC aux几何 | 占位符 1.0 | 真实面积和距离（用于并联加权） |
| FF聚类 | 整数 GlobalFFPoint ID | 量化弦KEY（浮点→整数，TOL=1e-4） |
| FF距离 | length × 0.5（段半长） | PointToSegmentDistance（质心到交线） |
| 裂缝段体积 | length × max(aperture, 1e-6) | max(area,1e-12) × max(aperture,1e-8) |
| FI连接存储 | FractureFace标量场 | FractureEdge标量场 |
| **通量组装循环** | **统一（与3D相同代码）** | **统一** |
| Jacobian结构 | 相同 | 相同 |

---

## 8. 典型执行序列（完整伪代码）

```
=== 初始化阶段（一次性）===
1. 初始化网格（MeshManager / MeshManager_3D）
2. 分配物性（RockSolidProps、FractureProps、CO₂Props、WaterProps）
3. 预计算导流能力（TransmissibilitySolver → T_MM/FI/NNC/FF）
4. 构建连接（FIM_TopologyBuilder → connMgr.PushConnection）
5. 并联聚合（connMgr.FinalizeAndAggregate）
6. 分配求解器（FIM_StateMap、FIM_BlockSparseMatrix）
7. FreezePattern → 冻结稀疏结构

=== 时间步循环 ===
for step = 1..max_steps:
    old_state = state
    select_profile (startup / long)
    apply_control_ramp (控制坡道)

    for iter = 0..max_newton_iter:
        global_mat.SetZero()

        // 1. 积累项：ADVar<N> 同时得到 R_acc + J_acc
        AssembleAccumulation(state, old_state, dt, phi, cr, ...)

        // 2. 通量项：统一循环（MM/FI/NNC/FF）
        SyncStateToFieldManager(state)  // P,T,Sw → ρ,μ,h,kr
        for conn in connMgr:
            evalFlux(seed_i) → f_wrt_i   // ∂F/∂u_i（ADVar双次seed）
            evalFlux(seed_j) → f_wrt_j   // ∂F/∂u_j
            AssembleFlux(i, j, f_wrt_i, f_wrt_j)

        // 3. 井源项（ADVar数值差分Jacobian）
        Assemble_Wells_2D/3D_FullJac(state)

        // 4. 边界项（Dirichlet/Neumann/Robin）
        Assemble_BoundaryConditions(state)

        A = ExportEigenSparseMatrix()
        b = ExportEigenResidual()
        conv_res = ||b||∞

        if conv_res < abs_res_tol: break (abs_res)
        if (iter>1) && rel_res < tol && rel_update < tol: break (rel_res_update)

        if enable_ptc: A_diag *= (1 + ptc_lambda)

        if enable_row_scaling: scale(A, b)
        dx = SolveLinearSystem(A, b)    // SparseLU / AMGCL / BiCGSTAB

        alpha = min(1, max_dP/|dx_P|, max_dT/|dx_T|, ...)  // 初始阻尼
        (state, alpha) = Armijo_LineSearch(state, dx, alpha)
        // 三级救援（Fallback→受控接受→PTC救援）

    if converged:
        t += dt
        adapt_dt (基于iter_used四段策略)
        SyncStateToFieldManager(state)
        OutputVTK (按间隔)
    else:
        state = old_state
        dt = max(dt × rollback_shrink, dt_min)
        step--   // 重试
```
