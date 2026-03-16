# FIM 隐式求解器完整流程说明

> 源文件根目录：`FIM_TransientEngine/`
> 主循环入口：`RunGeneric_impl.hpp` — `RunGenericFIMTransient<N, MeshMgrType, FieldMgrType>()`
> 模板参数 `N`：2 = 单相（水/CO₂）+ 热；3 = 三相（水+CO₂）+ 热

---

## 1. 概述：FIM 求解器定位与模块层次

```
RunGenericFIMTransient<N>()          ← 主循环入口（RunGeneric_impl.hpp:34）
 ├── 初始化阶段
 │    ├── TransmissibilitySolver          ← 导流能力预计算
 │    ├── FIM_TopologyBuilder             ← 连接关系构建
 │    ├── FIM_ConnectionManager           ← 连接聚合
 │    └── FIM_BlockSparseMatrix<N>        ← 冻结稀疏模式
 └── 时间步循环  for(step=1; step<=max_steps; ++step)
      └── Newton 迭代  for(iter=0; iter<max_newton_iter; ++iter)
           ├── 残差/Jacobian 组装
           │    ├── 积累项  AssembleAccumulation
           │    ├── 通量项  AssembleFlux
           │    ├── 井源项  Assemble_Wells_*D_FullJac
           │    └── 边界项  Assemble_*D
           ├── PTC 对角增强（可选）
           ├── 线性求解  SolveLinearSystem<N>
           │    ├── SparseLU（Eigen）
           │    ├── AMGCL（+SparseLU 回退）
           │    └── BiCGSTAB
           ├── 线搜索  Armijo + 非单调窗口 + 三级救援
           └── 收敛判断（多模式）
```

核心设计思想：**一次组装同时得到残差向量 b 与 Jacobian 矩阵 A**，借助 `ADVar<N>`（N 维对偶数）的自动微分，无需有限差分。

---

## 2. 状态管理与变量定义

### 2.1 FIM_StateMap<N>

文件：`FIM_StateMap.h:14`

```cpp
template<int N>
struct FIM_StateMap {
    std::vector<double> P;   // 压力 [Pa]
    std::vector<double> T;   // 温度 [K]
    std::vector<double> Sw;  // 水饱和度（仅 N=3 时有效）

    void InitSizes(size_t totalBlocks);
};
```

- `P[i]`、`T[i]` 对所有 N 均存在；`Sw[i]` 仅当 N=3 时分配内存（`if constexpr (N==3)`，`FIM_StateMap.h:27`）。
- 数组长度 = `totalBlocks` = 矩阵块数 + 裂缝段数。

### 2.2 三时间层

在每个时间步入口（`RunGeneric_impl.hpp:280`），同时维护三份状态：

| 变量 | 语义 |
|------|------|
| `old_state` | 时间步起始（已接受）状态，用于积累项中的 $u^n$ |
| `state` | 当前迭代状态 $u^{k}$，每次 Newton 更新后就地修改 |
| `best_state` | 步内最低残差快照，用于最佳迭代守卫 |

### 2.3 ADVar<N> 自动微分变量

在残差组装循环中，`ADVar<N>` 同时携带值和梯度：

```cpp
// RunGeneric_impl.hpp:612
ADVar<N> P(state.P[bi]); P.grad(0) = 1.0;      // ∂/∂P
ADVar<N> T(state.T[bi]); T.grad((N==2)?1:2) = 1.0; // ∂/∂T
// N=3 时还有：
ADVar<N> Sw(state.Sw[bi]); Sw.grad(1) = 1.0;    // ∂/∂Sw
```

Jacobian 的每个块由 `ADVar` 的 `.grad[j]` 分量直接给出，无需额外差分。

### 2.4 FIM_BlockSparseMatrix<N>

文件：`FIM_BlockSparseMatrix.h:29`

```
diag_blocks_[i]       : N×N 块（对角）
off_diag_blocks_[i][j]: N×N 块（非对角，仅连接存在时）
residual_[i]          : N×1 右端向量
```

`global_mat.FreezePattern()`（`RunGeneric_impl.hpp:158`）在初始化后冻结稀疏模式，后续迭代只填充数值，避免重复动态分配。

---

## 3. 时间步循环

### 3.1 主循环结构

```cpp
// RunGeneric_impl.hpp:278
for (int step = 1; step <= params.max_steps &&
     (target_end_time_s <= 0.0 || t < target_end_time_s - 1e-12); ++step)
```

每步开始时：
1. 如设有目标终止时间，截断 `dt = min(dt, t_remaining)`（`RunGeneric_impl.hpp:302`）。
2. 根据当前模拟时刻选择活动参数档位（startup/long）。
3. 重置步内缓存（线搜索历史、PTC 增强因子）。

### 3.2 两阶段参数配置

当 `params.enable_two_stage_profile == true` 时（`RunGeneric_impl.hpp:261`）：

```
t < params.startup_end_time_s  →  使用 startup_profile
t ≥ params.startup_end_time_s  →  使用 long_profile
```

- `startup_profile`：通常使用较小 `dt_max`、较低收敛容差，适合非稳态早期阶段。
- `long_profile`：放宽限制，允许更大时间步，提高效率。

切换时会打印 `[PROFILE]` 日志（`RunGeneric_impl.hpp:322`）。

### 3.3 控制坡道（Control Ramp）

当 `active_enable_control_ramp == true`（`RunGeneric_impl.hpp:360`）：

$$\text{ramp}(s) = r_0 + (1-r_0) \cdot s, \quad s = \frac{\text{step}-1}{\text{ramp\_steps}-1}$$

- 速率控制：`w.target_value *= control_ramp`
- BHP 控制：`w.target_value = p_anchor + ramp * (BHP_target - p_anchor)`
- 注入温度也按相同比例坡道（`RunGeneric_impl.hpp:401`）

目的是避免在仿真初期施加完整井约束导致残差爆炸。

---

## 4. Newton 迭代流程

### 4.1 外层迭代

```cpp
// RunGeneric_impl.hpp:576
for (int iter = 0; iter < active_max_newton_iter; ++iter) {
    global_mat.SetZero();
    // 1) 积累项
    // 2) 通量项
    // 3) 井源项
    // 4) 边界项
    // → 提取 A, b
    // → PTC 增强
    // → 线性求解
    // → 线搜索 + 状态更新
    // → 收敛判断
}
```

### 4.2 积累项组装

**N=2（单相+热）**（`RunGeneric_impl.hpp:627`）：

$$R_0^{\text{acc}} = \frac{V_i}{\Delta t}\bigl(\rho_w \phi - \rho_w^n \phi^n\bigr)$$
$$R_1^{\text{acc}} = \frac{V_i}{\Delta t}\bigl(e_w - e_w^n\bigr)$$

其中 $e_w = \rho_w(h_w - P/\rho_w) + (1-\phi)\rho_r c_{pr} T$（内能 + 岩石热容）。

岩石孔隙度含压缩性（`RunGeneric_impl.hpp:623`）：
$$\phi = \phi_{\text{ref}}\bigl(1 + c_r(P - P_{\text{ref}})\bigr)$$

**N=3（水+CO₂+热）**（`RunGeneric_impl.hpp:647`）：

$$R_0^{\text{acc}} = \frac{V_i}{\Delta t}(\rho_w \phi S_w - \rho_w^n \phi^n S_w^n)$$
$$R_1^{\text{acc}} = \frac{V_i}{\Delta t}(\rho_g \phi S_g - \rho_g^n \phi^n S_g^n)$$
$$R_2^{\text{acc}} = \frac{V_i}{\Delta t}\bigl[(e_{\text{fluid}}\phi + e_{\text{rock}}) - (e_{\text{fluid}}^n\phi^n + e_{\text{rock}}^n)\bigr]$$

组装调用：`FIM_GlobalAssembler<N,ADVar<N>>::AssembleAccumulation(bi, acc_eqs, global_mat)`（`RunGeneric_impl.hpp:656`）。

### 4.3 通量项组装

对每条连接 `conn(i,j)`（`RunGeneric_impl.hpp:677`）：

**N=2**：
```
dPhi = P_i - P_j - rho_avg * g * (z_i - z_j)   // 势差（含重力）
mob = rho/mu                                      // 迁移率（上风格式）
F[0] = T_Flow * mob_up * dPhi                    // 质量通量
F[1] = T_Heat * (T_i - T_j) + F[0] * h_up       // 热通量（传导+对流）
```

**N=3**（`RunGeneric_impl.hpp:701`）：
额外引入毛管压力 $P_c(S_w)$、相对渗透率 $k_{rw}$、$k_{rg}$（Mualem-van Genuchten）。

两次求导：`evalFlux(true)` 和 `evalFlux(false)` 分别对节点 i 和节点 j 进行 seed，得到 Jacobian 的两列贡献（`RunGeneric_impl.hpp:731`）。

### 4.4 Jacobian 提取

`ADVar<N>` 的 `.val` 直接给出残差分量，`.grad[j]` 给出 $\partial R / \partial u_j$。由 `FIM_GlobalAssembler` 将 N×N 子块填入 `FIM_BlockSparseMatrix`，最后通过：

```cpp
auto A = global_mat.ExportEigenSparseMatrix();  // RunGeneric_impl.hpp:908
auto b = global_mat.ExportEigenResidual();
```

导出 Eigen 稀疏矩阵和残差向量。

---

## 5. 线性求解

文件：`StepKernels_impl.hpp:26`

### 5.1 行缩放（可选）

```cpp
// StepKernels_impl.hpp:42
if (params.enable_row_scaling) {
    for (int r = 0; r < totalEq; ++r) {
        double denom = max(|D_acc[r]|, |A[r,r]|, row_scale_floor);
        double s = clamp(1/denom, row_scale_min, row_scale_max);
        A[r,:] *= s;  b[r] *= s;
    }
}
```

缩放分母选择积累项对角（`D_acc`）与矩阵对角中的最大值，保证缩放具有物理意义。

### 5.2 SparseLU（Eigen 直接法）

```cpp
// StepKernels_impl.hpp:111
// 模式未变则只重新分解（factorize），模式变化则完整 compute
if (!shape_changed) {
    cache.sparse_lu_solver.factorize(A_lu);   // 快速路径
} else {
    cache.sparse_lu_solver.compute(A_lu);     // 完整路径（含符号分析）
}
dx = cache.sparse_lu_solver.solve(-b_work);
```

模式哈希（FNV-1a，`StepKernels_impl.hpp:7`）用于检测稀疏结构是否改变，缓存符号分解，减少重复计算。

### 5.3 AMGCL（代数多重网格）

```cpp
// StepKernels_impl.hpp:66
if (!cache.amgcl_solver_ready) {
    cache.amgcl_solver = make_unique<AMGCLSolver<N>>(A_work, cache.amgcl_prm);
} else {
    cache.amgcl_solver->precond().rebuild(A_work);  // 仅重建预条件
}
auto [iters, error] = (*cache.amgcl_solver)(rhs, dx);
```

若 AMGCL 未收敛（`error > params.amgcl_tol`）且启用了 `amgcl_use_fallback_sparselu`，自动回退到 SparseLU（`StepKernels_impl.hpp:90`）。

---

## 6. 线搜索策略

### 6.1 初始阻尼 α

在进入 Armijo 循环之前，先根据最大允许变化量约束 α（`RunGeneric_impl.hpp:1345`）：

```cpp
alpha = min(alpha, max_dP_eff / (|dx[eqP]| + 1e-14));
alpha = min(alpha, max_dT_eff / (|dx[eqT]| + 1e-14));
// N=3:
alpha = min(alpha, max_dSw_eff / (|dx[eqSw]| + 1e-14));
```

有效上限含时间步感知缩放（`RunGeneric_impl.hpp:1340`）：
```
damp_scale = min(1, dt_eff / dt_ref)
max_dP_eff = max(1e-3, params.max_dP * damp_scale)
```

N=3 时还有安全两相阻尼（`enable_alpha_safe_two_phase`），防止 Sw 越过 [ε, 1-ε]（`RunGeneric_impl.hpp:1358`）。

### 6.2 Armijo 条件

```cpp
// RunGeneric_impl.hpp:1466
armijo_rhs = max(0, ref_res - c1 * alpha_try * conv_res_for_line_search);
// 试探点残差 ≤ armijo_rhs → 接受
```

参数：
- `c1`（Armijo 常数）：`clamp(armijo_c1, 1e-8, 1e-2)`
- 回退步长：`alpha_try *= bt_beta`，`bt_beta ∈ [0.1, 0.95]`
- 最大回退次数：`armijo_max_backtracks`

### 6.3 非单调窗口

```cpp
// RunGeneric_impl.hpp:1162
if (params.enable_nonmonotone_line_search) {
    line_search_hist.push_back(conv_res);
    // 保留最近 nonmonotone_window 个残差
    ref_res = max(line_search_hist);  // 参考残差取窗口最大值
}
```

允许残差在窗口内有小幅振荡，避免 Armijo 过于保守导致步长缩至极小。

### 6.4 三级救援

| 级别 | 触发条件 | 行为 |
|------|----------|------|
| **级别 1：Fallback 接受** | 所有 bt 均未满足 Armijo，但最优试探残差 ≤ 1.1 × ref | 接受最优试探点（`RunGeneric_impl.hpp:1556`）|
| **级别 2：受控接受（iter=1）** | iter=1 首次线搜索失败，且最优试探更新量在限制内 | 放宽 `controlled_accept_relax` 倍后接受（`RunGeneric_impl.hpp:1581`）|
| **级别 3：PTC 救援** | 连续 `ls_fail_rescue_threshold` 步在 iter=1 线搜索失败 | 不缩 dt，增大 PTC 对角增强因子 `ptc_rescue_boost`，重试当前步（`RunGeneric_impl.hpp:1680`）|

---

## 7. 收敛判断

每次 Newton 迭代结束后，依次检查（`RunGeneric_impl.hpp:1064`）：

### 7.1 绝对残差

```cpp
if (conv_res < params.abs_res_tol) {
    converged = true; converge_mode = "abs_res";
}
```

### 7.2 双重准则（相对残差 + 相对更新）

```cpp
// RunGeneric_impl.hpp:1069
double rel_res = conv_res / res_iter1;       // 与首次迭代残差之比
if (iter_used > 1
    && rel_res <= active_rel_res_tol
    && last_rel_update <= active_rel_update_tol) {
    converged = true; converge_mode = "rel_res_update";
}
```

- `rel_res_tol`（默认典型值 1e-3）：残差需降低至少 3 个量级
- `rel_update_tol`（默认典型值 1e-4）：相对更新量 `max(|dP|/P_ref, |dT|/T_ref)` 足够小

### 7.3 停滞接受（Stagnation Accept）

在最后一次 Newton 迭代时，若残差虽未达到容差但已显著下降（`RunGeneric_impl.hpp:1147`）：

```cpp
bool stagnation_ok = (growth <= stagnation_growth_tol)
    && (conv_res <= stagnation_abs_res_tol)
    && (conv_res < res_iter1)
    && (rel_drop >= stagnation_min_drop);
```

满足则以 `converge_mode = "stagnation"` 接受，并缩减下一步 `dt *= 0.85`。

### 7.4 最佳迭代守卫（Best-Iter Guard）

```cpp
// RunGeneric_impl.hpp:1079
if (enable_best_iter_guard && iter_used >= best_iter_guard_min_iter
    && best_quality_ok
    && grow_from_best > best_iter_growth_trigger) {
    state = best_state;
    converge_mode = "best_iter_guard";
}
```

当后续 Newton 步从已充分下降的"最优点"反弹时，回退到最优快照，防止"假收敛后爆炸"。

类似机制还有 `dt_floor_best`（在 dt 已达最小值时）和 `dt_floor_hold`（在 iter=1 时最优即为初始状态）（`RunGeneric_impl.hpp:1102`、`1127`）。

---

## 8. 自适应时间步

时间步调整在每步成功后执行（`RunGeneric_impl.hpp:1722`）：

### 8.1 绝对残差收敛时（`abs_res`）

| 使用迭代次数 | 操作 |
|-------------|------|
| `iter_used ≤ 3` | `dt = min(dt * 1.2, dt_max)` |
| `iter_used ≤ 5` | `dt = min(dt * 1.05, dt_max)` |
| `iter_used > 5` | `dt = max(dt * 0.8, dt_min)` |

### 8.2 相对残差收敛时（`rel_res_update`）

由四段阈值控制（`RunGeneric_impl.hpp:1728`）：

```
iter_used ≤ dt_relres_iter_grow_hi           → dt *= dt_relres_grow_factor   (扩大)
iter_used ≤ dt_relres_iter_neutral_hi        → dt *= dt_relres_neutral_factor (保持)
iter_used ≤ dt_relres_iter_soft_shrink_hi    → dt *= dt_relres_soft_shrink_factor (软缩)
iter_used >  dt_relres_iter_soft_shrink_hi   → dt *= dt_relres_hard_shrink_factor (硬缩)
```

典型默认值示例：

| 参数 | 典型含义 | 代码位置 |
|------|---------|---------|
| `dt_relres_iter_grow_hi` | ≤3 次迭代视为"快速收敛" | `RunGeneric_impl.hpp:1728` |
| `dt_relres_grow_factor` | 快速收敛放大因子，如 1.5 | `RunGeneric_impl.hpp:1731` |
| `dt_relres_neutral_factor` | 中等收敛保持系数，如 1.0 | `RunGeneric_impl.hpp:1732` |
| `dt_relres_soft_shrink_factor` | 慢速收敛轻微缩短，如 0.8 | `RunGeneric_impl.hpp:1733` |
| `dt_relres_hard_shrink_factor` | 接近发散强制缩短，如 0.5 | `RunGeneric_impl.hpp:1734` |

### 8.3 其他收敛模式

`stagnation / dt_floor_best / best_iter_guard / dt_floor_hold`：`dt *= 0.85`（`RunGeneric_impl.hpp:1744`）。

### 8.4 回滚时

```cpp
// RunGeneric_impl.hpp:1699
dt = max(dt * rollback_shrink_factor, dt_min);
state = old_state;
step--;  // 重试当前时间步
```

---

## 9. 关键参数速查表（TransientSolverParams / TransientStageProfile）

### 9.1 时间步控制

| 参数 | 作用 |
|------|------|
| `dt_init` | 初始时间步 [s] |
| `dt_min` | 最小时间步（回滚下限）|
| `dt_max` / `long_profile.dt_max` | 最大时间步上限 |
| `max_steps` | 最大时间步数 |
| `target_end_time_s` | 目标模拟终止时刻（≤0 表示不限）|
| `rollback_shrink_factor` | 回滚时 dt 缩减因子（0.05~1.0）|

### 9.2 Newton 迭代

| 参数 | 作用 |
|------|------|
| `max_newton_iter` | 每步最大 Newton 迭代次数 |
| `abs_res_tol` | 绝对残差收敛容差 |
| `rel_res_tol` | 相对残差收敛容差（相对于首次迭代）|
| `rel_update_tol` | 相对状态更新量容差 |
| `max_dP` / `max_dT` / `max_dSw` | 单步最大允许更新量（阻尼上限）|
| `min_alpha` | 最小允许线搜索步长 |

### 9.3 PTC（伪瞬态连续）

| 参数 | 作用 |
|------|------|
| `enable_ptc` | 开关 |
| `ptc_lambda_init` | 初始增强系数 |
| `ptc_lambda_decay` | 每次迭代衰减率（0~1）|
| `ptc_lambda_min` | 最小增强系数下限 |
| `ls_fail_rescue_threshold` | 触发 PTC 救援的连续 iter=1 线搜索失败次数 |
| `ls_fail_rescue_max` | 每步最多触发 PTC 救援次数 |
| `ptc_rescue_boost` | PTC 救援时增强因子放大倍率 |

### 9.4 线搜索

| 参数 | 作用 |
|------|------|
| `enable_armijo_line_search` | 开启 Armijo 线搜索 |
| `armijo_c1` | Armijo 常数（足够下降系数）|
| `armijo_beta` | 回退步长收缩因子 |
| `armijo_max_backtracks` | 最大回退次数 |
| `enable_nonmonotone_line_search` | 非单调模式开关 |
| `nonmonotone_window` | 非单调历史窗口大小 |
| `enable_controlled_accept_iter1` | iter=1 受控接受开关 |
| `controlled_accept_relax` | 受控接受放宽倍率 |

### 9.5 收敛 / 停滞

| 参数 | 作用 |
|------|------|
| `enable_stagnation_accept` | 停滞接受开关 |
| `stagnation_growth_tol` | 最大允许残差增长比 |
| `stagnation_abs_res_tol` | 停滞接受的绝对残差上限 |
| `stagnation_min_drop` | 最小要求残差下降比 |
| `enable_best_iter_guard` | 最佳迭代守卫开关 |
| `best_iter_guard_min_iter` | 触发守卫的最低迭代次数 |
| `best_iter_growth_trigger` | 从最优点反弹倍率触发阈值 |

### 9.6 线性求解

| 参数 | 作用 |
|------|------|
| `lin_solver` | `SparseLU` / `AMGCL` / `BiCGSTAB` |
| `enable_row_scaling` | 行缩放开关 |
| `row_scale_floor` | 缩放分母最小值（防除零）|
| `row_scale_min` / `row_scale_max` | 缩放系数范围 |
| `amgcl_tol` | AMGCL 收敛容差 |
| `amgcl_use_fallback_sparselu` | AMGCL 失败时回退 SparseLU |

### 9.7 两阶段配置

| 参数 | 作用 |
|------|------|
| `enable_two_stage_profile` | 开关 |
| `startup_end_time_s` | startup 阶段结束时刻 |
| `startup_profile` | `TransientStageProfile`（startup 参数包）|
| `long_profile` | `TransientStageProfile`（long 参数包）|

---

## 10. 典型执行序列（伪代码摘要）

```
初始化 state = ic, old_state = state
预计算 T_Flow, T_Heat（导流能力）
冻结稀疏模式

for step = 1..max_steps:
    old_state = state           // 保存本步起始状态
    select active_profile       // startup vs long
    apply control_ramp          // 坡道处理井约束

    for iter = 0..max_newton_iter:
        global_mat.SetZero()
        accumulation(state, old_state, dt)   // ADVar → R + J
        flux(state, connMgr)                  // ADVar → R + J
        wells(state)                          // 数值差分 Jacobian
        boundaries(state)                     // 数值 bc Jacobian

        A = ExportEigenSparseMatrix()
        b = ExportEigenResidual()
        conv_res = ||b||_inf

        if conv_res → converged: break

        if enable_ptc: augment A diagonal
        dx = SolveLinearSystem(A, b)
        alpha = compute_initial_damping(dx)   // max_dP/max_dT 约束
        (state, alpha) = armijo_line_search(state, dx, alpha)
        last_rel_update = max(|alpha*dx|/u_ref)

    if converged:
        t += dt; adapt dt (iter-based scaling)
    else:
        state = old_state; dt *= rollback_shrink; step--
```

---

*文档基于代码版本（git HEAD `54e2ec9`），生成于 2026-03-16。*
