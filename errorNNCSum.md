# EDFM 框架问题报告：NNC / FF 通量计算与矩阵组装

> **日期**：2026-03-16
> **性质**：纯分析报告，不含任何代码修改
> **范围**：NNC（基岩-裂缝）、FF（裂缝-裂缝）导流能力计算、拓扑装载、矩阵组装

---

## 目录

1. [问题总览](#1-问题总览)
2. [BUG-01（严重）：2D FF Star-Delta 同裂缝元对双重计数](#2-bug-01严重2d-ff-star-delta-同裂缝元对双重计数)
3. [BUG-02（显著）：3D FF 存在相同的双重计数结构问题](#3-bug-02显著3d-ff-存在相同的双重计数结构问题)
4. [BUG-03（显著）：NNC 裂缝侧渗透率 k_n → k_t 无声回退](#4-bug-03显著nnc-裂缝侧渗透率-k_n--k_t-无声回退)
5. [设计问题-01：2D NNC aux 几何占位符](#5-设计问题-012d-nnc-aux-几何占位符)
6. [设计问题-02：3D FF PointToSegmentDistance = 0 边界情况](#6-设计问题-023d-ff-pointtosegmentdistance--0-边界情况)
7. [精度问题：3D NNC distMatrixToFracPlane 的近似精度](#7-精度问题3d-nnc-distmatrixtofracplane-的近似精度)
8. [确认正确的模块](#8-确认正确的模块)
9. [修复优先级建议](#9-修复优先级建议)

---

## 1. 问题总览

| # | 级别 | 所在文件 | 核心问题 | 定量影响 |
|---|------|---------|---------|---------|
| BUG-01 | **严重** | `TransmissibilitySolver_2D.cpp:412` | 2D FF Star-Delta 包含同裂缝元对，与 FI 连接叠加 | 穿通交叉节点处同裂缝导流能力高估 **50%~100%** |
| BUG-02 | **显著** | `TransmissibilitySolver_3D.cpp:356` | 3D FF 同一结构问题 | 量级视交线几何而定，结构与 BUG-01 完全一致 |
| BUG-03 | **显著** | `TransmissibilitySolver_2D.cpp:269`<br>`TransmissibilitySolver_3D.cpp:268` | NNC 裂缝侧渗透率静默回退 $k_n \to k_t$ | 充填裂缝场景 $T_\text{NNC}$ 高估；所有场景无警告 |
| 设计-01 | 轻微 | `FIM_TopologyBuilder2D.h:108` | 2D NNC aux_area/aux_dist 用占位符 1.0 | 当前无功能 Bug，但语义错误 |
| 设计-02 | 轻微 | `TransmissibilitySolver_3D.cpp:510` | 3D FF 距离为零时强制截断至 1e-6 | 极端几何情形产生非物理大传导率 |
| 精度 | 中等 | `TransmissibilitySolver_3D.cpp:314` | 3D NNC `distMatrixToFracPlane` 可能为点距而非体积加权均值 | 大/不规则网格精度低于 2D 的 `avgDistance` |

---

## 2. BUG-01（严重）：2D FF Star-Delta 同裂缝元对双重计数

### 2.1 问题位置

- **文件**：`TransmissibilitySolver_2D.cpp:412`，函数 `Calculate_Transmissibility_FF`
- **关联文件**：`FIM_TopologyBuilder2D.h:51`（`_loadFI`），`FIM_TopologyBuilder2D.h:116`（`_loadFF`）

### 2.2 问题描述

Star-Delta 循环将同一 Junction（`GlobalFFPoint ID`）中的**所有**裂缝段做全配对，不区分这些段是否属于**同一条父裂缝**：

```cpp
// TransmissibilitySolver_2D.cpp:533-554
for (size_t i = 0; i < nElems; ++i) {
    for (size_t j = i + 1; j < nElems; ++j) {
        T_FF_Flow[ffIdx] = (half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow;
        // ↑ 包含了同一裂缝的相邻元对 (A1, A2)，type = Fracture_Fracture
        fieldMgr.ff_topology.emplace_back(elemSolverIndices[i], elemSolverIndices[j]);
        ffIdx++;
    }
}
```

与此同时，`_loadFI` 已经为同一裂缝上相邻的元对（A1, A2）生成了 `type = Fracture_Internal` 的连接：

```cpp
// FIM_TopologyBuilder2D.h:72-83
for (size_t i = 0; i + 1 < frac.elements.size(); ++i) {
    int nodeI = frac.elements[i].solverIndex;
    int nodeJ = frac.elements[i + 1].solverIndex;
    connMgr.PushConnection(nodeI, nodeJ, ..., ConnectionType::Fracture_Internal);
}
```

由于 `FinalizeAndAggregate` 以 `{type, nodeI, nodeJ}` 为键，`Fracture_Internal` 和 `Fracture_Fracture` 是**两个独立键**，二者不会合并，而是**同时进入通量主循环**，等效为并联叠加：

```
T_total(A1, A2) = T_FI(A1, A2)  +  T_FF_StarDelta(A1, A2)   ← 双重计数
```

### 2.3 定量误差推导

设交叉节点处有 $n$ 条等参数裂缝段（$K$、$W_f$、$h$、段长 $L$ 均相等），则：

$$G_k = \frac{K_t W_f h}{L/2} = \frac{2K_t W_f h}{L}$$

$$T_\text{FI}(A1, A2) = \frac{h}{d_{A1}/\text{cond}_{A1} + d_{A2}/\text{cond}_{A2}} = \frac{K_t W_f h}{L/2 + L/2} = \frac{K_t W_f h}{L}$$

$$T_\text{FF\_StarDelta}(A1, A2) = \frac{G_{A1} \cdot G_{A2}}{\sum_{k=1}^{n} G_k} = \frac{G^2}{nG} = \frac{G}{n} = \frac{2K_t W_f h}{nL}$$

$$\text{误差} = \frac{T_\text{FF\_StarDelta}}{T_\text{FI}} = \frac{2}{n}$$

| 交叉类型 | n（总支路数） | $T_\text{FF}/T_\text{FI}$ | 等效高估比例 |
|---------|-------------|--------------------------|-------------|
| T 型（1 条穿通 + 1 条终止） | 3 | 2/3 ≈ 67% | **+67%** |
| X 型（2 条裂缝均穿通） | 4 | 2/4 = 50% | **+50%** |
| 多裂缝汇聚（n 条） | n | 2/n | 随 n 增大递减 |

对于一条**仅有两个元**的纯线段（无其他裂缝介入，n=2）：

$$T_\text{FF\_StarDelta}(A1, A2) = \frac{G_{A1} G_{A2}}{G_{A1}+G_{A2}} = T_\text{FI}$$

此时 $T_\text{total} = 2 \cdot T_\text{FI}$，误差 **100%**。

### 2.4 物理解释与正确做法

**错误根因**：Star-Delta 变换的物理假设是"将若干支路汇聚到同一无体积节点"。对于裂缝网络：
- **跨裂缝对**（A1-B1 等）：必须通过 Junction 节点传递，Star-Delta 是唯一连接路径，**应该生成**
- **同裂缝对**（A1-A2）：它们在物理上通过裂缝本体直接相连，FI 连接已经覆盖，Star-Delta **不应再生成**

标准 EDFM 文献（Li & Lee 2008；Moinfar et al. 2014）中，FF Star-Delta **仅作用于不同裂缝之间的元对**。

**修复方向**（仅分析，不修改代码）：在 Star-Delta 内层循环的配对阶段，加入父裂缝 ID 的过滤条件，仅对 `parentFractureID_i ≠ parentFractureID_j` 的元对生成 FF 连接。

---

## 3. BUG-02（显著）：3D FF 存在相同的双重计数结构问题

### 3.1 问题位置

- **文件**：`TransmissibilitySolver_3D.cpp:417–438`（聚类阶段），`534–543`（Star-Delta 展开阶段）

### 3.2 问题描述

在 3D 中，两裂缝面的交叉是一条**线段**（交线）。当交线穿越同一裂缝（设为裂缝 A）的**多个单元格**时，多个属于裂缝 A 的元会出现在同一交线簇中：

```cpp
// TransmissibilitySolver_3D.cpp:417-438
for (const auto& inter : ffIntersections) {
    for (const auto& seg : inter.segments) {
        std::string key = getUndirectedKey(seg.start, seg.end);
        auto& cluster = clusterMap[key];
        if (seg.solverIndex_1 >= 0) cluster.solverIndices.push_back(seg.solverIndex_1);
        if (seg.solverIndex_2 >= 0) cluster.solverIndices.push_back(seg.solverIndex_2);
        // ↑ 来自同一裂缝的多个段（A1, A2...）同样被收入同一簇
    }
}
```

后续 Star-Delta 展开时，A1-A2（同属裂缝 A）会被配对生成 FF 连接，而它们之间已有 FI（`FractureEdge`）连接，形成与 BUG-01 完全相同的双重计数。

### 3.3 严重程度评估

3D 中的误差通常小于 2D，原因：
- 交线通常只穿越较少的同裂缝元格（受网格细化程度影响）
- 每簇中跨裂缝元对数量相对更多，稀释了同裂缝 Star-Delta 值

但**结构上完全一致**，属于同源问题，修复时需同步处理。

---

## 4. BUG-03（显著）：NNC 裂缝侧渗透率 k_n → k_t 无声回退

### 4.1 问题位置

- **文件**：`TransmissibilitySolver_2D.cpp:268–269`，`TransmissibilitySolver_3D.cpp:267–268`

### 4.2 问题代码

```cpp
// 2D（3D 相同）
auto p_Kf = fieldMgr.getFractureScalar(frac_str.k_n_tag);
if (!p_Kf) p_Kf = fieldMgr.getFractureScalar(frac_str.k_t_tag); // Fallback
```

NNC 导流能力公式中，裂缝侧的物理路径方向**垂直于裂缝面**：

$$T_\text{NNC} = \frac{A}{d_m / K_{mn} + (W_f/2) / K_f}$$

此处 $K_f$ 应为裂缝的**法向渗透率** $K_n$（抵抗垂直穿越的阻力），而非**切向渗透率** $K_t$（沿裂缝面流动的立方定律渗透率）。

### 4.3 两种渗透率的物理差异

| 渗透率 | 物理意义 | 典型量级（清洁裂缝） | 适用连接 |
|--------|---------|-------------------|---------|
| $K_t$（切向） | 沿裂缝面的流动能力，立方定律 $W_f^2/12$ | 极高（$10^{-10}$–$10^{-7}$ m²） | FI / FF |
| $K_n$（法向） | 穿越裂缝面的流动能力，取决于充填物 | 可低至 $10^{-18}$ m²（充填）| **NNC** |

### 4.4 误差影响分析

对于清洁开放裂缝：$(W_f/2)/K_t \to 0$，矩阵侧阻力 $d_m/K_{mn}$ 占主导，回退误差**可忽略**。

对于充填/部分充填裂缝：$K_n \ll K_t$，$(W_f/2)/K_n$ 不可忽略，回退使 $T_\text{NNC}$ **被高估**，导致裂缝-基岩交换速率偏快。

**更关键的风险**：回退是**静默的**——若用户未配置 `k_n_tag`（常见配置遗漏），程序不报错、不警告，使用 $K_t$ 继续计算。这在生产代码中是危险的，容易掩盖配置错误。

---

## 5. 设计问题-01：2D NNC aux 几何占位符

### 5.1 问题位置

`FIM_TopologyBuilder2D.h:105–112`，函数 `_loadNNC`

### 5.2 问题描述

```cpp
connMgr.PushConnection(
    nodeI, nodeJ,
    pFlow->data[idx], pHeat->data[idx],
    1.0, 1.0,                          // ← aux_area, aux_dist 均为占位符
    ConnectionType::Matrix_Fracture
);
```

而 3D 传递的是真实几何量：

```cpp
// FIM_TopologyBuilder3D.h:70-73
connMgr.PushConnection(
    pairs[i].matrixSolverIndex, pairs[i].fracCellSolverIndex,
    pFlow->data[i], pHeat->data[i],
    pairs[i].intersectionArea,          // ← 真实交互多边形面积
    pairs[i].distMatrixToFracPlane,     // ← 真实法向距离
    ConnectionType::Matrix_Fracture
);
```

### 5.3 影响

**现状**：由于一个（matrixCell, fracSeg）对在几何上唯一，NNC 的 `FinalizeAndAggregate` 实际上从不触发聚合，1.0 占位符不影响结果。

**潜在风险**：若未来网格生成允许同一 NNC 对出现多次（如裂缝被裂成多段穿过同一基岩单元），聚合时 `aux_dist` 的加权平均将给出无意义的 `1.0`，而不是真实的平均距离。

传递真实的 `pElem->length`（面积）和 `pElem->avgDistance`（距离）成本几乎为零，应与 3D 保持一致。

---

## 6. 设计问题-02：3D FF PointToSegmentDistance = 0 边界情况

### 6.1 问题位置

`TransmissibilitySolver_3D.cpp:510`

### 6.2 问题代码

```cpp
double d = std::max(
    Geometry_3D::PointToSegmentDistance(pElem->centroid, cluster.start, cluster.end),
    1e-6  // ← 强制截断，但未记录
);
double cond_flow = Kt[fLoc] * Wf[fLoc];
double t_f = (cond_flow * cluster.length) / d;   // d→1e-6 时 t_f 极大
```

### 6.3 影响

当裂缝元质心恰好落在交线段上（`d = 0`），截断至 `1e-6` 后产生 `G_k → ∞`。

在 Star-Delta 中：
$$T_\text{ij} = \frac{G_i G_j}{\sum G_k}$$

若某支路 $G_k \to \infty$，则其与其余所有支路的 $T_\text{ij}$ 趋近于有限值（收敛到 $G_{\text{other}}$），但自身对 $T_\text{ii}$（已被过滤）无意义，整体计算不会崩溃。然而：

- 数值稳定性依赖于截断值 `1e-6` 的选取；
- 物理上，质心恰在交线说明该元被交线对称分割，用单一质心距离代表流动路径本身就有模糊性；
- 该极端情形未在代码中记录，可能导致难以排查的静默误差。

---

## 7. 精度问题：3D NNC distMatrixToFracPlane 的近似精度

### 7.1 背景

EDFM 要求的正确 $d_m$ 是基岩单元内各点到裂缝面距离的**空间加权平均**：

$$d_m = \frac{1}{A_\text{int}} \iint_{\Omega_\text{int}} \left| \mathbf{n} \cdot (\mathbf{x} - \mathbf{x}_\text{frac}) \right| \, dA$$

### 7.2 可能的问题

`distMatrixToFracPlane` 的命名暗示这是**单元中心**到裂缝平面的点距：

$$d_m \approx |\mathbf{n} \cdot (\mathbf{x}_\text{center} - \mathbf{x}_\text{frac})|$$

而 2D 中的 `avgDistance` 根据其命名和用途（预计算的"平均距离"）更接近积分均值。

### 7.3 影响

两者的差异在以下情况下显著：
- **大网格**：基岩单元尺度与裂缝位置偏心量相当时，点距与积分均值差异可达 30%+
- **裂缝接近单元边界**（非穿心）：点距会低估平均距离，导致 $T_\text{NNC}$ 被高估
- **不规则多面体网格**：单元中心可能不是几何质心

若 3D 确实仅用点距而非积分均值，其 NNC 精度低于 2D，属于算法层面的不对等。建议核查 `InteractionPair::distMatrixToFracPlane` 的实际计算逻辑（位于 3D 几何模块）。

---

## 8. 确认正确的模块

以下模块经逐行比对，逻辑与物理公式均正确：

| 模块 | 文件 | 结论 |
|------|------|------|
| 串联阻力 / 两点传导率算子 | `FVM_Ops.h:82,100` | ✓ 公式 $T = A/(d_1/K_1 + d_2/K_2)$ 正确 |
| 2D / 3D FI 导流能力（调和平均） | `TransmissibilitySolver_2D.cpp:137`，`3D.cpp` | ✓ 标准 1D FVM 串联公式 |
| 2D NNC 法向渗透率投影 | `TransmissibilitySolver_2D.cpp:364-367` | ✓ $K_{mn} = n_x^2 K_{xx} + n_y^2 K_{yy}$，法向量从裂缝切向旋转 |
| 3D NNC 渗透率投影（含 z） | `TransmissibilitySolver_3D.cpp:311` | ✓ $K_{mn} = n_x^2 K_{xx} + n_y^2 K_{yy} + n_z^2 K_{zz}$ |
| FinalizeAndAggregate 聚合规则 | `FIM_ConnectionManager.h:105` | ✓ NNC/FF 累加，MM/FI 异常；排序确定性 |
| AssembleFlux 四块 Jacobian 填充 | `FIM_GlobalAssembler.h:42` | ✓ 对角/非对角符号与 FIM Jacobian 推导一致 |
| ADVar 双次 seed 机制 | `RunGeneric_impl.hpp:702-756` | ✓ 正确分离 $\partial F/\partial u_i$、$\partial F/\partial u_j$ |
| 统一通量主循环（无 type 分支） | `RunGeneric_impl.hpp:696` | ✓ 所有连接类型对称处理，设计简洁 |
| 热通量 = 传导 + 对流 | `RunGeneric_impl.hpp:718,750` | ✓ 标准公式 $F_h = T_H \Delta T + F_m h_\text{up}$ |
| PushConnection 输入校验 | `FIM_ConnectionManager.h:85` | ✓ NaN/Inf/负值检查，小负值截零 |
| 3D FF 量化 KEY 确定性排序 | `TransmissibilitySolver_3D.cpp:443-449` | ✓ 消除 unordered_map 迭代随机性 |

---

## 9. 修复优先级建议

```
优先级 1（严重，影响核心物理结果）
─────────────────────────────────────
BUG-01：2D FF Star-Delta 过滤同裂缝元对
BUG-02：3D FF 同步修复（同源问题）

优先级 2（显著，可能静默引入错误）
─────────────────────────────────────
BUG-03：NNC k_n 回退时至少输出 WARNING，
        或在文档中明确标注适用条件

优先级 3（设计改进，不影响当前正确性）
─────────────────────────────────────
设计-01：2D NNC 传递真实 aux_area / aux_dist
精度问题：核查 3D distMatrixToFracPlane 是否为积分均值

优先级 4（防御性处理）
─────────────────────────────────────
设计-02：3D FF PointToSegmentDistance=0 时
         记录 WARNING 或在几何层面预处理
```

---

*报告生成：2026-03-16*
*分析依据：源码逐行比对 + EDFM 标准文献（Li & Lee 2008；Moinfar et al. 2014）*
