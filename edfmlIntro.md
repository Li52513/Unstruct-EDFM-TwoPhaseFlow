# EDFM 框架技术文档：通量计算、矩阵组装与求解流程

> 本文系统梳理 EDFM（嵌入式离散裂缝模型）框架中基岩-裂缝（NNC）、裂缝-裂缝（FF）通量计算、拓扑构建、并联聚合与 FIM 求解的完整流程，涵盖 2D 与 3D 两种实现。
>
> **日期**：2026-03-16
> **代码根目录**：`2D-Unstr-Quadrilateral-EDFM/`

---

## 目录

1. [总体架构](#1-总体架构)
2. [连接类型定义](#2-连接类型定义)
3. [底层数学算子](#3-底层数学算子)
4. [2D_EDFM：FI 裂缝内部导流能力](#4-2d_edfmfi-裂缝内部导流能力)
5. [2D_EDFM：NNC 导流能力（裂缝-基岩）](#5-2d_edfmnnc-导流能力裂缝-基岩)
6. [2D_EDFM：FF 导流能力（裂缝-裂缝）](#6-2d_edfmff-导流能力裂缝-裂缝)
7. [3D_EDFM：NNC 导流能力](#7-3d_edfmnnc-导流能力)
8. [3D_EDFM：FF 导流能力](#8-3d_edfmff-导流能力)
9. [拓扑装载](#9-拓扑装载)
10. [并联聚合](#10-并联聚合)
11. [通量组装](#11-通量组装)
12. [Jacobian 块贡献](#12-jacobian-块贡献)
13. [2D vs 3D 差异汇总表](#13-2d-vs-3d-差异汇总表)

---

## 1. 总体架构

```
┌─────────────────────────────────────────────────────────────────────┐
│  静态预计算阶段（初始化一次）                                         │
│                                                                     │
│  TransmissibilitySolver                                             │
│  ┌─────────────────────────────────────────────────────┐           │
│  │  Matrix-Matrix (MM)  → T_Matrix_Flow / T_Matrix_Heat │           │
│  │  Fracture-Internal (FI) → T_FI_Flow / T_FI_Heat      │           │
│  │  Matrix-Fracture (NNC) → T_NNC_Flow / T_NNC_Heat     │           │
│  │  Fracture-Fracture (FF) → T_FF_Flow / T_FF_Heat      │           │
│  └─────────────────────────────────────────────────────┘           │
│            │                                                        │
│            ▼                                                        │
│  FIM_TopologyBuilder2D / 3D                                         │
│  ┌─────────────────────────────────────────────────────┐           │
│  │  _loadMatrix() → _loadFI() → _loadNNC() → _loadFF() │           │
│  │  调用 FIM_ConnectionManager::PushConnection(...)      │           │
│  └─────────────────────────────────────────────────────┘           │
│            │                                                        │
│            ▼                                                        │
│  FIM_ConnectionManager::FinalizeAndAggregate()                      │
│  ┌─────────────────────────────────────────────────────┐           │
│  │  NNC/FF 并联累加 T_total = ΣT_k                      │           │
│  │  MM/FI  重复则抛出异常                               │           │
│  │  按 type→nodeI→nodeJ 字典序排序                      │           │
│  └─────────────────────────────────────────────────────┘           │
└─────────────────────────────────────────────────────────────────────┘
              │
              ▼  Newton 迭代内（每步重复）
┌─────────────────────────────────────────────────────────────────────┐
│  RunGeneric_impl.hpp：统一通量组装循环                               │
│  for (const auto& conn : connMgr.GetConnections())                  │
│    evalFlux(wrt_i) → evalFlux(wrt_j)                                │
│    FIM_GlobalAssembler::AssembleFlux(i, j, f_wrt_i, f_wrt_j, mat)  │
└─────────────────────────────────────────────────────────────────────┘
              │
              ▼
  AMGCL / 直接求解器 → 更新状态 → 收敛判断
```

**关键设计原则**：
所有连接类型（MM / FI / NNC / FF）在组装阶段**共用同一段通量代码**，差异仅在预计算的 `T_Flow` / `T_Heat` 值。

---

## 2. 连接类型定义

**文件**：`FIM_ConnectionManager.h:17`

```cpp
enum class ConnectionType {
    Matrix_Matrix     = 0,   // MM：基岩-基岩（相邻网格）
    Fracture_Internal = 1,   // FI：裂缝段内部相邻
    Matrix_Fracture   = 2,   // NNC：裂缝-基岩非邻连接（EDFM 核心）
    Fracture_Fracture = 3    // FF：裂缝-裂缝交叉点
};
```

**Connection 结构体**（`FIM_ConnectionManager.h:24`）：

```cpp
struct Connection {
    int nodeI;        // 全局求解器索引（较小值）
    int nodeJ;        // 全局求解器索引（较大值）
    double T_Flow;    // 流动导流能力 [m³/(Pa·s)]
    double T_Heat;    // 热传导导流能力 [W/K]
    double aux_area;  // 辅助几何：面积（供聚合加权）
    double aux_dist;  // 辅助几何：距离
    ConnectionType type;
};
```

构造函数自动将 `nodeI = min(i,j)`，`nodeJ = max(i,j)`，确保连接方向一致性。

**全局求解器索引约定**：
- `0 .. nMat-1`：基岩单元
- `nMat .. nMat+nFrac-1`：裂缝段

---

## 3. 底层数学算子

**文件**：`FVM_Ops.h:82`

### 3.1 串联电阻

```cpp
// FVM_Ops.h:82
inline double Op_Math_SeriesResistance(double d1, double K1, double d2, double K2)
{
    double term1 = d1 / (std::max(K1, 1e-30) + kEpsilon);
    double term2 = d2 / (std::max(K2, 1e-30) + kEpsilon);
    return term1 + term2;
}
// R = d1/K1 + d2/K2
```

### 3.2 两点传导率

```cpp
// FVM_Ops.h:100
inline double Op_Math_Transmissibility(double d1, double K1, double d2, double K2, double area)
{
    double resistance = Op_Math_SeriesResistance(d1, K1, d2, K2);
    return area / (resistance + kEpsilon);
}
// T = Area / (d1/K1 + d2/K2)
```

这两个函数是所有导流能力计算（MM/FI/NNC/FF）的共同入口。`area` 参数对 NNC 来说是"裂缝段在基岩中的有效截面积"，对 FF 是"裂缝厚度×交线长度"。

---

## 4. 2D_EDFM：FI 裂缝内部导流能力

**文件**：`TransmissibilitySolver_2D.cpp:137`
**函数**：`Calculate_Transmissibility_FractureInternal`

**物理模型**：将裂缝段视为 1D FVM 单元，相邻段之间的通量由段中心到公共节点的半长 d1、d2 决定。

```
T_FI = (K_t · W_f · h) / (d1 + d2)
     ≡ Op_Math_Transmissibility(d1, K_t·W_f, d2, K_t·W_f, h)
```

| 符号 | 含义 | 代码变量 |
|------|------|---------|
| `K_t` | 裂缝切向渗透率 | `Kt[s1]` |
| `W_f` | 裂缝张开度（孔径） | `Wf[s1]` |
| `h` | 2D 厚度（= 1.0 m） | `thickness = 1.0` |
| `d1` | 段 1 半长 = `length × 0.5` | `e1.length * 0.5` |
| `d2` | 段 2 半长 = `length × 0.5` | `e2.length * 0.5` |

**关键代码**（行 216–229）：

```cpp
double cond1 = Kt[s1] * Wf[s1];
double cond2 = Kt[s2] * Wf[s2];
T_FI_Flow[currentConnIdx] = FVM_Ops::Op_Math_Transmissibility(d1, cond1, d2, cond2, thickness);
```

有效热导率使用体积混合规则：
```cpp
double lam_eff = phi * lam_fluid + (1.0 - phi) * lam_solid;
```

结果写入场：`"T_FI_Flow"`，`"T_FI_Heat"`（FractureFace 标量场）。

---

## 5. 2D_EDFM：NNC 导流能力（裂缝-基岩）

**文件**：`TransmissibilitySolver_2D.cpp:247`
**函数**：`Calculate_Transmissibility_NNC`

### 5.1 物理公式

```
T_NNC = (L_seg × h) / ( d_m/K_mn  +  W_f/2 / K_f )
```

| 符号 | 含义 | 代码来源 |
|------|------|---------|
| `L_seg` | 裂缝段在基岩内的长度 | `pElem->length` |
| `h` | 2D 单位厚度 = 1.0 m | `thickness = 1.0` |
| `d_m` | 基岩单元中心到裂缝线的**平均垂直距离** | `pElem->avgDistance` |
| `K_mn` | 基岩法向渗透率投影 | `nx²·Kxx + ny²·Kyy` |
| `W_f/2` | 裂缝半张开度（裂缝侧物理路径长） | `Wf[fLocIdx] / 2.0` |
| `K_f` | 裂缝法向渗透率（字段 `k_n_tag`，回退 `k_t_tag`） | `Kf[fLocIdx]` |

### 5.2 法向渗透率投影

```cpp
// TransmissibilitySolver_2D.cpp:358-367
Vector tangent = p2 - p1;                     // 裂缝切向量
Vector normal(-tangent.m_y, tangent.m_x, 0);  // 法向量（旋转90°）
normal = normalize(normal);
double nx = normal.m_x, ny = normal.m_y;
double k_m_dir = nx*nx * Kxx[mIdx] + ny*ny * Kyy[mIdx];
```

### 5.3 关键代码（行 369–386）

```cpp
double d_m  = std::max(pElem->avgDistance, 1e-6);
double area = pElem->length * thickness;

// Flow
T_Flow[job.nncIdx] = FVM_Ops::Op_Math_Transmissibility(
    d_m, k_m_dir, Wf[fLocIdx]/2.0, Kf[fLocIdx], area);

// Heat
double lam_f_eff = phi * lam_fluid + (1.0 - phi) * lam_solid;
T_Heat[job.nncIdx] = FVM_Ops::Op_Math_Transmissibility(
    d_m, lam_m, Wf[fLocIdx]/2.0, lam_f_eff, area);
```

使用 `#pragma omp parallel for` 并行化（行 327）。
结果写入场：`"T_NNC_Flow"`，`"T_NNC_Heat"`（NNC 标量场）。

---

## 6. 2D_EDFM：FF 导流能力（裂缝-裂缝）

**文件**：`TransmissibilitySolver_2D.cpp:412`
**函数**：`Calculate_Transmissibility_FF`

### 6.1 Star-Delta（星角）变换

对汇聚到同一交叉节点（Junction）的 n 条裂缝段支路，先计算每条支路的半导流能力，再做 Star-Delta 展开：

**Step 1：支路半导流能力**
```
G_k = (K_t · W_f · h) / d_k
```
其中 `d_k = pElem->length * 0.5`（段中心到节点的距离，行 514）。

**Step 2：Star-Delta 展开**
```
T_ij = G_i × G_j / Σ G_k
```

**关键代码**（行 516–537）：

```cpp
// 支路半导流能力
double cond = KtF[fLocIdx] * Wf[fLocIdx];
double t_f  = (cond * thickness) / d;     // d = pElem->length * 0.5
half_T_Flow[i] = t_f;
sum_T_Flow += t_f;

// Star-Delta 展开
T_FF_Flow[ffIdx] = (half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow;
```

### 6.2 聚类策略

以 `GlobalFFPoint ID`（整数）为键聚类，枚举 `elem.gIDstart` 和 `elem.gIDend`（行 451–455）。
排序 junctionIDs 以保证确定性（行 488）。

结果写入场：`"T_FF_Flow"`，`"T_FF_Heat"`（FF 标量场）；
同时填充 `fieldMgr.ff_topology`（`std::vector<std::pair<int,int>>`，记录每对 FF 连接的 solverIndex）。

---

## 7. 3D_EDFM：NNC 导流能力

**文件**：`TransmissibilitySolver_3D.cpp:248`
**函数**：`Calculate_Transmissibility_NNC`

### 7.1 物理公式

```
T_NNC = A_int / ( d_m/K_mn  +  W_f/2 / K_f )
```

与 2D 公式形式相同，但几何量含义完全不同：

| 符号 | 2D 含义 | 3D 含义 |
|------|--------|--------|
| 面积 | `L_seg × h`（长度×单位厚度） | `A_int`：基岩单元与裂缝平面的**真实交互多边形面积** |
| 距离 `d_m` | `avgDistance`（预计算平均值） | `distMatrixToFracPlane`：精确法向距离 |
| 渗透率投影 | `nx²Kxx + ny²Kyy` | `nx²Kxx + ny²Kyy + nz²Kzz` |

### 7.2 关键代码（行 307–316）

```cpp
// 3D 法向渗透率投影（含 z 分量）
double k_m_dir = nx*nx*Kxx[mIdx] + ny*ny*Kyy[mIdx] + nz*nz*Kzz[mIdx];

// 面积来自交互多边形，距离来自精确几何计算
double dist = std::max(pair.distMatrixToFracPlane, 1e-6);
T_Flow[i] = FVM_Ops::Op_Math_Transmissibility(
    dist, k_m_dir, Wf[fLocIdx]/2.0, Kf[fLocIdx],
    pair.intersectionArea);          // ← 真实多边形面积
```

法向量 `pair.polygonNormal` 是裂缝平面的单位法向量，由交互多边形几何决定。

---

## 8. 3D_EDFM：FF 导流能力

**文件**：`TransmissibilitySolver_3D.cpp:356`
**函数**：`Calculate_Transmissibility_FF`

### 8.1 量化弦 KEY 聚类（行 391–415）

3D FF 中，两条裂缝平面的交叉是一条**线段**而非点，需以线段为键聚类：

```cpp
const double TOL = 1e-4;
auto quantize = [TOL](const Vector& v) {
    return std::make_tuple(
        (long long)(floor(v.m_x / TOL + 0.5)),
        (long long)(floor(v.m_y / TOL + 0.5)),
        (long long)(floor(v.m_z / TOL + 0.5))
    );
};
// 无向键（端点较小者排前）
if (q1 > q2) std::swap(q1, q2);
key = "x1_y1_z1-x2_y2_z2";
```

相同交线的所有裂缝段归入同一簇，取最长线段的几何保留代表性坐标。

### 8.2 确定性排序（行 443–449）

```cpp
// 提取所有 KEY 后字典序排序，消除 unordered_map 迭代顺序的随机性
std::sort(sortedKeys.begin(), sortedKeys.end());
```

保证 CSV 报表和矩阵装配在不同编译器/运行平台上的完全可复现性。

### 8.3 PointToSegmentDistance（行 510）

```cpp
double d = std::max(
    Geometry_3D::PointToSegmentDistance(pElem->centroid, cluster.start, cluster.end),
    1e-6
);
```

取裂缝段质心到交线段的距离，代替 2D 中的 `length * 0.5`，更准确地反映 3D 几何。

### 8.4 Star-Delta（行 512–537）

```cpp
double cond_flow = Kt[fLoc] * Wf[fLoc];
double t_f = (cond_flow * cluster.length) / d;   // cluster.length = 交线长度

// 展开
T_FF_Flow[ffIdx] = (half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow;
```

Flow 与 Heat **独立判定**有效性（行 519）：某支路 Heat 数据缺失时，该支路 Heat 贡献置 0，不阻断 Flow 计算。

---

## 9. 拓扑装载

### 9.1 FIM_TopologyBuilder2D（`FIM_TopologyBuilder2D.h:11`）

```cpp
static void LoadAllConnections(...) {
    _loadMatrix(connMgr, meshMgr, fieldMgr);   // MM
    _loadFI    (connMgr, meshMgr, fieldMgr);   // FI
    _loadNNC   (connMgr, meshMgr, fieldMgr);   // NNC
    _loadFF    (connMgr, meshMgr, fieldMgr);   // FF
}
```

**_loadNNC（行 87–113）**：

```cpp
const auto& pairsMap = meshMgr.getNNCTopologyMap();  // {matIdx → [fracSolverIdx,...]}
size_t idx = 0;
for (const auto& kv : pairsMap) {
    int nodeI = kv.first;                     // 基岩单元全局索引
    for (int nodeJ : kv.second) {             // 裂缝段全局求解器索引
        connMgr.PushConnection(
            nodeI, nodeJ,
            pFlow->data[idx], pHeat->data[idx],
            1.0, 1.0,                         // ← 2D 不传递真实几何（占位符）
            ConnectionType::Matrix_Fracture
        );
        ++idx;
    }
}
```

注意：2D NNC 的 `aux_area=1.0`, `aux_dist=1.0` 是占位值，不参与实际聚合加权。

### 9.2 FIM_TopologyBuilder3D（`FIM_TopologyBuilder3D.h:59`）

```cpp
const auto& pairs = meshMgr.getInteractionPairs();   // InteractionPair vector
for (size_t i = 0; i < pairs.size(); ++i) {
    connMgr.PushConnection(
        pairs[i].matrixSolverIndex,
        pairs[i].fracCellSolverIndex,
        pFlow->data[i], pHeat->data[i],
        pairs[i].intersectionArea,             // ← 3D 传递真实面积
        pairs[i].distMatrixToFracPlane,        // ← 3D 传递精确距离
        ConnectionType::Matrix_Fracture
    );
}
```

3D NNC 传递真实的交互多边形面积和法向距离，用于后续的并联聚合加权。

### 9.3 FF 拓扑装载差异

| 维度 | aux_area | aux_dist |
|------|---------|---------|
| 2D | `0.5*(L_i*W_i + L_j*W_j)`（近似有效截面） | `0.5*(L_i + L_j)` |
| 3D | `0.5*(A_i + A_j)`（裂缝单元面积均值） | `centroid_i - centroid_j` 距离 |

---

## 10. 并联聚合

**文件**：`FIM_ConnectionManager.h:105`
**函数**：`FinalizeAndAggregate`

### 10.1 NNC / FF 的累加规则

相同 `(nodeI, nodeJ, type)` 的多条路径表示**物理并联**（同一基岩单元-裂缝对，来自多个几何片段）：

```cpp
// NNC 或 FF 出现重复时累加
it->second.T_Flow += raw.T_Flow;
it->second.T_Heat += raw.T_Heat;

// aux_dist 按面积加权平均
double totalArea = it->second.aux_area + raw.aux_area;
it->second.aux_dist =
    (it->second.aux_dist * it->second.aux_area + raw.aux_dist * raw.aux_area) / totalArea;
it->second.aux_area = totalArea;
```

物理含义：并联导流能力相加：$T_\text{total} = \sum_k T_k$。

### 10.2 MM / FI 的重复检测

```cpp
// 非 NNC/FF 出现重复时抛出异常
throw std::logic_error(
    "[FIM Topology Error] Duplicate Matrix/FI connection at (" +
    std::to_string(raw.nodeI) + "," + std::to_string(raw.nodeJ) + ")"
);
```

基岩-基岩面和裂缝内部连接在物理上不允许并联（每条连接唯一对应一个界面/边）。

### 10.3 几何重复 vs 物理并联

聚合时会区分：
- **几何重复**（`sameGeomAndTrans == true`）：完全相同的几何和导流能力，视为代码层面的重复插入，仅计数不累加
- **物理并联**（`sameGeom == false`）：不同几何片段产生的独立路径，累加

### 10.4 排序输出

```cpp
std::sort(globalConnections_.begin(), globalConnections_.end(),
    [](const Connection& a, const Connection& b) {
        if (a.type != b.type) return int(a.type) < int(b.type);
        if (a.nodeI != b.nodeI) return a.nodeI < b.nodeI;
        return a.nodeJ < b.nodeJ;
    });
```

排序顺序：MM → FI → NNC → FF，每类内按 `(nodeI, nodeJ)` 字典序。

---

## 11. 通量组装

### 11.1 统一循环（无类型分支）

**文件**：`FIM_TransientEngine/RunGeneric_impl.hpp:696`

```cpp
for (const auto& conn : connMgr.GetConnections()) {
    // conn.type 可以是 MM / FI / NNC / FF，通量公式完全相同
    int i = conn.nodeI, j = conn.nodeJ;

    auto evalFlux = [&](bool wrt_i) -> std::vector<ADVar<N>> {
        std::vector<ADVar<N>> F(N);
        ADVar<N> P_i(state.P[i]), T_i(state.T[i]);
        ADVar<N> P_j(state.P[j]), T_j(state.T[j]);

        if (wrt_i) { P_i.grad(0) = 1.0; T_i.grad(1) = 1.0; }  // seed i 侧
        else       { P_j.grad(0) = 1.0; T_j.grad(1) = 1.0; }  // seed j 侧
        // ... 物理通量计算（见下）...
        return F;
    };

    auto f_wrt_i = evalFlux(true);   // ∂F/∂u_i
    auto f_wrt_j = evalFlux(false);  // ∂F/∂u_j
    FIM_GlobalAssembler<N, ADVar<N>>::AssembleFlux(i, j, f_wrt_i, f_wrt_j, global_mat);
}
```

**ADVar 双次 seed 机制**：同一个 `evalFlux` lambda 被调用两次，分别以 i 侧或 j 侧变量为 AD 自变量，得到两列偏导数，填入 Jacobian 的两个方向。

### 11.2 通量物理公式（N=2 单相热流）

```cpp
// 行 712-718（N==2 分支）
ADVar<N> rho_avg_w = 0.5 * (pW_i.rho + pW_j.rho);
ADVar<N> dPhi = Compute_Potential_Diff(P_i, P_j, rho_avg_w, x_i, x_j, g);
                                   // dPhi = (P_i - P_j) - rho_avg * g * dz
ADVar<N> mob_i = pW_i.rho / pW_i.mu;
ADVar<N> mob_j = pW_j.rho / pW_j.mu;
ADVar<N> up_mob = Op_Upwind_AD(dPhi, mob_i, mob_j);    // 上风格式迁移率

F[0] = conn.T_Flow * up_mob * dPhi;                    // 质量通量
F[1] = conn.T_Heat * (T_i - T_j) + F[0] * up_h;       // 热通量（传导+对流）
```

对 NNC 连接：`conn.T_Flow = T_NNC`，物理上代表裂缝-基岩界面的**等效串联导流能力**。

### 11.3 两相三方程（N=3，CO₂+水+热）

```cpp
// 行 741-750（N==3 分支）
F[0] = T_Flow * up_mobW * dPhi_w;          // 水相质量通量
F[1] = T_Flow * up_mobG * dPhi_g;          // CO₂ 相质量通量
F[2] = Compute_Heat_Flux(T_Heat, T_i, T_j,
        F[0], F[1], up_h_w, up_h_g);        // 热通量（双相对流 + 传导）
```

---

## 12. Jacobian 块贡献

**文件**：`FIM_GlobalAssembler.h:42`
**函数**：`AssembleFlux`

```cpp
// FIM_GlobalAssembler.h:60-68
global_mat.AddResidual(block_i, eq,  f_i);    // R[i] += F_ij
global_mat.AddResidual(block_j, eq, -f_i);    // R[j] -= F_ij

for (int var = 0; var < N; ++var) {
    // 对角块（自耦）
    global_mat.AddDiagJacobian   (block_i, eq, var,  flux_wrt_i[eq].grad(var));
    global_mat.AddDiagJacobian   (block_j, eq, var, -flux_wrt_j[eq].grad(var));
    // 非对角块（互耦）
    global_mat.AddOffDiagJacobian(block_i, block_j, eq, var,  flux_wrt_j[eq].grad(var));
    global_mat.AddOffDiagJacobian(block_j, block_i, eq, var, -flux_wrt_i[eq].grad(var));
}
```

**四块贡献汇总表**：

| 矩阵位置 | 值 | 来源 | 物理含义 |
|---------|-----|------|---------|
| 残差 `R[i]` | `+F_ij.val` | `flux_wrt_i[eq].val` | 节点 i 流出 |
| 残差 `R[j]` | `-F_ij.val` | `flux_wrt_i[eq].val` | 节点 j 流入 |
| 对角块 `J[i,i]` | `+∂F/∂u_i` | `flux_wrt_i[eq].grad(var)` | i 对自身的贡献 |
| 对角块 `J[j,j]` | `-∂F/∂u_j` | `flux_wrt_j[eq].grad(var)` | j 对自身的贡献 |
| 非对角块 `J[i,j]` | `+∂F/∂u_j` | `flux_wrt_j[eq].grad(var)` | **j→i 耦合（NNC 独特贡献）** |
| 非对角块 `J[j,i]` | `-∂F/∂u_i` | `flux_wrt_i[eq].grad(var)` | **i→j 耦合（NNC 独特贡献）** |

**NNC 的矩阵耦合意义**：
对一条 NNC 连接 `(matrixCell_m, fracCell_f)`，非对角块 `J[m,f]` 和 `J[f,m]` 将基岩块和裂缝块**直接耦合**到同一 Newton 线性系统中，这是 FIM（全隐式法）相对 IMPES 的核心优势——裂缝-基岩交换完全隐式处理。

---

## 13. 2D vs 3D 差异汇总表

| 维度 | 2D_EDFM | 3D_EDFM |
|------|---------|---------|
| **NNC 面积** | `L_seg × h`（长度×单位厚度） | `A_int`（真实交互多边形面积） |
| **NNC 距离** | `avgDistance`（预计算均值，`FractureElement` 字段） | `distMatrixToFracPlane`（精确法向距离，`InteractionPair` 字段） |
| **渗透率投影** | `nx²·Kxx + ny²·Kyy` | `nx²·Kxx + ny²·Kyy + nz²·Kzz` |
| **裂缝段体积** | `length × max(aperture, 1e-6)` | `max(area, 1e-12) × max(aperture, 1e-8)` |
| **NNC 拓扑源** | `getNNCTopologyMap()`（`map<int, vector<int>>`） | `getInteractionPairs()`（`vector<InteractionPair>`） |
| **NNC aux 几何** | `aux_area=1.0, aux_dist=1.0`（占位） | 传递真实 `intersectionArea` 和 `distMatrixToFracPlane` |
| **FF 聚类策略** | 按 `GlobalFFPoint` 整数 ID | 量化弦 KEY（浮点坐标→整数，TOL=1e-4） |
| **FF 距离计算** | `pElem->length × 0.5`（段半长） | `PointToSegmentDistance(centroid, start, end)` |
| **FF 交线长度** | 裂缝段长度（1D 交叉点） | `cluster.length`（3D 交线长度） |
| **Heat 缺失处理** | 整体判断一次 | 每条支路独立判定（`branchHasHeat` flag） |
| **FI 连接存储** | `FractureFace` 标量场（面拓扑） | `FractureEdge` 标量场（边拓扑） |
| **通量组装循环** | 统一（与 MM 共用代码） | 统一（与 MM 共用代码） |
| **Jacobian 结构** | 相同 | 相同 |

---

## 附录 A：关键文件索引

| 文件 | 主要内容 |
|------|---------|
| `FIM_ConnectionManager.h` | `ConnectionType`、`Connection`、`PushConnection`、`FinalizeAndAggregate` |
| `FVM_Ops.h` | `Op_Math_SeriesResistance`、`Op_Math_Transmissibility` |
| `TransmissibilitySolver_2D.cpp` | 2D MM/FI/NNC/FF 导流能力计算 |
| `TransmissibilitySolver_3D.cpp` | 3D MM/FI/NNC/FF 导流能力计算 |
| `FIM_TopologyBuilder2D.h` | 2D 连接装载入 ConnectionManager |
| `FIM_TopologyBuilder3D.h` | 3D 连接装载入 ConnectionManager |
| `FIM_GlobalAssembler.h` | `AssembleAccumulation`、`AssembleFlux`、`AssembleSource` |
| `FIM_BlockSparseMatrix.h` | 块稀疏矩阵（对角块/非对角块存储） |
| `ADVar.hpp` | 自动微分变量（值+梯度向量） |
| `FIM_TransientEngine/RunGeneric_impl.hpp` | Newton 迭代主循环（行 696：统一通量组装） |

---

## 附录 B：NNC 串联阻力模型示意

```
基岩单元中心                    裂缝段中心
     ●─────────────────────────────●
     |        |          |         |
     |<─ d_m ─>|< W_f/2 >|         |
     |        裂缝平面              |

     阻力 R_matrix = d_m / K_mn
     阻力 R_frac   = (W_f/2) / K_f
     总导流能力 T_NNC = A / (R_matrix + R_frac)

     2D: A = L_seg × h          （线×厚）
     3D: A = A_int               （多边形面积）
```

---

*文档生成时间：2026-03-16*
*作者：Claude Code（基于源码自动提取）*
