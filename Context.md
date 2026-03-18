# 2D-Unstr-Quadrilateral-EDFM 工程全景上下文

> 本文档供 AI 辅助开发时快速建立工程全貌。每次新建对话时优先阅读本文件。

---

## 1. 工程定位与研究目标

**工程名称**：2D/3D 非结构化四边形网格嵌入式离散裂缝模型（EDFM）数值模拟器

**物理问题**：超临界 CO₂ 在裂缝性多孔介质中的热-流耦合渗流模拟，面向增强型地热系统（EGS）和CO₂地质封存应用。

**求解物理模式**：
- **N=2（单相+热）**：压力 P、温度 T —— CO₂ 或水单相流动与传热
- **N=3（两相+热）**：水压 P、水饱和度 Sw、温度 T —— CO₂-水两相渗流与传热

**论文章节映射**：
| 章节 | 状态 | 内容 |
|------|------|------|
| 第一章 | 已实现 | 非结构化网格 EDFM 快速生成与耦合（2D + 3D） |
| 第二章 | 在研 | 超临界 CO₂ 单相热-流耦合数值模型 |
| 第三章 | 在研 | CO₂-水两相渗流-传热耦合模型 |
| 第四章 | 待开展 | CO₂-H₂O 交替注采取热过程与运行策略优化 |

---

## 2. 工程目录树

```
2D-Unstr-Quadrilateral-EDFM/
├── Context.md                          ← 本文件
├── PROJECT_CONTEXT.md                  ← 实施进度与快速接手清单
├── REGRESSION_CASES.md                 ← 最小回归测试集
├── PLAN_7D_CHECKLIST.md                ← 7天实施计划验收清单
│
├── help/
│   ├── Intro_EDFM.md                   ← EDFM 架构、连接类型、导流能力公式
│   ├── Intro_FIM.md                    ← FIM 求解器架构、ADVar、Newton迭代
│   ├── Intro_FVMandSolverProcess.md    ← FVM-EDFM-FIM 全景框架
│   └── Error_NNCandFF.md               ← NNC/FF 错误诊断手册
│
├── FIM_TransientEngine/
│   ├── RunGeneric.hpp                  ← 公开入口（薄包装）
│   ├── RunGeneric_impl.hpp             ← 主求解循环（119K，核心！）
│   ├── StepKernels.hpp                 ← 时间步调度头
│   ├── StepKernels_impl.hpp            ← Newton迭代实现（通量/井/边界装配）
│   ├── FIM_GlobalAssembler.h           ← 积累/通量/源项装配接口
│   ├── FIM_BlockSparseMatrix.h         ← N×N 块稀疏矩阵（Pattern Freeze）
│   ├── FIM_ConnectionManager.h         ← 连接聚合管理（MM/FI/NNC/FF）
│   ├── StateSync.hpp                   ← 状态同步（old/current/best）
│   ├── Diagnostics.hpp                 ← 诊断工具
│   └── Types.hpp                       ← FIM 内部类型定义
│
├── Test/                               ← 测试输出目录（vtk, dat, csv等）
│
└── [源代码文件 ~153个，见第3节]
```

---

## 3. 模块架构与数据流

```
主调度 (main.cpp, --case dispatcher)
  ↓
网格初始化 (MeshManager / 3D_MeshManager)
  ├── 基岩网格生成：Gmsh API → Cell, Face, Node
  └── 裂缝网络生成、求交、拓扑构建
      ├── 3D EDFM（非匹配嵌入式，裂缝=三维域内的二维平面通道）:
      │       Fracture_2D → FractureEdge_2D → FractureNetwork_2D
      │       注：文件名前缀 `2D_` 指裂缝被降维为二维平面对象，非"2D模拟"
      └── 2D EDFM（传统匹配嵌入式，裂缝=二维域内的一维线段）:
              Fracture → FractureElement → FractureNetwork
  ↓
静态传导率计算 (TransmissibilitySolver_2D/3D)
  ├── MM: T = A·K/(d1+d2)（面积矢量E，各向异性K）
  ├── FI: T = Kt·Wf·h / (L/2)（裂缝段内相邻）
  ├── NNC: T = L_seg·h·K_proj / avgDistance（EDFM核心）
  └── FF: Star-Delta 等效导流（裂缝交叉点）
  ↓
场数据注册 (FieldManager_2D/3D → FieldRegistry)
  ├── matrixFields: P, T, Sw, rho, mu, ... (nMatrix 个单元)
  ├── fractureFields: 同上 (nFrac 个裂缝单元)
  └── nncFields: T_NNC_Flow, T_NNC_Heat, ... (nNNC 条连接)
  ↓
物性初始化 (PhysicalPropertiesManager_2D/3D)
  ├── 岩石: Kxx, Kyy, phi, cp_rock, rho_rock
  ├── CO2: Span-Wagner → CO2PropertyTable → AD_FluidEvaluator
  └── Water: IAPWS-95 → WaterPropertyTable → AD_FluidEvaluator
  ↓
FIM 拓扑装载 (FIM_TopologyBuilder2D/3D → FIM_ConnectionManager)
  └── 转换为统一 Connection 列表（MM/FI/NNC/FF 四类）
  ↓
矩阵冻结 (FIM_BlockSparseMatrix::FreezePattern)
  └── 固化稀疏结构，后续迭代仅 O(nnz) 填值
  ↓
FIM主循环 (RunGenericFIMTransient)
  ├── 时间步循环（startup/long 两阶段参数）
  └── Newton迭代
      ├── 积累项：V/Δt·(ρφ - ρⁿφⁿ) + 岩石压缩项
      ├── 通量项：统一循环 MM/FI/NNC/FF
      │   ├── 势差：ΔΦ = ΔP - ρ_avg·g·Δz（含重力）
      │   ├── 两相：ΔΦ += ΔPc（毛管压差）
      │   ├── 迎风：Op_Upwind_AD(ΔΦ, mob_i, mob_j)
      │   └── 通量：q = T·mob_upwind·ΔΦ
      ├── 井源项：WI*(P_cell - P_bhp) 或 定流量
      ├── 边界项：Dirichlet/Neumann/Robin
      └── 线性求解 → 线搜索(Armijo) → 收敛判断
  ↓
后处理 (PostProcess_2D/3D)
  ├── VTK Legacy ASCII (.vtk)
  └── Tecplot (.dat)
```

---

## 4. 连接类型（EDFM 核心）

```cpp
enum ConnectionType {
    Matrix_Matrix     = 0,  // MM：基岩-基岩相邻面（FVM 标准）
    Fracture_Internal = 1,  // FI：裂缝段内部相邻
    Matrix_Fracture   = 2,  // NNC：裂缝-基岩非邻连接（EDFM核心）
    Fracture_Fracture = 3   // FF：裂缝-裂缝交叉点
};
```

| 类型 | 面积 A | 距离 d | 渗透率 K |
|------|--------|--------|---------|
| MM | `face.vectorE.Mag() × thickness` | 法向投影距离 | `nx²Kxx + ny²Kyy` |
| FI | `thickness × Wf` | 段半长 L/2 | `Kt × Wf` |
| NNC | `L_seg × h` | `avgDistance`（几何平均距离） | `nx²Kxx + ny²Kyy` |
| FF | Star-Delta：`Gi×Gj/ΣGk` | `L/2` | `Kt × Wf` |

---

## 5. 关键数据结构

### ADVar<N>（自动微分变量，核心！）
```cpp
// 文件：ADVar.hpp
template<int N>
struct ADVar {
    double val;           // 值
    Eigen::Vector<double, N> grad;  // 梯度（N=2: [∂/∂P, ∂/∂T]; N=3: +∂/∂Sw）
    // 完整运算符重载：+, -, *, /, pow, sqrt, exp, log, abs, max, min
};
// 使用：ADVar<2> P; P.val=1e7; P.grad[0]=1.0;（seeding）
```

### FIM_BlockSparseMatrix<N>（块稀疏矩阵）
```cpp
// 文件：FIM_TransientEngine/FIM_BlockSparseMatrix.h
// 内部存储：Eigen::SparseMatrix<double> 导出给 AMGCL
// Pattern Freeze 优化：
//   第一次：Triplet 路径建立 CSR 缓存（O(nnz log nnz)）
//   后续：直接写入 valuePtr（O(nnz)，零分配）
// 关键方法：
AddResidual(block_idx, eq, val);               // 填充残差向量
AddDiagJacobian(block_idx, eq, var, jac_val);  // 填充对角块
AddOffDiagJacobian(i, j, eq, var, jac_val);    // 填充非对角块
```

### FieldManager_2D（三域场管理）
```cpp
// 文件：2D_FieldManager.h
// 三个域：matrixFields / fractureFields / nncFields
// 使用：
auto& P = fm.matrixFields.get<volScalarField>("Pressure");   // 标量场
auto& P_ad = fm.matrixFields.get<volADField<N>>("Pressure"); // AD场
```

### Connection（连接结构体）
```cpp
// 文件：FIM_ConnectionManager.h
struct Connection {
    ConnectionType type;
    int cell_i, cell_j;  // 两端求解器索引
    double T_flow;        // 流动传导率 [m³·s/kg 或类似]
    double T_heat;        // 热传导率 [W/K 或类似]
    // 对于 NNC/FF：还有面积、方向等辅助信息
};
```

---

## 6. 自动微分与流体物性

### AD_FluidEvaluator（EOS 桥接）
```
文件：AD_FluidEvaluator.h
架构：PropertyTable查找表 + 4层鲁棒数值差分 + AssembleADVar
  Tier1: 中心差分（Central）
  Tier2: 前向差分（Forward）
  Tier3: 后向差分（Backward）
  Tier4: 微小导数地板（DerivativeFloor）

关键接口：
  ADFluidProperties<N> evaluateCO2(P_ADVar, T_ADVar)
  ADFluidProperties<N> evaluateWater(P_ADVar, T_ADVar)
  void AssembleADVar(target, val, d_dP, d_dT, P, T)
```

### CO2PropertyTable（查找表）
```
文件：CO2PropertyTable.h/cpp
  - P 轴：若干离散压力点（Pa）
  - T 轴：若干离散温度点（K）
  - 2D Catmull-Rom 插值（双三次样条）
  - 属性：rho, mu, cp, cv, h, k（密度、粘度、比热、焓、热导率）
  - 物理范围：P ∈ [1MPa, 60MPa]，T ∈ [270K, 650K]
  - 注意：std::filesystem（已修复 experimental/filesystem 废弃问题）
```

### CapRelPerm_HD（毛管压力+相对渗透率）
```
文件：CapRelPerm_HD.h，CapRelPerm_HD_AD.h（AD版本）
模型：Mualem-van Genuchten
  krw = Se^L · (1-(1-Se^(1/m))^m)²
  krg = (1-Se)^L · (1-Se^(1/m))^(2m)
  Pc = (Se^(-1/m)-1)^(1/n) / α
参数结构：VGParams{alpha, n, Swr, Sgr}, RelPermParams{L}
```

---

## 7. 非正交修正（当前状态）

面 E/T 向量分解（`face.cpp::computeFaceVectors`）：
```
Aj = 面积矢量（normal × length）
ej = owner → neighbor 单位向量
MinimumCorrection：  vectorE = (Aj·ej) ej,   vectorT = Aj - vectorE
OrthogonalCorrection: vectorE = |Aj| ej,       vectorT = Aj - vectorE
OverRelaxed：        vectorE = |Aj|²/(Aj·ej) ej, vectorT = Aj - vectorE
```

**当前问题**：`TransmissibilitySolver` 仅用 `vectorE.Mag()` 作为传导率有效面积，
`vectorT`（非正交修正项）**未在通量循环中施加**。
非正交网格（如斜四边形）存在一阶截断误差，需在 `StepKernels_impl.hpp` 通量循环添加延迟修正项。

---

## 8. 井管理（当前状态）

井当前实现为**源项**（非独立 DOF）：
- BHP 控制：`q = WI × (P_cell - P_bhp)`（`FVM_Ops_AD.h::Op_Well_BHP_Source_AD`）
- Rate 控制：定值质量流率（`Op_Well_Rate_Source_AD`）
- 能量源项：`q_E = q_w·h_w + q_g·h_g`（`Op_Well_Energy_Source_AD`）
- 装配：`FIM_GlobalAssembler::AssembleSource(block_idx, source_vec, mat)`

**已实现（Step 3，2026-03-18）**：独立 Well DOF（`WellDOFManager.h`）— 每口井新增一个额外的 FIM 矩阵块：
- BHP 控制：`R_well = P_wbh - P_target = 0`（Dirichlet 约束）
- Rate 控制：`R_well = WI_mob·(P_res - P_wbh) - q_spec = 0`（Peaceman 模型）
- 储层侧额外 off-diagonal：`∂R_res/∂P_wbh = -WI_mob`
- `WI_mob` 从 BoundaryAssembler 输出的 `w_jac3[g_eq_p][0]` 提取

---

## 9. 外部依赖

| 库 | 版本 | 路径 | 用途 |
|----|------|------|------|
| **Eigen** | 5.0.0 | `D:\1Eigen\eigen-5.0.0\eigen-5.0.0` | 线性代数、ADVar::grad 存储 |
| **Gmsh** | 最新 | `D:\1gmesh\gmsh-source\` | 网格生成（API调用） |
| **AMGCL** | master | `D:\1AMGCL\amgcl-master\amgcl-master` | 代数多重网格线性求解器 |
| **CoolProp** | 6.x | `D:\1CoolProp\CoolPropStaticLibrary\` | 流体物性（已链接，待激活） |

---

## 10. 构建说明

**IDE**：Visual Studio 2022（v143 toolset），C++17，x64

**编译**：
```powershell
msbuild .\2D-Unstr-Quadrilateral-EDFM.sln /p:Configuration=Debug /p:Platform=x64
```

**运行**：
```bash
./exe --list                                 # 列出所有 case
./exe --case=day6_transient_2d_sp_injprod    # 2D 单相瞬态注采
./exe --case=day6_transient_2d_tp_injprod    # 2D 两相瞬态注采
./exe --case=day5_global_jac_2d              # Jacobian 验证
./exe --case=day4_well_viz                   # 生成可视化
```

**已知构建问题（已修复）**：
- ✅ `CO2PropertyTable.cpp`、`WaterPropertyTable.cpp`：`experimental/filesystem` → `std::filesystem`
- ✅ `vcxproj` Debug：链接 `CoolPropd.lib`（Debug）而非 `CoolProp.lib`（Release）
- ✅`vcxproj` Release：已配置链接 `bulid\Release` Gmsh，**需用户执行 `cmake --build . --config Release` 生成 Release Gmsh**

---

## 11. 完整模块索引

### 基础网格数据结构
| 文件 | 主要类/结构 | 功能 |
|------|-----------|------|
| `node.h/cpp` | `Node` | 顶点坐标与 ID |
| `face.h/cpp` | `Face` | 面几何、E/T向量分解、FVM参数 |
| `cell.h/cpp` | `Cell` | 单元几何（中心、体积、AABB） |
| `mesh.h/cpp` | `Mesh` | 网格拓扑管理（2D/3D） |
| `MeshDefinitions.h` | `DistanceMetric`, `IntersectionSearchStrategy_2D` | 枚举定义 |

### 3D EDFM 裂缝模块（核心原创）
> **非匹配嵌入式**：在三维多面体非结构化基岩网格中，将裂缝降维为二维平面流动通道独立剖分，
> 通过八叉树空间索引与前沿推进式求交算法构建 NNC 数据集。
> **文件名前缀 `2D_` 表示"裂缝为二维平面对象"，与"2D模拟"无关。**

| 文件 | 主要类 | 功能 |
|------|--------|------|
| `2D_Fracture.h/cpp` | `Fracture_2D` | 单条裂缝面：独立非结构化剖分、与基岩的求交 |
| `2D_FractureMeshElement.h` | `FractureElement_2D` | 裂缝面微元（tri/quad，含 solverIndex） |
| `2D_FractureEdge.h/cpp` | `FractureEdge_2D` | 裂缝面内部边（FVM 向量、FI 连接） |
| `2D_FractureNetwork.h/cpp` | `FractureNetwork_2D` | 全局裂缝网络管理、全局索引分配 |
| `2D_FracIndex.h/cpp` | `FracElemIndex_2D` | 全局求解器索引偏移表 |
| `2D_FractureIntersection.h` | `FracFracIntersectionObject` | 裂缝面-裂缝面交集多边形数据结构 |

### 2D EDFM 裂缝模块（核心原创）
> **传统匹配嵌入式**：在二维非结构化基岩网格（tri/quad）中，将裂缝嵌入为一维线段单元，
> 网格与裂缝轨迹共形剖分，直接建立裂缝-基岩 NNC 连接关系。

| 文件 | 主要类 | 功能 |
|------|--------|------|
| `Fracture.h/cpp` | `Fracture` | 单条裂缝线段：匹配网格剖分、与基岩面求交 |
| `FractureCommon.h` | — | 裂缝公共数据类型与常量定义 |
| `FractureElement.h/cpp` | `FractureElement` | 裂缝一维单元（含 solverIndex） |
| `FractureIntersectionPoint.h/cpp` | `FractureIntersectionPoint` | 裂缝-基岩面交点几何 |
| `FractureNetwork.h/cpp` | `FractureNetwork` | 全局裂缝网络管理与索引分配 |
| `FracIndex.h/cpp` | `FracElemIndex` | 全局求解器索引偏移表 |
| `EDFM_Geometry_3D.h/cpp` | — | 裂缝-基岩交集几何计算（3D 坐标变换辅助） |
| `MatrixEdge.h` | `MatrixEdge` | 基岩边数据结构（求交加速） |
| `GeometryCalculate.h` | — | 通用几何计算工具（面积、距离、法向量等） |

### 空间索引加速
| 文件 | 功能 |
|------|------|
| `AABB.h` | 轴对齐包围盒（构建、相交测试） |
| `FaceIndexedOctree.h` | 面索引八叉树（加速裂缝-面求交） |
| `8_DOP.h/cpp` | 8 方向离散多面体 |
| `14_DOP.h/cpp` | 14 方向离散多面体 |

### 网格管理器
| 文件 | 主要类 | 关键方法 |
|------|--------|---------|
| `MeshManager.h/cpp` | `MeshManager` | `BuildMesh2D_EDFM()`, `distributeSolverIndices()` |
| `3D_MeshManager.h/cpp` | `MeshManager_3D` | 3D 版本 |

### 传导率求解器
| 文件 | 主要类 | 四类传导率 |
|------|--------|-----------|
| `TransmissibilitySolver_2D.h/cpp` | `TransmissibilitySolver_2D` | `Calculate_Transmissibility_Matrix/FI/NNC/FF()` |
| `TransmissibilitySolver_3D.h/cpp` | `TransmissibilitySolver_3D` | 3D 版本 |

### FIM 拓扑与连接
| 文件 | 功能 |
|------|------|
| `FIM_TopologyBuilder2D.h` | `LoadAllConnections()` → 装载 MM/FI/NNC/FF |
| `FIM_TopologyBuilder3D.h` | 3D 版本 |
| `FIM_ConnectionManager.h` | `FinalizeAndAggregate()` → 并联聚合、排序 |

### 场数据管理
| 文件 | 功能 |
|------|------|
| `FieldRegistry.h` | `create<T>(name, size)`, `get<T>(name)` |
| `VolField.h` | `volScalarField`, `volADField<N>`, `faceScalarField` |
| `2D_FieldManager.h` | `InitSizes(nM, nF, nNNC, nFF, ...)` |
| `3D_FieldManager.h` | 3D 版本 |
| `FaceFieldRegistry.h` | 面字段注册 |
| `FieldAccess.h` | 字段统一访问接口 |

### 物性管理
| 文件 | 功能 |
|------|------|
| `2D_PhysicalPropertiesManager.h/cpp` | `UpdateFluidEOS_CO2()`, `UpdateAllFluidEOS()` |
| `2D_CO2_Properties.h/cpp` | CO2 EOS 更新（查找表接口） |
| `2D_Water_Properties.h/cpp` | 水 EOS 更新 |
| `CO2PropertyTable.h/cpp` | Catmull-Rom 2D 插值查找表（P-T 网格） |
| `WaterPropertyTable.h/cpp` | 同上（水） |
| `2D_RockSolidProperties.h/cpp` | 岩石物性（Kxx, Kyy, phi, cp, rho） |
| `2D_FractureProperties.h/cpp` | 裂缝物性（Wf, Fcd, phi_f） |
| `CapRelPerm_HD.h` | Van Genuchten 毛管压力+相对渗透率 |
| `CapRelPerm_HD_AD.h` | AD 版本（供 FIM 用） |
| `AD_FluidEvaluator.h` | EOS 桥接：PropertyTable → ADVar<N> |

### FVM 算子
| 文件 | 功能 |
|------|------|
| `FVM_Ops.h` | 传导率、散度等基础算子 |
| `FVM_Ops_AD.h` | AD 版本：`Op_Upwind_AD`, `Compute_Potential_Diff`, `Compute_Mass_Flux`, `Op_Well_BHP_Source_AD` |
| `FVM_Grad.h/cpp` | Green-Gauss / 最小二乘梯度重建 |

### FIM 求解器核心
| 文件 | 功能 |
|------|------|
| `FIM_BlockSparseMatrix.h` | N×N 块稀疏矩阵，Pattern Freeze 优化 |
| `FIM_GlobalAssembler.h` | `AssembleAccumulation/Flux/Source()` |
| `FIM_StateMap.h` | `FIM_StateMap<N>`：old/current/best 三层状态 |
| `FIM_JacobianVerifier.h` | Jacobian 有限差分验证工具 |
| `ADVar.hpp` | `ADVar<N>` 自动微分变量（全运算符重载） |

### 变量初始化
| 文件 | 功能 |
|------|------|
| `2D_VariableInitializer.h/cpp` | `InitFIMState<N>()`, `InitSinglePhaseState()` |
| `3D_VariableInitializer.h/cpp` | 3D 版本 |
| `InitParams.h` | `LinearInitParams`：线性初始场参数 |

### 边界与井
| 文件 | 功能 |
|------|------|
| `BoundaryConditionManager.h/cpp` | BC 类型管理（Dirichlet/Neumann/Robin） |
| `BoundaryAssembler.h/cpp` | BC 项装配到方程组 |
| `BoundaryFaceClassify.h` | 按 Gmsh PhysicalTag 分类边界面 |
| `Well_WellControlTypes.h` | `WellScheduleStep` 结构体，BHP/Rate/WAG |
| `Well_WellScheduleManager.h/cpp` | 井时间表加载与查询 |

### 后处理
| 文件 | 功能 |
|------|------|
| `2D_PostProcess.h/cpp` | `ExportVTK()`, `ExportVTU()`, `ExportPVD()`, `ExportTecplot()`, `SyncADFieldToScalar<N>()` |
| `3D_PostProcess.h/cpp` | 3D 版本 |
| `OutPutMesh.h` | 网格数据输出接口 |

### 辅助工具
| 文件 | 功能 |
|------|------|
| `UserDefineVarType.h` | `Vector`（= `Variable3D<double>`） |
| `Variable3D.hpp` | 3D 向量类（Mag, 运算符重载） |
| `VariableTensor.hpp` | 2D 张量类 |
| `GeneralMethods.h/cpp` | 通用辅助函数 |
| `GeometryCalculate.h/cpp` | 几何计算工具 |
| `VectorHash.h` | 向量哈希（用于 unordered_map） |
| `SolverContrlStrName_op.h` | 场名字符串常量 |

---

## 12. 测试框架索引

| 文件 | 对应 case 名称 | 功能 |
|------|--------------|------|
| `Test_Day6_TransientSolver.h/cpp` | `day6_transient_2d_sp_injprod` 等 | 完整瞬态 FIM 测试套件 |
| `Test_FluidEvaluator.h` | — | EOS 物性验证 |
| `Test_ADVar.h` | — | ADVar 单元测试 |
| `Test_FIM_Topology.h` | — | FIM 拓扑验证 |
| `Test_Day5_GlobalAssembly_Jacobian.h` | `day5_global_jac_2d/3d` | Jacobian 矩阵验证 |
| `test_Transmissibility_2D.h` | — | 2D 传导率基准 |
| `2D_EDFM_MeshTest_Benchmark.h` | — | 2D 网格性能基准 |
| `3D_EDFM_MeshTest_Benchmark.h` | — | 3D 网格性能基准 |

---

## 13. 在研与待实现项目

### 已实现（近期完成）
1. **CoolProp EOS 替换** ✅ `AD_FluidEvaluator.h` — `USE_COOLPROP_EOS` 宏已激活，实现：
   - 优先路径：`first_partial_deriv()` 解析导数（非两相区）
   - 降级路径：CoolProp 数值差分（两相/近临界区）
   - 线程安全：`thread_local AbstractState*` 支持 OpenMP

2. **vcxproj 运行时冲突修复** ✅ Debug 链接 `CoolPropd.lib`（/MDd），Release 链接 `CoolProp.lib` + Gmsh `bulid\Release`

3. **非正交修正项** ✅（Step 2，延迟修正策略）：
   - `Connection::vectorT` 携带每个 MM 连接的非正交修正向量
   - `FIM_TopologyBuilder2D::_loadMatrix()` 从 `face.vectorT` 填充
   - Newton 循环前用 Green-Gauss 重建单元中心梯度
   - AssembleFlux 后显式施加修正（无 Jacobian 贡献）
   - 启用：`params.enable_non_orthogonal_correction = true`

4. **VTU/PVD 后处理** ✅（Step 4）：
   - `PostProcess_2D::ExportVTU(filename, time)` — ParaView XML 单时刻快照
   - `PostProcess_2D::ExportPVD(pvd, vtu_files, times)` — 时间序列动画集合

5. **井独立 DOF** ✅（Step 3）：`WellDOFManager.h` 已创建，每口井新增独立矩阵块
   - `FIM_BlockSparseMatrix` 扩展至 `nMatrix + nFrac + nWells` 块
   - `StepKernels_impl.hpp::ApplyTrialUpdate` 新增 `num_well_blocks` 参数
   - `RunGeneric_impl.hpp` 集成：7 处修改（pattern 注册、状态初始化、装配、求解、更新）

### 尚未实现
（全部 Step 0-4 已完成）

### 已确认
4. **两相流迎风**：`Op_Upwind_AD` 已实现，确认 `StepKernels_impl.hpp` 两相通量循环使用基于势差的迎风（potential-based），而非简单的 owner-based

### 已修复（2026-03-18）
5. ✅ `CO2PropertyTable.cpp`、`WaterPropertyTable.cpp`：`experimental/filesystem` → `std::filesystem`
6. ✅ `vcxproj` Debug：`CoolProp.lib` → `CoolPropd.lib`
7. ✅ `vcxproj` Release：Gmsh lib 路径 `bulid\Debug` → `bulid\Release`（需用户编译 Release Gmsh）

---

## 14. Solver Index 规则

```
全局 DOF 顺序（solverIndex 规则）：
  [0, nMatrix)         → 基岩网格单元（按 MeshManager 分配）
  [nMatrix, nMatrix+nFrac) → 裂缝单元（跨所有裂缝，按 FractureNetwork 分配）
  nNNC = 裂缝-基岩 NNC 连接数（仅在 FieldManager nncFields 中记录）
  nFF  = 裂缝-裂缝 FF 连接数（同上）

关键变量（MeshManager 持有）：
  int nMatrixCells    // 基岩单元总数
  int nFracElems      // 裂缝单元总数（= FractureNetwork.total）
  int nNNC            // NNC 连接数
  int nFF             // FF 连接数
  int nTotalDOF = nMatrixCells + nFracElems
```

---

## 15. 关键算法参考

### EDFM NNC 平均距离
$$d_{avg} = \frac{\int_\Omega d(\mathbf{x}, \text{frac}) \, dV}{\int_\Omega dV}$$
对矩形单元有解析公式，对一般四边形单元用数值积分（`Fracture_2D::computeGeometryCouplingCoefficient`）。

### FIM Newton 迭代（隐式求解）
$$J \cdot \delta x = -R$$
- $R$：残差向量（积累项 + 通量项 + 井/边界源项）
- $J$：Jacobian（由 ADVar 自动提取，`AddDiagJacobian/AddOffDiagJacobian`）
- 线搜索：Armijo 准则 + 非单调窗口 + 三级救援

### Star-Delta FF 传导率
$$T_{FF}^{i-j} = \frac{G_i \cdot G_j}{\sum_k G_k}, \quad G_k = \frac{T_k \cdot W_k}{L_k/2}$$
