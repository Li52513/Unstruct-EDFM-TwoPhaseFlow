# 非结构化网格嵌入式离散裂缝模型 C++ 代码编写规则与注释规范

版本: v0.1 (草案)

目标: 基于《路书》确定的模块与算法，规范 C++ 代码的命名、组织、接口、注释与调用关系，确保四类问题的独立、准确与可维护实现。

----------------------------------------
## 1. 语言/依赖/基础类型

- 标准: C++17 或以上（与现有工程一致）
- 线性代数: Eigen（参考 `LinearSolver_Eigen.h`）
- 基础类型:
  - `using index_t = int32_t or int64_t;`（统一索引）
  - `using real_t = double;`（统一实数）
  - 枚举使用 `enum class`，避免隐式转换
  - 向量/矩阵首选 Eigen 的 `VectorXd/SpMat` 或轻量自定义结构

----------------------------------------
## 2. 命名/命名空间/文件

- 命名空间: `namespace edfm { ... }`，子空间按模块分：`edfm::mesh`, `edfm::frac`, `edfm::fvm`, `edfm::edfm`, `edfm::prop`, `edfm::bi`, `edfm::solver`, `edfm::poster`
- 文件名: `Module_Component.{h,cpp}`，如 `FVM_Assemble.h`，`EDFM_Coupling.cpp`
- 类/结构体: 首字母大写（`Mesh`, `FracPiece`, `Fields`）
- 函数: 小写加下划线（`build_face_trans`, `assemble_coupling`）或驼峰（保持全工程一致）
- 变量: 小写下划线（`cell_count`, `face_area`）；常量全大写（`MAX_ITERS`）
- 前缀: 模块前缀在类型/文件中体现；函数落在对应命名空间而不强制加字符串前缀

----------------------------------------
## 3. 目录结构建议

- `src/mesh`：网格拓扑与几何
- `src/frac`：裂缝读入、相交与 NNC 构建
- `src/fields`：场数据定义与存取
- `src/prop`：物性与本构
- `src/fvm`：MM 面通量与单元装配
- `src/edfm`：MF/FF 耦合与装配
- `src/bi`：边界与初始
- `src/solver`：时间循环与外迭代
- `src/poster`：导出与可视化
- `docs/`：文档（本文件与路书）

保持与现有工程兼容，可渐进式迁移。

----------------------------------------
## 4. 接口与数据结构

4.1 Mesh/FracMesh
- `struct Mesh { ... }; // cells, faces, neighbors, geom`
- `struct FracPiece { index_t cell; real_t area_or_len; Eigen::Vector3d n; Eigen::Vector3d c; real_t b; }`
- `struct NNC { index_t a, b; real_t A; real_t d_a, d_b; Eigen::Vector3d n; }; // MF/FF 通用`

4.2 Fields（SoA 首选）
- `struct Fields { std::vector<real_t> p, T, Sw; /* fracture DOFs if any */ };`
- 按需为裂缝片建立自由度容器：`std::vector<real_t> p_f, T_f, Sw_f;`

4.3 Properties 参数
- `struct FluidProps { /* tables or params */ };`
- `struct RockProps { /* K tensor, phi law, cp, lambda */ };`
- `struct RelPermParams { /* Corey params */ };`

----------------------------------------
## 5. 函数设计原则

- 纯函数优先：入参 `const&`，出参返回值或输出结构体；避免隐藏状态
- 明确单位/坐标系：注释标明 SI 单位与 z 方向定义
- 不抛异常的热路径：返回 `bool` 与错误码；调试使用 `assert`
- 避免重复分配：预分配容器，传入可写 `span`/引用
- 稀疏结构重用：稀疏模式（图）与数值分离；时间步内多次复用模式

----------------------------------------
## 6. 注释与文档（Doxygen 风格）

头注模板：
/**
 * brief: 计算基岩-裂缝 NNC 的等效通量系数
 * eq: T = A / (d_i/κ_i + d_f/κ_f)
 * @param mesh         基岩网格几何/拓扑
 * @param frac_piece   单元内裂缝片几何与属性（b, n, center）
 * @param rock_prop    单元各向异性 K 及相关参数
 * @param frac_prop    裂缝等效 K_f, b 等
 * @return real_t      T_{i,f}
 */

实现内注释：仅对复杂几何与边界约定进行说明；循环、简单代数不赘述。

----------------------------------------
## 7. 模块 API 规范（与《路书》对应）

7.1 mesh
- build_topology(const RawMesh&, Mesh& out)
- compute_geometry(const Mesh&, Geom& out)

7.2 frac
- intersect_cell(const Mesh&, index_t cell, const FracInput&, std::vector<FracPiece>& out)
- build_nnc(const Mesh&, const std::vector<FracPiece>&, std::vector<NNC>& mf, std::vector<NNC>& ff)

7.3 prop
- fluid_co2(real_t p, real_t T, const FluidProps&, FluidState& out)
- fluid_h2o(real_t p, real_t T, const FluidProps&, FluidState& out)
- relperm_corey(real_t Sw, const RelPermParams&, real_t& krw, real_t& krn)
- capillary_pc(real_t Sw, const PcParams&, real_t& pc)

7.4 fvm
- build_face_trans(const Mesh&, const RockProps&, std::vector<real_t>& T_face)
- mass_flux_face(face_id, states, T_face, gravity, out_flux)
- energy_flux_face(face_id, states, lambda_face, out_flux)
- assemble_cell(const Mesh&, cell_id, fluxes, sources, SpMat& A, Vector& b)

7.5 edfm
- mf_trans(const Mesh&, const FracPiece&, const RockProps&, const FracProps&) -> real_t
- ff_trans(const FracPiece& a, const FracPiece& b, const FracProps&) -> real_t
- assemble_nnc(const EqSelector&, const std::vector<NNC>&, const States&, SpMat& A, Vector& b)

7.6 bi
- apply_dirichlet(SpMat& A, Vector& b, const Set&, real_t value)
- apply_neumann(Vector& b, const Set&, real_t flux)
- init_fields(Fields&, const IC&)

7.7 solver
- step_fully_implicit(State&, Δt, Ctrl&, Stats&)
- step_impers(State&, Δt, Ctrl&, Stats&)
- time_loop(DriverConfig&, Callbacks&)

7.8 poster
- write_tecplot(const Mesh&, const std::vector<FracPiece>&, const Fields&, const std::string& path)

----------------------------------------
## 8. 多物理问题的装配与选择

- 通过 `ProblemConfig` 选择方程集：
  - 单相热: enable_mass=true, enable_energy=true, phases=1
  - 两相热: enable_mass=true, enable_energy=true, phases=2
- `EqSelector` 控制装配行列（p, S_w, T）以及裂缝自由度
- `AssembleContext` 缓存物性、上风侧、T_{pq}，减少重复计算

----------------------------------------
## 9. 线性求解与收敛

- 稀疏矩阵: Eigen::SparseMatrix<real_t>；预留非零模式
- 线性解算器: BiCGSTAB/GMRES + ILU(0) 或外部库占位接口
- 非线性: Newton 带阻尼；残差范数与变量相对增量双判据

----------------------------------------
## 10. 性能与鲁棒性

- SoA/紧致索引；避免临时对象；热点分支消除
- 面/连接遍历顺序与内存局部性匹配
- 几何与 T_{pq} 可在时间步外缓存（随物性非线性更新时再校正）
- 上风选择统一函数，确保质量/能量一致性

----------------------------------------
## 11. 测试与验证（建议）

- 单元测试：几何相交、T_{ij}/T_{i,f} 正确性、守恒性
- 基础算例：
  - 单相稳态与导热稳态解析解
  - 含单裂缝压力/温度降（对比等效渗透率）
  - 两相 Buckley–Leverett 式前沿（无裂缝）
  - 含裂缝两相驱替的体积分数守恒检查

----------------------------------------
## 12. 编码示例（片段）

// 计算 MF 连接通量系数（示意）
real_t mf_trans(const Mesh& mesh, const FracPiece& fp, const RockProps& rock, const FracProps& frac) {
    const auto& K = rock.K_cell(fp.cell);
    Eigen::Vector3d n = fp.n.normalized();
    real_t kappa_i = (n.transpose() * K * n).value();
    real_t kappa_f = frac.k_normal(fp); // b^2/12 等价
    real_t d_i = geom_distance_cell_to_piece(mesh, fp.cell, fp);
    real_t d_f = 0.5 * std::max(real_t(1e-9), fp.b);
    real_t A   = fp.area_or_len;
    return A / (d_i / kappa_i + d_f / kappa_f);
}

//----------------------------------------
## 13. 代码评审清单

- 接口是否与《路书》一致（命名/输入/输出）
- 守恒性：质量与能量方程上风一致；NNC 与 MM 装配符号正确
- 物性调用集中、无重复计算；异常输入防护
- 可读性：头注齐全；复杂几何有注释；命名统一
- 性能：热点处无多余分配/拷贝；缓存命中率考虑

----------------------------------------
## 14. 入门实操模板与检查表（逐步完成）

14.1 目录与命名空间（先用 VS 过滤器分组即可）
- [ ] 命名空间：`edfm` 顶层；子空间 `edfm::mesh, ::frac, ::fvm, ::edfm, ::prop, ::bi, ::solver, ::poster`
- [ ] 头文件骨架：
  - [ ] `include/FVM/Fields.h`（`real_t/index_t`, `Field`, `TensorField`）
  - [ ] `include/FVM/Schemes.h`（`TimeScheme`, `DiffScheme`, `NonOrthoCorrection`）
  - [ ] `include/FVM/Context.h`（`AssembleContext`, `PropertiesAccessor`, `BoundaryAndInitial` 前置声明）
  - [ ] `include/FVM/Terms.h`（`Term` 抽象与 `ddt/diff/source` 工厂）
  - [ ] `include/FVM/Equation.h`（收集项与 `assemble`）

14.2 注释模板（Doxygen）
/**
 * brief: 装配扩散项 ∇·(K/μ ∇u) 的 TPFA 系数
 * formula: T = A / (d_i/κ_i + d_j/κ_j), F = -T(Δu - ρ g Δz)
 * note: 非正交修正采用最小修正法，可关闭
 * @param ctx   网格/物性/时间步等上下文
 * @param asmb  稀疏装配器（add_to_A, add_to_b）
 */

14.3 函数签名约定（与《路书》一致）
- `ddt(const Field& u, const TimeScheme& ts)` → 对角系数 + 旧时刻 RHS
- `diff(const TensorField& KoverMu, const Field& u, const DiffScheme&, const NonOrthoCorrection&)`
- `source(const Field& q)` → 仅 RHS
- `Equation::assemble(const AssembleContext&, SparseAssembler&)` → 先 LHS 后 RHS 再边界

14.4 SparseAssembler 约定
- [ ] 接口：`add_to_A(i,j,val)`, `add_to_b(i,val)`, `apply_dirichlet(i,val)`, `apply_neumann(face,flux)`
- [ ] 后端：Eigen Triplet + `SparseMatrix<real_t>`
- [ ] 边界顺序：Dirichlet 强制后再加 Neumann 到 b

14.5 适配旧实现（先打通闭环）
- [ ] `diff.apply` 内部可调用现有 `Conv_FirstOrder_ComputerBoundaryFace.h` 的面通量例程（做参数映射）
- [ ] 求解保持现有 `LinearSolver_Eigen.h`，新框架仅提供 A,b
- [ ] `Solver_TimeLoopDriver.h` 中以新 `Equation` 门面替换旧装配调用

14.6 最小闭环任务清单
- [ ] 实现 `source.apply` 与 `ddt.apply`（仅对角）
- [ ] 实现 `diff.apply` 的 MM TPFA（可暂不含非正交修正）
- [ ] `Equation::assemble` 与边界强制逻辑
- [ ] 小算例验证（单相稳态）：总流入=总流出，<1e-8×规模

14.7 扩展顺序（建议）
- [ ] 非正交修正开关（最小修正法）
- [ ] 能量方程通量一致性（对流/扩散同上风）
- [ ] EDFM：MF → FF 连接为独立 Term 并装配到矩阵
- [ ] 两相 IMPERS：压力步/饱和度步/温度步三段式
- [ ] 3D 几何量与裂缝多边形处理

14.8 示例占位（仅接口调用）
```cpp
using namespace FVM;
Equation eq;
eq += ddt(p, TimeScheme{TimeScheme::Type::BackwardEuler});
eq += diff(K_over_mu, p, DiffScheme{DiffScheme::Type::TPFA_CDS},
           NonOrthoCorrection{NonOrthoCorrection::Type::MinCorrection});
eq == source(qm);
eq.assemble(ctx, asmb); // 产生 A, b 后交给现有 Solver
```

