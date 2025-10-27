# 非结构化网格嵌入式离散裂隙模型搭建与求解“路书”

版本: v0.1 (草案)

目标: 给出从网格/裂缝建模到方程离散、EDFM耦合、时间推进与外迭代、以及输出的完整模块化流程与公式-伪代码，覆盖四类问题：
- 2D 基岩+复杂缝网 CO2 单相渗流-换热
- 2D 基岩+复杂缝网 CO2-H2O 两相渗流-换热
- 3D 基岩+复杂缝网 CO2 单相渗流-换热
- 3D 基岩+复杂缝网 CO2-H2O 两相渗流-换热

----------------------------------------
## 1. 模块全景与命名前缀

统一模块及前缀（目录与命名建议）：
- Mesh_*: 基岩网格读入、几何量、拓扑、面信息与邻接
- FracMesh_*: 裂缝网络读入与几何处理（线段/多边形与单元相交、分段）
- Fields: 主变量与物性字段的存储与访问（网格单元/面/裂缝段）
- Properties_*: 流体与岩石物性计算（ρ, μ, c_p, λ 等；相对渗透率、毛管等）
- FVM_*: 有限体积离散与矩阵/源项组装（面通量、扩散通量、源项）
- EDFM_*: 裂缝-基岩、裂缝-裂缝非邻接连接(NNC)与耦合通量/源项组装
- BoundaryandInitial_*: 初边值条件施加与更新
- Solver_*: 时间推进与外迭代框架（Fully-Implicit、IMPES/IMPERS、非线性求解）
- Poster_*: 输出（Tecplot/VTK等），与后处理辅助

数据命名（建议）：
- 索引: `index_t`；标量: `real_t`
- 单元/面/裂缝字段采用结构化数组（SoA）以利缓存
- 字段命名：`Fields::p`, `Fields::T`, `Fields::Sw` 等；面量如 `Fields::face_flux`

----------------------------------------
## 2. 四类物理模型与主变量

符号约定：
- ϕ 孔隙度；K 绝对渗透率张量；k_rα 相对渗透率；μ 动黏度；ρ 密度；
- c_p 比热容；λ 有效导热系数（含流体与岩石）；g 重力；z 高度；
- p 压力；u 速度；T 温度；S_w 水相饱和度；S_n 非润湿相（CO2）饱和度；p_c 毛管压；
- 下标 α∈{w, n} 表示相；在裂缝域用下标 f；

2.1 单相（2D/3D）
- 质量守恒: ∂(ϕρ)/∂t + ∇·(ρu) = q_m
- Darcy 定律: u = - (K/μ) (∇p - ρ g ∇z)
- 能量方程（体积平均）:
  ∂[(ϕρ c_p + (1-ϕ)ρ_r c_{p,r})T]/∂t + ∇·(ρ c_p u T - λ ∇T) = q_h

主变量: Fully-Implicit 可取 {p, T}；分裂类框架可先压强场后温度。

2.2 两相 CO2-H2O（2D/3D，非混相近似）
- 相质量守恒（各相）: ∂(ϕ ρ_α S_α)/∂t + ∇·(ρ_α u_α) = q_α
- 相速率: u_α = - (K k_{rα}(S_w)/μ_α) (∇p_α - ρ_α g ∇z)
- 毛管关系: p_n - p_w = p_c(S_w)
- 约束: S_w + S_n = 1
- 能量方程（相对流-导热耦合）:
  ∂[C_T T]/∂t + ∇·(h_w ρ_w u_w + h_n ρ_n u_n - λ ∇T) = q_h
  其中 C_T = ϕ(ρ_w c_{p,w} S_w + ρ_n c_{p,n} S_n) + (1-ϕ)ρ_r c_{p,r}
  h_α 可取 c_{p,α} T 或更一般焓模型。

主变量选择：
- IMPES/IMPERS: {p, S_w, T}；压力以油水平均压或某相压为主变量，配合 p_c(S_w)
- Fully-Implicit: {p, S_w, T} 同步 Newton

----------------------------------------
## 3. 时空离散与装配框架（FVM_*)

3.1 时间离散
- 后向欧拉（Fully-Implicit）: (X^{n+1}-X^n)/Δt + R(X^{n+1}) = 0
- IMPES/IMPERS: 压力隐式，饱和度/温度显式或半隐式（可对 T 采用隐式）
- 自适应时间步：基于非线性迭代计数、残差衰减、CFL（对显式分量）

3.2 空间离散（非结构网格，2D/3D）
- 有限体积积分：
  对单元 i，质量方程离散为：
  (V_i (ϕρ)_i^{n+1} - V_i (ϕρ)_i^{n})/Δt + Σ_{f∈∂i} (ρ u_n A)_f^{n+1} = Q_{m,i}
  其中面法向流量 u_n A 采用两点通量近似（TPFA），非正交修正可选：
  T_{ij} = A_f / (d_i/κ_i^n + d_j/κ_j^n)，κ_i^n = n_f·K_i·n_f
  通量：F_{ij} = - T_{ij} (p_j - p_i - ρ̄ g (z_j - z_i))，上风 ρ̄
- 能量方程：对流项使用与质量守恒一致的上风；扩散项用面导热系数 λ_f 的调和平均。

3.3 非线性与雅可比
- Fully-Implicit: 组装残量 R(x) 与稀疏雅可比 J = ∂R/∂x；线性解算器用 Eigen（见 LinearSolver_Eigen.h）
- 阻尼/线搜索；收敛准则（相对与绝对残差、变量增量）

----------------------------------------
## 4. EDFM_*：嵌入式离散裂缝耦合

4.1 几何嵌入与 NNC
- 2D: 裂缝为线段；3D: 裂缝为多边形。对每个网格单元 i，计算裂缝与单元多边形交线（2D 为线段，3D 为多边形片）
- 构建三类连接：
  - MM: 基岩-基岩常规面连接（来自 Mesh_*）
  - MF: 基岩-裂缝 非邻接连接（NNC）
  - FF: 裂缝-裂缝 NNC（同单元内相邻裂缝片）

4.2 裂缝通量与等效参数
- 裂缝内流动近似平行板：k_f ≈ b^2/12（b 为孔隙开度/缝宽），各向异性时取 n_f·K_f·n_f
- 面向 MF 连接的“半通量系数”与调和平均：
  T_{i,f} = A_{i,f} / (d_i/κ_i + d_f/κ_f)
  - A_{i,f}: 单元 i 中裂缝片与单元的交线长度(2D)/面积(3D)
  - d_i: 单元体心到裂缝片等效连线的几何距离；d_f: 裂缝半厚度（2D 取 b/2；3D 取法向半厚度）
  - κ_i = n_{i,f}·K_i·n_{i,f}；κ_f = n_{i,f}·K_f·n_{i,f}
- FF 连接：对同一单元内两裂缝片 a,b，
  T_{a,b} = A_{a,b} / (d_a/κ_a + d_b/κ_b)
  A_{a,b} 为两片公共边长(2D 为交点等效，常用局部插值长度；3D 为公共线段长度 or 邻接面积近似)；d_a,b 为片心到公共界面的距离
- 重力修正：对任一连接 pq，通量 F_{pq} = - T_{pq} (p_q - p_p - ρ̄ g (z_q - z_p))
- 两相：替换 κ → κ λ_α，其中 λ_α = k_{rα}(S)/μ_α 为相动率；多相混配采用上风选择

4.3 方程耦合与装配要点
- MF/FF 连接产生额外的系数与源项，等价于在稀疏矩阵中添加非邻接条目
- 单相质量：对单元 i 的残量，加入 Σ_f F_{i,f}；对裂缝自由度 f 的残量加入 -Σ_i F_{i,f}
- 能量一致：对流项采用相同上风，扩散项在裂缝内使用 λ_f（可含导热增强）
- 两相：对每相通量分别装配，或采用压力-饱和度分裂形式（见 6 节）

----------------------------------------
## 5. Properties_*：物性与本构

建议抽象：
- 流体：ρ(p,T), μ(p,T), c_p(p,T), h(T)；相对渗透率 k_{rα}(S)；p_c(S)
- 岩石：ϕ(p,T) 可压缩；K 可各向异性；ρ_r, c_{p,r}, λ_r
- 裂缝：k_f(b), ϕ_f≈1（或给定）；热导 λ_f；相对渗透率可与基岩不同

示例（可替换更高保真模型）：
- k_{rα}(S): Corey 型；p_c(S): Brooks–Corey or van Genuchten；
- λ 有效：λ_eff = ϕ Σ_α S_α λ_α + (1-ϕ) λ_r

----------------------------------------
## 6. Solver_*：时间推进与外迭代框架

6.1 Fully-Implicit（单相/两相/热）
伪代码：
1) t=0 初始化 Fields；
2) while t < t_end:
   - 预测 x^{n+1,0}（延用上步值）
   - for k=0..k_max:
       a) 评估 Properties_* at x^{n+1,k}
       b) FVM_* 组装矩阵 A 与右端 b（含 EDFM_* NNC）
       c) 线性求解 A δx = -R
       d) x^{n+1,k+1} = x^{n+1,k} + ω δx（线搜索）
       e) 检查收敛（残差与增量）
   - 接受步长，t+=Δt；自适应调整 Δt

变量选择：
- 单相热：x=[p, T]
- 两相热：x=[p, S_w, T]；采用 p_c(S_w) 与 λ_α(S_w) 装配耦合雅可比

6.2 IMPES/IMPERS（两相热）
- 压力步（隐式）：
  - 以总动率 λ_t = Σ_α λ_α 与加权密度上风，求解压力或平均压 p
  - 求 u_t 与相速度 u_α = f_α u_t - K λ_α ρ_α g ∇z，其中 f_α = λ_α/λ_t
- 饱和度步（显式或半隐式）：
  - 采用迎风/限定器：S^{n+1} = S^n - (Δt/V) Σ_f (F_w)_f
- 温度步（可隐式 IMPERS）：
  - 若与压力弱耦合，可对 T 做隐式解（扩散主导），或半隐式保守对流项
- 迭代外循环：若热-流强耦合，每时间步内进行 P→S→T→校正 的弱耦合循环至收敛

----------------------------------------
## 7. 边界/初始条件（BoundaryandInitial_*）

类型：
- 压力/温度 Dirichlet；通量/热通量 Neumann；混合 Robin（换热）
- 两相边界：指定相压/总压与饱和度/分数流；温度同理
实现：将 Dirichlet 以强制法装入矩阵（对角置 1，右端置边值），或弱式转移为源项

----------------------------------------
## 8. Poster_*：输出

- Tecplot/VTK 导出：单元标量/向量场，裂缝片（线段/多边形）场
- 建议与现有 `TecplotExporter.h` 对齐，提供分区导出（基岩、裂缝）

----------------------------------------
## 9. 子函数清单（名称/输入/输出/要点）

9.1 Mesh_*
- Mesh_BuildTopology(mesh) → {cells, faces, neighbors, geom}
  - 输入：原始网格点与单元连接
  - 输出：面法向、面积/长度、中心、邻接(i,j)、单元体积/面积

9.2 FracMesh_*
- FracMesh_IntersectCell(cell_poly, frac_entities) → {segments or polygons in cell}
  - 输出：每单元内裂缝片列表，几何量（长度/面积、法向、中心）
- FracMesh_BuildNNC(mesh, frac_in_cells) → {MF_list, FF_list}
  - 输出：NNC 连接对与其几何量 A, d, n

9.3 Properties_*
- Prop_Fluid_CO2(p,T) → {ρ, μ, c_p, λ}
- Prop_Fluid_H2O(p,T) → {ρ, μ, c_p, λ}
- Prop_Rock(params) → {ρ_r, c_{p,r}, λ_r, ϕ(p)}
- RelPerm_Corey(Sw, params) → {k_rw, k_rn}
- Capillary_BrooksCorey(Sw, params) → p_c

9.4 FVM_*（基岩 MM 通量）
- FVM_BuildFaceTrans(mesh, K) → {T_ij per face}
- FVM_MassFlux(face, states) → F_{ij} = -T_{ij}(Δp - ρ̄ g Δz)（上风 ρ̄）
- FVM_EnergyFlux(face, states) → (ρ c u T)_upwind - λ_f ∂T/∂n A
- FVM_AssembleCell(i, faces, fluxes, sources) → 填充 A,b

9.5 EDFM_*
- EDFM_MatrixFractureT(cell_i, frac_piece) → T_{i,f}（含 κ、A、d 计算）
- EDFM_FractureFractureT(frac_a, frac_b) → T_{a,b}
- EDFM_AssembleCoupling(eqn, connections, states) → 在稀疏矩阵中添加 NNC 条目

9.6 BoundaryandInitial_*
- BI_ApplyDirichlet(A,b,set, value)
- BI_ApplyNeumann(b,set, flux)
- BI_InitFields(fields, IC)

9.7 Solver_*
- Solver_FullyImplicit_Step(fields, Δt) → {converged, stats}
- Solver_IMPERS_Step(fields, Δt) → 同上
- TimeStepControl(update_rules)

9.8 Poster_*
- Poster_WriteTecplot(mesh, frac, fields, path)

----------------------------------------
## 10. 关键公式到可嵌入表达式（示例）

10.1 面通量（单相）：
- κ_i = n·K_i·n, κ_j = n·K_j·n
- T_{ij} = A / (d_i/κ_i + d_j/κ_j)
- ρ̄ = upwind(ρ_i, ρ_j, sign(p_j - p_i - ρ̄ g Δz)) 迭代中可用 Picard 固定上风
- F_{ij} = - T_{ij} ( (p_j - p_i) - ρ̄ g (z_j - z_i) )

10.2 MF/FF（单相）：
- T_{i,f} 与 T_{a,b} 形式同上，参数替换为裂缝等效值
- 质量项装配：R_i += Σ_f F_{i,f}；R_f -= Σ_i F_{i,f}

10.3 两相动率与分裂：
- λ_α = k_{rα}(S_w)/μ_α；λ_t = Σ_α λ_α；f_α = λ_α/λ_t
- 压力步通量：F_t = - T (Δp - ρ̄_t g Δz)，ρ̄_t 上风总密度；相通量 F_α = f_α F_t - T λ_α (ρ_α g Δz)
- 饱和度更新（显式）：S_w^{n+1} = S_w^n - (Δt/V) Σ_f (F_w)_f

10.4 能量通量：
- 对流： (ρ c_p u T)_f 采用压力/速度上风侧的 T 与 c_p
- 扩散： q_cond = - λ_f (ΔT/Δn) A，λ_f 调和平均

----------------------------------------
## 11. 主流程伪代码（四问题通用骨架）

Setup:
- Mesh_* 构建基岩网格拓扑与几何
- FracMesh_* 读入裂缝，并与单元相交生成片，构建 NNC（MF/FF）
- Fields 分配 {p, T, [S_w], [裂缝变量]} 并初始化
- 准备 Properties_* 参数与表；BoundaryandInitial_* 施加边界与初值

Time Loop:
- for t in [0, t_end):
  - 选择 Solver_*（FI 或 IMPES/IMPERS）推进一步
  - 动态调步；周期性 Poster_* 输出

Fully-Implicit Assemble (per Newton iter):
- 评估物性于单元/面/裂缝片（上风需要面向）
- FVM_* 计算 MM 面通量系数与能量扩散
- EDFM_* 计算 MF/FF 通量系数
- 组装质量与能量残量和雅可比（包含 NNC）
- 线性求解与更新

----------------------------------------
## 12. 维度扩展与实现注意

- 2D→3D：面积/体积、法向与几何距离计算替换；裂缝从线段→多边形片
- 非正交/各向异性：如需精度提升，引入非正交修正项与全张量 K 处理
- 稳定性：两相与热的上风、限制器、相容性（能量与质量伴随）
- 可扩展性：SoA 数据布局、稀疏模式缓存、跨步重用 T_{pq}

----------------------------------------
## 13. 与现有代码对齐（参考）

- 现有文件：
  - `Solver_TimeLoopDriver.h`: 可映射至 Solver_* 外层驱动
  - `TimeIterm_Euler_SinglePhase_PressureEq.h`, `TimeIterm_Euler_SinglePhase_TemperatureEQ.h`: 对应单相全隐/半隐步进组件
  - `Conv_FirstOrder_ComputerBoundaryFace.h`: 一阶迎风对流计算可复用于面通量
  - `TecplotExporter.h`: Poster_* 输出接口

----------------------------------------
## 14. 清单检查（交付验收）

- 四类问题均有：控制方程、未知量定义、离散表达、NNC 耦合表达
- 模块与子函数命名、输入/输出、装配伪代码齐备
- 求解框架：FI 与 IMPES/IMPERS 两套均可落地
 - 与现有代码接口对齐，并可逐步替换/增强

----------------------------------------
## 15. 入门填充模板与检查表（逐步完成）

本节提供四类问题按模块分解的“填空式模板”。建议一次只填完一个问题的一到两个模块，并以算例验证后再前进。

15.1 选择问题与主变量
- [ ] 问题选择：
  - [ ] 2D 单相热（p, T）
  - [ ] 2D 两相热（p, S_w, T）
  - [ ] 3D 单相热（p, T）
  - [ ] 3D 两相热（p, S_w, T）
- [ ] 主变量：填写本问题用到的变量集合（如 p/T 或 p/Sw/T）
- [ ] 时间推进：Fully-Implicit / IMPERS（压力隐式、饱和度显/半隐、温度隐式）

15.2 Mesh_/FracMesh_（几何与 NNC）
- [ ] Mesh_BuildTopology：输入（点、连通）、输出（cell/face/邻接/几何量）：
- [ ] FracMesh_IntersectCell：输入（cell 多边形、裂缝对象）、输出（片段/法向/长度或面积）：
- [ ] FracMesh_BuildNNC：输出（MF/FF 列表，A/d/法向）：

15.3 Properties_（物性与本构）
- [ ] 流体模型：ρ(p,T), μ(p,T), c_p(p,T), λ(p,T)
- [ ] 相对渗透率：k_rw(Sw), k_rn(Sw)（Corey 参数）：
- [ ] 毛管压：p_c(Sw)（Brooks–Corey/vG 参数）：
- [ ] 岩石：ϕ(p), K 张量, λ_r, ρ_r, c_{p,r}

15.4 控制方程与离散（FVM_）
- [ ] 方程清单：
  - 质量：
  - 能量：
  - 两相附加：p_c(Sw)、λ_α(Sw)
- [ ] 算子表达（按“ddt+diff(+div)+gravity=source”形式写）：
- [ ] 时间离散：后向欧拉/步长控制：
- [ ] 空间离散：TPFA/中心差分；非正交修正（最小修正法）：
- [ ] 面量上风：密度/动率的上风侧定义：
- [ ] 雅可比结构（FI 情况下的变量顺序与耦合项）：

15.5 EDFM_（裂缝耦合）
- [ ] 几何：单元内裂缝片（中心、法向、长度/面积、缝宽 b）：
- [ ] MF 连接：T_{i,f} = A / (d_i/κ_i + d_f/κ_f)；κ_f = n·K_f·n（b^2/12）：
- [ ] FF 连接：T_{a,b} 定义（公共界面长度/距离）：
- [ ] 重力修正：F_{pq} = -T_{pq}(Δp - ρ̄ g Δz)：
- [ ] 两相：κ→κ λ_α，上风选择：
- [ ] 组装：矩阵非邻接条目与 RHS 贡献：

15.6 BoundaryandInitial_（定解条件）
- [ ] Dirichlet 集合和值：
- [ ] Neumann 集合与通量：
- [ ] Robin（若有）：
- [ ] 初值设定：
- [ ] 施加策略：强制法（对角置 1）/转源项：

15.7 Solver_（时间步与外迭代）
- [ ] Fully-Implicit 框架：
  - 变量：
  - Newton：阻尼、收敛准则：
  - 线性解：BiCGSTAB/GMRES + 预处理：
- [ ] IMPERS 框架：
  - 压力步：λ_t/ρ̄_t、u_t 计算：
  - 饱和度步：迎风/限定器：
  - 温度步：隐式/半隐式：
  - 弱耦合外循环与停止准则：

15.8 Poster_（输出）
- [ ] 输出变量：p/T/[Sw]（单元/面/裂缝片/节点后处理）：
- [ ] Tecplot/VTK 路径与频率：
- [ ] 守恒检查与日志：

15.9 验证算例（每引入一个模块即验证）
- [ ] 单相稳态压差（无裂缝）：总流入=总流出（<1e-8×规模）
- [ ] 单裂缝等效导流（MF）：压降曲线对比参考：
- [ ] 能量稳态导热：解析核对：
- [ ] 两相驱替（无裂缝）：Buckley–Leverett 前沿：
- [ ] 含裂缝两相：体积分数守恒：

15.10 开放问题与记录
- [ ] 上风一致性（质量/能量用同一上风侧）
- [ ] 非正交修正敏感性
- [ ] 时间步自适应策略与阈值
- [ ] 线性求解鲁棒性与预处理选择

