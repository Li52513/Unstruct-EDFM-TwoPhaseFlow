## 快速说明 — 给 AI 编码代理的操作要点

本仓库是一个基于非结构化网格 + EDFM 的二维/准三维（2D/3D 扩展）多物理场单相/两相流求解器，核心用 C++ 实现（Visual Studio/msbuild 工程）。下面是能让 Agent 立即上手的可执行、代码级别约定和路径。

### 快速运行（开发迭代）
- 打开 `2D-Unstr-Quadrilateral-EDFM.sln` 用 Visual Studio 运行，或在 Developer Command Prompt 下使用：
  - `msbuild 2D-Unstr-Quadrilateral-EDFM.sln /p:Configuration=Debug`（Windows）
- 可在 `main.cpp` 中把网格划分数量调小（例：sectionNumX/Y = 3）用于快速本地迭代。

### 项目结构与“把握全局”要点
- 入口：`main.cpp` — 启动顺序与典型工作流位于此。重要步骤：构造 `MeshManager` -> 添加断层 -> `DetectAndSubdivideFractures()` -> 初始化物性`PhysicalPropertiesManager` -> `Initializer` 填充场 -> 组装并求解时间推进（`runTransient_singlePhase`）。
- 网格与裂缝管理：`MeshManager.h/.cpp`, `Mesh.h/.cpp`, `FractureNetwork.h/.cpp`。
- 求解器汇编（矩阵/稀疏格式）：`Solver_AssemblerCOO.h`（COO 构建、compress、toCSR、buildUnknownMap）。Assembly 函数示例：`assemblePressure_singlePhase_COO`、`assembleTemperature_singlePhase_COO`。
- 后处理/导出：`PostProcessor.cpp` 会将矩阵/裂缝状态导出为 CSV（目录默认 `out\matrix`、`out\fracture`），并包含 `times.csv` 索引。

### 重要命名与约定（必须保持一致）
- Face fields（面场）命名示例：`a_f_Diff_p_w`, `s_f_Diff_p_w`, `a_f_Diff_T`。Assembly 代码直接通过这些字符串访问面场。
- Cell（体/分段）时间项：`aC_time_p`, `bC_time_p`, `aC_time_T`, `bC_time_T`。
- 裂缝场（段级）示例：`pf_w`, `Sf_w`, `Tf`（参见 `PostProcessor::exportFracture_` 中的 preferred 列顺序）。
- 索引/ID 规则：face arrays 在代码中通常用 `face.id - 1` 作索引；cell 的 `id < 0` 表示 ghost cell；若要构造未知量索引，参见 `buildUnknownMap`。

### 常见改动点与注意事项（对 Agent 很重要）
- 修改线性代数/组装时，优先改 `Solver_AssemblerCOO.h` 的汇编函数并添加或保持现有 field 名称；很多调用通过字符串约定耦合。
- 当添加新场（field）时：
  1) 在合适位置填充 `FieldRegistry` 或 `FaceFieldRegistry`（见 `main.cpp` 的初始化流程）；
  2) 确保名称长度/大小与目标（cell-count 或 segment-count）一致，否则 `PostProcessor` 自动检测列会忽略它们。
- 网格/断层检测：`MeshManager::DetectAndSubdivideFractures()` 会调用 `FractureNetwork` 的一系列步骤（交点检测、去重、细分），修改时应用现有 distanceMetric 枚举与 `subdivide` 策略保持兼容。

### 调试与输出习惯
- 为快速调试，保持 `sectionNumX/Y` 小并在 `main.cpp` 使用 `mgr.exportMesh("mesh")` / `mgr.exportFractures("fractures")` 输出网格用于可视化。
- 后处理脚本：仓库包含 MATLAB 绘图脚本（`plot_matrix_state.m`, `plot_fracture_state.m`, `plot_CO2_results.m`），CSV 格式与脚本字段名应保持一致。

### 小样例（汇编函数约束）
- `assemblePressure_singlePhase_COO` 要求注册有：face fields `a_face`、`s_face` 与 cell fields `a_time`、`b_time`。示例调用见 `assemblePressure_CO2_singlePhase_COO`。

### 搜索点（快速定位修改位置）
- 初始化与运行流程：`main.cpp`。
- 网格/断层逻辑：`MeshManager.*`, `Mesh.*`, `FractureNetwork.*`。
- 组装器与稀疏矩阵：`Solver_AssemblerCOO.h`。
- 后处理与导出：`PostProcessor.cpp`。

如果你需要我把 Agent 的策略进一步细化（例如，如何在修改汇编后自动添加回归验证脚本、或生成小规模运行的 CI 任务样例），告诉我你想要的输出格式与优先级。你希望我现在把它合并到仓库并运行一次简单构建验证吗？
