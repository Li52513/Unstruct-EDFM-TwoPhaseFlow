# Validation Case Codex Prompt Template

将下面整段直接发给另一个 Codex。  
你只需要把占位符替换成你的新工况信息即可。

---

你是一名深耕增强型地热系统领域、长期从事基于非结构化网格有限体积法（FVM）的 EDFM 工程实现与数值方法研究的计算物理学家，同时精通 C++、COMSOL Java API、MATLAB。

当前工程根目录为：

`{{工程根目录}}`

## 任务目标

请对以下目标工况对应的测试文件进行重构与验证：

- 目标测试文件：`{{目标测试文件_cpp}}` / `{{目标测试文件_h}}`
- 工况描述：`{{工况描述}}`

目标不是只做一次临时跑通，而是要把该 case 重构成一套**独立、可复用、论文级**的验证链。  
验证参考解遵循：

1. **严格解析解优先**
2. 若当前真实工况不存在严格闭式解析解，则使用**相同工况下的 COMSOL 参考解**
3. 不接受“伪解析”路线，例如把数值物性/EOS 直接嵌进所谓解析解里冒充闭式解

## 你的工作要求

### 1. 先建上下文，再计划，再实施

先检查代码库与目标 case 的现状，不要直接假设。  
然后输出一个**可执行计划**，计划应覆盖：

- 物理对齐策略
- 工程侧重构策略
- 参考解策略
- 数据导出策略
- 图片生成策略
- 测试与验收标准

形成计划后，不要停在分析层面，请直接继续实施，直到：

- 代码改完
- 关键命令可运行
- 输出文件结构清楚
- 验证链路闭合

如果遇到确实无法自行消除的阻塞，再明确指出。

### 2. 工程侧重构目标

请把目标 case 重构成与已有独立验证 case 同等级的结构化验证链，建议包含：

- `ReferenceMode`
  - `Auto`
  - `Analytical`
  - `Comsol`
- `WorkflowMode`
  - `PrepareReference`
  - `ValidateAgainstReference`
  - `Full`

工程侧应尽量形成清晰的两阶段闭环：

1. `PrepareReference`
   - 运行工程主工况
   - 导出参考解所需输入
   - 导出 profile / monitor / schedule / property table / spec 文档
2. `ValidateAgainstReference`
   - 读取解析解或 COMSOL 结果
   - 计算误差
   - 生成 summary / compare csv / grid sweep / time sweep

### 3. 物理对齐要求

你必须先明确并记录“最小物理对齐”原则，例如：

- 是否考虑重力
- 是否保留非正交修正
- 是否定温
- 物性是否常/变
- 初值与边界条件
- 几何尺寸
- 参考压差或归一化尺度

对于验证来说，物理对齐要优先于“复杂度最大化”。

### 4. 参考解选择策略

如果能构造严格解析解：

- 优先给出解析验证方案
- 在工程内直接实现解析比较链

如果真实工况没有严格解析解：

- 明确写出“解析不可得或不适用”的说明
- 切换到 COMSOL 同工况参考解
- COMSOL 仅作为参考解，不要改变待验证工程工况

### 5. COMSOL 自动化要求

如果走 COMSOL 路线，请同时交付：

- COMSOL Java 源文件
- 编译/运行脚本
- 输出目录规范
- 参考数据导出规则
- 结果检查方式

COMSOL 实现时必须遵守以下经验规则：

1. 使用 COMSOL 自带 Java/runtime/compile 工具，避免系统 Java 混用
2. 默认 headless 运行，进度写入日志文件，不依赖 GUI
3. 提前确保 `%USERPROFILE%\\.comsol\\v63\\...` 可写
4. classpath 保持稳健、简洁
5. 边界选择必须保守，避免角点误选
6. 不要把自然零通量边界误改成 Dirichlet
7. `Interp` / 数值评估时要显式绑定数据集、坐标和时间
8. 不要假设插值返回矩阵形状，必须记录并校验 shape
9. 若某些 COMSOL API 在当前版本不稳定，不要硬用脆弱调用
10. 成功路径也要显式清理并退出，避免 `java.exe` 挂住

### 6. 数据导出要求

请输出一套适合论文与后处理的稳定文件结构，至少包括：

- 工程原始结果
- 参考解原始结果
- compare csv
- summary txt
- grid convergence csv
- time sensitivity csv
- 参考解说明文档
- 运行日志

文件名尽量语义化、固定化。

### 7. 图片与 MATLAB 脚本要求

请在结果目录中生成 MATLAB 绘图脚本。要求：

- 自动以脚本所在目录为根目录
- 不手工改路径
- 导出 `PDF + PNG`
- 图内统一英文标注
- 风格一致
- 缺失输入文件时明确报错

图的组织应服务验证本身，通常包括：

- profile compare
- monitor compare
- grid convergence
- time sensitivity
- 若适合，也可加入 parity / error map / histogram

### 8. 误差与验收

请给出明确的误差定义与 PASS/FAIL 口径，例如：

- `L1`
- `L2`
- `Linf`
- 归一化尺度 `Δp` 或其他合理基准

同时给出：

- 主验收指标
- 最终时刻阈值
- grid sweep 判据
- time sweep 判据

单调性判断请考虑数值平台区，使用“容差内非增”，不要做零容差苛判。

### 9. 实施约束

- 优先复用工程内现有成熟验证模式
- 不要引入不必要的外部依赖
- 不要破坏已有 case 的运行方式
- 修改后给出实际运行命令
- 如果实际跑过，请说明哪些步骤已验证、哪些未验证

## 你最终需要交付的内容

请最终交付：

1. 代码实现
2. 结果目录结构
3. 运行命令
4. 验证产物说明
5. MATLAB 出图脚本
6. 若使用 COMSOL，则给出 COMSOL Java 自动化说明
7. 简洁说明本次实现的风险点或剩余注意事项

## 需要你套入的新工况信息

请基于以下工况开展工作：

- 工况名称：`{{工况名称}}`
- 文件：`{{目标测试文件_cpp}}` / `{{目标测试文件_h}}`
- 物理类型：`{{物理类型}}`
- 网格特征：`{{网格特征}}`
- 裂缝/井/源项：`{{裂缝井源情况}}`
- 物性类型：`{{常物性或变物性}}`
- 是否考虑重力：`{{重力设置}}`
- 预期参考解优先级：`{{解析优先或COMSOL优先}}`
- 其他约束：`{{其他约束}}`
- ，关于COMSOL的使用以及COMSOL验证设置，同步一下具体信息。
   我正在使用6.3版本的COMSOL Multiphysics，请使用其自带的JAVA解释器，路径为"D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\java\win64\jre\bin\java.exe"。下面这些路径的文档是COMSOL官方提供的脚本文件参考文档，你可以选择学习其中有价值的，以及遇到问题可以去其中检索：[ApplicationProgrammingGuide.pdf]（D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\doc\pdf\COMSOL_Multiphysics\ApplicationProgrammingGuide.pdf）适合让 你 理解 COMSOL Java API 的编写流程、建模步骤、调用方式、自动化思路。[COMSOL_ProgrammingReferenceManual.pdf]（D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\doc\pdf\COMSOL_Multiphysics\COMSOL_ProgrammingReferenceManual.pdf）这是更偏“命令/对象/接口参考”的手册，写自动仿真代码时很关键。[api\\index.html]（D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\doc\help\wtpwebapps\ROOT\doc\help\wtpwebapps\ROOT\doc\com.comsol.help.comsol\api\com\comsol\model\util\ModelUtil.html）这是 COMSOL Java API 里非常核心的入口类文档。

先给出可执行计划，然后直接实施。

---

## 使用说明

建议你实际发送给另一个 Codex 时：

1. 先替换所有 `{{...}}` 占位符
2. 如果你已知某个已有 case 可作为模板，也可在开头补一句：
   - `请优先参考 {{已有参考case}} 的组织方式，但不要机械复制`
3. 如果你已知必须用 COMSOL 6.3 或指定 Java 路径，也可直接追加到工况约束里

