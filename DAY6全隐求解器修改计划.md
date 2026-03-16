### Day6 全隐式求解器「一次性全清单」代码修改计划

**Summary**
- 目标：一次性覆盖《终极清单》1-28 全项，并把当前“能跑但强烈锯齿/回滚”的状态修到“稳定收敛 + 可推进时间步 + 可扩展到 3D 两相”。
- 策略：按 12 个改动包执行，每个改动包都明确“改哪里、怎么改、解决哪条清单问题、完成判据”。

### 关键文件索引
1. F1: [RunGeneric_impl.hpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientEngine/RunGeneric_impl.hpp)  
2. F2: [StepKernels_impl.hpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientEngine/StepKernels_impl.hpp)  
3. F3: [FIM_BlockSparseMatrix.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_BlockSparseMatrix.h)  
4. F4: [StateSync.hpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientEngine/StateSync.hpp)  
5. F5: [BoundaryAssembler.cpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/BoundaryAssembler.cpp)  
6. F6: [FVM_Ops_AD.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FVM_Ops_AD.h)  
7. F7: [FIM_TransientCaseKit.hpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientCaseKit.hpp)  
8. F8: [Types.hpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientEngine/Types.hpp)  
9. F9: [AD_FluidEvaluator.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/AD_FluidEvaluator.h)  
10. F10: [TransmissibilitySolver_3D.cpp](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/TransmissibilitySolver_3D.cpp)  
11. F11: [FIM_TopologyBuilder3D.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TopologyBuilder3D.h) / [FIM_ConnectionManager.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_ConnectionManager.h)  
12. F12: [Well_WellControlTypes.h](D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Well_WellControlTypes.h)

### 修改清单（按改动包执行）
| 改动包 | 修改文件 | 怎么改 | 对应终极清单 |
|---|---|---|---|
| 01 时间步与牛顿策略基线 | F7, F8 | `long_profile.max_newton_iter=24`，`rel_update_tol=1e-3`，`dt_relres_grow_factor=1.05`，`dt_relres_soft_shrink_factor=0.90`，并默认开启 `enable_best_iter_guard=true`、`enable_stagnation_accept=true`。 | 1,2,3,4,5 |
| 02 收敛判据去“假收敛” | F1 | 分离 `raw_res` 与 `scaled_res`；`abs_res_tol` 用 `raw_res` 判据；行缩放仅用于线性求解稳健和诊断，不直接作为物理收敛唯一判据。 | 1-5（控制策略失效根因） |
| 03 岩石累积项去硬编码+压缩性 | F1, F4 | 去掉 `phi=0.2/cp_r=1000/rho_r=2600` 硬编码；从场读取 `phi/rho_r/cp_r/c_r`；引入 `phi_ad=phi_ref*(1+c_r*(P-P_ref))`，`P_ref` 用初始场；累积与 probe 两处同步改。 | Bug1, 21, 22 |
| 04 Limiter 封顶 + PTC 符号修正 | F1 | 将 `damp_scale` 改为 `min(1.0,max(1e-12,dt_eff/dt_ref))`；PTC 对角增强使用 `sign(d_acc)` 保号叠加。 | Bug2, Bug3 |
| 05 两相热力学一致性 | F1, F4, F9 | `Sw` 进入 `kr/Pc` 前强制安全截断；CO2 物性改用 `P_gas=P+Pc`；`AD_FluidEvaluator` 最终兜底不再返回全零梯度，给小导数地板。 | 雷区1,6,12 |
| 06 井模型隐式耦合修复 | F5, F12 | `Total` 注入分相支持总流度/注入组分逻辑；Rate 控制下 `fw` 用 AD 形式参与 Jacobian；注入焓计算使用常量井口压力（去掉 `P_cell` 梯度污染）；裂缝井 `WI` 缺失改为抛错；`Rate` 增加标准态体积转质量可选路径。 | 雷区2,4,5,9,24 |
| 07 边界项全耦合装配 | F5, F6, F1 | `Assemble_2D/3D` 从“标量对角”升级到“方程块 Jacobian”；边界通量显式乘物理系数（`k, mobility, lambda_eff`）；补全裂缝边界组装分支。 | 雷区7,15/20,19,11 |
| 08 裂缝域相渗/毛管模型分流 | F1, F4, F5 | 按 block 域判断：裂缝单元强制 `krw=Sw, krg=1-Sw, Pc=0`；基质保持 vG。 | 雷区8 |
| 09 跨网格导数与 EOS 缓存 | F1 | 保持双播种（wrt_i/wrt_j）并固化测试；新增单迭代物性缓存，避免同一 cell 多次 EOS 重算。 | 雷区16,17 |
| 10 NNC/拓扑与各向异性稳健化 | F11, F10 | 连接加载前 `reserve`；NNC 法向量归一后再做 `n^T K n`；对 pair/field 尺寸做强校验。 | 雷区10,14 |
| 11 稀疏模式固定化与 O(N) 更新 | F3, F1 | `ExportEigenSparseMatrix` 不再按阈值删元素；冻结后构建固定 pattern + value 索引，迭代仅更新数值。 | 雷区25,26 |
| 12 线性求解器内存与 CPR-AMG | F2, F8, F1 | `LinearSolverCache` 提升到跨 time-step 复用；减少重复深拷贝；新增 `AMGCL_CPR` 路径（压力块 AMG + 全耦合 ILU 校正）并保留 SparseLU 回退。 | 雷区27,28 |

### Test Plan（每个改动包完成后都执行）
1. 结构正确性：`day6_matrix_audit_2d_edfm`、`day6_matrix_audit_3d_edfm` 必须通过。  
2. 收敛行为：`day6_transient_2d_tp_injprod` 需显著优于当前基线（基线：`rollbacks=40/50`，`t_end=7.824e-02s`）。  
3. 覆盖回归：`day6_transient_2d_sp_injprod`、`day6_transient_3d_sp_injprod`、`day6_transient_3d_tp_injprod`。  
4. 稀疏稳定：连续牛顿迭代中矩阵 `nnz` 与 pattern hash 恒定。  
5. 线性稳健：AMGCL 失败回退 SparseLU 路径可达，且不出现 `dx_nan_inf`、`factorize` 崩溃。  

### Assumptions
- `WellControlMode::Rate` 现有输入默认仍按“质量流量”解释；标准态体积流量通过新增开关显式启用。  
- 若 `c_r` 场不存在，默认 `c_r=0` 并打印告警，不阻塞运行。  
- 井模型在过渡期可继续使用 `ADVar<3>`，但接口按模板化方向改造，保证 `N=2/3` 一致。  
- 你按上述 12 包逐步提交修改后，我按包做逐项 review（正确性、回归风险、性能副作用）。  
