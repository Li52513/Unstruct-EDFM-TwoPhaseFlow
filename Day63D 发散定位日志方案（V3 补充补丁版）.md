# Day6/3D 发散定位日志方案（V3 补充补丁版）

## Summary
- 保持你确认的 V3 计划 **完全不变**。
- 仅补充 2 项：`[SOLVER]` 常驻健康摘要、`[INCIDENT-CONST]` 两相本构快照。
- 其余模块、触发条件、分层日志机制全部沿用 V3。

## Key Additions
1. 在常驻摘要中补回 `[SOLVER]`（每迭代输出）
- `SparseLU` 输出：`solver=SparseLU compute_ok solve_ok nnzA info`
- `BiCGSTAB` 输出：`solver=BiCGSTAB compute_ok solve_ok iters error info`
- 若某字段不可读，输出 `NA`，但键名不变。
- 在线性失败分支（`linear_solve_fail`）前，强制再打印一次 `[SOLVER]` 最终状态，保证故障可追溯。

2. 在事故包中增加 `[INCIDENT-CONST]`（仅 `N==3`）
- 对当前 hot block 输出：`Sw krw krg Pc dkrw_dSw dkrg_dSw dPc_dSw`
- 计算口径与主方程保持一致：使用与通量装配同一套 VG/相渗参数（同一 `vg/rp`）。
- 若出现 `NaN/Inf` 或异常值，附加 `flag=invalid_constitutive`。
- 单相 (`N==2`) 不输出该段，避免噪声。

## Test Plan Additions
1. `SOLVER` 可见性
- 跑默认 `SparseLU` case，确认每迭代都有 `[SOLVER]`。
- 切换 `BiCGSTAB` 跑同一 case，确认 `iters/error` 正常打印。
- 人为构造线性失败，确认失败前有 `[SOLVER]` 终态行。

2. 本构快照可用性
- 在 `N==3` case 触发 incident，确认出现 `[INCIDENT-CONST]` 且 7 个量齐全。
- 在 `N==2` case 触发 incident，确认不输出 `[INCIDENT-CONST]`。

## Assumptions
- 除上述两项外，V3 原计划逐条不改。
- `kr/Pc` 及其对 `Sw` 导数以当前 AD 路径计算，不引入新物理模型。
