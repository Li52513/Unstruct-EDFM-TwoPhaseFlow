# PLAN_7D_CHECKLIST

## 0. Rules (Must Follow)

- This file is the daily execution tracker: target -> strict acceptance -> PASS sign-off.
- A day can be checked only when all implementation items and all acceptance gates are PASS.
- End of day must include evidence: commit id, log path, result summary.
- If a change touches `--case`, regression commands, or pass criteria, update in the same change:
  - `REGRESSION_CASES.md`
  - `PROJECT_CONTEXT.md` (at least `Last Updated` and `Git Commit`)

---

## D0 Baseline (Current State)

- [x] `main.cpp` has dispatcher mode (`--case / --list / --help`)
- [x] `PROJECT_CONTEXT.md` is available for new-session handover
- [x] `REGRESSION_CASES.md` defines minimum suite `R1-R5`

Evidence:
- Commit: `1b74f9c (dirty working tree baseline)`
- Date: `2026-03-04`

---

## Day 1 - Architecture Freeze and Connection Consistency

Status: `[x]`

Implementation checklist:
- [x] Freeze DOF mapping conventions for single-phase and two-phase paths
- [x] Unify MM/FI/NNC/FF connection key and deterministic ordering
- [x] Freeze connection stats output format (CI-friendly)

Strict acceptance gates:
- [x] R4 PASS (2D transmissibility + FIM conservation)
- [x] R5 PASS (3D transmissibility + FIM conservation)
- [x] Connection stats are reproducible for same input

Evidence (fill):
- Commit: dirty working tree (uncommitted)
- Logs: --case=day1_arch_conn, --case=day1_arch_conn_repro
- Verdict: PASS

---

## Day 2 - Potential, Upwind, and Phase-Flux Kernel

Status: `[x]`

Implementation checklist:
- [x] Finalize potential-difference interface (including gravity)
- [x] Finalize upwind-by-potential logic with correct AD derivative direction
- [x] Integrate two-phase phase-flux assembly on shared kernel

Strict acceptance gates:
- [x] R1 PASS (FVM-AD operator tests)
- [x] `--case=grad_all` PASS
- [x] No non-physical net flow in hydrostatic/gravity equilibrium case

Evidence (fill):
- Commit: dirty working tree (manual run, uncommitted)
- Logs: --case=day2_fvm_ad, --case=grad_all
- Verdict: PASS

---

## Day 3 - Boundary Conditions and Leakoff Discretization

Status: `[x]`

Implementation checklist:
- [x] Integrate Dirichlet/Neumann into residual and Jacobian
- [x] Add parameterized leakoff source/sink with switch
- [x] Ensure boundary derivatives enter correct blocks

Strict acceptance gates:
- [x] `--case=day3_bc_patch` PASS
- [x] `--case=day3_leakoff_switch` PASS
- [x] Linear pressure-drop and adiabatic patch tests PASS
- [x] 2D/3D grid-level boundary assembly PASS with row-offset check

Evidence (fill):
- Commit: dirty working tree (uncommitted)
- Logs: `--case=day3_bc_patch`, `--case=day3_leakoff_switch`, `--case=day3_viz`
- Verdict: PASS (including 2D/3D VTK visualization export and field sanity)

---

## Day 4 - Well Model and Injection/Production Control (WAG Base)

Status: `[x]`

Implementation checklist:
- [x] Integrate well terms (BHP/Rate) consistently into grid equations
- [x] Verify coupled response with leakoff on/off
- [x] Build switchable schedule control skeleton for WAG extension

Strict acceptance gates:
- [x] Single-well pressure response is physically consistent
- [x] Leakoff on/off trend is correct and explainable
- [x] At least one of R4/R5 PASS (depending on impact scope)

Evidence (fill):
- Commit: d9ddb9e (dirty working tree)
- Logs: `--case=day4_well_patch`, `--case=day4_well_viz`, `--case=trans_2d`, `--case=trans_3d`
- Verdict: PASS

---

## Day 5 - Global Assembly and Jacobian Verification

Status: `[x]`

Implementation checklist:
- [x] Finish key block assembly for single-phase and two-phase equations
- [x] Ensure accumulation and flux derivatives are assembled consistently
- [x] Provide FD vs AD Jacobian comparison case/script

Strict acceptance gates:
- [x] FD vs AD max relative error `< 1e-6`
- [x] R1 PASS
- [x] At least one of R4/R5 PASS for impacted dimension

Evidence (fill):
- Commit: d9ddb9e (dirty working tree)
- Logs: `--case=day5_block_matrix_robust`, `--case=day5_global_jac_2d`, `--case=day5_global_jac_3d`, `--case=day2_fvm_ad`, `--case=trans_2d`, `--case=trans_3d`
- Verdict: PASS

---

## Day 6 - Transient Nonlinear Solver Robustness

Status: `[ ]`

Implementation checklist:
- [x] Add runnable Day6 dispatcher entries for full scenario matrix:
  - `day6_transient_2d_sp_injprod`
  - `day6_transient_2d_tp_injprod`
  - `day6_transient_2d_tp_multiwell`
  - `day6_transient_3d_sp_injprod`
  - `day6_transient_3d_tp_injprod`
  - `day6_transient_3d_tp_multiwell`
- [x] Add Newton/FIM or IMPES iteration logs per step: `step, iter, residual, dt, rollback`
- [x] Add limiter / adaptive-dt / rollback triggers and explicit trigger counters in logs
- [x] Enforce physical guards (`S in [0,1]`, pressure/temperature lower bounds, no NaN/Inf)
- [x] Integrate PostProcess VTK export in Day6 transient cases:
  - 2D: use `PostProcess_2D::ExportVTK(...)`
  - 3D: use `PostProcess_3D::ExportVTK(...)`
  - output fields include at least `P`, `T`, (`S_w` for two-phase), optional `DomainID`/well source markers

Strict acceptance gates:
- [ ] Gate A (Transient Stability, all Day6 scenarios)
      Commands:
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_sp_injprod`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_tp_injprod`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_tp_multiwell`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_sp_injprod`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_tp_injprod`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_tp_multiwell`
      Pass: each case completes at least `50` continuous steps; no `Error/Exception/Fatal`
- [ ] Gate B (Nonlinear Convergence Traceability)
      Pass: every step prints `iter/residual/dt`; non-converged step has rollback record; final step is converged
- [ ] Gate C (Physical Validity)
      Pass: no `NaN/Inf`; saturation bounded in `[0,1]`; guarded state variables do not violate lower bounds
- [ ] Gate D (VTK Export Integrity)
      Pass:
      - each Day6 case exports at least one `.vtk` file
      - target files exist under `Test/Transient/Day6/...`
      - files are non-empty and include expected scalar arrays (`P`,`T`,`S_w` when applicable)
      - log contains explicit export PASS line per case
- [ ] Gate E (ParaView Visualization Acceptance, manual)
      Pass:
      - `.vtk` loads without reader error in ParaView
      - scalar fields are visible and non-uniform where expected (`P`,`T`,`S_w`)
      - well influence region and injection/production trend are visually consistent
      - multi-well case shows distinguishable interference/plume fronts
- [ ] Gate F (Regression Guard)
      Command: `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day2_fvm_ad`
      Pass: R1 PASS
- [ ] Gate G (Impacted-Dimension Transmissibility Guard)
      2D path: `--case=trans_2d`; 3D path: `--case=trans_3d`
      Pass: at least impacted dimension gate PASS

Evidence (fill):
- Commit:
- Logs:
  - Day6 transient:
    - `--case=day6_transient_2d_sp_injprod`
    - `--case=day6_transient_2d_tp_injprod`
    - `--case=day6_transient_2d_tp_multiwell`
    - `--case=day6_transient_3d_sp_injprod`
    - `--case=day6_transient_3d_tp_injprod`
    - `--case=day6_transient_3d_tp_multiwell`
  - Guard rails: `--case=day2_fvm_ad` + `--case=trans_2d|trans_3d`
  - Key log lines: `step/iter/residual/dt/rollback/limiter_count/vtk_export`
  - VTK outputs (examples):
    - `Test/Transient/Day6/2D_sp_injprod/final.vtk`
    - `Test/Transient/Day6/2D_tp_multiwell/final.vtk`
    - `Test/Transient/Day6/3D_sp_injprod/final.vtk`
    - `Test/Transient/Day6/3D_tp_multiwell/final.vtk`
  - ParaView screenshots path (optional but recommended):
    - `Test/Transient/Day6/ParaViewShots/...`
- Verdict: PASS / FAIL
---

## Day 7 - Two-Phase Thermal Closed Loop and WAG Regression

Status: `[ ]`

Implementation checklist:
- [ ] Add runnable Day7 dispatcher entries for closed-loop + WAG full matrix:
  - `day7_closedloop_2d_sp_baseline`
  - `day7_closedloop_2d_tp_wag_singlewell`
  - `day7_closedloop_2d_tp_wag_multiwell`
  - `day7_closedloop_3d_sp_baseline`
  - `day7_closedloop_3d_tp_wag_singlewell`
  - `day7_closedloop_3d_tp_wag_multiwell`
- [ ] Implement reproducible WAG schedule driver:
  - fixed cycle definition (`W-inj -> soak/transition -> CO2-inj -> prod`)
  - configurable cycle count and stage duration
  - explicit schedule log per step (active phase, active wells, control mode)
- [ ] Integrate non-isothermal two-phase closed-loop exchange (matrix/fracture/NNC/FF + wells)
- [ ] Export paper-grade metrics time series:
  - `P, T, S_w, q_w, q_g, q_h, cumulative_heat_extraction, injected_mass, produced_mass`
  - per-well and field-total summaries
- [ ] Integrate PostProcess output for Day7:
  - VTK snapshots (`initial/mid/final`) for each scenario
  - optional CSV summaries for plotting and reproducibility

Strict acceptance gates:
- [ ] Gate A (Scenario Coverage)
      Commands:
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_sp_baseline`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_tp_wag_singlewell`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_tp_wag_multiwell`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_sp_baseline`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_tp_wag_singlewell`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_tp_wag_multiwell`
      Pass: all required scenarios run to completion with no `Error/Exception/Fatal`
- [ ] Gate B (WAG Schedule Correctness)
      Pass:
      - logged stage sequence matches configured WAG cycle definition
      - phase switching occurs at configured step/time boundaries
      - injection/production control mode transitions are explicit and traceable
- [ ] Gate C (Closed-Loop Physical Consistency)
      Pass:
      - no `NaN/Inf`
      - two-phase bounds valid (`S_w in [0,1]`, implied `S_g in [0,1]`)
      - pressure/temperature remain within configured physical bounds
      - net heat extraction direction is physically explainable
- [ ] Gate D (Mass/Energy Accounting)
      Pass:
      - cumulative injected vs produced mass/energy are reported
      - imbalance statistics are printed and within configured tolerance
      - no unexplained monotonicity breaks in key cumulative metrics
- [ ] Gate E (VTK Export Integrity)
      Pass:
      - each Day7 scenario exports `.vtk` snapshots (`initial/mid/final`)
      - files exist under `Test/Transient/Day7/...` and are non-empty
      - fields include `P`, `T`, and `S_w` for two-phase scenarios
- [ ] Gate F (ParaView Visualization Acceptance, manual)
      Pass:
      - all exported `.vtk` files load without reader errors
      - plume/thermal fronts and well interference are visible and coherent
      - WAG stage transitions are visually identifiable in scalar fields
- [ ] Gate G (Core Regression Guard)
      Commands:
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day2_fvm_ad`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=trans_2d`
      - `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=trans_3d`
      Pass: R1/R4/R5 PASS
- [ ] Gate H (Full Baseline Guard)
      Pass: R1-R5 all PASS and Day6 impacted gates remain PASS after Day7 changes

Evidence (fill):
- Commit:
- Logs:
  - Day7 scenario logs:
    - `--case=day7_closedloop_2d_sp_baseline`
    - `--case=day7_closedloop_2d_tp_wag_singlewell`
    - `--case=day7_closedloop_2d_tp_wag_multiwell`
    - `--case=day7_closedloop_3d_sp_baseline`
    - `--case=day7_closedloop_3d_tp_wag_singlewell`
    - `--case=day7_closedloop_3d_tp_wag_multiwell`
  - Guard rails:
    - `--case=day2_fvm_ad`
    - `--case=trans_2d`
    - `--case=trans_3d`
  - Key log lines:
    - `step/iter/residual/dt`
    - `wag_stage/schedule_switch`
    - `mass_balance/energy_balance/cumulative_heat_extraction`
    - `vtk_export`
  - VTK outputs (examples):
    - `Test/Transient/Day7/2D_tp_wag_multiwell/initial.vtk`
    - `Test/Transient/Day7/2D_tp_wag_multiwell/final.vtk`
    - `Test/Transient/Day7/3D_tp_wag_multiwell/initial.vtk`
    - `Test/Transient/Day7/3D_tp_wag_multiwell/final.vtk`
  - ParaView screenshots path (optional but recommended):
    - `Test/Transient/Day7/ParaViewShots/...`
- Verdict: PASS / FAIL

---

## Daily Sign-off Log (Required)

| Date | Day | Owner | Commit | Required Checks | Result | Notes |
|---|---|---|---|---|---|---|
| 2026-03-04 | Day 2 | Yongwei | dirty working tree | R1 + grad_all + hydrostatic gate | PASS | Day2 FVM-AD suite all PASS; 2D/3D grad patch tests PASS |
| 2026-03-04 | Day 3 | Yongwei | dirty working tree | day3_bc_patch + day3_leakoff_switch + day3_viz | PASS | BC patch + Leakoff OFF/ON + two-phase AD + 2D/3D grid-level assembly PASS; VTK export PASS; area-weighted boundary residual pattern observed |
| 2026-03-05 | Day 4 | Yongwei | d9ddb9e (dirty working tree) | day4_well_patch + day4_well_viz + trans_2d + trans_3d | PASS | 6 Day4 PASS lines present; 2D/3D VTK exported to fixed paths; trans_2d/trans_3d PASS; no Error/Exception/Fatal |
| 2026-03-07 | Day 5 | Yongwei | d9ddb9e (dirty working tree) | day5_block_matrix_robust + day5_global_jac_2d + day5_global_jac_3d + day2_fvm_ad + trans_2d + trans_3d | PASS | Day5 infra gate PASS; 2D/3D FD-vs-AD rel err < 1e-6; R1/R4/R5 PASS; no Error/Exception/Fatal |
| YYYY-MM-DD | Day N |  |  | R? / case? | PASS/FAIL |  |
