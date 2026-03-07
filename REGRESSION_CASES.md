# REGRESSION_CASES

## 0. Purpose

This file defines the minimum regression suite for this EDFM-FVM codebase.
It must cover:

- 2D/3D EDFM geometry embedding pipeline
- NNC/FF transmissibility generation
- FIM topology aggregation and conservation checks
- FVM-AD operator baseline stability
- Day5 global assembly and FD-vs-AD Jacobian verification gates

---

## 1. Prerequisites

- Build completed: `x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe`
- Run commands from repository root
- Check case names first:

```powershell
.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --list
```

---

## 2. Required Minimum Regression Matrix

| ID | Scope | Command | Pass Criteria |
|---|---|---|---|
| R1 | FVM-AD operators | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day2_fvm_ad` | Exit code `0`, no unhandled exception |
| R2 | 2D geometry chain (DFN) | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=2d_geom_benchmark_dfn` | Contains `Benchmark Completed`; files updated under `Test/MeshTest/GeomIndexTest/2D_EDFM` |
| R3 | 3D geometry chain | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=3d_mesh_benchmark` | Contains `Benchmark Completed`; files updated under `Test/MeshTest/GeomIndexTest/3D_EDFM` |
| R4 | 2D transmissibility + FIM conservation | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=trans_2d` | Contains `[PASS] 2D Flow & Heat Transmissibility Conservation Validated.` |
| R5 | 3D transmissibility + FIM conservation | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=trans_3d` | Contains `[PASS] Flow & Heat Transmissibility Conservation Validated.` |

---

## 2.1 Day1 Explicit Dispatcher Cases (Recommended)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day1 Gate | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day1_arch_conn` | Run Day1 main gate in one command (R4 + R5) | Both 2D and 3D conservation PASS lines appear |
| Day1 Repro Gate | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day1_arch_conn_repro` | Run reproducibility sequence (`trans_2d x2 + trans_3d x2`) | Aggregation stats are consistent across repeated runs |

---

## 2.2 Day3 Explicit Dispatcher Cases (Mandatory for Boundary/Leakoff Edits)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day3 BC Patch | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day3_bc_patch` | Verify Dirichlet/Neumann boundary operators and AD derivatives | Contains `[PASS] 线性压降 patch PASS` and `[PASS] 绝热 patch PASS` |
| Day3 Leakoff Switch | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day3_leakoff_switch` | Verify single/two-phase Leakoff operators and grid-level 2D/3D assembly path | Contains `[PASS] Leakoff OFF=0 PASS`, `[PASS] Leakoff ON>0 PASS`, `[PASS] 两相 w/g sink 与导数 PASS`, `[PASS] 2D 网格装配与多块(P/T)分流行写入 PASS`, and `[PASS] 3D 网格装配与分流行写入 PASS` |

---
## 2.3 Day4 Explicit Dispatcher Case (Mandatory for Well-Control Edits)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day4 Well Patch | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day4_well_patch` | Verify 2D/3D well assembly paths, matrix+fracture completion, mass+energy coupling, and CSV-driven WAG schedule (including alternating water/gas injection/production) | Contains `[PASS] WI 2D geometric-eq`, `[PASS] WI 3D geometric-eq`, `[PASS] BHP/Rate mass coupling`, `[PASS] Well energy coupling`, `[PASS] WAG config switching`, and `[PASS] Matrix+Fracture completion` |
| Day4 Well Viz | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day4_well_viz` | Export fixed-path Day4 2D/3D well-source VTK for ParaView acceptance | Contains `[PASS] Day4 2D VTK exported`, `[PASS] Day4 3D VTK exported`, and `[PASS] Day4 well visualization VTK exported`; outputs `Test/BoundaryTest/day4_well_viz_2d.vtk` and `Test/BoundaryTest/day4_well_viz_3d.vtk` |

---
## 2.4 Day5 Explicit Dispatcher Cases (Mandatory for Global Assembly/Jacobian Edits)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day5 Infra Gate | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day5_block_matrix_robust` | Validate block sparse matrix container robustness (bounds/pattern/zeroing/export stability) | Contains `[PASS] FIM_BlockSparseMatrix Infrastructure is strictly PROD-ready.` |
| Day5 2D Jacobian Gate | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day5_global_jac_2d` | Validate 2D global assembly with FD vs AD Jacobian checks (single-phase + two-phase) | Contains `[PASS] Day5 2D Global Assembly Passed!` and per-scenario `FD vs AD max relative error < 1e-6` |
| Day5 3D Jacobian Gate | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day5_global_jac_3d` | Validate 3D global assembly with FD vs AD Jacobian checks (single-phase + two-phase) | Contains `[PASS] Day5 3D Global Assembly Passed!` and per-scenario `FD vs AD max relative error < 1e-6` |

---
## 2.5 Day6 Explicit Dispatcher Cases (Mandatory for Transient/VTK Visualization Edits)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day6 2D Single-Phase Inj/Prod | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_sp_injprod` | 2D single-phase transient stability + well response + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK exported |
| Day6 2D Two-Phase Inj/Prod | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_tp_injprod` | 2D two-phase transient stability + S/T coupling + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK with `P/T/S_w` |
| Day6 2D Two-Phase Multiwell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_2d_tp_multiwell` | 2D multi-well interference/plume transient + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK exported |
| Day6 3D Single-Phase Inj/Prod | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_sp_injprod` | 3D single-phase transient stability + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK exported |
| Day6 3D Two-Phase Inj/Prod | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_tp_injprod` | 3D two-phase transient stability + S/T coupling + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK with `P/T/S_w` |
| Day6 3D Two-Phase Multiwell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day6_transient_3d_tp_multiwell` | 3D multi-well interference/plume transient + VTK export | `>=50` steps stable, no `Error/Exception/Fatal`, VTK exported |

Note (ParaView acceptance):
- Load exported `.vtk` from `Test/Transient/Day6/...` in ParaView.
- Verify scalar fields render correctly (`P`,`T`,`S_w` when applicable) and flow/thermal trends are physically consistent.

---
## 2.6 Day7 Explicit Dispatcher Cases (Mandatory for Closed-Loop/WAG Regression Edits)

| Case | Command | Purpose | Pass Criteria |
|---|---|---|---|
| Day7 2D Single-Phase Baseline | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_sp_baseline` | 2D closed-loop baseline for thermal-flow trend reference | No `Error/Exception/Fatal`; VTK snapshots exported |
| Day7 2D Two-Phase WAG Singlewell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_tp_wag_singlewell` | 2D two-phase thermal WAG with single injector/producer | Stage switching logged; metrics and VTK exported |
| Day7 2D Two-Phase WAG Multiwell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_2d_tp_wag_multiwell` | 2D multiwell interference and plume response under WAG | Stage switching logged; metrics and VTK exported |
| Day7 3D Single-Phase Baseline | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_sp_baseline` | 3D closed-loop baseline for thermal-flow trend reference | No `Error/Exception/Fatal`; VTK snapshots exported |
| Day7 3D Two-Phase WAG Singlewell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_tp_wag_singlewell` | 3D two-phase thermal WAG with single injector/producer | Stage switching logged; metrics and VTK exported |
| Day7 3D Two-Phase WAG Multiwell | `.\x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe --case=day7_closedloop_3d_tp_wag_multiwell` | 3D multiwell interference and plume response under WAG | Stage switching logged; metrics and VTK exported |

Mandatory Day7 evidence:
- Exported VTK snapshots exist under `Test/Transient/Day7/...` (`initial/mid/final`).
- ParaView can load exported files and show coherent `P/T/S_w` evolution.
- Key metrics table includes cumulative heat extraction and mass/energy balance stats.

---
## 3. Optional Extended Regression

- `--case=3d_distance_accuracy` for geometry distance accuracy
- `--case=3d_nnc_ff_static` for static NNC/FF analytic consistency
- `--case=property_sweep` for CO2/water property sweep stability
- `--case=grad_all` for 2D/3D gradient operator checks
- `--case=day1_arch_conn` for one-shot Day1 gate
- `--case=day1_arch_conn_repro` for one-shot Day1 reproducibility check
- `--case=day4_well_patch` for Day4 well-control framework gate
- `--case=day4_well_viz` for Day4 well-source VTK visualization gate
- `--case=day5_block_matrix_robust` for Day5 block sparse matrix infra gate
- `--case=day5_global_jac_2d` for Day5 2D global assembly/Jacobian gate
- `--case=day5_global_jac_3d` for Day5 3D global assembly/Jacobian gate
- `--case=day6_transient_2d_sp_injprod` for Day6 2D single-phase transient + VTK gate
- `--case=day6_transient_2d_tp_multiwell` for Day6 2D multiwell transient + VTK gate
- `--case=day6_transient_3d_sp_injprod` for Day6 3D single-phase transient + VTK gate
- `--case=day6_transient_3d_tp_multiwell` for Day6 3D multiwell transient + VTK gate
- `--case=day7_closedloop_2d_tp_wag_singlewell` for Day7 2D two-phase WAG gate
- `--case=day7_closedloop_2d_tp_wag_multiwell` for Day7 2D two-phase multiwell WAG gate
- `--case=day7_closedloop_3d_tp_wag_singlewell` for Day7 3D two-phase WAG gate
- `--case=day7_closedloop_3d_tp_wag_multiwell` for Day7 3D two-phase multiwell WAG gate

---

## 4. Recommended Frequency

- Before merge to `main`: run at least `R1-R5`
- After geometry/topology edits: `R2-R5` mandatory
- After AD/FVM operator edits: `R1` mandatory, `day3_bc_patch` + `day3_leakoff_switch` mandatory when boundary/leakoff logic changes, `grad_all` recommended
- After well-control edits: `day4_well_patch` + `day4_well_viz` + at least one transmissibility gate (`R4` or `R5`) mandatory
- After Day5 global assembly edits: `day5_block_matrix_robust` + (`day5_global_jac_2d` or `day5_global_jac_3d`) + `R1` + impacted dimension transmissibility gate (`R4` for 2D / `R5` for 3D) mandatory
- After Day6 transient/visualization edits: impacted Day6 transient gates + `R1` + impacted transmissibility gate (`R4`/`R5`) mandatory; ParaView manual check required
- After Day7 closed-loop/WAG edits: required Day7 scenario gates + `R1-R5` all PASS mandatory; VTK + ParaView + metrics evidence mandatory

---

## 5. Result Log Template

```text
Date: YYYY-MM-DD
Commit: <short_sha or dirty>
Runner: <name>

R1 day2_fvm_ad: PASS/FAIL
R2 2d_geom_benchmark_dfn: PASS/FAIL
R3 3d_mesh_benchmark: PASS/FAIL
R4 trans_2d: PASS/FAIL
R5 trans_3d: PASS/FAIL
D5a day5_block_matrix_robust: PASS/FAIL
D5b day5_global_jac_2d: PASS/FAIL
D5c day5_global_jac_3d: PASS/FAIL
D6a day6_transient_2d_sp_injprod: PASS/FAIL
D6b day6_transient_2d_tp_multiwell: PASS/FAIL
D6c day6_transient_3d_sp_injprod: PASS/FAIL
D6d day6_transient_3d_tp_multiwell: PASS/FAIL
D7a day7_closedloop_2d_tp_wag_singlewell: PASS/FAIL
D7b day7_closedloop_2d_tp_wag_multiwell: PASS/FAIL
D7c day7_closedloop_3d_tp_wag_singlewell: PASS/FAIL
D7d day7_closedloop_3d_tp_wag_multiwell: PASS/FAIL

Notes:
- ...
```

---

## 6. Workflow Integration

- `PLAN_7D_CHECKLIST.md` is the daily execution tracker.
- A day can be marked PASS only after all required checks for that day are PASS.
- If new `--case` commands are introduced, or pass criteria are changed, update this file and `PROJECT_CONTEXT.md` in the same change set.

