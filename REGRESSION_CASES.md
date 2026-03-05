# REGRESSION_CASES

## 0. Purpose

This file defines the minimum regression suite for this EDFM-FVM codebase.
It must cover:

- 2D/3D EDFM geometry embedding pipeline
- NNC/FF transmissibility generation
- FIM topology aggregation and conservation checks
- FVM-AD operator baseline stability

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
## 3. Optional Extended Regression

- `--case=3d_distance_accuracy` for geometry distance accuracy
- `--case=3d_nnc_ff_static` for static NNC/FF analytic consistency
- `--case=property_sweep` for CO2/water property sweep stability
- `--case=grad_all` for 2D/3D gradient operator checks
- `--case=day1_arch_conn` for one-shot Day1 gate
- `--case=day1_arch_conn_repro` for one-shot Day1 reproducibility check
- `--case=day4_well_patch` for Day4 well-control framework gate
- `--case=day4_well_viz` for Day4 well-source VTK visualization gate

---

## 4. Recommended Frequency

- Before merge to `main`: run at least `R1-R5`
- After geometry/topology edits: `R2-R5` mandatory
- After AD/FVM operator edits: `R1` mandatory, `day3_bc_patch` + `day3_leakoff_switch` mandatory when boundary/leakoff logic changes, `grad_all` recommended
- After well-control edits: `day4_well_patch` + `day4_well_viz` + at least one transmissibility gate (`R4` or `R5`) mandatory

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

Notes:
- ...
```

---

## 6. Workflow Integration

- `PLAN_7D_CHECKLIST.md` is the daily execution tracker.
- A day can be marked PASS only after all required checks for that day are PASS.
- If new `--case` commands are introduced, or pass criteria are changed, update this file and `PROJECT_CONTEXT.md` in the same change set.

