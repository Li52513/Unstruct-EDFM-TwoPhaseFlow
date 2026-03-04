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

## 3. Optional Extended Regression

- `--case=3d_distance_accuracy` for geometry distance accuracy
- `--case=3d_nnc_ff_static` for static NNC/FF analytic consistency
- `--case=property_sweep` for CO2/water property sweep stability
- `--case=grad_all` for 2D/3D gradient operator checks
- `--case=day1_arch_conn` for one-shot Day1 gate
- `--case=day1_arch_conn_repro` for one-shot Day1 reproducibility check

---

## 4. Recommended Frequency

- Before merge to `main`: run at least `R1-R5`
- After geometry/topology edits: `R2-R5` mandatory
- After AD/FVM operator edits: `R1` mandatory, `grad_all` recommended

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
