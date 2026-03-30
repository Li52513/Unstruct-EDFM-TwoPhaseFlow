# COMSOL 6.3 Reference Workflow

This directory contains the COMSOL Java automation for:

- case: `Test/Transient/FullCaseTest/H_T_CO2_ConstPP/h_t_co2_constpp_singlefrac_nowell`
- physics: 2D single-phase CO2, constant properties, pressure-temperature coupling
- fracture route: explicit lower-dimensional fracture on an internal line segment, coupled to the surrounding matrix by 2D-1D exchange terms
- gravity: off

## Inputs

Run the engineering prepare phase first so these files exist:

- `engineering/profile_station_definitions.csv`
- `engineering/monitor_point_definitions.csv`
- `engineering/profile_report_schedule.csv`
- `engineering/monitor_sample_schedule.csv`
- `reference/comsol_input/property_table.csv`

## Commands

From the workspace root:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_SingleFrac_NoWell\run_comsol_reference.ps1 -Mode Compile
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_SingleFrac_NoWell\run_comsol_reference.ps1 -Mode Run
```

One-shot:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_SingleFrac_NoWell\run_comsol_reference.ps1 -Mode All
```

Skip the fine-mesh consistency check:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_SingleFrac_NoWell\run_comsol_reference.ps1 -Mode Run -SkipFineCheck
```

## Outputs

The workflow writes to `reference/comsol/`:

- `comsol_profile_matrix_horizontal_*.csv`
- `comsol_profile_fracture_tangent_*.csv`
- `comsol_profile_cross_normal_*.csv`
- `comsol_monitor_timeseries.csv`
- `comsol_model.mph`
- `comsol_model_refined.mph`
- `comsol_reference_mesh_check.txt`
- `comsol_run_summary.md`
- `comsol_progress.log`
- `comsol_java_stdout.log`
- `comsol_java_stderr.log`

## Notes

- The wrapper uses COMSOL's bundled `javac.exe` and builds a support JAR for auxiliary classes automatically.
- The default runtime is `DirectJava`, aligned with the existing `VaryPP_NoFrac_NoWell` COMSOL 6.3 workflow that has already produced reference outputs in this project.
- `ComsolBatch` is retained as an experimental fallback path.
- Headless COMSOL still needs a writable `%USERPROFILE%\.comsol\v63\...` tree for the default `DirectJava` route, but the wrapper also redirects recovery and temp files to a dedicated ASCII-only runtime directory.
- The wrapper treats output files as the primary success criterion and will fail if the expected CSV or MPH files are missing.
- If a previous run summary reports `thin_band_fallback`, the wrapper archives that legacy surrogate output tree into `reference/comsol_surrogate_legacy/<timestamp>/` before running the explicit-fracture route.
- Final paper-grade validation should be based on a run without `-SkipFineCheck`; the skip flag is only for early runtime debugging.
