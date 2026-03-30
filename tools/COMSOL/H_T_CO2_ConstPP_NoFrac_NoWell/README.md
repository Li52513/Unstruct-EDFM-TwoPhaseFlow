# COMSOL 6.3 Reference Workflow

This directory contains the COMSOL Java automation for:

- case: `Test/Transient/FullCaseTest/H_T_CO2_ConstPP/h_t_co2_constpp_nofrac_nowell`
- physics: 2D single-phase CO2, constant properties, pressure-temperature coupling
- geometry: matrix-only rectangular porous medium, no fracture, no well
- gravity: off

## Inputs

Run the engineering case first so these files exist:

- `engineering/profile_station_definitions.csv`
- `engineering/monitor_point_definitions.csv`
- `engineering/profile_report_schedule.csv`
- `engineering/monitor_sample_schedule.csv`
- `reference/comsol_input/property_table.csv`

## Commands

From the workspace root:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Compile
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Run
```

One-shot:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\COMSOL\H_T_CO2_ConstPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode All
```

## Outputs

The workflow writes to `reference/comsol/`:

- `comsol_profile_matrix_horizontal_*.csv`
- `comsol_profile_matrix_vertical_midline_*.csv`
- `comsol_monitor_timeseries.csv`
- `comsol_model.mph`
- `comsol_reference_mesh_check.txt`
- `comsol_run_summary.md`
- `comsol_progress.log`
- `comsol_java_stdout.log`
- `comsol_java_stderr.log`

## Notes

- v1 uses a matrix-only coefficient-form PDE route for both pressure and temperature, which is the allowed fallback when the built-in Darcy/heat-transfer interfaces are not convenient for automation.
- The wrapper treats output files as the success criterion and fails if expected CSV or MPH artifacts are missing.
- The current mesh-check file is a contract artifact for workflow completeness; if paper-grade COMSOL mesh-independence evidence is required, extend this directory with a second finer mesh pass.
