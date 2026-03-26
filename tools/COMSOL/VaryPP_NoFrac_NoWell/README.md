# COMSOL 6.3 Reference Workflow

This directory contains the COMSOL Java automation for the VaryPP validation case:

- case: `Test/Transient/FullCaseTest/H_CO2_VaryPP/h_co2_varypp_nofrac_nowell`
- physics: 2D `Coefficient Form PDE`
- equation: `S(p) * dp/dt - div(lambda(p) * grad(p)) = 0`
- temperature: `360 K`
- pressure window: `12 MPa -> 8 MPa`
- gravity: off

## Prerequisites

Run the engineering case once before COMSOL so these files exist in the case directory:

- `profile_station_definitions.csv`
- `monitor_point_definitions.csv`
- `profile_report_schedule.csv`
- `monitor_sample_schedule.csv`
- `reference/comsol_input/rho_fun.txt`
- `reference/comsol_input/mu_fun.txt`
- `reference/comsol_input/cf_fun.txt`
- `reference/comsol_input/lambda_fun.txt`
- `reference/comsol_input/S_fun.txt`

## Recommended Workflow

Run from the workspace root:

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1
```

The wrapper script uses the COMSOL 6.3 bundled Java toolchain by default:

- compile: `javac.exe`
- run: `java.exe`

This is the preferred route when the workspace path contains non-ASCII characters.

## Compile Only

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Compile
```

If you explicitly want the COMSOL native compiler wrapper:

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Compile -Compiler ComsolCompile
```

## Run Only

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Run
```

Skip the fine-mesh reference check:

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Run -SkipFineCheck
```

If you explicitly want the COMSOL batch launcher:

```powershell
.\Tools\COMSOL\VaryPP_NoFrac_NoWell\run_comsol_reference.ps1 -Mode Run -Runtime ComsolBatch
```

## Outputs

The Java script writes to:

- `reference/comsol/comsol_profile_t010pct.csv`
- `reference/comsol/comsol_profile_t050pct.csv`
- `reference/comsol/comsol_profile_t100pct.csv`
- `reference/comsol/comsol_monitor_timeseries.csv`
- `reference/comsol/comsol_model.mph`
- `reference/comsol/comsol_model_refined.mph`
- `reference/comsol/comsol_reference_mesh_check.txt`
- `reference/comsol/comsol_batch.log`
- `reference/comsol/comsol_java_stdout.log`
- `reference/comsol/comsol_java_stderr.log`
- `reference/comsol/comsol_batch_stdout.log`
- `reference/comsol/comsol_batch_stderr.log`

## Notes

- The direct Java runtime expects the current Windows user profile to be writable because COMSOL creates `.comsol\v63\...` runtime folders there.
- On a normal local Windows session this is usually automatic. If needed, create these folders once under your user profile and rerun:
  - `.comsol\v63\comsol.recoveries`
  - `.comsol\v63\configuration\comsol`
  - `.comsol\v63\workspace\comsol`
