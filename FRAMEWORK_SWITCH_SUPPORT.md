# Framework Switch Support Matrix

This document records the intended public switch surface for the transient EDFM framework.

## Public External Switches

- Geometry/topology:
  - `2D / 3D`
  - `NoFrac / SingleFrac / CrossFrac / DFN`
- Physics:
  - `N=1 / N=2 / N=3`
  - `Water / CO2`
  - `Const / EOS`
  - specialized `single_phase_no_convection`
- Controls:
  - `NoWell / BHP / Rate / MultiCompletion`
  - `FIM / IMPES`
  - `AMGCL / SparseLU / AMGCL_CPR / BiCGSTAB`

## Current Status

### Supported

- `FIM` route is the active transient route.
- `N=1` pressure-only cases without wells.
- `N=2` single-phase cases with unified external fluid config.
- `N=3` two-phase water/CO2 cases with unified external `EOS / Constant` switch.
- `N=2 / N=3` well model with explicit BHP DOF.
- `AMGCL / SparseLU / AMGCL_CPR / BiCGSTAB` linear solver selection through solver params.

### Reserved

- `SolverRoute::IMPES` remains a public switch, but is reserved and intentionally fails fast until implemented.

### Explicitly Not Supported

- `N=1` with wells.
- fracture well without explicit `wi_override`.
- `N=3` saturation boundary conditions as a fully independent Neumann/Robin control variable.

## Unified Fluid Entry

Preferred public entry is `FIM_Engine::UnifiedFluidModelConfig` through:

- `TransientOptionalModules::SetFluidModelConfig(...)`
- `FIM_CaseKit::PropertyPreset2D/3D::fluid_model`

This config now exposes:

- `N=1` pressure-only `Const / EOS`
- `N=2` single-phase `Water / CO2` and `Const / EOS`
- `N=3` two-phase `Water+CO2 EOS / Water+CO2 Constant`

Legacy fields remain for backward compatibility:

- `single_phase_fluid`
- `fluid_property_eval`
- `pressure_only_property_mode`
- `pressure_only_baseline_*`

The runtime resolves unified config first, then derives the legacy per-route fields used by existing internals.

## Dynamic Baseline Coverage

Current baseline runtime coverage is intentionally limited to:

- `test_h_co2_constpp_nofrac_nowell`
- `test_h_co2_constpp_singlefrac_nowell`
- `test_h_co2_varypp_nofrac_nowell`
- `test_h_co2_varypp_singlefrac_nowell`
- `test_h_t_co2_constpp_nofrac_nowell`

`Day6` and `ladder` cases are not part of this baseline contract.
