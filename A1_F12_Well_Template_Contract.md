# A1_F12 Well Template Contract

## 1. Scope
- This contract freezes the 2D injector-producer sample templates `A7`, `B7`, and `C7`.
- The later 2D well cases `A8-A12`, `B8-B12`, and `C8-C12` must inherit this contract unless an explicit later milestone updates both this document and its check script.
- Fracture completions and COMSOL automation payloads are excluded; those are frozen separately under `C1`.

## 2. Purpose
- Make `A7/B7/C7` the canonical source for thin well-template structure in 2D.
- Freeze default well placement, default control policy, stage semantics, engineering/reference/report placement, and well-output CSV schemas.
- Declare which outputs are common to all well templates and which are intentionally case-family-specific.

## 3. Stage Semantics
- `prepare_reference`
  - Materializes weak-coupling reference-input materials without invoking COMSOL directly.
  - Must write the canonical directory layout, `reference_contract.txt`, `template_status.md`, `property_table.csv`, `well_schedule.csv`, and `reference_spec.md`.
- `solve_only`
  - Runs the engineering solve and writes engineering-side outputs.
  - `A7` writes the pressure-only engineering payload.
  - `B7/C7` additionally write `engineering/well_timeseries.csv`.
- `validate_only`
  - Re-runs the case under the case root, consumes existing engineering/reference payloads, and writes validation outputs.
  - Must fall back to `missing_reference` cleanly when COMSOL data is absent.
- `full_workflow`
  - Chains `prepare_reference -> solve_only -> validate_only`.

## 4. Canonical Directory Contract
- Every well-template case root must contain:
  - `studies/`
  - `figures/`
  - `engineering/`
  - `reference/`
  - `report/`
- `report/scripts/` must always exist.
- Directory responsibilities are frozen:
  - `studies/`: later well-side studies, not yet populated by the sample templates.
  - `figures/`: later well-side figure exports.
  - `engineering/`: stage manifest, engineering CSV payloads, VTK outputs, and run summaries.
  - `reference/`: weak-coupling contracts plus future COMSOL payloads.
  - `report/`: status markdown, validation summary, and plotting scripts.

## 5. Frozen Default Well Template
- `well_mode = injprod`
- `control_policy = injector_rate_producer_bhp`
- `injector_name = INJ`
- `producer_name = PROD`
- `injector_location = { label=injector, x_fraction=0.25, y_fraction=0.5, z_fraction=0.5, domain=matrix, axis=None }`
- `producer_location = { label=producer, x_fraction=0.75, y_fraction=0.5, z_fraction=0.5, domain=matrix, axis=None }`
- `injector_control_mode = rate`
- `producer_control_mode = bhp`
- `injector_component_mode = gas`
- `producer_component_mode = total`
- `injector_target_value = -1.0`
- `producer_target_value = 9.5e6 Pa`
- `avoid_direct_fracture_completion = true`
- The target positions implied by the frozen `48x6`, `400 m x 40 m` sample mesh are:
  - injector target `(100 m, 20 m)`
  - producer target `(300 m, 20 m)`

## 6. Common Required Files
- `engineering/stage_manifest.txt`
- `engineering/property_table.csv`
- `engineering/well_schedule.csv`
- `engineering/reference_spec.md`
- `engineering/run_summary.txt`
- `reference/reference_contract.txt`
- `report/template_status.md`
- `report/scripts/plot_validation_results.m`

The manifest remains the execution truth source and must be able to expose states such as `prepared_reference_inputs`, `completed`, and `missing_reference`.

## 7. Case-Family Differences
- `A7`
  - Physics family: `N=1` pressure-only CO2 injector-producer.
  - Injection fluid text is frozen to `co2`.
  - `thermal_policy = none`
  - `injector_temperature = -1.0`
  - Mandatory validation variables are `pressure`, `well_bhp`, and `well_rate`.
  - `engineering/well_timeseries.csv` and `report/validation_summary.md` are not yet mandatory for the pressure-only sample template.
- `B7`
  - Physics family: `N=2` thermal CO2 injector-producer.
  - Injection fluid text is frozen to `co2`.
  - `thermal_policy = cold_injection`
  - `injector_temperature = 320 K`
  - Mandatory validation variables are `pressure`, `temperature`, `well_bhp`, `well_rate`, and `production_temperature`.
  - Must provide `engineering/well_timeseries.csv` and `report/validation_summary.md`.
- `C7`
  - Physics family: `N=3` two-phase thermal CO2/H2O injector-producer.
  - Injection fluid text is frozen to `co2 into brine`.
  - `thermal_policy = cold_injection`
  - `injector_temperature = 320 K`
  - Mandatory validation variables are `pressure`, `temperature`, `co2_saturation`, `well_bhp`, `well_rate`, and `production_temperature`.
  - Must provide `engineering/well_timeseries.csv` and `report/validation_summary.md`.

## 8. Frozen Required Artifacts By Case
- `A7` must provide:
  - `engineering/stage_manifest.txt`
  - `engineering/property_table.csv`
  - `engineering/well_schedule.csv`
  - `engineering/reference_spec.md`
  - `engineering/run_summary.txt`
  - `reference/reference_contract.txt`
  - `report/template_status.md`
  - `report/scripts/plot_validation_results.m`
- `B7` must additionally provide:
  - `engineering/well_timeseries.csv`
  - `report/validation_summary.md`
- `C7` must additionally provide:
  - `engineering/well_timeseries.csv`
  - `report/validation_summary.md`

## 9. Frozen CSV Schemas
- `A7` `engineering/well_schedule.csv`
```csv
role,well_name,domain,control_mode,target_value,component_mode,completion_id,target_x_m,target_y_m,actual_x_m,actual_y_m,injection_is_co2
```

- `B7/C7` `engineering/well_schedule.csv`
```csv
role,well_name,domain,control_mode,target_value,component_mode,completion_id,target_x_m,target_y_m,actual_x_m,actual_y_m,injection_is_co2,injection_temperature_k
```

- `B7/C7` `engineering/well_timeseries.csv`
```csv
sample_id,time_s,well_name,role,control_mode,component_mode,target_value,actual_bhp_pa,actual_rate_kg_per_s,completion_pressure_pa,completion_temperature_k,reported_temperature_k,target_x_m,target_y_m,actual_x_m,actual_y_m
```

## 10. Transitional Outputs
- Additional loose artifacts may still appear under the case root or under `engineering/`.
- Those extra files are not part of the frozen contract unless they are listed above.
- Future automation and consumers must rely only on the canonical directories and files defined here.

## 11. Change Rule
- Future well cases may only add case-specific differences on top of this contract.
- If a later milestone needs new well CSV columns, different default well placement, or additional mandatory report files, this contract and `check_m4_s5_well_contract.ps1` must be updated together.
