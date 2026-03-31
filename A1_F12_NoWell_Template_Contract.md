# A1_F12 No-Well Template Contract

## 1. Scope
- This contract freezes the 2D no-well baseline templates `A1`, `B1`, and `C1`.
- The later 2D no-well cases `A2-A6`, `B2-B6`, and `C2-C6` must inherit this contract.
- Well cases are excluded and will be frozen separately in `M4-S5`.

## 2. Purpose
- Make `A1/B1/C1` the only canonical source for 2D no-well thin templates.
- Freeze stage semantics, canonical directory layout, engineering/reference exports, report placement, and Matlab script placement.
- Declare that loose case-root artifacts are transitional and are not part of the contract.

## 3. Stage Semantics
- `prepare_reference`
  - Prepares reference-input materials, sampling definitions, schedules, property tables, and reference contracts.
  - Does not invoke COMSOL directly.
- `solve_only`
  - Runs the engineering solve and writes engineering-side outputs.
- `validate_only`
  - Consumes existing engineering/reference data and writes validation outputs.
  - Must fall back to `missing_reference` cleanly when reference data is absent.
- `full_workflow`
  - Chains `prepare_reference -> solve_only -> validate_only`.

## 4. Canonical Directory Contract
- Every no-well case root must contain:
  - `studies/`
  - `figures/`
  - `engineering/`
  - `reference/`
  - `report/`
- `report/scripts/` must always exist.
- Directory responsibilities are frozen:
  - `studies/`: grid and time-step studies.
  - `figures/`: exported figures.
  - `engineering/`: raw engineering CSV, schedules, station definitions, and stage manifest.
  - `reference/`: reference contracts, COMSOL inputs, and reference results.
  - `report/`: status markdown, validation summary, and Matlab scripts.

## 5. Common Required Files
- `engineering/stage_manifest.txt`
- `reference/reference_contract.txt`
- `report/template_status.md`
- `report/scripts/`

The manifest is the execution truth source and must be able to expose states such as `prepared_reference_inputs`, `completed`, and `missing_reference`.

## 6. A1 Contract
- `A1` is the analytical single-phase pressure-diffusion baseline.
- Required canonical files:
  - `studies/grid_convergence.csv`
  - `studies/time_sensitivity.csv`
  - `engineering/stage_manifest.txt`
  - `reference/reference_contract.txt`
  - `report/template_status.md`
- Any additional analytical summary file at the case root is transitional and not part of the frozen contract.

## 7. B1 and C1 Engineering Export Contract
- Both `B1` and `C1` must provide:
  - `engineering/property_table.csv`
  - `engineering/profile_station_definitions.csv`
  - `engineering/monitor_point_definitions.csv`
  - `engineering/profile_report_schedule.csv`
  - `engineering/monitor_sample_schedule.csv`
  - `engineering/eng_monitor_timeseries.csv`
- The profile families are frozen to:
  - `matrix_horizontal`
  - `matrix_vertical_midline`
- The snapshot tags are frozen to:
  - `t010pct`
  - `t050pct`
  - `t100pct`
- Therefore the following engineering files are mandatory for both `B1` and `C1`:
  - `engineering/eng_profile_matrix_horizontal_t010pct.csv`
  - `engineering/eng_profile_matrix_horizontal_t050pct.csv`
  - `engineering/eng_profile_matrix_horizontal_t100pct.csv`
  - `engineering/eng_profile_matrix_vertical_midline_t010pct.csv`
  - `engineering/eng_profile_matrix_vertical_midline_t050pct.csv`
  - `engineering/eng_profile_matrix_vertical_midline_t100pct.csv`

## 8. B1 Reference and Report Contract
- `B1` must additionally provide:
  - `reference/comsol_input/property_table.csv`
  - `report/validation_summary.md`
  - `report/validation_summary.csv`
- Matlab scripts may continue evolving, but they must live under `report/scripts/`.

## 9. C1 Reference and Report Contract
- `C1` must additionally provide:
  - `report/validation_summary.md`
  - `report/scripts/plot_validation_results.m`
- When reference data is missing, `validate_only/full_workflow` must end in `missing_reference` while still materializing the canonical five-directory layout.

## 10. Frozen CSV Schemas
- `profile_station_definitions.csv`
```csv
station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m
```

- `B1` `eng_profile_*.csv`
```csv
station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,p_num_pa,t_num_k
```

- `C1` `eng_profile_*.csv`
```csv
station_id,label,family,location,target_axis_m,target_x_m,target_y_m,block_id,actual_x_m,actual_y_m,target_time_s,actual_time_s,p_num_pa,t_num_k,sw_num
```

- `B1` `eng_monitor_timeseries.csv` header pattern
```text
sample_id,target_time_s,actual_time_s,p_num_<label>...,t_num_<label>...
```

- `C1` `eng_monitor_timeseries.csv` header pattern
```text
sample_id,target_time_s,actual_time_s,p_num_<label>...,t_num_<label>...,sw_num_<label>...
```

- `A1` `studies/grid_convergence.csv`
```csv
study,label,case_dir,nx,ny,dt_init,h_char,t_end,steps,l1_norm,l2_norm,linf_norm
```

## 11. Transitional Outputs
- Some legacy outputs may still appear directly under the case root.
- Those loose artifacts are not part of the frozen contract.
- Future consumers must rely only on the five canonical directories and the files listed in this contract.

## 12. Change Rule
- Future no-well cases may only add case-specific differences on top of this contract.
- If a later step needs new profile families, new snapshot tags, or additional CSV columns, this contract and its check script must be updated together.
