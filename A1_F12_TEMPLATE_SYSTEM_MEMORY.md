# A1-F12 Template System Memory

## Frozen Decisions
- System title is `A1-F12`; 2D cases are `A1-C12` and 3D cases are `D1-F12`.
- Keep descriptive `.h/.cpp` template filenames; `A1-F12` is exposed through a central catalog/dispatcher mapping.
- Every catalog case exposes four staged entry points: `RunSolveOnly`, `RunPrepareReference`, `RunValidateOnly`, `RunFullWorkflow`.
- All case artifacts follow the five-directory contract: `studies/`, `figures/`, `engineering/`, `reference/`, `report/`.
- COMSOL remains weakly coupled: templates emit and consume reference artifacts, while external Java/PowerShell automation executes COMSOL.
- Well defaults are frozen as injector-rate + producer-BHP; two-phase well cases inject CO2; thermal well cases use cold injection.
- `CaseMetadata` is frozen to carry `case_code`, `dispatcher_key`, `case_slug`, `description`, `reference_mode`, `implementation_status`, `output_root`, `well_control_policy`, `injection_fluid`, `thermal_injection_policy`, plus the six enum axes `dimension`, `equation_mode`, `property_mode`, `fracture_mode`, `well_mode`.
- `CaseArtifactPaths` is frozen to carry the five-directory contract plus the derived file paths `engineering_stage_manifest_path`, `reference_contract_path`, and `report_status_markdown_path`.
- `ValidationSpec2D/3D` is frozen as the shared contract for validation variables, profile snapshots, characteristic lines, observation points, grid-study tags, time-step study values, and required well-series outputs; 3D additionally carries characteristic slice definitions.
- `WellTemplateSpec` is frozen as the shared contract for injector/producer naming, control policy, injection fluid, thermal policy, default quarter/three-quarter placement, control modes, target values, injection temperature, and the default rule to avoid direct fracture completion.
- Family default acceptance-policy placeholders are frozen by equation family:
  - `A/D`: pressure profile + pressure monitor; when wells exist, require `well_bhp` and `well_rate`.
  - `B/E`: pressure/temperature profile + pressure/temperature monitor; when wells exist, require `well_bhp`, `well_rate`, and `production_temperature`.
  - `C/F`: pressure/temperature/CO2 saturation profile + monitor; when wells exist, require `well_bhp`, `well_rate`, `production_temperature`, `phase_fraction`, `water_cut`, and `co2_production_rate`.
  - All three families require reference contour overlays, profile overlays, and monitor overlays; numeric tolerances remain placeholders (`-1.0`, `is_placeholder=true`) until later tightening.
- Template/common-module boundary is frozen:
  - Case templates may only keep case metadata, geometry/mesh construction, boundary and well setup, property switches, discretization/assembly hooks, solve orchestration, case-specific postprocess hooks, and the four staged entry points.
  - Shared validation, reference I/O, grid/time-step studies, Matlab script generation, summary/report aggregation, and generic acceptance checks must not be reintroduced inside case templates.
  - The intended landing zones are `Case2D_Validation.*`, `Case2D_ReferenceIO.*`, `Case2D_Studies.*`, `Case2D_Matlab.*`, and the parallel `Case3D_*` modules.
- Documentation follows a dual-file workflow:
  - [A1_F12_TEMPLATE_SYSTEM_MEMORY.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_TEMPLATE_SYSTEM_MEMORY.md) stores long-lived decisions, the 72-case status snapshot, blockers, and next steps.
  - [A1_F12_进度记录.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_进度记录.md) is the execution log of milestones, ordered steps, completion evidence, and recent updates.
- Commit timing is governed by [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md); after each completed step, check this file before deciding whether to commit.

## Current Phase
- Phase: `M2` in progress; `M2-S1` through `M2-S4` are complete, and the next step is `M2-S5` donor-template thinning in isolated worktree `codex/a1-f12-exec`.

## Completed Items
- Added shared artifact helpers in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Artifacts.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.cpp).
- Added memory writer helper in [CaseCommon_Memory.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Memory.h) and [CaseCommon_Memory.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Memory.cpp).
- Added 72-case central catalog in [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h) and [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp).
- Added skeleton bootstrap helper in [CaseCommon_Skeleton.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.h) and [CaseCommon_Skeleton.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.cpp).
- Retrofitted staged APIs into the currently implemented 2D no-well templates: `A1`, `A2`, `A3`, `A4`, `B1`, `B3`.
- Disabled template-internal COMSOL autorun defaults for `B1` and `B3` to preserve weak coupling.
- Bootstrapped staged skeleton templates for `C1`, `A7`, `B7`, and `C7`.
- Updated [main.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp) to support `--case=A1` and `--stage=...`.
- Created the root-level execution tracker [A1_F12_进度记录.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_进度记录.md) in the isolated worktree and seeded M0 progress.
- Verified the isolated worktree can build in `Debug|x64` and that `--list` enumerates the full A1-F12 catalog plus legacy and auxiliary cases.
- Added [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md) to govern when the branch should commit and which local support files should stay uncommitted.
- Added 2D donor extraction helpers [Case2D_Validation.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Validation.h), [Case2D_Validation.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Validation.cpp), [Case2D_Studies.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Studies.h), and [Case2D_Studies.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Studies.cpp), and rewired A1 donor validation/study logic through them.
- Verified the `M2-S1` A1 donor extraction with fresh `Debug|x64` MSBuild and a full `--case=A1 --stage=solve_only` smoke run that completed successfully and regenerated the analytical summary plus grid/time study CSVs under `Test/Transient/FullCaseTest/H_CO2_ConstPP/h_co2_constpp_nofrac_nowell`.
- Added shared B1 donor extraction helpers [Case2D_ReferenceIO.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_ReferenceIO.h), [Case2D_ReferenceIO.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_ReferenceIO.cpp), [Case2D_Matlab.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.h), and [Case2D_Matlab.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.cpp), then rewired B1 donor `engineering/reference` CSV export, validation summary output, and no-fracture PT Matlab script generation through them.
- Hardened the shared `Case2DReferenceIO` output path against worktree-length regressions by adding the short-ASCII staging-file workflow and Win32 extended-path (`\\\\?\\`) commit fallback needed under `.worktrees\\codex\\a1-f12-exec`.
- Verified the `M2-S2` B1 donor extraction with fresh `Debug|x64` MSBuild and a full `--case=B1 --stage=validate_only` smoke run that no longer fails during engineering export; it now reaches the expected `COMSOL temperature reference files are missing` branch after completing the legacy engineering solve.
- Extended [Case2D_Matlab.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.h) and [Case2D_Matlab.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.cpp) with a fracture-aware single-fracture PT validation script writer that reuses the shared path-safe ASCII staging workflow for Matlab output under long worktree paths.
- Rewired the B3 donor [Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp) so fracture-aware profile/monitor engineering exports, station and schedule exports, profile reference loading, monitor reference loading, and Matlab script generation now route through `Case2DReferenceIO` and `Case2DMatlab` instead of donor-local file writers.
- Verified the `M2-S3` B3 donor extraction with fresh `Debug|x64` MSBuild and a worktree-local `--case=B3 --stage=prepare_reference` smoke run that completed to `validation_status=prepared_reference_inputs`, emitted fracture-aware engineering outputs plus validation summary artifacts under `Test/Transient/FullCaseTest/H_T_CO2_ConstPP/h_t_co2_constpp_singlefrac_nowell`, and correctly reported missing COMSOL profile/monitor reference files instead of failing during engineering/profile export.
- Consolidated the four 2D shared modules into a clearer public surface by documenting their roles in the headers and moving `Case2D_Validation.*` plus `Case2D_Studies.*` onto the same path-safe `Case2DReferenceIO::WriteAsciiFile(...)` output contract already used by `Case2D_ReferenceIO.*` and `Case2D_Matlab.*`.
- Verified the `M2-S4` module-surface freeze with fresh `Debug|x64` MSBuild and a worktree-local `--case=A1 --stage=solve_only` smoke run that regenerated `analytical_summary.txt`, `grid_convergence.csv`, and `time_sensitivity.csv` through the now-unified path-safe module outputs.
- Completed `M1-S1` by freezing the `CaseMetadata` field set in [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h) and [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp), adding catalog-side defaults and metadata contract validation, then re-verifying `Debug|x64` build and `--list`.
- Completed `M1-S2` by freezing `CaseArtifactPaths`, `ValidationSpec2D/3D`, and `WellTemplateSpec` in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Artifacts.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.cpp), then routing [CaseCommon_Skeleton.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.cpp) through the unified artifact paths and re-verifying a full `Debug|x64` build.
- Completed `M1-S3` by freezing family-level acceptance-policy placeholders in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h), wiring `A/D`, `B/E`, and `C/F` defaults in [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp), and re-verifying `Debug|x64` build while intentionally keeping numeric tolerances as placeholders.
- Completed `M1-S4` by freezing the template/common-module boundary in the execution documents: case templates must remain thin and may no longer absorb generic validation, reference I/O, study loops, Matlab generation, or report aggregation logic; those responsibilities are reserved for the planned `Case2D_*` and `Case3D_*` modules.
- Completed `M1-S5` by auditing `A1/A2/A3/A4/B1/B3/A7/B7/C1/C7` against the staged-entry and five-directory contracts: all ten cases are catalog-routable with four staged entry points; the four skeleton cases honor the unified artifact contract via `CaseCommon_Skeleton`; the implemented templates still contain the expected legacy drift that M2 must remove, including collapsed stage semantics in `A1/A2/A3/A4`, retained template-internal validation/reference logic in `B1/B3`, and the known `B1 validate_only` / `B3 RunSolveOnly` semantic deviations.

## Donor Templates
- `A1`: donor for analytic validation, feature-line profiles, and grid/dt studies.
- `B1`: donor for reference I/O, observation monitors, validation summaries, and Matlab generation.
- `B3`: donor for fracture-aware validation, fracture sampling, and fracture figure generation.

## Replaced Legacy Entries
- Legacy standalone entries for `A1/B1/A3/A2/A4/B3` now route through `RunFullWorkflow()` instead of raw `RunTestCase()`.
- Legacy descriptive `test_h_*` entries are still kept for backward compatibility during migration.

## How To Read Case Codes
### Family Prefix Meaning
- `A`: 2D, single-phase CO2, `N=1`, pressure diffusion.
- `B`: 2D, single-phase CO2, `N=2`, flow and heat transport.
- `C`: 2D, two-phase CO2/H2O, `N=3`, flow and heat transport.
- `D`: 3D, single-phase CO2, `N=1`, pressure diffusion.
- `E`: 3D, single-phase CO2, `N=2`, flow and heat transport.
- `F`: 3D, two-phase CO2/H2O, `N=3`, flow and heat transport.

### Index Meaning
- `1`: constant properties, no fracture, no well.
- `2`: variable properties, no fracture, no well.
- `3`: constant properties, single fracture, no well.
- `4`: variable properties, single fracture, no well.
- `5`: constant properties, complex fracture, no well.
- `6`: variable properties, complex fracture, no well.
- `7`: constant properties, no fracture, one injector + one producer.
- `8`: variable properties, no fracture, one injector + one producer.
- `9`: constant properties, single fracture, one injector + one producer.
- `10`: variable properties, single fracture, one injector + one producer.
- `11`: constant properties, complex fracture, one injector + one producer.
- `12`: variable properties, complex fracture, one injector + one producer.

### Status Legend
- `implemented`: catalog-routable and backed by a callable implementation.
- `skeleton`: has dedicated files, staged entry points, and artifact contract, but not yet a real solver/validation implementation.
- `planned`: registered in the catalog, but no dedicated implementation yet.

## Case Status Matrix
### 2D: A1-C12
- `A1` implemented
- `A2` implemented
- `A3` implemented
- `A4` implemented
- `A5` planned
- `A6` planned
- `A7` skeleton
- `A8` planned
- `A9` planned
- `A10` planned
- `A11` planned
- `A12` planned
- `B1` implemented
- `B2` planned
- `B3` implemented
- `B4` planned
- `B5` planned
- `B6` planned
- `B7` skeleton
- `B8` planned
- `B9` planned
- `B10` planned
- `B11` planned
- `B12` planned
- `C1` skeleton
- `C2` planned
- `C3` planned
- `C4` planned
- `C5` planned
- `C6` planned
- `C7` skeleton
- `C8` planned
- `C9` planned
- `C10` planned
- `C11` planned
- `C12` planned

### 3D: D1-F12
- `D1` planned
- `D2` planned
- `D3` planned
- `D4` planned
- `D5` planned
- `D6` planned
- `D7` planned
- `D8` planned
- `D9` planned
- `D10` planned
- `D11` planned
- `D12` planned
- `E1` planned
- `E2` planned
- `E3` planned
- `E4` planned
- `E5` planned
- `E6` planned
- `E7` planned
- `E8` planned
- `E9` planned
- `E10` planned
- `E11` planned
- `E12` planned
- `F1` planned
- `F2` planned
- `F3` planned
- `F4` planned
- `F5` planned
- `F6` planned
- `F7` planned
- `F8` planned
- `F9` planned
- `F10` planned
- `F11` planned
- `F12` planned

## Blockers
- `N=1 + wells` is still blocked in [RunGeneric_impl.hpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/FIM_TransientEngine/RunGeneric_impl.hpp) because the pressure-only AD route explicitly rejects wells.
- `A1/A2/A3/A4` still collapse `solve_only/prepare_reference/validate_only/full_workflow` onto one legacy execution path and do not yet use the unified five-directory artifact split; `M2-S1` only extracted A1 donor validation/study helpers, not the staged-semantics split.
- `B1/B3` still keep large amounts of generic validation/reference/study/Matlab logic inside template files; `B1 validate_only` still replays a fresh engineering solve internally, and `B3 RunSolveOnly` currently routes through the prepare-reference plan.
- `B1` now routes much of its `engineering/reference` I/O and Matlab generation through shared modules, but it still contains transitional dead bodies behind wrapper returns (`WriteValidationSummary`, `WriteMatlabPlotScript`) and still needs `M2-S5` / `M3` cleanup to become a genuinely thin template.
- `B3` now routes fracture-aware engineering/reference I/O and Matlab generation through shared modules, but it still retains donor-local fracture geometry sampling, profile/monitor compare kernels, study aggregation, and summary/report bodies; `M2-S4/M2-S5` still need to extract or trim these remaining template-local validation pieces.
- The four 2D shared modules now have stable roles and a unified path-safe output contract, but `A1/B1/B3` still retain sizeable donor-local wrappers and dead bodies; the remaining work is template thinning, not another public-module reshuffle.
- Most of the 72 catalog cases are registered but still have `planned` or `skeleton` status.
- The VS/MSBuild build in this worktree currently depends on local `.vcxproj` edits to include the new `Case2D_*` source files; those project-file edits remain local build support and are not part of the tracked branch payload.

## Next Steps
- Start `M2-S5` by removing donor-local dead bodies and redundant wrappers from `A1/B1/B3`, while keeping the already frozen four-module public surface intact.
- Re-evaluate the branch against [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md), because `M2-S4` now forms a coherent module-surface-freeze checkpoint.
- Continue separating staged semantics from legacy execution paths in `A1/B1/B3`, not just moving helper functions.
- Lift the N=1 well restriction and promote `A7` from `skeleton` to `implemented`.

## Key Files
- [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md)
- [A1_F12_进度记录.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_进度记录.md)
- [A1_F12_TEMPLATE_SYSTEM_MEMORY.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_TEMPLATE_SYSTEM_MEMORY.md)
- [main.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp)
- [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h)
- [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp)
- [Case2D_Validation.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Validation.h)
- [Case2D_Validation.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Validation.cpp)
- [Case2D_Studies.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Studies.h)
- [Case2D_Studies.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Studies.cpp)
- [Case2D_ReferenceIO.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_ReferenceIO.h)
- [Case2D_ReferenceIO.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_ReferenceIO.cpp)
- [Case2D_Matlab.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.h)
- [Case2D_Matlab.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Case2D_Matlab.cpp)
- [Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp)
- [Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp)
- [Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp)
- [Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp)
