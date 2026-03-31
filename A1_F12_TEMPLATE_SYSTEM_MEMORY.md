# A1-F12 Template System Memory

## Frozen Decisions
- System title is `A1-F12`; 2D cases are `A1-C12` and 3D cases are `D1-F12`.
- Keep descriptive `.h/.cpp` template filenames; `A1-F12` is provided through a central catalog/dispatcher mapping.
- Every catalog case exposes four staged entry points: `RunSolveOnly`, `RunPrepareReference`, `RunValidateOnly`, `RunFullWorkflow`.
- All template artifacts follow the five-directory contract: `studies/`, `figures/`, `engineering/`, `reference/`, `report/`.
- COMSOL remains weakly coupled: templates emit/consume reference artifacts, while external Java/PowerShell automation executes COMSOL.
- Well defaults are frozen as injector-rate + producer-BHP; two-phase well cases inject CO2; thermal well cases use cold injection.
- `CaseMetadata` is frozen to carry `case_code`, `dispatcher_key`, `case_slug`, `description`, `reference_mode`, `implementation_status`, `output_root`, `well_control_policy`, `injection_fluid`, `thermal_injection_policy`, plus the six enum axes `dimension/equation_mode/property_mode/fracture_mode/well_mode`.
- `CaseArtifactPaths` is frozen to carry the five-directory contract plus the derived file paths `engineering_stage_manifest_path`, `reference_contract_path`, and `report_status_markdown_path`.
- `ValidationSpec2D/3D` is frozen as the shared contract for validation variables, profile snapshots, characteristic lines, observation points, grid-study tags, time-step study values, and required well-series outputs; 3D additionally carries characteristic slice definitions.
- `WellTemplateSpec` is frozen as the shared contract for injector/producer naming, control policy, injection fluid, thermal policy, default quarter/three-quarter placement, control modes, target values, injection temperature, and the “avoid direct fracture completion” default.
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
- Phase: M1 completed; next step is M2-S1 extraction of the A1 analytical validation/study chain in isolated worktree `codex/a1-f12-exec`

## Completed Items
- Added shared artifact helpers in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Artifacts.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.cpp).
- Added memory writer helper in [CaseCommon_Memory.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Memory.h) and [CaseCommon_Memory.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Memory.cpp).
- Added 72-case central catalog in [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h) and [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp).
- Added skeleton bootstrap helper in [CaseCommon_Skeleton.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.h) and [CaseCommon_Skeleton.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.cpp).
- Retrofitted staged APIs into the currently implemented 2D no-well templates: A1, A2, A3, A4, B1, B3.
- Disabled template-internal COMSOL autorun defaults for B1 and B3 to preserve weak coupling.
- Bootstrapped staged skeleton templates for C1, A7, B7 and C7.
- Updated [main.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp) to support `--case=A1` and `--stage=...`.
- Created the root-level execution tracker [A1_F12_进度记录.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_进度记录.md) in the isolated worktree and seeded M0 progress.
- Verified the isolated worktree can build in `Debug|x64` and that `--list` enumerates the full A1-F12 catalog plus legacy/auxiliary cases.
- Added [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md) to govern when the branch should commit and which local support files should stay uncommitted.
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
- Legacy standalone entries for A1/B1/A3/A2/A4/B3 now route through `RunFullWorkflow()` instead of raw `RunTestCase()`.
- Legacy descriptive `test_h_*` entries are still kept for backward compatibility during migration.

## How To Read Case Codes
### Family Prefix Meaning
- `A`: 2D, 单相 CO2, N=1, 压力扩散
- `B`: 2D, 单相 CO2, N=2, 渗流传热
- `C`: 2D, 两相 CO2/H2O, N=3, 渗流传热
- `D`: 3D, 单相 CO2, N=1, 压力扩散
- `E`: 3D, 单相 CO2, N=2, 渗流传热
- `F`: 3D, 两相 CO2/H2O, N=3, 渗流传热

### Index Meaning
- `1`: 常物性, 无裂缝, 无井
- `2`: 变物性, 无裂缝, 无井
- `3`: 常物性, 单裂缝, 无井
- `4`: 变物性, 单裂缝, 无井
- `5`: 常物性, 复杂裂缝, 无井
- `6`: 变物性, 复杂裂缝, 无井
- `7`: 常物性, 无裂缝, 一注一采
- `8`: 变物性, 无裂缝, 一注一采
- `9`: 常物性, 单裂缝, 一注一采
- `10`: 变物性, 单裂缝, 一注一采
- `11`: 常物性, 复杂裂缝, 一注一采
- `12`: 变物性, 复杂裂缝, 一注一采

### Examples
- `A1` = 2D, 单相 CO2, 常物性, 无裂缝, 无井, 压力扩散
- `B7` = 2D, 单相 CO2, 常物性, 无裂缝, 一注一采, 渗流传热
- `C10` = 2D, 两相 CO2/H2O, 变物性, 单裂缝, 一注一采, 渗流传热
- `E3` = 3D, 单相 CO2, 常物性, 单裂缝, 无井, 渗流传热
- `F12` = 3D, 两相 CO2/H2O, 变物性, 复杂裂缝, 一注一采, 渗流传热

### Status Legend
- `implemented`: 已接入 catalog 且存在可调用实现
- `skeleton`: 已有独立模板文件、阶段入口和 artifact 合同，但尚未完成真实求解/验证
- `planned`: 已在 catalog 登记，但尚未建立对应模板实现

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
- `A1/A2/A3/A4` still collapse `solve_only/prepare_reference/validate_only/full_workflow` onto one legacy execution path and do not yet use the unified five-directory artifact split.
- `B1/B3` still keep large amounts of generic validation/reference/study/Matlab logic inside template files; `B1 validate_only` still replays a fresh engineering solve internally, and `B3 RunSolveOnly` currently routes through the prepare-reference plan.
- Most of the 72 catalog cases are registered but still have `planned` or `skeleton` status.

## Next Steps
- Start `M2-S1` by extracting the A1 analytical validation chain, profile export, and grid/dt study helpers into the planned 2D common modules.
- Re-evaluate the branch against [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md), because M1 is now fully closed and the branch has a coherent audit-and-governance checkpoint.
- Extract B1 validation/reference logic into shared 2D modules so staged semantics are physically separated, not just API-separated.
- Lift the N=1 well restriction and promote A7 from `skeleton` to `implemented`.

## Key Files
- [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md)
- [A1_F12_进度记录.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_进度记录.md)
- [A1_F12_TEMPLATE_SYSTEM_MEMORY.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_TEMPLATE_SYSTEM_MEMORY.md)
- [main.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/main.cpp)
- [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h)
- [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp)
- [Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp)
- [Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp)
- [Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp)
- [Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp)
