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
- Phase: `M5-S1` is complete in isolated worktree `codex/a1-f12-exec`; `A2-A6` are now all `implemented`, with `A5/A6` promoted by dedicated complex-fracture no-well wrappers on top of the frozen no-well and well-template contracts. `M4` remains numerically in progress only because `C1` has not yet closed well-level reference automation, `M5` is now in progress, and the next mainline step is `M5-S2`. In parallel `C0`, `B1` true COMSOL reference closure and the `C1` carrier-model checkpoint are now complete, but both `M3` and `M4` stay numerically open until `C1` is backed by a true COMSOL physics/study/solver reference solve.

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
- Completed `M2-S5` donor-template thinning by adding [tools/checks/check_m2_s5_template_thinning.ps1](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/tools/checks/check_m2_s5_template_thinning.ps1) as a contract check, removing A1/B1/B3 donor-local pure wrappers, duplicate CSV helpers, and unreachable fallback bodies, and rewiring the remaining call sites to `Case2DValidation`, `Case2DReferenceIO`, and `Case2DMatlab` directly.
- Verified `M2-S5` with fresh `Debug|x64` MSBuild, a passing `check_m2_s5_template_thinning.ps1`, and staged smoke runs for `A1 --stage=solve_only`, `B1 --stage=validate_only`, and `B3 --stage=prepare_reference`.
- Completed `M1-S1` by freezing the `CaseMetadata` field set in [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h) and [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp), adding catalog-side defaults and metadata contract validation, then re-verifying `Debug|x64` build and `--list`.
- Completed `M1-S2` by freezing `CaseArtifactPaths`, `ValidationSpec2D/3D`, and `WellTemplateSpec` in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Artifacts.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.cpp), then routing [CaseCommon_Skeleton.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Skeleton.cpp) through the unified artifact paths and re-verifying a full `Debug|x64` build.
- Completed `M1-S3` by freezing family-level acceptance-policy placeholders in [CaseCommon_Artifacts.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Artifacts.h) and [CaseCommon_Catalog.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.h), wiring `A/D`, `B/E`, and `C/F` defaults in [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/CaseCommon_Catalog.cpp), and re-verifying `Debug|x64` build while intentionally keeping numeric tolerances as placeholders.
- Completed `M1-S4` by freezing the template/common-module boundary in the execution documents: case templates must remain thin and may no longer absorb generic validation, reference I/O, study loops, Matlab generation, or report aggregation logic; those responsibilities are reserved for the planned `Case2D_*` and `Case3D_*` modules.
- Completed `M1-S5` by auditing `A1/A2/A3/A4/B1/B3/A7/B7/C1/C7` against the staged-entry and five-directory contracts: all ten cases are catalog-routable with four staged entry points; the four skeleton cases honor the unified artifact contract via `CaseCommon_Skeleton`; the implemented templates still contain the expected legacy drift that M2 must remove, including collapsed stage semantics in `A1/A2/A3/A4`, retained template-internal validation/reference logic in `B1/B3`, and the known `B1 validate_only` / `B3 RunSolveOnly` semantic deviations.
- Completed `M3-S1` by splitting A1 into a real staged thin-template runner: added `RunStageByKeyImpl(...)`, analytical artifact-only `prepare_reference`, `solve_only` output under `engineering/`, and `validate_only/full_workflow` paths that now materialize `studies/grid_convergence.csv` plus `studies/time_sensitivity.csv` under the template-system case root; verified with a passing `check_m3_s1_a1_stage_split.ps1`, fresh `Debug|x64` MSBuild, and four staged smoke runs for `A1 --stage=prepare_reference`, `solve_only`, `validate_only`, and `full_workflow`.
- Completed `M3-S2` by turning B1 into a real staged thin-template donor: added `RunStageByKeyImpl(...)`, `ConfigureSummaryPaths(...)`, template-system five-directory path wiring (`studies/`, `figures/`, `engineering/`, `reference/`, `report/`, plus `report/scripts/`), stage manifest/reference-contract/status writers, and a real `prepare_reference` path that emits property tables plus reference/COMSOL specs without running the transient solve; verified with a passing `check_m3_s2_b1_thin_template.ps1`, fresh `Debug|x64` MSBuild, a successful `B1 --stage=prepare_reference`, a successful `B1 --stage=solve_only`, and a `B1 --stage=validate_only` smoke run that now exits through `missing_reference` while correctly backfilling `stage_manifest.txt` and `template_status.md`.
- Completed `M3-S3` by promoting C1 from `skeleton` to a real staged `implemented` no-well template: replaced the skeleton-only stage stubs with a real `RunStageByKeyImpl(...)`, wired `prepare_reference` to emit the C1 property/reference/profile/monitor contracts, wired `solve_only` to run the real 2D `N=3` CO2/H2O no-fracture/no-well transient solve, and wired `validate_only/full_workflow` to complete engineering output then stop through the expected weak-coupling `missing_reference` path when COMSOL payloads are absent; verified with a passing `check_m3_s3_c1_thin_template.ps1`, fresh `Debug|x64` MSBuild, a successful `C1 --stage=prepare_reference`, a successful `C1 --stage=solve_only`, and a `C1 --stage=validate_only` smoke run that now backfills `stage_manifest.txt`, `template_status.md`, `validation_summary.md`, and the engineering profile/monitor CSV payload under the template-system case root.
- Completed `M3-S4` by freezing the no-well sample/output contract across `A1/B1/C1`: added [A1_F12_NoWell_Template_Contract.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_NoWell_Template_Contract.md) to lock five-directory responsibilities, stage semantics, mandatory profile/monitor payloads, frozen snapshot tags, and CSV schemas; added `tools/checks/check_m3_s4_nowell_contract.ps1` with long-path-safe existence checks; verified the contract with a passing worktree run of `check_m3_s4_nowell_contract.ps1`.
- Completed `M3-S5` by replacing the remaining legacy primary no-well entry points in [main.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/main.cpp): introduced `CatalogAliasEntry`, `legacyCatalogAliases`, `PrintCatalogAliases(...)`, and `ResolveCatalogCaseAlias(...)`, removed the direct `A1/B1` primary case lambdas from the legacy case table, and redirected the compatibility names `test_h_co2_constpp_nofrac_nowell`, `test_h_t_co2_constpp_nofrac_nowell`, and `test_h_tp_co2h2o_constpp_nofrac_nowell` through the central A1/B1/C1 catalog path; verified with a passing `check_m3_s5_catalog_aliases.ps1`, fresh `Debug|x64` MSBuild, `--list`, and three alias-based `prepare_reference` smoke runs.
- Completed `M4-S1` by enabling `N=1 + wells` in the pressure-only transient engine: [RunGeneric_impl.hpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/FIM_TransientEngine/RunGeneric_impl.hpp) now accepts normalized well schedules, instantiates `WellDOFManager<1>`, assembles well equations, and updates well-block pressures; [WellDOFManager.h](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/WellDOFManager.h) now exposes the explicit `N==1` gradient mapping; verified with `check_m4_s1_n1_well_support.ps1` and fresh `Debug|x64` MSBuild.
- Completed `M4-S2` by promoting [Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.cpp) from `skeleton` to a real staged `implemented` A7 template: the file now owns a real `RunStageByKeyImpl(...)`, weak-coupling `prepare_reference` artifact generation, a real `solve_only` N=1 inj/prod run, and `validate_only/full_workflow` paths that stop through `missing_reference` when COMSOL payloads are absent; verified with `check_m4_s2_a7_thin_template.ps1`, fresh `Debug|x64` MSBuild, and a successful `--case=A7 --stage=solve_only` smoke run that completed 146 steps with two active well blocks and wrote engineering VTK plus run-summary artifacts under the template-system case root.
- Fixed the template-system output-root drift uncovered during `M4-S2`: [CaseCommon_Catalog.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/CaseCommon_Catalog.cpp) now anchors `BuildOutputRoot()` to the real worktree/repo root via `weakly_canonical(__FILE__.parent_path())`, so staged runs no longer depend on the launch working directory; verified with `check_m4_s2_output_root_anchor.ps1`, fresh `Debug|x64` MSBuild, and a repeat `A7 --stage=solve_only` launched from the main repository root that now reports and writes the output directory under the real `a1-f12-exec` worktree path.
- Completed `M4-S3` by promoting [Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.cpp) from `skeleton` to a real staged `implemented` B7 template: the file now owns a real `RunStageByKeyImpl(...)`, weak-coupling `prepare_reference` artifact generation, a real `solve_only` N=2 thermal inj/prod run, and `validate_only/full_workflow` paths that stop through `missing_reference` when COMSOL profile/monitor/well payloads are absent; verified with `check_m4_s3_b7_thin_template.ps1`, fresh `Debug|x64` MSBuild, and successful `--case=B7 --stage=prepare_reference`, `solve_only`, and `validate_only` smoke runs.
- Completed `M4-S4` by promoting [Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.cpp) from `skeleton` to a real staged `implemented` C7 template: the file now owns a real `RunStageByKeyImpl(...)`, weak-coupling `prepare_reference` artifact generation, a real `solve_only` N=3 two-phase thermal inj/prod run, and `validate_only/full_workflow` paths that stop through `missing_reference` when COMSOL profile/monitor/well payloads are absent; startup instability at the injector thermal front was stabilized by extending the control-ramp, and verification passed with `check_m4_s4_c7_thin_template.ps1`, fresh `Debug|x64` MSBuild, and successful `--case=C7 --stage=prepare_reference`, `solve_only`, and `validate_only` smoke runs.
- Completed `M4-S5` by freezing the injector-producer sample contract in [A1_F12_Well_Template_Contract.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_Well_Template_Contract.md): the document now freezes the canonical well-mode layout, default injector/producer placement and controls, family-specific thermal deltas, required `well_schedule` and `well_timeseries` schemas, and `report/scripts` obligations for `A7/B7/C7`. Verification passed with [tools/checks/check_m4_s5_well_contract.ps1](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/tools/checks/check_m4_s5_well_contract.ps1), which checks the frozen contract against the current `A7/B7/C7` template-system case roots.
- Completed `M5-S1` by extending the remaining 2D `N=1` no-well A-family to `A2-A6`: [Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.cpp) and [Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.cpp) now expose deterministic complex-fracture plan keys, while [Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.cpp) and [Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.cpp](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.cpp) provide dedicated A5/A6 wrappers that are routable from the catalog and legacy `main.cpp` names. Verification passed with `check_m5_s1_a_family_extension.ps1`, fresh `Debug|x64` MSBuild, `--list`, `A5 --stage=solve_only`, and a fresh post-tuning `A6 --stage=solve_only` that completed successfully on a smoke-scale `24x3`, `t_end=2.0e4 s` complex-fracture plan.
- In parallel `codex/a1-f12-c0`, completed the real `B1` COMSOL reference closure: shared automation now resolves worktree-relative paths safely, the donor wrapper is long-path safe, real `reference/comsol/*` payloads plus a real `.mph` can be produced, and `B1 --stage=validate_only` reaches `validation_status=passed`.
- In parallel `codex/a1-f12-c0`, completed the current `C1` external automation checkpoint: `prepare_reference` now exports a non-empty `wrapper_relpath`, the external Java/PowerShell builder emits a real COMSOL carrier-model `.mph`, the `sg_ref_*` saturation-carrying reference schema is frozen, and `C1 --stage=validate_only` now reaches `passed` against that carrier payload.

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
- `A5` implemented
- `A6` implemented
- `A7` implemented
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
- `B7` implemented
- `B8` planned
- `B9` planned
- `B10` planned
- `B11` planned
- `B12` planned
- `C1` implemented
- `C2` planned
- `C3` planned
- `C4` planned
- `C5` planned
- `C6` planned
- `C7` implemented
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
- `A7` is now a real staged template, but `validate_only/full_workflow` still rerun the engineering solve and then stop through `missing_reference` until COMSOL pressure reference payloads and persistent engineering snapshots are available.
- `B7` is now a real staged template, but `validate_only/full_workflow` still rerun the engineering solve and then stop through `missing_reference` until COMSOL thermal profile/monitor/well reference payloads and persistent engineering snapshots are available.
- `A2/A3/A4/A5/A6` still collapse `solve_only/prepare_reference/validate_only/full_workflow` onto one legacy execution path and do not yet use the unified five-directory artifact split; `A5/A6` are now dedicated complex-fracture wrappers, but they still delegate to the legacy A3/A4 execution path under the hood. Only `A1` has completed the staged-semantics split so far.
- `A1` still reruns the engineering solve during `validate_only` until engineering snapshot persistence is implemented.
- Parallel `C0` work has already shown that `B1` can consume a true COMSOL reference payload and reach `validation_status=passed`, but that closure is not yet merged back into `codex/a1-f12-exec`.
- Parallel `C0` work has already replaced `C1 missing_reference` with a carrier-model `.mph` plus `sg_ref_*` payload contract and a passing template-side acceptance run, but `C1` still lacks a true COMSOL physics/study/solver reference solve and therefore does not yet satisfy the stricter "independent COMSOL reference" completion line.
- `B3` still keeps large amounts of generic validation/reference/study/Matlab logic inside the template file, and `RunSolveOnly` still routes through the prepare-reference plan.
- `B3` now routes fracture-aware engineering/reference I/O and Matlab generation through shared modules, but it still retains donor-local fracture geometry sampling, profile/monitor compare kernels, study aggregation, and summary/report bodies; those remaining validation kernels are deferred to `M3`.
- The four 2D shared modules now have stable roles and a unified path-safe output contract; the remaining 2D no-well work is template thinning and staged-semantics separation, not another public-module reshuffle.
- Most of the 72 catalog cases are registered but still have `planned` or `skeleton` status.
- The VS/MSBuild build in this worktree currently depends on local `.vcxproj` edits to include the new `Case2D_*` source files; those project-file edits remain local build support and are not part of the tracked branch payload.

## Next Steps
- Re-evaluate the branch against [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md), because `M5-S1` now forms a coherent A-family extension checkpoint.
- Start `M5-S2` by extending the remaining 2D `N=2` no-well `B` family (`B2-B6`) on top of the frozen no-well/well-template contracts and the now-complete A-family `A2-A6` coverage.
- The two bullets immediately below are legacy carry-over notes from the previous checkpoint and are superseded by the `M5-S1` / `M5-S2` pointers above.
- Re-evaluate the branch against [A1_F12_提交策略.md](/D:/Yongwei/博士生涯/100-Research/110Code/111-2D_EDFM_FVM_CO2PlumingSystem/B-Code/2D-Unstr-Quadrilateral-EDFM/.worktrees/codex/a1-f12-exec/A1_F12_提交策略.md), because `M4-S5` now forms a coherent well-contract freeze checkpoint.
- Start `M5-S1` by extending the remaining 2D no-well `A` family (`A2-A6`) on top of the frozen `A1/B1/C1` and `A7/B7/C7` template contracts.
- Sync the documented `B1` true-COMSOL closure and `C1` carrier-model checkpoint from the parallel `C0` branch back into `codex/a1-f12-exec` when the merge window is clear.
- Keep `C0` open until `C1` replaces the carrier-model route with a true COMSOL physics/study/solver solve that exports independent `pressure/temperature/co2_saturation` reference payloads.

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
