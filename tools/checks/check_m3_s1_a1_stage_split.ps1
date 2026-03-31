$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$target = Join-Path $repoRoot "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp"
$content = Get-Content $target -Raw

$legacyPatterns = @(
    'void\s+RunSolveOnly\(\)\s*\{\s*ExecutePlanByKeyImpl\("h_co2_constpp_nofrac_nowell"\);\s*\}',
    'void\s+RunPrepareReference\(\)\s*\{\s*ExecutePlanByKeyImpl\("h_co2_constpp_nofrac_nowell"\);\s*\}',
    'void\s+RunValidateOnly\(\)\s*\{\s*ExecutePlanByKeyImpl\("h_co2_constpp_nofrac_nowell"\);\s*\}',
    'void\s+RunFullWorkflow\(\)\s*\{\s*ExecutePlanByKeyImpl\("h_co2_constpp_nofrac_nowell"\);\s*\}'
)

$requiredPatterns = @(
    'RunStageByKeyImpl\(',
    'CaseStage::SolveOnly',
    'CaseStage::PrepareReference',
    'CaseStage::ValidateOnly',
    'CaseStage::FullWorkflow'
)

$violations = @()

foreach ($pattern in $legacyPatterns) {
    if ($content -match $pattern) {
        $violations += "legacy_stage_wrapper::$pattern"
    }
}

foreach ($pattern in $requiredPatterns) {
    if ($content -notmatch $pattern) {
        $violations += "missing_stage_split_marker::$pattern"
    }
}

if ($violations.Count -gt 0) {
    Write-Host "M3-S1 A1 stage-split contract FAILED." -ForegroundColor Red
    foreach ($violation in $violations) {
        Write-Host ("  {0}" -f $violation)
    }
    exit 1
}

Write-Host "M3-S1 A1 stage-split contract PASSED." -ForegroundColor Green
exit 0
