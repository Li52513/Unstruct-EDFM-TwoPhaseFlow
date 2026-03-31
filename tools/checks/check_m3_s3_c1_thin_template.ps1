$ErrorActionPreference = "Stop"

$target = Join-Path $PSScriptRoot "..\..\Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.cpp"
$target = [System.IO.Path]::GetFullPath($target)

if (-not (Test-Path $target)) {
    throw "[M3-S3] Missing target file: $target"
}

$text = Get-Content -Raw -Path $target
$failures = New-Object System.Collections.Generic.List[string]

if ($text -notmatch 'RunStageByKeyImpl\s*\(') {
    $failures.Add("Missing RunStageByKeyImpl(...) staged dispatcher.")
}

if ($text -notmatch 'RunPrepareReference\s*\(\)\s*\{[\s\S]*RunStageByKeyImpl\s*\(\s*"h_tp_co2h2o_constpp_nofrac_nowell"\s*,\s*CaseCommon::CaseStage::PrepareReference\s*\)') {
    $failures.Add("RunPrepareReference() must dispatch through RunStageByKeyImpl(..., PrepareReference).")
}

if ($text -notmatch 'RunSolveOnly\s*\(\)\s*\{[\s\S]*RunStageByKeyImpl\s*\(\s*"h_tp_co2h2o_constpp_nofrac_nowell"\s*,\s*CaseCommon::CaseStage::SolveOnly\s*\)') {
    $failures.Add("RunSolveOnly() must dispatch through RunStageByKeyImpl(..., SolveOnly).")
}

if ($text -notmatch 'RunValidateOnly\s*\(\)\s*\{[\s\S]*RunStageByKeyImpl\s*\(\s*"h_tp_co2h2o_constpp_nofrac_nowell"\s*,\s*CaseCommon::CaseStage::ValidateOnly\s*\)') {
    $failures.Add("RunValidateOnly() must dispatch through RunStageByKeyImpl(..., ValidateOnly).")
}

if ($text -notmatch 'RunFullWorkflow\s*\(\)\s*\{[\s\S]*RunStageByKeyImpl\s*\(\s*"h_tp_co2h2o_constpp_nofrac_nowell"\s*,\s*CaseCommon::CaseStage::FullWorkflow\s*\)') {
    $failures.Add("RunFullWorkflow() must dispatch through RunStageByKeyImpl(..., FullWorkflow).")
}

if ($text -match 'CaseCommon::ThrowStageNotImplemented') {
    $failures.Add("C1 must no longer call CaseCommon::ThrowStageNotImplemented after promotion from skeleton.")
}

if ($text -match 'CaseCommon::PrepareSkeletonArtifacts') {
    $failures.Add("C1 must stop using CaseCommon::PrepareSkeletonArtifacts once it becomes implemented.")
}

if ($failures.Count -gt 0) {
    Write-Host "[M3-S3][FAIL] C1 thin-template contract is not satisfied."
    foreach ($failure in $failures) {
        Write-Host " - $failure"
    }
    exit 1
}

Write-Host "[M3-S3][PASS] C1 thin-template contract is satisfied."
