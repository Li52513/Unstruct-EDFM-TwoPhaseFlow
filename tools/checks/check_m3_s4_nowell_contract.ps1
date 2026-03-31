$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$contractDoc = Join-Path $repoRoot "A1_F12_NoWell_Template_Contract.md"
$caseRoot = Join-Path $repoRoot "Test\Transient\A1_F12_TemplateSystem"
$failures = New-Object System.Collections.Generic.List[string]

function Convert-ToLongPath {
    param(
        [string]$Path
    )

    if ([string]::IsNullOrWhiteSpace($Path)) {
        return $Path
    }

    if ($Path.StartsWith("\\?\")) {
        return $Path
    }

    if ($Path.StartsWith("\\")) {
        return "\\?\UNC\" + $Path.TrimStart("\")
    }

    return "\\?\" + $Path
}

function Test-ExistingPath {
    param(
        [string]$Path
    )

    $fullPath = [System.IO.Path]::GetFullPath($Path)
    $longPath = Convert-ToLongPath $fullPath

    return [System.IO.File]::Exists($longPath) -or [System.IO.Directory]::Exists($longPath)
}

function Require-Path {
    param(
        [string]$Path,
        [string]$Message
    )

    if (-not (Test-ExistingPath $Path)) {
        $failures.Add($Message + " [$Path]")
    }
}

Require-Path $contractDoc "Missing frozen no-well contract document."

$caseMap = @{
    "A1" = @{
        Root = Join-Path $caseRoot "A1_2d_n1_co2_constpp_nofrac_nowell"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "reference\reference_contract.txt",
            "report\template_status.md",
            "studies\grid_convergence.csv",
            "studies\time_sensitivity.csv"
        )
    }
    "B1" = @{
        Root = Join-Path $caseRoot "B1_2d_n2_co2_constpp_nofrac_nowell"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "engineering\property_table.csv",
            "engineering\profile_station_definitions.csv",
            "engineering\monitor_point_definitions.csv",
            "engineering\profile_report_schedule.csv",
            "engineering\monitor_sample_schedule.csv",
            "engineering\eng_monitor_timeseries.csv",
            "engineering\eng_profile_matrix_horizontal_t010pct.csv",
            "engineering\eng_profile_matrix_horizontal_t050pct.csv",
            "engineering\eng_profile_matrix_horizontal_t100pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t010pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t050pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t100pct.csv",
            "reference\reference_contract.txt",
            "reference\comsol_input\property_table.csv",
            "report\template_status.md",
            "report\validation_summary.md",
            "report\validation_summary.csv"
        )
    }
    "C1" = @{
        Root = Join-Path $caseRoot "C1_2d_n3_co2h2o_constpp_nofrac_nowell"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "engineering\property_table.csv",
            "engineering\profile_station_definitions.csv",
            "engineering\monitor_point_definitions.csv",
            "engineering\profile_report_schedule.csv",
            "engineering\monitor_sample_schedule.csv",
            "engineering\eng_monitor_timeseries.csv",
            "engineering\eng_profile_matrix_horizontal_t010pct.csv",
            "engineering\eng_profile_matrix_horizontal_t050pct.csv",
            "engineering\eng_profile_matrix_horizontal_t100pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t010pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t050pct.csv",
            "engineering\eng_profile_matrix_vertical_midline_t100pct.csv",
            "reference\reference_contract.txt",
            "report\template_status.md",
            "report\validation_summary.md",
            "report\scripts\plot_validation_results.m"
        )
    }
}

foreach ($caseCode in @("A1", "B1", "C1")) {
    $root = $caseMap[$caseCode].Root
    Require-Path $root "$caseCode case root is missing."
    foreach ($relative in $caseMap[$caseCode].Required) {
        Require-Path (Join-Path $root $relative) "$caseCode missing required no-well contract artifact: $relative"
    }
}

if ($failures.Count -gt 0) {
    Write-Host "[M3-S4][FAIL] No-well contract freeze is not satisfied." -ForegroundColor Red
    foreach ($failure in $failures) {
        Write-Host " - $failure"
    }
    exit 1
}

Write-Host "[M3-S4][PASS] No-well contract freeze is satisfied." -ForegroundColor Green
exit 0
