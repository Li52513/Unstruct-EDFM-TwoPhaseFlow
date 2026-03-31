$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$contractDoc = Join-Path $repoRoot "A1_F12_Well_Template_Contract.md"
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

function Read-Text {
    param(
        [string]$Path
    )

    $fullPath = [System.IO.Path]::GetFullPath($Path)
    $longPath = Convert-ToLongPath $fullPath
    if (-not [System.IO.File]::Exists($longPath)) {
        return $null
    }

    return [System.IO.File]::ReadAllText($longPath)
}

function Read-Lines {
    param(
        [string]$Path
    )

    $text = Read-Text $Path
    if ($null -eq $text) {
        return @()
    }

    return ($text -split "`r?`n")
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

function Require-TextContains {
    param(
        [string]$Path,
        [string]$Needle,
        [string]$Message
    )

    $text = Read-Text $Path
    if ($null -eq $text -or -not $text.Contains($Needle)) {
        $failures.Add($Message + " [$Path]")
    }
}

function Require-LineEquals {
    param(
        [string]$Path,
        [int]$Index,
        [string]$Expected,
        [string]$Message
    )

    $lines = Read-Lines $Path
    if ($lines.Count -le $Index -or $lines[$Index] -ne $Expected) {
        $failures.Add($Message + " [$Path]")
    }
}

function Require-LinePrefix {
    param(
        [string]$Path,
        [int]$Index,
        [string]$ExpectedPrefix,
        [string]$Message
    )

    $lines = Read-Lines $Path
    if ($lines.Count -le $Index -or -not $lines[$Index].StartsWith($ExpectedPrefix)) {
        $failures.Add($Message + " [$Path]")
    }
}

function Require-LineContains {
    param(
        [string]$Path,
        [int]$Index,
        [string]$Needle,
        [string]$Message
    )

    $lines = Read-Lines $Path
    if ($lines.Count -le $Index -or -not $lines[$Index].Contains($Needle)) {
        $failures.Add($Message + " [$Path]")
    }
}

Require-Path $contractDoc "Missing frozen well contract document."

$cases = @{
    "A7" = @{
        Root = Join-Path $caseRoot "A7_2d_n1_co2_constpp_nofrac_injprod"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "engineering\property_table.csv",
            "engineering\well_schedule.csv",
            "engineering\reference_spec.md",
            "engineering\run_summary.txt",
            "reference\reference_contract.txt",
            "report\template_status.md",
            "report\scripts\plot_validation_results.m"
        )
        StatusNeedles = @(
            '- Well policy: `injector rate + producer BHP`',
            '- Injection fluid: `co2`'
        )
    }
    "B7" = @{
        Root = Join-Path $caseRoot "B7_2d_n2_co2_constpp_nofrac_injprod"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "engineering\property_table.csv",
            "engineering\well_schedule.csv",
            "engineering\well_timeseries.csv",
            "engineering\reference_spec.md",
            "engineering\run_summary.txt",
            "reference\reference_contract.txt",
            "report\template_status.md",
            "report\validation_summary.md",
            "report\scripts\plot_validation_results.m"
        )
        StatusNeedles = @(
            '- Well policy: `injector rate + producer BHP`',
            '- Injection fluid: `co2`',
            '- Thermal policy: `cold_injection`',
            '- Injector temperature: `320 K`'
        )
    }
    "C7" = @{
        Root = Join-Path $caseRoot "C7_2d_n3_co2h2o_constpp_nofrac_injprod"
        Required = @(
            "studies",
            "figures",
            "engineering",
            "reference",
            "report",
            "report\scripts",
            "engineering\stage_manifest.txt",
            "engineering\property_table.csv",
            "engineering\well_schedule.csv",
            "engineering\well_timeseries.csv",
            "engineering\reference_spec.md",
            "engineering\run_summary.txt",
            "reference\reference_contract.txt",
            "report\template_status.md",
            "report\validation_summary.md",
            "report\scripts\plot_validation_results.m"
        )
        StatusNeedles = @(
            '- Well policy: `injector rate + producer BHP`',
            '- Injection fluid: `co2 into brine`',
            '- Thermal policy: `cold_injection`',
            '- Injector temperature: `320 K`'
        )
    }
}

$a7ScheduleHeader = "role,well_name,domain,control_mode,target_value,component_mode,completion_id,target_x_m,target_y_m,actual_x_m,actual_y_m,injection_is_co2"
$thermalScheduleHeader = "role,well_name,domain,control_mode,target_value,component_mode,completion_id,target_x_m,target_y_m,actual_x_m,actual_y_m,injection_is_co2,injection_temperature_k"
$wellTimeseriesHeader = "sample_id,time_s,well_name,role,control_mode,component_mode,target_value,actual_bhp_pa,actual_rate_kg_per_s,completion_pressure_pa,completion_temperature_k,reported_temperature_k,target_x_m,target_y_m,actual_x_m,actual_y_m"

foreach ($caseCode in @("A7", "B7", "C7")) {
    $root = $cases[$caseCode].Root
    Require-Path $root "$caseCode case root is missing."
    foreach ($relative in $cases[$caseCode].Required) {
        Require-Path (Join-Path $root $relative) "$caseCode missing required well contract artifact: $relative"
    }

    $statusPath = Join-Path $root "report\template_status.md"
    foreach ($needle in $cases[$caseCode].StatusNeedles) {
        Require-TextContains $statusPath $needle "$caseCode template status is missing frozen contract text: $needle"
    }
}

$a7SchedulePath = Join-Path $cases["A7"].Root "engineering\well_schedule.csv"
Require-LineEquals $a7SchedulePath 0 $a7ScheduleHeader "A7 well_schedule header drifted."
Require-LinePrefix $a7SchedulePath 1 "injector,INJ,matrix,rate,-1,gas," "A7 injector row drifted."
Require-LineContains $a7SchedulePath 1 ",100,20," "A7 injector target position drifted."
Require-LineContains $a7SchedulePath 1 ",true" "A7 injector CO2 flag drifted."
Require-LinePrefix $a7SchedulePath 2 "producer,PROD,matrix,bhp,9500000,total," "A7 producer row drifted."
Require-LineContains $a7SchedulePath 2 ",300,20," "A7 producer target position drifted."
Require-LineContains $a7SchedulePath 2 ",false" "A7 producer CO2 flag drifted."

foreach ($caseCode in @("B7", "C7")) {
    $schedulePath = Join-Path $cases[$caseCode].Root "engineering\well_schedule.csv"
    $timeseriesPath = Join-Path $cases[$caseCode].Root "engineering\well_timeseries.csv"

    Require-LineEquals $schedulePath 0 $thermalScheduleHeader "$caseCode well_schedule header drifted."
    Require-LinePrefix $schedulePath 1 "injector,INJ,matrix,rate,-1,gas," "$caseCode injector row drifted."
    Require-LineContains $schedulePath 1 ",100,20," "$caseCode injector target position drifted."
    Require-LineContains $schedulePath 1 ",true,320" "$caseCode injector thermal contract drifted."
    Require-LinePrefix $schedulePath 2 "producer,PROD,matrix,bhp,9500000,total," "$caseCode producer row drifted."
    Require-LineContains $schedulePath 2 ",300,20," "$caseCode producer target position drifted."
    Require-LineContains $schedulePath 2 ",false,-1" "$caseCode producer thermal contract drifted."

    Require-LineEquals $timeseriesPath 0 $wellTimeseriesHeader "$caseCode well_timeseries header drifted."
    Require-LinePrefix $timeseriesPath 1 "0,0,INJ,injector,rate,gas,-1," "$caseCode well_timeseries injector t0 row drifted."
    Require-LinePrefix $timeseriesPath 2 "0,0,PROD,producer,bhp,total,9500000," "$caseCode well_timeseries producer t0 row drifted."
}

if ($failures.Count -gt 0) {
    Write-Host "[M4-S5][FAIL] Well template contract freeze is not satisfied." -ForegroundColor Red
    foreach ($failure in $failures) {
        Write-Host " - $failure"
    }
    exit 1
}

Write-Host "[M4-S5][PASS] Well template contract freeze is satisfied." -ForegroundColor Green
exit 0
