param(
    [string]$CaseRoot = '.\Test\Transient\FullCaseTest\H_T_CO2_ConstPP\h_t_co2_constpp_singlefrac_nowell'
)

$ErrorActionPreference = 'Stop'

function Require-File {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Missing required file: $Path"
    }
}

function Get-Span {
    param(
        [object[]]$Rows,
        [string]$Column
    )
    $measure = $Rows | Measure-Object -Property $Column -Minimum -Maximum
    return ([double]$measure.Maximum - [double]$measure.Minimum)
}

$referenceDir = Join-Path $CaseRoot 'reference\comsol'
$fractureProfilePath = Join-Path $referenceDir 'comsol_profile_fracture_tangent_t100pct.csv'
$monitorPath = Join-Path $referenceDir 'comsol_monitor_timeseries.csv'

Require-File $fractureProfilePath
Require-File $monitorPath

$fractureProfile = Import-Csv -LiteralPath $fractureProfilePath
$monitorRows = Import-Csv -LiteralPath $monitorPath

$profilePressureSpan = Get-Span -Rows $fractureProfile -Column 'p_ref_pa'
$profileTemperatureSpan = Get-Span -Rows $fractureProfile -Column 't_ref_k'

$fracturePressureColumns = @($monitorRows[0].PSObject.Properties.Name | Where-Object { $_ -like 'p_ref_fr*' })
$fractureTemperatureColumns = @($monitorRows[0].PSObject.Properties.Name | Where-Object { $_ -like 't_ref_fr*' })

if ($fracturePressureColumns.Count -eq 0 -or $fractureTemperatureColumns.Count -eq 0) {
    throw "Missing fracture monitor columns in $monitorPath"
}

$finalMonitor = $monitorRows | Sort-Object { [double]$_.target_time_s } | Select-Object -Last 1
$finalPressureValues = @($fracturePressureColumns | ForEach-Object { [double]$finalMonitor.$_ })
$finalTemperatureValues = @($fractureTemperatureColumns | ForEach-Object { [double]$finalMonitor.$_ })
$monitorPressureSpan = (($finalPressureValues | Measure-Object -Minimum -Maximum).Maximum - ($finalPressureValues | Measure-Object -Minimum -Maximum).Minimum)
$monitorTemperatureSpan = (($finalTemperatureValues | Measure-Object -Minimum -Maximum).Maximum - ($finalTemperatureValues | Measure-Object -Minimum -Maximum).Minimum)

$pressureOk = $profilePressureSpan -gt 1.0e5
$temperatureOk = $profileTemperatureSpan -gt 1.0
$monitorPressureOk = $monitorPressureSpan -gt 1.0e5
$monitorTemperatureOk = $monitorTemperatureSpan -gt 1.0

Write-Host ("profile_pressure_span_pa={0}" -f $profilePressureSpan)
Write-Host ("profile_temperature_span_k={0}" -f $profileTemperatureSpan)
Write-Host ("monitor_pressure_span_pa={0}" -f $monitorPressureSpan)
Write-Host ("monitor_temperature_span_k={0}" -f $monitorTemperatureSpan)

if (-not ($pressureOk -and $temperatureOk -and $monitorPressureOk -and $monitorTemperatureOk)) {
    throw (
        "COMSOL fracture response check failed: " +
        "profile_pressure_ok=$pressureOk, profile_temperature_ok=$temperatureOk, " +
        "monitor_pressure_ok=$monitorPressureOk, monitor_temperature_ok=$monitorTemperatureOk"
    )
}
