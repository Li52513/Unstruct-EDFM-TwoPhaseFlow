$ErrorActionPreference = 'Stop'

$workspaceRoot = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path))
$exePath = Join-Path $workspaceRoot 'x64\Debug\Bin\2D-Unstr-Quadrilateral-EDFM.exe'

if (-not (Test-Path $exePath)) {
    throw "Missing test executable: $exePath"
}

$expectedCases = @(
    'test_h_t_co2_constpp_nofrac_nowell',
    'test_h_t_co2_constpp_nofrac_nowell_grid',
    'test_h_t_co2_constpp_nofrac_nowell_dt',
    'test_h_t_co2_constpp_nofrac_nowell_all'
)

$listOutput = & $exePath --list
foreach ($caseName in $expectedCases) {
    if (-not ($listOutput | Select-String -SimpleMatch $caseName)) {
        throw "Missing case in --list output: $caseName"
    }
}

$caseDir = Join-Path $workspaceRoot 'Test\Transient\FullCaseTest\H_T_CO2_ConstPP\h_t_co2_constpp_nofrac_nowell'
$engineeringDir = Join-Path $caseDir 'engineering'

if (-not (Test-Path $engineeringDir)) {
    try {
        & $exePath --case test_h_t_co2_constpp_nofrac_nowell | Out-Null
    } catch {
        Write-Host "Baseline case returned non-zero while preparing artifacts: $($_.Exception.Message)"
    }
}

$requiredPaths = @(
    $engineeringDir,
    (Join-Path $caseDir 'reference\analytic'),
    (Join-Path $caseDir 'reference\comsol_input'),
    (Join-Path $caseDir 'reference\comsol'),
    (Join-Path $caseDir 'report\validation_summary.md'),
    (Join-Path $caseDir 'metrics.csv'),
    (Join-Path $caseDir 'validation_summary.csv')
)

foreach ($path in $requiredPaths) {
    if (-not (Test-Path $path)) {
        throw "Missing validation artifact: $path"
    }
}

Write-Host 'HT const-property no-fracture validation contract check passed.'
