param(
    [string]$CaseRoot = '.\Test\Transient\FullCaseTest\H_T_CO2_ConstPP\h_t_co2_constpp_singlefrac_nowell',
    [double]$L2Threshold = 0.05,
    [double]$LinfThreshold = 0.10
)

$ErrorActionPreference = 'Stop'

function Require-File {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Missing required file: $Path"
    }
}

function Import-KeyValueCsv {
    param([string]$Path)
    $map = @{}
    foreach ($row in (Import-Csv -LiteralPath $Path)) {
        $map[$row.key] = [double]$row.value
    }
    return $map
}

function Get-Norms {
    param(
        [double[]]$Errors,
        [double]$Scale
    )
    if ($Errors.Count -eq 0) {
        throw 'No error samples available.'
    }

    $absErrors = @($Errors | ForEach-Object { [Math]::Abs($_) })
    $l1 = ($absErrors | Measure-Object -Average).Average
    $l2 = [Math]::Sqrt((($Errors | ForEach-Object { $_ * $_ }) | Measure-Object -Average).Average)
    $linf = ($absErrors | Measure-Object -Maximum).Maximum

    [pscustomobject]@{
        L1Norm = $l1 / $Scale
        L2Norm = $l2 / $Scale
        LinfNorm = $linf / $Scale
    }
}

$propertyTable = Join-Path $CaseRoot 'engineering\property_table.csv'
$matrixOnlyProfile = Join-Path $CaseRoot 'reference\comsol\comsol_matrix_only_profile_matrix_horizontal_t100pct.csv'

Require-File $propertyTable
Require-File $matrixOnlyProfile

$properties = Import-KeyValueCsv -Path $propertyTable
$profileRows = Import-Csv -LiteralPath $matrixOnlyProfile

$lx = $properties['lx']
$pLeft = $properties['p_left']
$pRight = $properties['p_right']
$deltaP = [Math]::Abs($pLeft - $pRight)

if ($deltaP -le 0.0) {
    throw 'Invalid pressure scale: |p_left - p_right| must be positive.'
}

$errors = New-Object System.Collections.Generic.List[double]
foreach ($row in $profileRows) {
    $x = [double]$row.target_x_m
    $pRef = [double]$row.p_ref_pa
    $pExpected = $pLeft + ($pRight - $pLeft) * ($x / $lx)
    $errors.Add($pRef - $pExpected)
}

$norms = Get-Norms -Errors $errors.ToArray() -Scale $deltaP
$pass = ($norms.L2Norm -le $L2Threshold) -and ($norms.LinfNorm -le $LinfThreshold)

Write-Host ("matrix_only_pressure_l1_norm={0}" -f $norms.L1Norm)
Write-Host ("matrix_only_pressure_l2_norm={0}" -f $norms.L2Norm)
Write-Host ("matrix_only_pressure_linf_norm={0}" -f $norms.LinfNorm)
Write-Host ("matrix_only_pressure_thresholds=l2<={0}, linf<={1}" -f $L2Threshold, $LinfThreshold)

if (-not $pass) {
    throw (
        "Matrix-only pressure alignment check failed: " +
        "l2_norm=$($norms.L2Norm), linf_norm=$($norms.LinfNorm)"
    )
}
