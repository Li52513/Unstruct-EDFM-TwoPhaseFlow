$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$runGenericPath = Join-Path $repoRoot 'FIM_TransientEngine\RunGeneric_impl.hpp'
$wellMgrPath = Join-Path $repoRoot 'WellDOFManager.h'

if (!(Test-Path $runGenericPath)) {
    throw "[check_m4_s1] missing file: $runGenericPath"
}
if (!(Test-Path $wellMgrPath)) {
    throw "[check_m4_s1] missing file: $wellMgrPath"
}

$runGenericLines = Get-Content $runGenericPath
$wellMgrText = Get-Content $wellMgrPath -Raw

$start = -1
$end = -1
for ($i = 0; $i -lt $runGenericLines.Count; $i++) {
    if ($runGenericLines[$i] -match 'inline void RunGenericFIMTransientPressureOnlyN1\(') {
        $start = $i
        continue
    }
    if ($start -ge 0 -and $runGenericLines[$i] -match '^\s*template <int N, typename MeshMgrType, typename FieldMgrType>$') {
        $end = $i - 1
        break
    }
}

if ($start -lt 0) {
    throw '[check_m4_s1] failed to locate RunGenericFIMTransientPressureOnlyN1.'
}
if ($end -lt $start) {
    $end = $runGenericLines.Count - 1
}

$n1Region = ($runGenericLines[$start..$end] -join "`n")

$failures = New-Object System.Collections.Generic.List[string]

if ($n1Region.Contains('[N=1] unsupported: wells are not enabled in pressure-only AD route.')) {
    $failures.Add('N=1 route still hard-throws on non-empty wells.')
}

$requiredN1Patterns = @(
    'WellDOFManager<1> well_mgr;',
    'SelectActiveAndNormalizeWells(wells, mgr, 0.0, false)',
    'well_mgr.Setup(',
    'well_mgr.RegisterPatternConnections(global_mat);',
    'well_mgr.AssembleWellEquations(',
    'ApplyTrialUpdate<1>('
)

foreach ($pattern in $requiredN1Patterns) {
    if (-not $n1Region.Contains($pattern)) {
        $failures.Add("N=1 route is missing required well-support pattern: $pattern")
    }
}

$expectedWellMgrPattern = @'
if constexpr (N == 1) {
            if (grad_idx == 0) return 0; // P
            return -1;
        }
'@

if (-not $wellMgrText.Contains($expectedWellMgrPattern)) {
    $failures.Add('WellDOFManager::GradIndexToReservoirDOF is not specialized for N=1 pressure-only wells.')
}

if ($failures.Count -gt 0) {
    Write-Host '[check_m4_s1] FAIL'
    foreach ($failure in $failures) {
        Write-Host " - $failure"
    }
    exit 1
}

Write-Host '[check_m4_s1] PASS'
exit 0
