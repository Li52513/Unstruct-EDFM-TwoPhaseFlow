$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$target = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.cpp'

if (-not (Test-Path $target)) {
    throw "B7 template file is missing: $target"
}

$content = Get-Content $target -Raw

$requiredPatterns = @(
    'RunStageByKeyImpl',
    'RunGenericFIMTransient<2>',
    'BuildWellSchedule',
    'WriteB7WellTimeSeries',
    'Case2DReferenceIO::WriteAsciiFile'
)

foreach ($pattern in $requiredPatterns) {
    if ($content -notmatch [regex]::Escape($pattern)) {
        throw "B7 template is missing required pattern: $pattern"
    }
}

$forbiddenPatterns = @(
    'CaseCommon_Skeleton.h',
    'PrepareSkeletonArtifacts',
    'ThrowStageNotImplemented'
)

foreach ($pattern in $forbiddenPatterns) {
    if ($content -match [regex]::Escape($pattern)) {
        throw "B7 template still contains skeleton marker: $pattern"
    }
}

Write-Host 'PASS: B7 template is staged and no longer a skeleton.'
