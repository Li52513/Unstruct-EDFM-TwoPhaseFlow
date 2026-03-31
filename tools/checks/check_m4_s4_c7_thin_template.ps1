$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$target = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.cpp'

if (-not (Test-Path $target)) {
    throw "C7 template file is missing: $target"
}

$content = Get-Content $target -Raw

$requiredPatterns = @(
    'RunStageByKeyImpl',
    'RunGenericFIMTransient<3>',
    'BuildWellSchedule',
    'WriteC7WellTimeSeries',
    'Case2DReferenceIO::WriteAsciiFile'
)

foreach ($pattern in $requiredPatterns) {
    if ($content -notmatch [regex]::Escape($pattern)) {
        throw "C7 template is missing required pattern: $pattern"
    }
}

$forbiddenPatterns = @(
    'CaseCommon_Skeleton.h',
    'PrepareSkeletonArtifacts',
    'ThrowStageNotImplemented'
)

foreach ($pattern in $forbiddenPatterns) {
    if ($content -match [regex]::Escape($pattern)) {
        throw "C7 template still contains skeleton marker: $pattern"
    }
}

Write-Host 'PASS: C7 template is staged and no longer a skeleton.'
