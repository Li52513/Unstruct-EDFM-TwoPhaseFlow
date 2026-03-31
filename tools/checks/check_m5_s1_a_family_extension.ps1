$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$catalog = Join-Path $repoRoot 'CaseCommon_Catalog.cpp'
$main = Join-Path $repoRoot 'main.cpp'

$requiredFiles = @(
    'Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.cpp',
    'Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.cpp'
)

foreach ($relativePath in $requiredFiles) {
    $fullPath = Join-Path $repoRoot $relativePath
    if (-not (Test-Path $fullPath)) {
        throw "M5-S1 missing required file: $relativePath"
    }
}

$catalogContent = Get-Content $catalog -Raw
$mainContent = Get-Content $main -Raw

$catalogPatterns = @(
    'Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.h',
    'int RunA5(CaseStage stage)',
    'int RunA6(CaseStage stage)',
    'if (caseCode == "A5") return BindingInfo{&RunA5, "implemented"};',
    'if (caseCode == "A6") return BindingInfo{&RunA6, "implemented"};'
)

foreach ($pattern in $catalogPatterns) {
    if ($catalogContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S1 catalog is missing required pattern: $pattern"
    }
}

$mainPatterns = @(
    'test_h_co2_constpp_complexfrac_nowell',
    'test_h_co2_varypp_complexfrac_nowell'
)

foreach ($pattern in $mainPatterns) {
    if ($mainContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S1 main entry list is missing required pattern: $pattern"
    }
}

$templateChecks = @(
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_CO2_ConstPP_ComplexFrac',
            'h_co2_constpp_complexfrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_CO2_VaryPP_ComplexFrac',
            'h_co2_varypp_complexfrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    }
)

foreach ($check in $templateChecks) {
    $content = Get-Content $check.Path -Raw
    foreach ($pattern in $check.Patterns) {
        if ($content -notmatch [regex]::Escape($pattern)) {
            throw "M5-S1 template $($check.Path) is missing required pattern: $pattern"
        }
    }
}

Write-Host 'PASS: M5-S1 A-family extension is wired for A5/A6.'
