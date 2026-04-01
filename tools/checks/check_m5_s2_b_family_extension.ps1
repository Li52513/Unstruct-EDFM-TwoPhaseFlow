$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$catalog = Join-Path $repoRoot 'CaseCommon_Catalog.cpp'
$main = Join-Path $repoRoot 'main.cpp'
$project = Join-Path $repoRoot '2D-Unstr-Quadrilateral-EDFM.vcxproj'
$filters = Join-Path $repoRoot '2D-Unstr-Quadrilateral-EDFM.vcxproj.filters'
$b1 = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp'
$b3 = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp'

$requiredFiles = @(
    'Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.cpp',
    'Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.cpp',
    'Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.cpp',
    'Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.cpp'
)

foreach ($relativePath in $requiredFiles) {
    $fullPath = Join-Path $repoRoot $relativePath
    if (-not (Test-Path $fullPath)) {
        throw "M5-S2 missing required file: $relativePath"
    }
}

$catalogContent = Get-Content $catalog -Raw
$mainContent = Get-Content $main -Raw
$projectContent = Get-Content $project -Raw
$filtersContent = Get-Content $filters -Raw
$b1Content = Get-Content $b1 -Raw
$b3Content = Get-Content $b3 -Raw

$catalogPatterns = @(
    'Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.h',
    'int RunB2(CaseStage stage)',
    'int RunB4(CaseStage stage)',
    'int RunB5(CaseStage stage)',
    'int RunB6(CaseStage stage)',
    'if (caseCode == "B2") return BindingInfo{&RunB2, "implemented"};',
    'if (caseCode == "B4") return BindingInfo{&RunB4, "implemented"};',
    'if (caseCode == "B5") return BindingInfo{&RunB5, "implemented"};',
    'if (caseCode == "B6") return BindingInfo{&RunB6, "implemented"};'
)

foreach ($pattern in $catalogPatterns) {
    if ($catalogContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S2 catalog is missing required pattern: $pattern"
    }
}

$mainPatterns = @(
    'test_h_t_co2_varypp_nofrac_nowell',
    'test_h_t_co2_varypp_singlefrac_nowell',
    'test_h_t_co2_constpp_complexfrac_nowell',
    'test_h_t_co2_varypp_complexfrac_nowell'
)

foreach ($pattern in $mainPatterns) {
    if ($mainContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S2 main entry list is missing required pattern: $pattern"
    }
}

foreach ($relativePath in $requiredFiles) {
    if ($projectContent -notmatch [regex]::Escape($relativePath)) {
        throw "M5-S2 project file is missing compile/include entry: $relativePath"
    }
    if ($filtersContent -notmatch [regex]::Escape($relativePath)) {
        throw "M5-S2 project filters are missing entry: $relativePath"
    }
}

$b1Patterns = @(
    'h_t_co2_varypp_nofrac_nowell',
    'MakeSinglePhaseCO2EOS'
)

foreach ($pattern in $b1Patterns) {
    if ($b1Content -notmatch [regex]::Escape($pattern)) {
        throw "M5-S2 B1 donor is missing required pattern: $pattern"
    }
}

$b3Patterns = @(
    'h_t_co2_varypp_singlefrac_nowell',
    'h_t_co2_constpp_complexfrac_nowell',
    'h_t_co2_varypp_complexfrac_nowell',
    'MakeSinglePhaseCO2EOS'
)

foreach ($pattern in $b3Patterns) {
    if ($b3Content -notmatch [regex]::Escape($pattern)) {
        throw "M5-S2 B3 donor is missing required pattern: $pattern"
    }
}

$templateChecks = @(
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_T_CO2_VaryPP_NoFrac',
            'h_t_co2_varypp_nofrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_T_CO2_VaryPP_SingleFrac',
            'h_t_co2_varypp_singlefrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_T_CO2_ConstPP_ComplexFrac',
            'h_t_co2_constpp_complexfrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_T_CO2_VaryPP_ComplexFrac',
            'h_t_co2_varypp_complexfrac_nowell',
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
            throw "M5-S2 template $($check.Path) is missing required pattern: $pattern"
        }
    }
}

Write-Host 'PASS: M5-S2 B-family extension is wired for B2-B6.'
