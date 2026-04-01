$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$catalog = Join-Path $repoRoot 'CaseCommon_Catalog.cpp'
$main = Join-Path $repoRoot 'main.cpp'
$project = Join-Path $repoRoot '2D-Unstr-Quadrilateral-EDFM.vcxproj'
$filters = Join-Path $repoRoot '2D-Unstr-Quadrilateral-EDFM.vcxproj.filters'
$c1 = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.cpp'

$requiredFiles = @(
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.cpp',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.cpp',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.cpp',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.cpp',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.cpp'
)

foreach ($relativePath in $requiredFiles) {
    $fullPath = Join-Path $repoRoot $relativePath
    if (-not (Test-Path $fullPath)) {
        throw "M5-S3 missing required file: $relativePath"
    }
}

$catalogContent = Get-Content $catalog -Raw
$mainContent = Get-Content $main -Raw
$projectContent = Get-Content $project -Raw
$filtersContent = Get-Content $filters -Raw
$c1Content = Get-Content $c1 -Raw

$catalogPatterns = @(
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.h',
    'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.h',
    'int RunC2(CaseStage stage)',
    'int RunC3(CaseStage stage)',
    'int RunC4(CaseStage stage)',
    'int RunC5(CaseStage stage)',
    'int RunC6(CaseStage stage)',
    'if (caseCode == "C2") return BindingInfo{&RunC2, "implemented"};',
    'if (caseCode == "C3") return BindingInfo{&RunC3, "implemented"};',
    'if (caseCode == "C4") return BindingInfo{&RunC4, "implemented"};',
    'if (caseCode == "C5") return BindingInfo{&RunC5, "implemented"};',
    'if (caseCode == "C6") return BindingInfo{&RunC6, "implemented"};'
)

foreach ($pattern in $catalogPatterns) {
    if ($catalogContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S3 catalog is missing required pattern: $pattern"
    }
}

$mainPatterns = @(
    'test_h_tp_co2h2o_varypp_nofrac_nowell',
    'test_h_tp_co2h2o_constpp_singlefrac_nowell',
    'test_h_tp_co2h2o_varypp_singlefrac_nowell',
    'test_h_tp_co2h2o_constpp_complexfrac_nowell',
    'test_h_tp_co2h2o_varypp_complexfrac_nowell'
)

foreach ($pattern in $mainPatterns) {
    if ($mainContent -notmatch [regex]::Escape($pattern)) {
        throw "M5-S3 main entry list is missing required pattern: $pattern"
    }
}

foreach ($relativePath in $requiredFiles) {
    if ($projectContent -notmatch [regex]::Escape($relativePath)) {
        throw "M5-S3 project file is missing compile/include entry: $relativePath"
    }
    if ($filtersContent -notmatch [regex]::Escape($relativePath)) {
        throw "M5-S3 project filters are missing entry: $relativePath"
    }
}

$c1Patterns = @(
    'h_tp_co2h2o_varypp_nofrac_nowell',
    'h_tp_co2h2o_constpp_singlefrac_nowell',
    'h_tp_co2h2o_varypp_singlefrac_nowell',
    'h_tp_co2h2o_constpp_complexfrac_nowell',
    'h_tp_co2h2o_varypp_complexfrac_nowell',
    'MakeTwoPhaseWaterCO2EOS',
    'addFracture'
)

foreach ($pattern in $c1Patterns) {
    if ($c1Content -notmatch [regex]::Escape($pattern)) {
        throw "M5-S3 C1 donor is missing required pattern: $pattern"
    }
}

$templateChecks = @(
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_TP_CO2H2O_VaryPP_NoFrac',
            'h_tp_co2h2o_varypp_nofrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_TP_CO2H2O_ConstPP_SingleFrac',
            'h_tp_co2h2o_constpp_singlefrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_TP_CO2H2O_VaryPP_SingleFrac',
            'h_tp_co2h2o_varypp_singlefrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_TP_CO2H2O_ConstPP_ComplexFrac',
            'h_tp_co2h2o_constpp_complexfrac_nowell',
            'RunSolveOnly',
            'RunPrepareReference',
            'RunValidateOnly',
            'RunFullWorkflow'
        )
    },
    @{
        Path = Join-Path $repoRoot 'Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.cpp'
        Patterns = @(
            'namespace Test_H_TP_CO2H2O_VaryPP_ComplexFrac',
            'h_tp_co2h2o_varypp_complexfrac_nowell',
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
            throw "M5-S3 template $($check.Path) is missing required pattern: $pattern"
        }
    }
}

Write-Host 'PASS: M5-S3 C-family extension is wired for C2-C6.'
