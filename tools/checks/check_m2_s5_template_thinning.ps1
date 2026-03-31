$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)

$checks = @(
    @{
        Path = Join-Path $repoRoot "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.cpp"
        Patterns = @(
            "double ClampFraction\(",
            "std::string MakeReportTag\(",
            "std::vector<double> BuildSortedFractions\("
        )
    },
    @{
        Path = Join-Path $repoRoot "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.cpp"
        Patterns = @(
            "std::vector<std::string> SplitCsvLine\(",
            "CsvTable ReadCsvTable\(",
            "double CsvGetDouble\(",
            "void WriteValidationSummary\(",
            "void WriteMatlabPlotScript\("
        )
    },
    @{
        Path = Join-Path $repoRoot "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.cpp"
        Patterns = @(
            "std::vector<std::string> SplitCsvLine\(",
            "struct CsvTable\s*\{",
            "CsvTable ReadCsvTable\(",
            "double CsvGetDouble\(",
            "double ClampFraction\(",
            "std::string MakeReportTag\(",
            "std::vector<double> BuildSortedFractions\("
        )
    }
)

$violations = @()

foreach ($check in $checks) {
    $content = Get-Content $check.Path -Raw
    foreach ($pattern in $check.Patterns) {
        if ($content -match $pattern) {
            $violations += [PSCustomObject]@{
                Path = $check.Path
                Pattern = $pattern
            }
        }
    }
}

if ($violations.Count -gt 0) {
    Write-Host "M2-S5 template-thinning contract FAILED." -ForegroundColor Red
    foreach ($violation in $violations) {
        Write-Host ("  {0} :: {1}" -f $violation.Path, $violation.Pattern)
    }
    exit 1
}

Write-Host "M2-S5 template-thinning contract PASSED." -ForegroundColor Green
exit 0
