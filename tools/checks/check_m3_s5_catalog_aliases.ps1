$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$mainPath = Join-Path $repoRoot "main.cpp"
$failures = New-Object System.Collections.Generic.List[string]

if (-not (Test-Path $mainPath)) {
    Write-Host "[M3-S5][FAIL] main.cpp is missing." -ForegroundColor Red
    exit 1
}

$content = Get-Content $mainPath -Raw

function Require-Contains {
    param(
        [string]$Needle,
        [string]$Message
    )

    if (-not $content.Contains($Needle)) {
        $failures.Add($Message)
    }
}

function Require-NotContains {
    param(
        [string]$Needle,
        [string]$Message
    )

    if ($content.Contains($Needle)) {
        $failures.Add($Message)
    }
}

Require-Contains "struct CatalogAliasEntry" "main.cpp must define a catalog-alias entry struct for legacy primary case names."
Require-Contains "legacyCatalogAliases" "main.cpp must maintain a legacy catalog-alias table."
Require-Contains '"test_h_co2_constpp_nofrac_nowell"' "A1 legacy primary alias must be present."
Require-Contains '"A1"' "A1 catalog target must be present in the alias mapping."
Require-Contains '"test_h_t_co2_constpp_nofrac_nowell"' "B1 legacy primary alias must be present."
Require-Contains '"B1"' "B1 catalog target must be present in the alias mapping."
Require-Contains '"test_h_tp_co2h2o_constpp_nofrac_nowell"' "C1 legacy alias must be present."
Require-Contains '"C1"' "C1 catalog target must be present in the alias mapping."
Require-Contains "PrintCatalogAliases" "main.cpp must print the catalog alias section for discoverability."
Require-Contains "ResolveCatalogCaseAlias" "main.cpp must resolve legacy aliases before catalog dispatch."

Require-NotContains '{"test_h_co2_constpp_nofrac_nowell", "Standalone test: 2D single-phase CO2 const-property no-fracture no-well (independent template)", []() { Test_H_CO2_ConstPP::RunFullWorkflow(); return 0; }}' "A1 legacy primary entry must not directly invoke the template anymore."
Require-NotContains '{"test_h_t_co2_constpp_nofrac_nowell", "Standalone test: 2D single-phase CO2 const-property P-T coupled no-fracture no-well", []() { Test_H_T_CO2_ConstPP_NoFrac::RunFullWorkflow(); return 0; }}' "B1 legacy primary entry must not directly invoke the template anymore."

if ($failures.Count -gt 0) {
    Write-Host "[M3-S5][FAIL] Legacy primary entry replacement is incomplete." -ForegroundColor Red
    foreach ($failure in $failures) {
        Write-Host " - $failure"
    }
    exit 1
}

Write-Host "[M3-S5][PASS] Legacy primary no-well entry points route through catalog aliases." -ForegroundColor Green
exit 0
