param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [switch]$IncludeNoisyPatterns
)

function ShouldSkipDiagnosticContext {
    param(
        [string]$RuleName,
        [string]$Line
    )

    if ($RuleName -notin @("VizFieldLiteral", "LegacyAliasLiteral")) {
        return $false
    }

    return $Line -match '(?:\b(?:const_)?snap\s*\[)|(?:std::c(?:out|err))|(?:nlohmann::json)'
}

$targetRelPaths = @(
    "TransmissibilitySolver_2D.cpp",
    "TransmissibilitySolver_3D.cpp",
    "FIM_TopologyBuilder2D.h",
    "FIM_TopologyBuilder3D.h",
    "BoundaryAssembler.cpp",
    "FIM_TransientEngine\\StateSync.hpp",
    "FIM_TransientEngine\\RunGeneric_impl.hpp"
)

$patterns = @(
    @{ Name = "TransientOutputRoot"; Regex = '"Test/Transient/Day6"' },
    @{ Name = "TransmissibilityFieldLiteral"; Regex = '"T_(Matrix|FI|NNC|FF)_(Flow|Heat)"' },
    @{ Name = "WellAuxFieldLiteral"; Regex = '"mob_density_[wg]"' },
    @{ Name = "NonOrthTempFieldLiteral"; Regex = '"n23_(probe|main)_(p|t|sw)_tmp"' }
)

if ($IncludeNoisyPatterns) {
    $patterns += @(
        @{ Name = "VizFieldLiteral"; Regex = '"(P|S_w)"' },
        @{ Name = "LegacyAliasLiteral"; Regex = '"(Pressure|Temperature|T_res|Saturation|Sw|s_w|sw|p_w|p)"' }
    )
}

$issues = New-Object System.Collections.Generic.List[object]

foreach ($relPath in $targetRelPaths) {
    $path = Join-Path $RepoRoot $relPath
    if (-not (Test-Path $path)) {
        Write-Warning "Skip missing file: $path"
        continue
    }

    $lineNo = 0
    Get-Content $path | ForEach-Object {
        $lineNo++
        $line = $_
        foreach ($pattern in $patterns) {
            if ($line -match $pattern.Regex) {
                if (ShouldSkipDiagnosticContext -RuleName $pattern.Name -Line $line) {
                    continue
                }
                $issues.Add([pscustomobject]@{
                    File = $relPath
                    Line = $lineNo
                    Rule = $pattern.Name
                    Text = $line.Trim()
                })
            }
        }
    }
}

if ($issues.Count -eq 0) {
    Write-Host "[PASS] No hardcoded field/path literals found in target modules."
    exit 0
}

$issues |
    Sort-Object File, Line, Rule |
    Format-Table File, Line, Rule, Text -AutoSize

Write-Host ""
Write-Host ("[FAIL] Found {0} potential hardcoded literal hit(s)." -f $issues.Count)
Write-Host "If you need legacy alias/viz scans too, rerun with -IncludeNoisyPatterns."
exit 1
