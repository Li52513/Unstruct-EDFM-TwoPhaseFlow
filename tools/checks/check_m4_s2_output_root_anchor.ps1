$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$catalogPath = Join-Path $repoRoot 'CaseCommon_Catalog.cpp'

if (-not (Test-Path -LiteralPath $catalogPath)) {
    throw "Missing CaseCommon_Catalog.cpp at: $catalogPath"
}

$content = Get-Content -Path $catalogPath -Raw -Encoding UTF8

if ($content -match 'return\s+"Test/Transient/A1_F12_TemplateSystem"\s*;') {
    throw 'FAIL: BuildOutputRoot() is still hard-coded as a relative path.'
}

if ($content -notmatch '__FILE__') {
    throw 'FAIL: BuildOutputRoot() is not anchored to the compiled source tree.'
}

Write-Host 'PASS: BuildOutputRoot() is anchored to the worktree/repo root.'
