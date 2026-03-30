param(
    [ValidateSet('All', 'Compile', 'Run')]
    [string]$Mode = 'All',
    [string]$CaseDir = 'Test\Transient\FullCaseTest\H_T_CO2_ConstPP\h_t_co2_constpp_nofrac_nowell',
    [ValidateSet('Javac', 'ComsolCompile')]
    [string]$Compiler = 'Javac',
    [ValidateSet('DirectJava', 'ComsolBatch')]
    [string]$Runtime = 'DirectJava',
    [switch]$SkipFineCheck
)

$ErrorActionPreference = 'Stop'

function Normalize-ProcessPathVariable {
    $pathValue = [Environment]::GetEnvironmentVariable('Path', 'Process')
    if ([string]::IsNullOrWhiteSpace($pathValue)) {
        $pathValue = [Environment]::GetEnvironmentVariable('PATH', 'Process')
    }
    [Environment]::SetEnvironmentVariable('Path', $pathValue, 'Process')
    [Environment]::SetEnvironmentVariable('PATH', $null, 'Process')
}

Normalize-ProcessPathVariable

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$workspaceRoot = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $scriptDir))
$javaSource = Join-Path $scriptDir 'ComsolHTCO2ConstPPNoFracNoWell.java'
$javaClass = Join-Path $scriptDir 'ComsolHTCO2ConstPPNoFracNoWell.class'

$comsolRoot = 'D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics'
$comsolJava = Join-Path $comsolRoot 'java\win64\jre\bin\java.exe'
$comsolJavac = Join-Path $comsolRoot 'java\win64\jre\bin\javac.exe'
$comsolCompile = Join-Path $comsolRoot 'bin\win64\comsolcompile.exe'
$comsolBatch = Join-Path $comsolRoot 'bin\win64\comsolbatch.exe'
$comsolPluginDir = Join-Path $comsolRoot 'plugins'
$comsolBinWin64 = Join-Path $comsolRoot 'bin\win64'

$resolvedCaseDir = if ([System.IO.Path]::IsPathRooted($CaseDir)) {
    $CaseDir
} else {
    Join-Path $workspaceRoot $CaseDir
}
$resolvedCaseDir = (Resolve-Path $resolvedCaseDir).Path
$referenceDir = Join-Path $resolvedCaseDir 'reference\comsol'
New-Item -ItemType Directory -Force -Path $referenceDir | Out-Null

function Get-ComsolClasspath {
    $pluginWildcard = Join-Path $comsolPluginDir '*'
    return ($pluginWildcard + ';' + $scriptDir)
}

function Get-ProfileTags {
    param([string]$CaseDirPath)
    $schedulePath = Join-Path $CaseDirPath 'engineering\profile_report_schedule.csv'
    if (-not (Test-Path $schedulePath)) {
        throw "Missing profile schedule: $schedulePath"
    }
    return (Import-Csv -Path $schedulePath | ForEach-Object { $_.tag })
}

function Assert-ExpectedOutputs {
    param([string]$CaseDirPath)
    $tags = Get-ProfileTags -CaseDirPath $CaseDirPath
    $families = @('matrix_horizontal', 'matrix_vertical_midline')
    $required = @(
        (Join-Path $CaseDirPath 'reference\comsol\comsol_monitor_timeseries.csv'),
        (Join-Path $CaseDirPath 'reference\comsol\comsol_model.mph'),
        (Join-Path $CaseDirPath 'reference\comsol\comsol_progress.log'),
        (Join-Path $CaseDirPath 'reference\comsol\comsol_run_summary.md'),
        (Join-Path $CaseDirPath 'reference\comsol\comsol_reference_mesh_check.txt')
    )
    foreach ($family in $families) {
        foreach ($tag in $tags) {
            $required += Join-Path $CaseDirPath ("reference\\comsol\\comsol_profile_{0}_{1}.csv" -f $family, $tag)
        }
    }
    $missing = @($required | Where-Object { -not (Test-Path $_) })
    if ($missing.Count -gt 0) {
        throw ("COMSOL run finished but expected outputs are missing:`n" + ($missing -join "`n"))
    }
}

function Invoke-LoggedProcess {
    param(
        [string]$FilePath,
        [string[]]$Arguments,
        [string]$StdoutPath,
        [string]$StderrPath
    )

    Remove-Item $StdoutPath -ErrorAction SilentlyContinue
    Remove-Item $StderrPath -ErrorAction SilentlyContinue

    $proc = Start-Process `
        -FilePath $FilePath `
        -ArgumentList $Arguments `
        -WorkingDirectory $workspaceRoot `
        -RedirectStandardOutput $StdoutPath `
        -RedirectStandardError $StderrPath `
        -NoNewWindow `
        -PassThru `
        -Wait

    return $proc.ExitCode
}

function Invoke-LoggedDirectJava {
    param(
        [string]$FilePath,
        [string[]]$Arguments,
        [string]$StdoutPath,
        [string]$StderrPath
    )

    Remove-Item $StdoutPath -ErrorAction SilentlyContinue
    Remove-Item $StderrPath -ErrorAction SilentlyContinue

    $launcher = 'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe'
    $quotedArgs = $Arguments | ForEach-Object { "'" + ($_ -replace "'", "''") + "'" }
    $command = "& '$FilePath' $($quotedArgs -join ' ') 1> '$StdoutPath' 2> '$StderrPath'; exit `$LASTEXITCODE"
    $proc = Start-Process `
        -FilePath $launcher `
        -ArgumentList @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-Command', $command) `
        -WorkingDirectory $workspaceRoot `
        -NoNewWindow `
        -PassThru `
        -Wait
    return $proc.ExitCode
}

function Ensure-ComsolUserFolders {
    param(
        [Parameter(Mandatory = $true)][string]$RecoveryDir,
        [Parameter(Mandatory = $true)][string]$PrefsDir
    )

    $configArea = Join-Path $PrefsDir 'configuration\comsol'
    $workspaceArea = Join-Path $PrefsDir 'workspace\comsol'

    New-Item -ItemType Directory -Force -Path $RecoveryDir | Out-Null
    New-Item -ItemType Directory -Force -Path $configArea | Out-Null
    New-Item -ItemType Directory -Force -Path $workspaceArea | Out-Null
}

Push-Location $workspaceRoot
try {
    if ($Mode -eq 'All' -or $Mode -eq 'Compile') {
        Remove-Item (Join-Path $scriptDir '*.class') -ErrorAction SilentlyContinue
        if ($Compiler -eq 'Javac') {
            & $comsolJavac `
                -cp (Join-Path $comsolPluginDir 'com.comsol.api_1.0.0.jar') `
                $javaSource
        } else {
            & $comsolCompile $javaSource
        }

        if ($LASTEXITCODE -ne 0) {
            throw "COMSOL Java compile failed using compiler mode '$Compiler'."
        }
    }

    if ($Mode -eq 'All' -or $Mode -eq 'Run') {
        $env:COMSOL_CASE_DIR = $resolvedCaseDir
        if ($SkipFineCheck) {
            $env:COMSOL_SKIP_FINE_CHECK = '1'
        } else {
            Remove-Item Env:COMSOL_SKIP_FINE_CHECK -ErrorAction SilentlyContinue
        }
        Remove-Item (Join-Path $referenceDir 'comsol_java_error.log') -ErrorAction SilentlyContinue
        Remove-Item (Join-Path $referenceDir 'comsol_runtime') -Recurse -Force -ErrorAction SilentlyContinue
        $runtimeBase = Join-Path $env:TEMP ("codex_comsol_" + (Split-Path $resolvedCaseDir -Leaf))
        $runtimeRoot = Join-Path $runtimeBase 'runtime'
        $recoveryDir = Join-Path $runtimeRoot 'recoveries'
        $tempDir = Join-Path $runtimeRoot 'tmp'
        $prefsDir = Join-Path $runtimeBase 'prefs'
        New-Item -ItemType Directory -Force -Path $recoveryDir | Out-Null
        New-Item -ItemType Directory -Force -Path $tempDir | Out-Null
        New-Item -ItemType Directory -Force -Path $prefsDir | Out-Null

        if ($Runtime -eq 'DirectJava') {
            Ensure-ComsolUserFolders -RecoveryDir $recoveryDir -PrefsDir $prefsDir
            $stdoutLog = Join-Path $referenceDir 'comsol_java_stdout.log'
            $stderrLog = Join-Path $referenceDir 'comsol_java_stderr.log'
            $classpath = Get-ComsolClasspath
            $exitCode = Invoke-LoggedDirectJava `
                -FilePath $comsolJava `
                -Arguments @(
                    "-Djava.library.path=$comsolBinWin64",
                    "-Dcodex.comsol.runtimeRoot=$runtimeRoot",
                    "-Dcs.recoverydir=$recoveryDir",
                    "-Dcs.tmpdir=$tempDir",
                    "-Dcs.prefsdir=$prefsDir",
                    "-Djava.io.tmpdir=$tempDir",
                    '-cp',
                    $classpath,
                    'ComsolHTCO2ConstPPNoFracNoWell'
                ) `
                -StdoutPath $stdoutLog `
                -StderrPath $stderrLog
        } else {
            $stdoutLog = Join-Path $referenceDir 'comsol_batch_stdout.log'
            $stderrLog = Join-Path $referenceDir 'comsol_batch_stderr.log'
            $exitCode = Invoke-LoggedProcess `
                -FilePath $comsolBatch `
                -Arguments @(
                    '-inputfile',
                    $javaClass,
                    '-batchlog',
                    (Join-Path $referenceDir 'comsol_batch.log')
                ) `
                -StdoutPath $stdoutLog `
                -StderrPath $stderrLog
        }

        if ($exitCode -ne 0) {
            throw "COMSOL runtime failed using runtime mode '$Runtime'. See '$stderrLog' and '$stdoutLog'."
        }

        Assert-ExpectedOutputs -CaseDirPath $resolvedCaseDir
    }
}
finally {
    Remove-Item Env:COMSOL_CASE_DIR -ErrorAction SilentlyContinue
    Remove-Item Env:COMSOL_SKIP_FINE_CHECK -ErrorAction SilentlyContinue
    Pop-Location
}
