param(
    [ValidateSet('All', 'Compile', 'Run')]
    [string]$Mode = 'All',
    [string]$CaseDir = 'Test\Transient\FullCaseTest\H_T_CO2_ConstPP\h_t_co2_constpp_singlefrac_nowell',
    [ValidateSet('Javac', 'ComsolCompile')]
    [string]$Compiler = 'Javac',
    [ValidateSet('DirectJava', 'ComsolBatch')]
    [string]$Runtime = 'DirectJava',
    [switch]$SkipFineCheck
)

$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$workspaceRoot = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $scriptDir))
$javaSource = Join-Path $scriptDir 'ComsolHTCO2ConstPPSingleFracNoWell.java'
$javaClass = Join-Path $scriptDir 'ComsolHTCO2ConstPPSingleFracNoWell.class'
$javaSupportJar = Join-Path $scriptDir 'ComsolHTCO2ConstPPSingleFracNoWell-support.jar'

$comsolRoot = 'D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics'
$comsolJava = Join-Path $comsolRoot 'java\win64\jre\bin\java.exe'
$comsolJavac = Join-Path $comsolRoot 'java\win64\jre\bin\javac.exe'
$comsolCompile = 'D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\bin\win64\comsolcompile.exe'
$comsolBatch = 'D:\1Comsol6.3\Comosol6.3\COMSOL63\Multiphysics\bin\win64\comsolbatch.exe'
$comsolPathList = Join-Path $comsolRoot 'bin\comsolpath.txt'
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
    $families = @('matrix_horizontal', 'fracture_tangent', 'cross_normal')
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

function Get-ReferenceRepresentation {
    param([string]$SummaryPath)
    if (-not (Test-Path $SummaryPath)) {
        return $null
    }
    $match = Select-String -Path $SummaryPath -Pattern 'Fracture representation:\s*(.+)$' | Select-Object -First 1
    if ($null -eq $match) {
        return $null
    }
    return $match.Matches[0].Groups[1].Value.Trim()
}

function Archive-LegacySurrogateOutputs {
    param([string]$ReferenceDirPath)

    $summaryPath = Join-Path $ReferenceDirPath 'comsol_run_summary.md'
    $representation = Get-ReferenceRepresentation -SummaryPath $summaryPath
    if ([string]::IsNullOrWhiteSpace($representation) -or $representation -ne 'thin_band_fallback') {
        return
    }

    $timestamp = Get-Date -Format 'yyyyMMdd_HHmmss'
    $legacyRoot = Join-Path (Split-Path -Parent $ReferenceDirPath) 'comsol_surrogate_legacy'
    $archiveDir = Join-Path $legacyRoot $timestamp
    New-Item -ItemType Directory -Force -Path $archiveDir | Out-Null

    $entries = Get-ChildItem -Path $ReferenceDirPath -Force
    foreach ($entry in $entries) {
        Move-Item -Force -Path $entry.FullName -Destination (Join-Path $archiveDir $entry.Name)
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

    $pathValue = [System.Environment]::GetEnvironmentVariable('Path')
    [System.Environment]::SetEnvironmentVariable('PATH', $null, 'Process')
    if (-not [string]::IsNullOrWhiteSpace($pathValue)) {
        [System.Environment]::SetEnvironmentVariable('Path', $pathValue, 'Process')
    }

    $sanitizedArguments = @()
    foreach ($arg in $Arguments) {
        if ($null -ne $arg -and -not [string]::IsNullOrWhiteSpace([string]$arg)) {
            $sanitizedArguments += [string]$arg
        }
    }
    if ($sanitizedArguments.Count -eq 0) {
        throw "Refusing to launch '$FilePath' with an empty argument list."
    }

    $proc = Start-Process `
        -FilePath $FilePath `
        -ArgumentList $sanitizedArguments `
        -WorkingDirectory $workspaceRoot `
        -RedirectStandardOutput $StdoutPath `
        -RedirectStandardError $StderrPath `
        -NoNewWindow `
        -PassThru `
        -Wait

    return $proc.ExitCode
}

function Ensure-ComsolUserFolders {
    $userHome = [Environment]::GetFolderPath('UserProfile')
    if ([string]::IsNullOrWhiteSpace($userHome)) {
        throw 'Failed to resolve Windows user profile path for COMSOL runtime.'
    }

    $recoveries = Join-Path $userHome '.comsol\v63\comsol.recoveries'
    $configArea = Join-Path $userHome '.comsol\v63\configuration\comsol'
    $workspaceArea = Join-Path $userHome '.comsol\v63\workspace\comsol'

    try {
        New-Item -ItemType Directory -Force -Path $recoveries | Out-Null
        New-Item -ItemType Directory -Force -Path $configArea | Out-Null
        New-Item -ItemType Directory -Force -Path $workspaceArea | Out-Null
    } catch {
        throw "COMSOL runtime needs a writable Windows user profile folder. Failed to initialize '$userHome\\.comsol\\v63\\...'. $($_.Exception.Message)"
    }
}

function Resolve-ComsolRuntimeLayout {
    $homePath = Join-Path $workspaceRoot '.comsol_userprofile'
    New-Item -ItemType Directory -Force -Path $homePath | Out-Null
    $layout = [ordered]@{
        Home         = $homePath
        ConfigDir    = Join-Path $homePath '.comsol\v63\configuration\comsol_batch_validation'
        DataDir      = Join-Path $homePath '.comsol\v63\data\comsol_batch_validation'
        PrefsDir     = Join-Path $homePath '.comsol\v63\prefs\comsol_batch_validation'
        RecoveryDir  = Join-Path $homePath '.comsol\v63\recoveries\comsol_batch_validation'
        TempDir      = Join-Path $homePath '.comsol\v63\tmp\comsol_batch_validation'
    }

    foreach ($path in $layout.Values) {
        New-Item -ItemType Directory -Force -Path $path | Out-Null
    }

    return $layout
}

function Resolve-DirectJavaRuntimeRoot {
    $baseTemp = $env:TEMP
    if ([string]::IsNullOrWhiteSpace($baseTemp)) {
        $baseTemp = [System.IO.Path]::GetTempPath()
    }
    if ([string]::IsNullOrWhiteSpace($baseTemp)) {
        throw 'Failed to resolve a writable TEMP directory for COMSOL DirectJava runtime.'
    }

    $runtimeRoot = Join-Path $baseTemp 'codexcomsolhtsinglefrac'
    New-Item -ItemType Directory -Force -Path $runtimeRoot | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runtimeRoot 'recoveries') | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runtimeRoot 'prefs') | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $runtimeRoot 'tmp') | Out-Null
    return $runtimeRoot
}

function Build-JavaSupportJar {
    Remove-Item $javaSupportJar -ErrorAction SilentlyContinue
    $tempZip = [System.IO.Path]::ChangeExtension($javaSupportJar, '.zip')
    Remove-Item $tempZip -ErrorAction SilentlyContinue
    $classFiles = Get-ChildItem -Path $scriptDir -Filter '*.class' | Select-Object -ExpandProperty Name
    if (-not $classFiles -or $classFiles.Count -eq 0) {
        throw "Missing compiled class files under '$scriptDir'."
    }
    Push-Location $scriptDir
    try {
        Compress-Archive -Path $classFiles -DestinationPath $tempZip -Force
    }
    finally {
        Pop-Location
    }
    Move-Item -Force $tempZip $javaSupportJar
}

Push-Location $workspaceRoot
try {
    if ($Mode -eq 'All' -or $Mode -eq 'Compile') {
        Remove-Item (Join-Path $scriptDir '*.class') -ErrorAction SilentlyContinue
        Remove-Item $javaSupportJar -ErrorAction SilentlyContinue
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

        Build-JavaSupportJar
    }

    if ($Mode -eq 'All' -or $Mode -eq 'Run') {
        if (-not (Test-Path $javaSupportJar)) {
            Build-JavaSupportJar
        }
        Archive-LegacySurrogateOutputs -ReferenceDirPath $referenceDir
        Remove-Item (Join-Path $referenceDir 'comsol_java_stdout.log') -ErrorAction SilentlyContinue
        Remove-Item (Join-Path $referenceDir 'comsol_java_stderr.log') -ErrorAction SilentlyContinue
        Remove-Item (Join-Path $referenceDir 'comsol_java_error.log') -ErrorAction SilentlyContinue
        Remove-Item (Join-Path $referenceDir 'comsol_batch_stdout.log') -ErrorAction SilentlyContinue
        Remove-Item (Join-Path $referenceDir 'comsol_batch_stderr.log') -ErrorAction SilentlyContinue
        $runtimeLayout = Resolve-ComsolRuntimeLayout

        if ($Runtime -eq 'DirectJava') {
            Ensure-ComsolUserFolders

            $stdoutLog = Join-Path $referenceDir 'comsol_java_stdout.log'
            $stderrLog = Join-Path $referenceDir 'comsol_java_stderr.log'
            $classpath = Get-ComsolClasspath
            $directRuntimeRoot = Resolve-DirectJavaRuntimeRoot
            $directRecoveryDir = Join-Path $directRuntimeRoot 'recoveries'
            $directPrefsDir = Join-Path $directRuntimeRoot 'prefs'
            $directTmpDir = Join-Path $directRuntimeRoot 'tmp'
            $runArgs = @(
                "-Djava.library.path=$comsolBinWin64",
                "-Dcs.recoverydir=$directRecoveryDir",
                "-Dcs.prefsdir=$directPrefsDir",
                "-Dcs.tmpdir=$directTmpDir",
                "-Dcodex.comsol.runtimeRoot=$directRuntimeRoot",
                '-cp',
                $classpath,
                'ComsolHTCO2ConstPPSingleFracNoWell',
                $resolvedCaseDir
            )
            if ($SkipFineCheck) {
                $runArgs += '--skip-fine-check'
            }
            $exitCode = Invoke-LoggedProcess `
                -FilePath $comsolJava `
                -Arguments $runArgs `
                -StdoutPath $stdoutLog `
                -StderrPath $stderrLog
        } else {
            $stdoutLog = Join-Path $referenceDir 'comsol_batch_stdout.log'
            $stderrLog = Join-Path $referenceDir 'comsol_batch_stderr.log'
            $runArgs = @(
                '-autosave',
                'off',
                '-configuration',
                $runtimeLayout.ConfigDir,
                '-data',
                $runtimeLayout.DataDir,
                '-prefsdir',
                $runtimeLayout.PrefsDir,
                '-recoverydir',
                $runtimeLayout.RecoveryDir,
                '-tmpdir',
                $runtimeLayout.TempDir,
                '-classpathadd',
                $javaSupportJar,
                '-inputfile',
                $javaClass,
                '-dev',
                $javaSupportJar,
                '-batchlog',
                (Join-Path $referenceDir 'comsol_batch.log'),
                '-prodargs',
                $resolvedCaseDir
            )
            if ($SkipFineCheck) {
                $runArgs += '--skip-fine-check'
            }
            $exitCode = Invoke-LoggedProcess `
                -FilePath $comsolBatch `
                -Arguments $runArgs `
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
    Pop-Location
}
