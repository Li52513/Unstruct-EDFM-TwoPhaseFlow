param(
    [ValidateSet('All', 'Compile', 'Run')]
    [string]$Mode = 'All',
    [string]$CaseDir = 'Test\Transient\FullCaseTest\H_CO2_VaryPP\h_co2_varypp_nofrac_nowell',
    [ValidateSet('Javac', 'ComsolCompile')]
    [string]$Compiler = 'Javac',
    [ValidateSet('DirectJava', 'ComsolBatch')]
    [string]$Runtime = 'DirectJava',
    [switch]$SkipFineCheck
)

$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$workspaceRoot = Split-Path -Parent (Split-Path -Parent (Split-Path -Parent $scriptDir))
$javaSource = Join-Path $scriptDir 'ComsolVaryPPNoFracNoWell.java'
$javaClass = Join-Path $scriptDir 'ComsolVaryPPNoFracNoWell.class'

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

Push-Location $workspaceRoot
try {
    if ($Mode -eq 'All' -or $Mode -eq 'Compile') {
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

        if ($Runtime -eq 'DirectJava') {
            Ensure-ComsolUserFolders

            $stdoutLog = Join-Path $referenceDir 'comsol_java_stdout.log'
            $stderrLog = Join-Path $referenceDir 'comsol_java_stderr.log'
            $classpath = Get-ComsolClasspath
            $exitCode = Invoke-LoggedProcess `
                -FilePath $comsolJava `
                -Arguments @(
                    "-Djava.library.path=$comsolBinWin64",
                    '-cp',
                    $classpath,
                    'ComsolVaryPPNoFracNoWell'
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
    }
}
finally {
    Remove-Item Env:COMSOL_CASE_DIR -ErrorAction SilentlyContinue
    Remove-Item Env:COMSOL_SKIP_FINE_CHECK -ErrorAction SilentlyContinue
    Pop-Location
}
