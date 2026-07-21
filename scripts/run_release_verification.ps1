# run_release_verification.ps1
#
# Milestone 7: end-to-end release pipeline orchestrator -- source-tree
# tests, build (or verify an existing) release artifact, optionally install
# it into a GAUSS package directory, then run the installed-package public
# API gate (tests/package_public_api.e) against the INSTALLED copy via
# `library quaids;` (not the source tree via #include, which every other
# test in this repo uses). Adapted from gauss-qardl's
# scripts/run_release_verification.ps1, scaled down: this repo does not yet
# have gauss-qardl's separate new-model-benchmark/validation-benchmark/
# examples-smoke scripts -- the equivalent validation already lives inside
# run_source_tests.ps1's 7 tgauss test files (in particular
# quaids_synthetic_validation_test.e and quaids_published_validation_test.e).
#
# -InstallRoot defaults to <GaussHome>/pkgs, i.e. installing for real
# alongside every other package on the machine (c:\gauss26\pkgs\quaids) --
# only happens when -InstallArtifact is passed explicitly. Point
# -InstallRoot elsewhere (e.g. a scratch directory) to verify the pipeline
# without touching the shared GAUSS installation.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [string]$GaussHome = "",
    [string]$InstallRoot = "",
    [switch]$BuildArtifact,
    [switch]$ForceArtifact,
    [switch]$InstallArtifact,
    [switch]$SkipInstalledPackageTest
)

function Invoke-Step {
    param(
        [string]$Name,
        [scriptblock]$Command
    )

    Write-Host ""
    Write-Host "==> $Name"
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed"
    }
}

function Invoke-GaussBatch {
    param(
        [string]$Exe,
        [string[]]$Arguments
    )

    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $Exe
    $psi.Arguments = (($Arguments | ForEach-Object {
        if ($_ -match '[\s"]') {
            '"' + ($_ -replace '"', '\"') + '"'
        } else {
            $_
        }
    }) -join " ")
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true

    $proc = [System.Diagnostics.Process]::Start($psi)
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()

    [pscustomobject]@{
        ExitCode = $proc.ExitCode
        Output = ($stdout + $stderr)
    }
}

$testsDir = Join-Path $RepoRoot "tests"
$scriptsDir = Join-Path $RepoRoot "scripts"
$packagePath = Join-Path $RepoRoot "package.json"
$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$artifactPath = Join-Path $RepoRoot "$($pkg.name) $($pkg.version).zip"

if ([string]::IsNullOrWhiteSpace($GaussHome)) {
    $GaussHome = Split-Path -Parent $GaussExe
}

if ([string]::IsNullOrWhiteSpace($InstallRoot)) {
    $InstallRoot = Join-Path $GaussHome "pkgs"
}

Invoke-Step "Source tests" {
    & powershell -ExecutionPolicy Bypass -File (Join-Path $testsDir "run_source_tests.ps1") -RepoRoot $RepoRoot -GaussExe $GaussExe
}

if ($BuildArtifact) {
    Invoke-Step "Build release artifact" {
        $buildArgs = @("-ExecutionPolicy", "Bypass", "-File", (Join-Path $scriptsDir "build_package.ps1"), "-RepoRoot", $RepoRoot)
        if ($ForceArtifact) {
            $buildArgs += "-Force"
        }
        & powershell @buildArgs
    }
} else {
    Invoke-Step "Verify release artifact" {
        & powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "verify_release_artifact.ps1") -RepoRoot $RepoRoot -ArtifactPath $artifactPath
    }
}

if ($InstallArtifact) {
    Invoke-Step "Install release artifact into GAUSS package directory" {
        $installDir = Join-Path $InstallRoot $pkg.name
        if (Test-Path -LiteralPath $installDir) {
            Remove-Item -LiteralPath $installDir -Recurse -Force
        }
        New-Item -ItemType Directory -Path $installDir | Out-Null
        Expand-Archive -LiteralPath $artifactPath -DestinationPath $installDir -Force
        & powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "build_lcg.ps1") -PackageRoot $installDir -PackageName $pkg.name

        $installedPackagePath = Join-Path $installDir "package.json"
        if (-not (Test-Path -LiteralPath $installedPackagePath)) {
            throw "installed artifact did not create package.json at $installedPackagePath"
        }

        $installedPkg = Get-Content -LiteralPath $installedPackagePath -Raw | ConvertFrom-Json
        if ([string]$installedPkg.version -ne [string]$pkg.version) {
            throw "installed package version '$($installedPkg.version)' does not match release version '$($pkg.version)'"
        }
    }
}

if (-not $SkipInstalledPackageTest) {
    Invoke-Step "Installed-package public API test" {
        $wrapper = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_pkg_" + [System.Guid]::NewGuid().ToString("N") + ".e")
        $gaussTestsDir = $testsDir -replace "\\", "/"
        Set-Content -Path $wrapper -Value @(
            "new;",
            "chdir `"$gaussTestsDir`";",
            "run `"$gaussTestsDir/package_public_api.e`";"
        )
        # `library quaids;` inside package_public_api.e resolves against
        # GAUSS's own configured package directory (normally
        # <GaussHome>/pkgs), not -InstallRoot -- this step is only a real
        # end-to-end check when -InstallRoot matches that directory (the
        # default). A non-default -InstallRoot still exercises
        # build/install staging correctly; it just means this specific step
        # will not find the package unless GAUSS is separately configured
        # to look there.

        try {
            $result = Invoke-GaussBatch -Exe $GaussExe -Arguments @("-nb", "-b", "-x", $wrapper)
            $result.Output
            if ($result.ExitCode -ne 0 -or ($result.Output -match "Program execute failed|error G[0-9]+|Program compile failed")) {
                exit 1
            }
        } finally {
            Remove-Item -LiteralPath $wrapper -ErrorAction SilentlyContinue
        }
    }
}

Write-Host ""
Write-Host "run_release_verification.ps1: PASS"
