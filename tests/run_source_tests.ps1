# run_source_tests.ps1
#
# Milestone 7: runs the package-manifest consistency check, then every
# source-tree test (#include-based, not library-based) in tests/, in one
# shot. Adapted from gauss-qardl's tests/run_source_tests.ps1.
#
# This repo's tests print their own "PASS"/"FAIL" line per check and a
# final "...: ALL N CHECKS PASSED" (or "N CHECKS FAILED") summary line --
# CLAUDE.md documents that tgauss's process exit code is NOT a reliable
# pass/fail signal for this harness, so this runner checks the printed
# summary line (and any GAUSS-level compile/execute error) rather than
# relying on exit code alone.
#
# quaids_pubtable_test.e requires the pubtable package to be installed
# (this machine has it at c:\gauss26\pkgs\pubtable) -- pass -SkipPubtable
# to skip it on a machine without pubtable.
#
# quaids_curvature_test.e requires the optmt package to be installed
# (this machine has it at c:\gauss26\pkgs\optmt, and it is now a real
# package.json dependency -- see "Milestone 10" in CLAUDE.md) -- pass
# -SkipCurvature to skip it on a machine without optmt.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$SkipPubtable,
    [switch]$SkipCurvature
)

$testsDir = Join-Path $RepoRoot "tests"

& powershell -ExecutionPolicy Bypass -File (Join-Path $testsDir "verify_package_manifest.ps1") -RepoRoot $RepoRoot
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

$gaussTests = @(
    "quaids_schema_test.e",
    "quaids_formula_parity_test.e",
    "quaids_synthetic_validation_test.e",
    "quaids_published_validation_test.e",
    "quaids_hypothesis_tests_test.e",
    "quaids_elasticities_test.e",
    "quaids_welfare_test.e"
)

if (-not $SkipPubtable) {
    $gaussTests += "quaids_pubtable_test.e"
}

if (-not $SkipCurvature) {
    $gaussTests += "quaids_curvature_test.e"
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
    $psi.WorkingDirectory = $testsDir

    $proc = [System.Diagnostics.Process]::Start($psi)
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()

    [pscustomobject]@{
        ExitCode = $proc.ExitCode
        Output = ($stdout + $stderr)
    }
}

$failed = @()

foreach ($test in $gaussTests) {
    Write-Host ""
    Write-Host "==> $test"
    $result = Invoke-GaussBatch -Exe $GaussExe -Arguments @("-b", "-x", $test)
    $output = $result.Output
    $output

    $hasGaussError = $output -match "Program execute failed|error G[0-9]+|Program compile failed"
    $hasPassSummary = $output -match "ALL \d+ CHECKS PASSED"
    $hasFailSummary = $output -match "\d+ CHECKS FAILED"

    if ($hasGaussError -or $hasFailSummary -or -not $hasPassSummary) {
        $failed += $test
    }
}

if ($failed.Count -gt 0) {
    Write-Host ""
    Write-Host "run_source_tests.ps1: FAIL -- $($failed -join ', ')"
    exit 1
}

Write-Host ""
Write-Host "run_source_tests.ps1: PASS"
