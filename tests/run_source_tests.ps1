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
#
# quaids_curvature_bootstrap_test.e (Milestone 15) refits the whole
# quaidsFit()+quaidsCurvatureFit() pipeline B times per model, which adds
# roughly 45-50s to this script's runtime even at the small B values used
# there -- close to doubling this file's own ~30s baseline. Opt-in via
# -SkipBootstrap's ABSENCE is the default (i.e. it runs unless skipped),
# but the automatic push-triggered CI workflow (.github/workflows/tests.yml)
# passes -SkipBootstrap so routine pushes stay fast; run without the flag
# locally, or as part of release verification, to exercise it.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$SkipPubtable,
    [switch]$SkipCurvature,
    [switch]$SkipBootstrap
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
    "quaids_welfare_test.e",
    "quaids_reliability_regression_test.e"
)

if (-not $SkipPubtable) {
    $gaussTests += "quaids_pubtable_test.e"
}

if (-not $SkipCurvature) {
    $gaussTests += "quaids_curvature_test.e"
}

if (-not $SkipBootstrap) {
    $gaussTests += "quaids_curvature_bootstrap_test.e"
}

function Invoke-GaussBatch {
    param(
        [string]$Exe,
        [string[]]$Arguments
    )

    # Milestone 15 finding: reading stdout fully (ReadToEnd()) before
    # touching stderr is a classic .NET Process deadlock -- if the child
    # writes enough to BOTH streams to fill their OS pipe buffers before
    # either is drained, the child blocks mid-write while this script
    # blocks reading the other stream, and neither side ever proceeds.
    # quaids_curvature_bootstrap_test.e's QUAIDS block routinely prints
    # dozens to hundreds of "Optmt: function evaluation failed" lines to
    # stderr (a normal, expected side effect of optmt hitting a bad
    # bootstrap resample -- see that test's own header), enough volume to
    # hit exactly this deadlock; a run through this function hung for
    # hours where running the same file directly (tgauss -b -x ...) never
    # did, since a direct console run has no pipe buffer to fill. Fixed by
    # draining both streams asynchronously via events instead of
    # sequential ReadToEnd() calls.
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

    $proc = [System.Diagnostics.Process]::new()
    $proc.StartInfo = $psi

    $outputBuilder = [System.Text.StringBuilder]::new()
    $errorBuilder = [System.Text.StringBuilder]::new()

    $outputEvent = Register-ObjectEvent -InputObject $proc -EventName OutputDataReceived -Action {
        if ($null -ne $EventArgs.Data) { [void]$Event.MessageData.AppendLine($EventArgs.Data) }
    } -MessageData $outputBuilder
    $errorEvent = Register-ObjectEvent -InputObject $proc -EventName ErrorDataReceived -Action {
        if ($null -ne $EventArgs.Data) { [void]$Event.MessageData.AppendLine($EventArgs.Data) }
    } -MessageData $errorBuilder

    try {
        [void]$proc.Start()
        $proc.BeginOutputReadLine()
        $proc.BeginErrorReadLine()
        $proc.WaitForExit()
    } finally {
        Unregister-Event -SourceIdentifier $outputEvent.Name -ErrorAction SilentlyContinue
        Unregister-Event -SourceIdentifier $errorEvent.Name -ErrorAction SilentlyContinue
    }

    [pscustomobject]@{
        ExitCode = $proc.ExitCode
        Output = ($outputBuilder.ToString() + $errorBuilder.ToString())
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
