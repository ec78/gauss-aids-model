# run_examples_smoke.ps1
#
# Public release roadmap PR-405: smoke-tests every example in examples/
# by actually running it, checking for a GAUSS compile/execute error
# (examples have no "ALL N CHECKS PASSED" assertions of their own --
# Milestone 31's own design is manual/eyeball-comparison demonstration
# scripts, not a test suite -- so "did it run to completion without
# erroring" is the whole signal, matching gauss-qardl's own
# tests/run_examples_smoke.ps1, the proven pattern this file adapts).
#
# Location independence (PR-405's "make includes resolve relative to the
# program/artifact rather than the caller's current working directory"):
# empirically confirmed (not assumed) that GAUSS's bare `#include
# filename` resolves via gauss.cfg's src_path wildcard search
# ($(PACKAGEDIR)\*\src;$(PACKAGEDIR)\*\examples across every installed
# package) regardless of the process's actual working directory, as long
# as the referenced package (here, quaids) is installed -- confirmed by
# running an unmodified example from a completely unrelated directory,
# both with a relative filename (cwd set via this script's own
# WorkingDirectory) and via an absolute script path from elsewhere. This
# is why examples/10_curvature_imposition.e and
# examples/13_pubtable_reporting.e were changed from `#include
# ../src/quaidscurvature.src`/`#include ../src/pubtable_quaids.src`
# (a source-tree-relative path, which only resolves when cwd happens to
# already be examples/) to bare `#include quaidscurvature.src`/`#include
# pubtable_quaids.src` (resolves via the installed package, from
# anywhere). This script itself still sets each child process's own
# working directory to examples/ (via ProcessStartInfo.WorkingDirectory,
# not a GAUSS-level `chdir` -- confirmed empirically that a GAUSS `chdir`
# statement placed before an `#include` in the same file has NO effect on
# that #include's resolution, since GAUSS resolves every #include in a
# full compile pass before any runtime statement -- including an earlier
# `chdir` -- ever executes), so this script itself can be invoked from
# any directory and still correctly test every example, including
# examples/00_real_data_quickstart.e's own loadd() call, which -- unlike
# #include -- has no package-search fallback and genuinely needs
# examples/ as the process cwd (a real, confirmed GAUSS/OS file-I/O
# limitation, not a bug in that example -- see its own header comment).
#
# Skips mirror tests/run_source_tests.ps1's existing -SkipCurvature/
# -SkipPubtable convention for the two optional-package examples.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$SkipCurvature,
    [switch]$SkipPubtable
)

$examplesDir = Join-Path $RepoRoot "examples"

$examples = @(
    "00_real_data_quickstart.e",
    "01_basic_estimation.e",
    "02_dataframe_input.e",
    "03_preflight_diagnostics.e",
    "04_hypothesis_tests.e",
    "05_elasticities_shares_slutzky.e",
    "06_welfare_analysis.e",
    "07_zero_share_correction.e",
    "08_robust_standard_errors.e",
    "09_replicate_weights.e"
)

if (-not $SkipCurvature) {
    $examples += "10_curvature_imposition.e"
}

$examples += "11_survey_weighted_estimation.e"
$examples += "12_applied_workflow.e"

if (-not $SkipPubtable) {
    $examples += "13_pubtable_reporting.e"
}

# Same async-stream-drain fix as tests/run_source_tests.ps1's own
# Invoke-GaussBatch (Milestone 15 finding): sequential ReadToEnd() calls
# on stdout then stderr can deadlock once a child writes enough to both
# streams to fill their OS pipe buffers before either is drained --
# 10_curvature_imposition.e's QUAIDS curvature bootstrap block routinely
# prints many "Optmt: function evaluation failed" lines to stderr, real
# volume in the same class that caused this deadlock originally.
function Invoke-GaussBatch {
    param(
        [string]$Exe,
        [string[]]$Arguments,
        [string]$WorkingDirectory
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
    $psi.WorkingDirectory = $WorkingDirectory

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

# Every file a running example can generate in examples/ -- matches
# scripts/build_package.ps1's own $generatedTestFiles list and
# scripts/verify_release_artifact.ps1's own forbidden-entries pattern
# (kept in sync by hand across all three; see that script's own comment).
$generatedFiles = @(
    "blanciforti_results.txt",
    "quaids_coefficients.tex",
    "quaids_coefficients.md",
    "quaids_coefficients.csv",
    "quaids_income_elasticities.md",
    "quaids_uncompensated_elasticities.tex",
    "quaids_compensated_elasticities.csv",
    "quaids_workflow"
)

function Remove-GeneratedExampleFiles {
    foreach ($name in $generatedFiles) {
        $path = Join-Path $examplesDir $name
        if (Test-Path -LiteralPath $path) {
            Remove-Item -LiteralPath $path -Force -ErrorAction SilentlyContinue
        }
    }
}

$failed = @()

foreach ($example in $examples) {
    Write-Host ""
    Write-Host "==> $example"
    $result = Invoke-GaussBatch -Exe $GaussExe -Arguments @("-b", "-x", $example) -WorkingDirectory $examplesDir
    $output = $result.Output
    $output

    $hasGaussError = $output -match "Program execute failed|error G[0-9]+|Program compile failed"
    if ($hasGaussError) {
        Write-Host "FAIL  $example"
        $failed += $example
    } else {
        Write-Host "PASS  $example"
    }

    Remove-GeneratedExampleFiles
}

if ($failed.Count -gt 0) {
    Write-Host ""
    Write-Host "run_examples_smoke.ps1: FAIL -- $($failed -join ', ')"
    exit 1
}

Write-Host ""
Write-Host "run_examples_smoke.ps1: PASS ($($examples.Count) examples)"
