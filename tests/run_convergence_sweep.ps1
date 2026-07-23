# run_convergence_sweep.ps1
#
# Milestone 12: runs tests/quaids_convergence_sweep.e (the seed-sweep
# convergence-reliability diagnostic) and captures its full output to
# tests/convergence_sweep_report.txt so results are re-runnable and
# diffable across code changes -- addressing the fact that the original
# informal 8-seed probe behind CLAUDE.md's "Seed sensitivity" note never
# survived as a committed artifact anywhere in this repo's history.
#
# This is a DIAGNOSTIC runner, not a pass/fail gate -- unlike
# run_source_tests.ps1, it does not check for an "ALL N CHECKS PASSED"
# line (quaids_convergence_sweep.e deliberately never prints one) and is
# not invoked by run_source_tests.ps1. Run it manually whenever you want
# to re-measure the iterated estimator's convergence-failure rate (e.g.
# before/after a change to src/quaids.src's iteration loop).

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [string]$OutFile = (Join-Path $PSScriptRoot "convergence_sweep_report.txt")
)

$testsDir = Join-Path $RepoRoot "tests"

$psi = [System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName = $GaussExe
$psi.Arguments = "-b -x quaids_convergence_sweep.e"
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.WorkingDirectory = $testsDir

Write-Host "Running quaids_convergence_sweep.e -- this may take a minute..."
$proc = [System.Diagnostics.Process]::Start($psi)
$stdout = $proc.StandardOutput.ReadToEnd()
$stderr = $proc.StandardError.ReadToEnd()
$proc.WaitForExit()

$output = $stdout + $stderr
$output | Out-File -Encoding utf8 $OutFile

if ($output -match "Program execute failed|error G[0-9]+|Program compile failed") {
    Write-Host ""
    Write-Host "run_convergence_sweep.ps1: FAIL -- GAUSS reported an error, see $OutFile"
    exit 1
}

Write-Host ""
Write-Host "run_convergence_sweep.ps1: report written to $OutFile"
