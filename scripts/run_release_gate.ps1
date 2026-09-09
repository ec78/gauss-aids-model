# run_release_gate.ps1
#
# PR-601: a single go/no-go command for a release candidate. Thin
# orchestration, not a reimplementation -- scripts/run_release_verification.ps1
# (called with -BuildArtifact -ForceArtifact -InstallArtifact and no
# -SkipInstalledPackageTest) already chains, in order: the full manifest/
# public-API/documentation-consistency/documentation-quality checks
# (inside tests/run_source_tests.ps1's own first four steps), the full
# source-tree test suite INCLUDING bootstrap tests (no -SkipBootstrap is
# passed anywhere in this file), the release artifact build and its own
# self-verification (name/version/CHANGELOG entry/required entries/
# forbidden generated artifacts/archive-level link check --
# scripts/verify_release_artifact.ps1), a clean install into a real GAUSS
# package directory, the installed-package public API test against that
# exact artifact, and all 14 example smoke tests against that exact
# install. That one call covers PR-601's Work-list items 2 through 8.
# This script adds only the pieces nothing else does: item 1 (clean
# worktree), item 9 (the numerical convergence sweep, captured not
# gated -- see below for why), and item 10 (artifact checksum and a
# structured release record).
#
# Exits nonzero on the FIRST failed gate (Invoke-GateStep re-throws),
# matching the acceptance evidence -- this is a real gate, not a
# best-effort report.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [string]$GaussHome = "",
    [string]$InstallRoot = "",
    [switch]$AllowDirtyWorktree,
    [switch]$SkipConvergenceSweep
)

$ErrorActionPreference = "Stop"

function Invoke-GateStep {
    param(
        [string]$Name,
        [scriptblock]$Command
    )

    Write-Host ""
    Write-Host "==> [$Name]"
    $stepStart = Get-Date
    try {
        & $Command
        if ($LASTEXITCODE -ne 0) {
            throw "$Name failed (exit code $LASTEXITCODE)"
        }
        $status = "pass"
    } catch {
        $status = "fail"
        $script:record.steps += [pscustomobject]@{
            name = $Name
            status = $status
            durationSeconds = [math]::Round(((Get-Date) - $stepStart).TotalSeconds, 1)
        }
        Write-Host ""
        Write-Host "run_release_gate.ps1: NO-GO -- '$Name' failed: $($_.Exception.Message)"
        Save-ReleaseRecord -Status "no-go"
        exit 1
    }
    $script:record.steps += [pscustomobject]@{
        name = $Name
        status = $status
        durationSeconds = [math]::Round(((Get-Date) - $stepStart).TotalSeconds, 1)
    }
}

function Save-ReleaseRecord {
    param([string]$Status)

    $script:record.status = $Status
    $script:record.finishedAt = (Get-Date).ToString("o")
    $script:record.totalDurationSeconds = [math]::Round(((Get-Date) - $script:gateStart).TotalSeconds, 1)

    $recordsDir = Join-Path $RepoRoot "release_records"
    if (-not (Test-Path -LiteralPath $recordsDir)) {
        New-Item -ItemType Directory -Path $recordsDir | Out-Null
    }
    $stamp = (Get-Date).ToString("yyyyMMdd-HHmmss")
    $recordPath = Join-Path $recordsDir "$($script:record.version)-$stamp.json"
    $script:record | ConvertTo-Json -Depth 6 | Out-File -Encoding utf8 $recordPath
    Write-Host ""
    Write-Host "Release record written to $recordPath"
}

$gateStart = Get-Date
$packagePath = Join-Path $RepoRoot "package.json"
$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$version = [string]$pkg.version
$packageName = [string]$pkg.name

$record = [pscustomobject]@{
    package = $packageName
    version = $version
    startedAt = $gateStart.ToString("o")
    finishedAt = $null
    totalDurationSeconds = $null
    status = "in-progress"
    sourceCommit = $null
    gitWorktreeClean = $null
    toolVersions = [pscustomobject]@{
        gauss = $null
        powershell = $PSVersionTable.PSVersion.ToString()
        os = [System.Environment]::OSVersion.VersionString
    }
    steps = @()
    convergenceSweep = $null
    artifact = $null
}

# --- Item 1: clean worktree and synchronized release metadata ---
# (metadata synchronization itself -- package.json/CITATION.cff/
# docs/public-api.json/CHANGELOG.md agreeing on $version -- is verified
# by scripts/verify_public_api.ps1, run as part of the pipeline below;
# this step covers the other half: no uncommitted changes sitting in the
# worktree that wouldn't be part of the tagged commit.)

Write-Host "==> [Clean worktree]"
$gitStatus = & git -C $RepoRoot status --porcelain
$record.gitWorktreeClean = [string]::IsNullOrWhiteSpace($gitStatus)
if (-not $record.gitWorktreeClean -and -not $AllowDirtyWorktree) {
    Write-Host "run_release_gate.ps1: NO-GO -- worktree is not clean:"
    Write-Host $gitStatus
    Write-Host "Commit or stash these changes first, or pass -AllowDirtyWorktree to override (not recommended for a real release)."
    Save-ReleaseRecord -Status "no-go"
    exit 1
}
$record.sourceCommit = (& git -C $RepoRoot rev-parse HEAD).Trim()
Write-Host "worktree clean: $($record.gitWorktreeClean); source commit: $($record.sourceCommit)"

# --- GAUSS version (for the record) ---
# Uses a real Process object rather than `& $GaussExe ... 2>&1`, since
# with $ErrorActionPreference = "Stop" (set above), PowerShell treats ANY
# stderr line from a native executable as a terminating error, even on a
# genuinely successful run -- confirmed directly (this exact call failed
# the gate's own first real run for exactly this reason before being
# fixed). The Process object reads stdout/stderr as plain data, not
# PowerShell error records.
$versionProbe = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_gate_ver_" + [System.Guid]::NewGuid().ToString("N") + ".e")
Set-Content -Path $versionProbe -Value @("new;", "print `"probe complete`";")
try {
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $GaussExe
    $psi.Arguments = "-b -x `"$versionProbe`""
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $verProc = [System.Diagnostics.Process]::Start($psi)
    $verOutput = $verProc.StandardOutput.ReadToEnd() + $verProc.StandardError.ReadToEnd()
    $verProc.WaitForExit()
    if ($verOutput -match "(GAUSS [\d.]+ \([^)]+\) \d+-bit)") {
        $record.toolVersions.gauss = $Matches[1]
    }
} finally {
    Remove-Item -LiteralPath $versionProbe -ErrorAction SilentlyContinue
}

# --- Items 2-8: manifest/API/docs checks, full source suite (including
# bootstrap), example smoke tests, artifact build+verify, clean install,
# installed-package public API test -- all via the existing pipeline. ---

# Hashtable splatting does not survive crossing a `powershell -File`
# subprocess boundary (confirmed directly -- it failed the gate's own
# first full run: -BuildArtifact arrived at the child as the literal
# string "True" instead of a switch, which PowerShell's own parameter
# binder then rejected). Build a plain string argument array instead.
$verificationArgList = @(
    "-ExecutionPolicy", "Bypass", "-File", (Join-Path $RepoRoot "scripts\run_release_verification.ps1"),
    "-RepoRoot", $RepoRoot,
    "-GaussExe", $GaussExe,
    "-BuildArtifact", "-ForceArtifact", "-InstallArtifact"
)
if (-not [string]::IsNullOrWhiteSpace($GaussHome)) { $verificationArgList += @("-GaussHome", $GaussHome) }
if (-not [string]::IsNullOrWhiteSpace($InstallRoot)) { $verificationArgList += @("-InstallRoot", $InstallRoot) }

Invoke-GateStep "Full release verification (manifest, API inventory, docs, full source suite with bootstrap, artifact build/verify, clean install, installed-package test, example smoke tests)" {
    & powershell @verificationArgList
}

# --- Item 9: numerical benchmark/sweep results ---
#
# Deliberately CAPTURED, not GATED: tests/quaids_convergence_sweep.e is,
# by this project's own established and documented design, a diagnostic
# report generator with no pass/fail threshold -- there is no known
# convergence guarantee for the iterated estimator to gate on, and the
# roadmap's own Production-Readiness (not alpha) Exit Criteria is what
# introduces a predeclared success threshold. For the public-alpha gate,
# PR-601's own acceptance evidence only requires the release record to
# CAPTURE benchmark results, not enforce a numeric bar -- so a bad sweep
# result here is recorded, visible in the release record, and worth a
# human's attention before tagging, but does not by itself flip this
# script's own go/no-go exit code the way every other step does.
if (-not $SkipConvergenceSweep) {
    Write-Host ""
    Write-Host "==> [Convergence sweep (captured, not gated -- see comment above)]"
    $sweepReport = Join-Path $RepoRoot "tests\convergence_sweep_report.txt"
    & powershell -ExecutionPolicy Bypass -File (Join-Path $RepoRoot "tests\run_convergence_sweep.ps1") -RepoRoot $RepoRoot -GaussExe $GaussExe -OutFile $sweepReport
    if (Test-Path -LiteralPath $sweepReport) {
        $sweepText = Get-Content -LiteralPath $sweepReport -Raw
        # Real captured format, confirmed by reading an actual report
        # (GAUSS's own print padding, not assumed): "SUMMARY:Iterated AIDS
        # (linear)  (       200.00000 seeds, tobs=       3000.0000 )" --
        # no space after "SUMMARY:", the model name itself can contain
        # spaces/parens ("Iterated AIDS (linear)"), so the model-name
        # capture stops at the first run of 2+ spaces before the seed-
        # count "(", not at the first "(" (which could be inside the name).
        $summaryBlocks = [regex]::Matches($sweepText, "SUMMARY:(.+?)\s{2,}\(\s*(\d+)(?:\.\d+)?\s*seeds,\s*tobs=\s*(\d+)(?:\.\d+)?\s*\)[\s\S]*?never-converged:\s*[\d.]+\s*\(\s*([\d.]+)\s*%\)[\s\S]*?converged-but-wrong:\s*[\d.]+\s*\(\s*([\d.]+)\s*%\)[\s\S]*?converged-correctly:\s*[\d.]+\s*\(\s*([\d.]+)\s*%\)")
        $sweepSummary = @()
        foreach ($m in $summaryBlocks) {
            $sweepSummary += [pscustomobject]@{
                model = $m.Groups[1].Value.Trim()
                seeds = [int]$m.Groups[2].Value
                tobs = [int]$m.Groups[3].Value
                neverConvergedPct = [double]$m.Groups[4].Value
                convergedWrongPct = [double]$m.Groups[5].Value
                convergedCorrectlyPct = [double]$m.Groups[6].Value
            }
        }
        $record.convergenceSweep = [pscustomobject]@{
            reportFile = "tests/convergence_sweep_report.txt"
            summary = $sweepSummary
        }
        Write-Host "convergence sweep summary:"
        $sweepSummary | Format-Table | Out-String | Write-Host
    } else {
        Write-Host "convergence sweep did not produce a report file -- recording as unavailable."
        $record.convergenceSweep = [pscustomobject]@{ reportFile = $null; summary = @() }
    }
} else {
    Write-Host ""
    Write-Host "==> [Convergence sweep skipped (-SkipConvergenceSweep)]"
}

# --- Item 10: artifact checksum and final release notes ---

$artifactPath = Join-Path $RepoRoot "$packageName $version.zip"
if (-not (Test-Path -LiteralPath $artifactPath)) {
    throw "expected release artifact not found at $artifactPath after the build step"
}
$hash = Get-FileHash -LiteralPath $artifactPath -Algorithm SHA256
$changelogPath = Join-Path $RepoRoot "CHANGELOG.md"
$changelogText = Get-Content -LiteralPath $changelogPath -Raw
$notesMatch = [regex]::Match($changelogText, "(?ms)^##\s+$([regex]::Escape($version))\s*.*?(?=^##\s+\S|\z)")
$releaseNotes = if ($notesMatch.Success) { $notesMatch.Value.Trim() } else { $null }

$record.artifact = [pscustomobject]@{
    path = "$packageName $version.zip"
    sha256 = $hash.Hash
    sizeBytes = (Get-Item -LiteralPath $artifactPath).Length
    releaseNotes = $releaseNotes
}

Write-Host ""
Write-Host "artifact: $($record.artifact.path)"
Write-Host "sha256:   $($record.artifact.sha256)"

Save-ReleaseRecord -Status "go"

Write-Host ""
Write-Host "run_release_gate.ps1: GO -- release candidate $packageName $version passed every gate."
