# verify_sslib_pin.ps1
#
# TVP-AIDS initiative, Stage 6: checks the shared sslib install
# (sslib.pin.json's own installPath) against the handful of files that
# have actually drifted or mattered in past incidents (Stage 3's
# init_diffTVP arity break; Stage 4's stale sstvp.src/sskalman.src -- see
# sslib.pin.json's own header comment). Prints a CLEAR warning up front
# instead of letting drift surface later as a confusing low-level GAUSS
# error deep inside a kalmanFilter*/ssFitTVP call -- the actual failure
# mode both past incidents took the form of.
#
# Exit code is 0 (even with a mismatch) unless -Strict is passed -- a
# newer sslib commit may still work fine, so this is advisory by default;
# use -Strict for a hard gate (e.g. before a release).
#
# Wired into tests/run_source_tests.ps1's -SkipTVPKalman-false branch,
# run once before the sslib-dependent tests themselves.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$PinFile = $null,
    [switch]$Strict
)

if ([string]::IsNullOrWhiteSpace($PinFile)) {
    $PinFile = Join-Path $RepoRoot "sslib.pin.json"
}
if (-not (Test-Path -LiteralPath $PinFile)) {
    throw "sslib.pin.json not found at $PinFile"
}

$pin = Get-Content -LiteralPath $PinFile -Raw | ConvertFrom-Json
$installPath = [string]$pin.installPath
# .PSObject.Properties.Name, wrapped in @(), yields a one-element array
# containing $null for a truly EMPTY PSCustomObject (a PowerShell 5.1
# collection-unwrapping quirk, confirmed directly) -- check .Count on
# .Properties itself first, and only enumerate Name via ForEach-Object
# (safe either way) once known non-empty.
$trackedFiles = @()
if ($pin.verifiedFiles.PSObject.Properties.Count -gt 0) {
    $trackedFiles = @($pin.verifiedFiles.PSObject.Properties | ForEach-Object { $_.Name })
}

if (-not (Test-Path -LiteralPath $installPath)) {
    Write-Warning "verify_sslib_pin.ps1: sslib install not found at $installPath (pinned commit: $($pin.pinnedCommit)) -- sslib-dependent tests will fail to compile."
    if ($Strict) { exit 1 }
    exit 0
}

if ($trackedFiles.Count -eq 0) {
    Write-Warning "verify_sslib_pin.ps1: sslib.pin.json has no verifiedFiles yet -- run scripts\sync_sslib.ps1 to populate it. Skipping the drift check for now."
    exit 0
}

$mismatches = @()
$missing = @()
foreach ($relPath in $trackedFiles) {
    $fullPath = Join-Path $installPath $relPath
    if (-not (Test-Path -LiteralPath $fullPath)) {
        $missing += $relPath
        continue
    }
    $actualHash = (Get-FileHash -LiteralPath $fullPath -Algorithm SHA256).Hash.ToLowerInvariant()
    $expectedHash = [string]$pin.verifiedFiles.$relPath
    if ($actualHash -ne $expectedHash) {
        $mismatches += $relPath
    }
}

if ($missing.Count -gt 0 -or $mismatches.Count -gt 0) {
    Write-Warning "verify_sslib_pin.ps1: installed sslib at $installPath does NOT match the pinned commit $($pin.pinnedCommit) (last verified $($pin.lastVerified))."
    if ($missing.Count -gt 0) {
        Write-Warning "  missing: $($missing -join ', ')"
    }
    if ($mismatches.Count -gt 0) {
        Write-Warning "  content differs: $($mismatches -join ', ')"
    }
    Write-Warning "  Run 'scripts\sync_sslib.ps1 -SourceRepoPath <path to a gauss-state-space checkout>' to resync, or -Commit <hash> to repin -- see that script's own header. A newer commit may still work fine; this is advisory, not necessarily broken."
    if ($Strict) {
        Write-Host "verify_sslib_pin.ps1: FAIL (-Strict)"
        exit 1
    }
    Write-Host "verify_sslib_pin.ps1: WARN -- drift detected, continuing (pass -Strict to make this a hard failure)"
    exit 0
}

Write-Host "verify_sslib_pin.ps1: PASS -- installed sslib at $installPath matches the pinned commit $($pin.pinnedCommit) ($($trackedFiles.Count) files verified)"
