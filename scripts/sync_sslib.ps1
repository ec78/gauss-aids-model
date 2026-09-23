# sync_sslib.ps1
#
# TVP-AIDS initiative, Stage 6: makes "resync the shared sslib
# (gauss-state-space) install to this repo's documented pin" a single
# reproducible command, replacing the ad hoc `git archive`-by-hand process
# used every TVP-AIDS stage so far -- the actual root cause of the two
# real drift incidents this initiative already hit (Stage 3's
# init_diffTVP arity break; Stage 4's stale sstvp.src/sskalman.src): no
# repeatable resync path meant a partial/stale state could persist
# unnoticed. See sslib.pin.json (repo root) for the current pin and
# CLAUDE.md's "sslib (gauss-state-space)" note for background.
#
# Two modes:
#   - Resync (default): rebuilds the shared install from sslib.pin.json's
#     OWN pinnedCommit -- use this to recover from drift or after a fresh
#     machine setup.
#   - Repin (-Commit <hash>): also MOVES the pin itself to a new commit
#     before syncing -- use this deliberately, after confirming with any
#     concurrent gauss-state-space session (see CLAUDE.md's shared-
#     package-directory collision-risk note) that the new commit is a
#     real, intended, finished push, not someone's in-progress work.
#
# Reads the pinned/target commit from a LOCAL git checkout of
# gauss-state-space (-SourceRepoPath) via `git archive` -- never from that
# repo's own live working tree state (only its committed git history),
# and never from the shared install itself (which this script overwrites,
# not reads from).
#
# After replacing the install's contents, rebuilds its lib/sslib.lcg
# catalog via this repo's own build_lcg.ps1 (confirmed package-agnostic),
# then recomputes and rewrites sslib.pin.json's verifiedFiles/
# lastVerified so scripts/verify_sslib_pin.ps1 has something current to
# check future runs against.

param(
    [Parameter(Mandatory = $true)]
    [string]$SourceRepoPath,
    [string]$Commit = $null,
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$PinFile = $null
)

if ([string]::IsNullOrWhiteSpace($PinFile)) {
    $PinFile = Join-Path $RepoRoot "sslib.pin.json"
}
if (-not (Test-Path -LiteralPath $PinFile)) {
    throw "sslib.pin.json not found at $PinFile"
}

$pin = Get-Content -LiteralPath $PinFile -Raw | ConvertFrom-Json

if ([string]::IsNullOrWhiteSpace($Commit)) {
    $targetCommit = [string]$pin.pinnedCommit
    Write-Host "sync_sslib.ps1: resyncing to the currently-pinned commit $targetCommit (pass -Commit to repin to a different one)"
} else {
    $targetCommit = $Commit
    Write-Host "sync_sslib.ps1: REPINNING from $($pin.pinnedCommit) to $targetCommit"
}

if ([string]::IsNullOrWhiteSpace($targetCommit)) {
    throw "no target commit -- sslib.pin.json's pinnedCommit is empty and -Commit was not given"
}

if (-not (Test-Path -LiteralPath (Join-Path $SourceRepoPath ".git"))) {
    throw "SourceRepoPath does not look like a git checkout (no .git found): $SourceRepoPath"
}

# Confirm the target commit actually exists in this checkout's history
# before doing anything destructive.
& git -C $SourceRepoPath cat-file -e "$targetCommit^{commit}" 2>$null
if ($LASTEXITCODE -ne 0) {
    throw "commit $targetCommit not found in $SourceRepoPath's git history -- fetch it first (git fetch) or check the commit hash"
}

$installPath = [string]$pin.installPath
if ([string]::IsNullOrWhiteSpace($installPath)) {
    throw "sslib.pin.json has no installPath"
}

$stagingDir = Join-Path ([System.IO.Path]::GetTempPath()) ("sslib_sync_" + [System.Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $stagingDir | Out-Null
$zipPath = Join-Path $stagingDir "sslib.zip"
$extractDir = Join-Path $stagingDir "extracted"

try {
    Write-Host "sync_sslib.ps1: archiving $targetCommit from $SourceRepoPath"
    & git -C $SourceRepoPath archive --format=zip --output=$zipPath $targetCommit
    if ($LASTEXITCODE -ne 0) {
        throw "git archive failed (exit $LASTEXITCODE)"
    }

    Expand-Archive -LiteralPath $zipPath -DestinationPath $extractDir -Force

    if (-not (Test-Path -LiteralPath (Join-Path $extractDir "package.json"))) {
        throw "archived commit $targetCommit has no package.json at its root -- not a valid sslib package tree"
    }

    Write-Host "sync_sslib.ps1: replacing $installPath with the archived tree (robocopy /MIR)"
    if (-not (Test-Path -LiteralPath $installPath)) {
        New-Item -ItemType Directory -Path $installPath | Out-Null
    }
    # /MIR mirrors the source into the destination, removing anything in
    # the destination that is no longer in the source -- a clean replace,
    # not an additive copy (stale files from an older commit must not
    # survive a resync). robocopy's own "success" exit codes are 0-7, not
    # just 0.
    $robocopyOutput = & robocopy $extractDir $installPath /MIR /NFL /NDL /NJH /NJS
    if ($LASTEXITCODE -ge 8) {
        throw "robocopy failed (exit $LASTEXITCODE): $robocopyOutput"
    }

    Write-Host "sync_sslib.ps1: rebuilding lib/sslib.lcg"
    & (Join-Path $PSScriptRoot "build_lcg.ps1") -PackageRoot $installPath -PackageName "sslib"

    # Recompute verifiedFiles hashes over the same handful of files that
    # have actually drifted or mattered in past incidents (see
    # sslib.pin.json's own header comment / this file's own header) --
    # not the whole tree, to keep this maintainable.
    # See verify_sslib_pin.ps1's matching comment -- @() around
    # .PSObject.Properties.Name yields a one-element array containing
    # $null for a truly empty PSCustomObject, not an empty array.
    $trackedFiles = @()
    if ($pin.verifiedFiles.PSObject.Properties.Count -gt 0) {
        $trackedFiles = @($pin.verifiedFiles.PSObject.Properties | ForEach-Object { $_.Name })
    }
    if ($trackedFiles.Count -eq 0) {
        $trackedFiles = @("src/sstvp.src", "src/sskalman.src", "src/ssmain.src", "src/ssstructural.src")
    }

    $newVerifiedFiles = [ordered]@{}
    foreach ($relPath in $trackedFiles) {
        $fullPath = Join-Path $installPath $relPath
        if (-not (Test-Path -LiteralPath $fullPath)) {
            Write-Warning "sync_sslib.ps1: tracked file not found after sync, skipping: $relPath"
            continue
        }
        $hash = (Get-FileHash -LiteralPath $fullPath -Algorithm SHA256).Hash.ToLowerInvariant()
        $newVerifiedFiles[$relPath] = $hash
    }

    $pin.pinnedCommit = $targetCommit
    $pin.verifiedFiles = $newVerifiedFiles
    $pin.lastVerified = (Get-Date -Format "yyyy-MM-dd")
    $pinJson = $pin | ConvertTo-Json -Depth 6
    $utf8NoBom = New-Object System.Text.UTF8Encoding $false
    [System.IO.File]::WriteAllText($PinFile, $pinJson, $utf8NoBom)

    Write-Host "sync_sslib.ps1: done -- $installPath now matches $targetCommit ($($newVerifiedFiles.Count) files hashed into $PinFile)"
} finally {
    Remove-Item -LiteralPath $stagingDir -Recurse -Force -ErrorAction SilentlyContinue
}
