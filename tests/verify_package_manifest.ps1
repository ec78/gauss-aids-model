# verify_package_manifest.ps1
#
# Milestone 7: sanity-checks package.json against the actual src/ directory
# -- no duplicate/missing/unlisted entries, quaids.sdf loads first (structs
# before procedures that reference them). Adapted from gauss-qardl's
# tests/verify_package_manifest.ps1.
#
# Does NOT yet cross-check package.json against docs/COMMAND_REFERENCE.md
# the way gauss-qardl's version does -- docs/ is Milestone 8 in this repo's
# roadmap and does not exist yet. Add that check here once
# docs/COMMAND_REFERENCE.md exists.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
)

$packagePath = Join-Path $RepoRoot "package.json"
$srcDir = Join-Path $RepoRoot "src"

if (-not (Test-Path -LiteralPath $packagePath)) {
    throw "package.json not found at $packagePath"
}

$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$srcEntries = @($pkg.src)

if ($srcEntries.Count -eq 0) {
    throw "package.json src array is empty"
}

$duplicates = $srcEntries | Group-Object | Where-Object { $_.Count -gt 1 }
if ($duplicates.Count -gt 0) {
    $names = ($duplicates | ForEach-Object { $_.Name }) -join ", "
    throw "package.json src has duplicate entries: $names"
}

if ($srcEntries[0] -ne "quaids.sdf") {
    throw "package.json src must list quaids.sdf first so structures are registered before procedures"
}

$missing = @()
foreach ($entry in $srcEntries) {
    $path = Join-Path $srcDir $entry
    if (-not (Test-Path -LiteralPath $path)) {
        $missing += $entry
    }
}

if ($missing.Count -gt 0) {
    throw "package.json src references missing files: $($missing -join ', ')"
}

# pubtable_quaids.src is deliberately excluded from package.json's src array
# -- it has a hard compile-time dependency on pubtable.sdf's struct types
# (see CLAUDE.md's "Milestone 6: reporting via pubtable" section), so
# listing it would make pubtable a hard dependency for the whole package to
# even compile. Any other .src/.sdf file added to src/ is expected to be a
# required part of the package and must be listed.
$intentionallyUnlisted = @("pubtable_quaids.src")

$actualSrc = Get-ChildItem -LiteralPath $srcDir -File |
    Where-Object { $_.Extension -in ".src", ".sdf" } |
    ForEach-Object { $_.Name }

$unlisted = $actualSrc | Where-Object { $srcEntries -notcontains $_ -and $intentionallyUnlisted -notcontains $_ }
if ($unlisted.Count -gt 0) {
    throw "src files not listed in package.json (and not in the intentionally-unlisted allowlist): $($unlisted -join ', ')"
}

$missingAllowlisted = $intentionallyUnlisted | Where-Object { $actualSrc -notcontains $_ }
if ($missingAllowlisted.Count -gt 0) {
    throw "verify_package_manifest.ps1's intentionally-unlisted allowlist references files that no longer exist: $($missingAllowlisted -join ', ')"
}

if ([string]::IsNullOrWhiteSpace([string]$pkg.name)) {
    throw "package.json name is empty"
}

if ([string]::IsNullOrWhiteSpace([string]$pkg.version)) {
    throw "package.json version is empty"
}

if ([string]::IsNullOrWhiteSpace([string]$pkg.license)) {
    throw "package.json license is empty"
}

Write-Host "verify_package_manifest.ps1: PASS"
