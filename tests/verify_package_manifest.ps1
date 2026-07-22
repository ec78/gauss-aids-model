# verify_package_manifest.ps1
#
# Milestone 7: sanity-checks package.json against the actual src/ directory
# -- no duplicate/missing/unlisted entries, quaids.sdf loads first (structs
# before procedures that reference them). Adapted from gauss-qardl's
# tests/verify_package_manifest.ps1.
#
# Milestone 8: also cross-checks docs/COMMAND_REFERENCE.md against the
# actual source -- every documented proc must actually be defined
# somewhere in src/ (including src/pubtable_quaids.src, since that file's
# procs are documented too, even though it's intentionally excluded from
# package.json's src array -- see the allowlist below), and every linked
# command-reference page must exist.

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

# --- Milestone 8: docs/COMMAND_REFERENCE.md cross-check ---

$commandRefPath = Join-Path $RepoRoot "docs\COMMAND_REFERENCE.md"
if (-not (Test-Path -LiteralPath $commandRefPath)) {
    throw "docs/COMMAND_REFERENCE.md not found at $commandRefPath"
}

$commandRef = Get-Content -LiteralPath $commandRefPath -Raw
$linkMatches = [regex]::Matches($commandRef, '\[([A-Za-z_][A-Za-z0-9_]*)\]\(command-reference/([^)]+\.md)\)')

$docCommands = New-Object System.Collections.Generic.List[string]
$linkedDocPaths = New-Object System.Collections.Generic.List[string]
foreach ($match in $linkMatches) {
    $docCommands.Add($match.Groups[1].Value)
    $linkedDocPaths.Add($match.Groups[2].Value)
}
$docCommands = @($docCommands | Sort-Object -Unique)

if ($docCommands.Count -eq 0) {
    throw "docs/COMMAND_REFERENCE.md does not list any public commands"
}

$missingDocPages = @()
foreach ($relPath in ($linkedDocPaths | Sort-Object -Unique)) {
    $fullPath = Join-Path (Join-Path $RepoRoot "docs\command-reference") $relPath
    if (-not (Test-Path -LiteralPath $fullPath)) {
        $missingDocPages += $relPath
    }
}
if ($missingDocPages.Count -gt 0) {
    throw "docs/COMMAND_REFERENCE.md references missing command pages: $($missingDocPages -join ', ')"
}

# Every .src file actually present in src/ -- both package.json's required
# src array and the intentionally-unlisted allowlist (pubtable_quaids.src)
# -- is fair game for documented procs, since this repo documents the
# optional pubtable adapter too.
$allSrcFiles = $srcEntries + $intentionallyUnlisted | Sort-Object -Unique
$sourceText = ""
foreach ($entry in $allSrcFiles) {
    if ([System.IO.Path]::GetExtension($entry) -ne ".src") {
        continue
    }
    $sourceText += "`n"
    $sourceText += Get-Content -LiteralPath (Join-Path $srcDir $entry) -Raw
}

# Matches all three GAUSS proc-declaration forms (see the identical fix in
# scripts/build_lcg.ps1's header comment for why: "proc (struct X) =
# name(...)" alone is not enough for this codebase).
$procMatches = [regex]::Matches($sourceText, '(?m)^\s*proc\s*(?:\(([^)]*)\)|\d+)?\s*(?:=\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\(')
$exportedProcs = @($procMatches | ForEach-Object { $_.Groups[2].Value } | Sort-Object -Unique)

$missingDocumentedProcs = $docCommands | Where-Object { $exportedProcs -notcontains $_ }
if ($missingDocumentedProcs.Count -gt 0) {
    throw "docs/COMMAND_REFERENCE.md documents procedures not found in src/: $($missingDocumentedProcs -join ', ')"
}

Write-Host "verify_package_manifest.ps1: PASS"
