# verify_public_api.ps1
#
# Public release contract check (PUBLIC_RELEASE_ROADMAP.md Phase 0):
#
# 1. PR-001 acceptance evidence -- "one automated metadata check confirms
#    the version across all release files": package.json, CITATION.cff,
#    docs/public-api.json, and the CHANGELOG.md's top-level entry must all
#    agree on the same version string.
# 2. PR-003 acceptance evidence -- "a machine-readable public API
#    inventory is committed and checked in CI": every procedure listed in
#    docs/public-api.json's supported_procedures/compatibility_procedures
#    must actually be defined in src/ (package.json's own src array, plus
#    the optional pubtable_quaids.src adapter for optional_adapter_procedures)
#    and must have a linked docs/command-reference/*.md page; every struct
#    listed in public_structs.names must actually be defined in
#    src/quaids.sdf.
#
# This intentionally does not require the reverse (every proc in src/ must
# appear in public-api.json) -- this codebase has many private `_quaids*`
# helpers by established convention, and enumerating all of them in a
# public-facing inventory would be noise, not signal.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
)

function Get-ProcNames {
    param([string]$Text)
    # Same three-proc-declaration-form regex as build_lcg.ps1/
    # verify_package_manifest.ps1 -- "proc (struct X) = name(...)",
    # "proc N = name(...)", and "proc name(...)".
    $matches_ = [regex]::Matches($Text, '(?m)^\s*proc\s*(?:\(([^)]*)\)|\d+)?\s*(?:=\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\(')
    return @($matches_ | ForEach-Object { $_.Groups[2].Value } | Sort-Object -Unique)
}

$packagePath = Join-Path $RepoRoot "package.json"
$citationPath = Join-Path $RepoRoot "CITATION.cff"
$changelogPath = Join-Path $RepoRoot "CHANGELOG.md"
$publicApiPath = Join-Path $RepoRoot "docs\public-api.json"
$commandRefPath = Join-Path $RepoRoot "docs\COMMAND_REFERENCE.md"
$srcDir = Join-Path $RepoRoot "src"

foreach ($p in @($packagePath, $citationPath, $changelogPath, $publicApiPath, $commandRefPath)) {
    if (-not (Test-Path -LiteralPath $p)) {
        throw "required file not found: $p"
    }
}

$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$publicApi = Get-Content -LiteralPath $publicApiPath -Raw | ConvertFrom-Json
$citationText = Get-Content -LiteralPath $citationPath -Raw
$changelogText = Get-Content -LiteralPath $changelogPath -Raw

# --- 1. Version consistency across release metadata ---

$version = [string]$pkg.version
if ([string]::IsNullOrWhiteSpace($version)) {
    throw "package.json version is empty"
}

$citationMatch = [regex]::Match($citationText, '(?m)^version:\s*"?([^"\r\n]+)"?\s*$')
if (-not $citationMatch.Success) {
    throw "CITATION.cff has no 'version:' entry"
}
$citationVersion = $citationMatch.Groups[1].Value.Trim()
if ($citationVersion -ne $version) {
    throw "version mismatch: package.json is '$version', CITATION.cff is '$citationVersion'"
}

$publicApiVersion = [string]$publicApi.version
if ($publicApiVersion -ne $version) {
    throw "version mismatch: package.json is '$version', docs/public-api.json is '$publicApiVersion'"
}

if ($changelogText -notmatch "(?m)^##\s+$([regex]::Escape($version))\b") {
    throw "CHANGELOG.md does not contain a top-level entry for version $version"
}

Write-Host "verify_public_api.ps1: release metadata version '$version' consistent across package.json, CITATION.cff, docs/public-api.json, CHANGELOG.md"

# --- 2. Public API inventory reconciliation ---

$srcEntries = @($pkg.src)
$sourceText = ""
foreach ($entry in $srcEntries) {
    $entryPath = Join-Path $srcDir $entry
    if (-not (Test-Path -LiteralPath $entryPath)) {
        throw "package.json src entry not found: $entryPath"
    }
    if ([System.IO.Path]::GetExtension($entry) -eq ".src") {
        $sourceText += "`n" + (Get-Content -LiteralPath $entryPath -Raw)
    }
}

# Optional modules (public release roadmap PR-101/PR-003, extended in
# TVP-AIDS Stage 6): files deliberately excluded from package.json's src
# array (a hard compile-time dependency on another package's struct types
# would otherwise force that package on every quaids user) --
# pubtable_quaids.src (needs pubtable), quaidscurvature.src (needs optmt),
# and six TVP-AIDS files (need sslib -- see sslib.pin.json/CLAUDE.md's
# "Optional modules" section). Each still ships in the installed
# package's src/ directory and is real, documented, supported API; a
# caller opts in via an explicit #include, per each module's own "setup"
# string in docs/public-api.json.
$corePocs = Get-ProcNames -Text $sourceText
# module.source is normally a single file path (curvature/pubtable, one file
# each), but a module's public procs can genuinely span several files (the
# TVP-AIDS optional module -- its public API is spread across the six files
# Stages 1-6 built it in). @() treats a single string as a one-element array,
# so this stays backward-compatible with the existing single-string entries.
# Keyed by module NAME, not by source, since source may no longer be a valid
# scalar hashtable key.
$moduleProcsByName = @{}
foreach ($module in @($publicApi.optional_modules)) {
    $moduleSources = @($module.source)
    $moduleText = ""
    foreach ($modulePath in $moduleSources) {
        $fullModulePath = Join-Path $RepoRoot ([string]$modulePath)
        if (-not (Test-Path -LiteralPath $fullModulePath)) {
            throw "docs/public-api.json optional_modules['$($module.name)'].source not found: $modulePath"
        }
        $moduleText += "`n" + (Get-Content -LiteralPath $fullModulePath -Raw)
    }
    $moduleProcsByName[[string]$module.name] = Get-ProcNames -Text $moduleText
}
$allModuleProcs = @($moduleProcsByName.Values | ForEach-Object { $_ } | Sort-Object -Unique)
$allKnownProcs = @($corePocs + $allModuleProcs | Sort-Object -Unique)

$commandRefText = Get-Content -LiteralPath $commandRefPath -Raw
$linkMatches = [regex]::Matches($commandRefText, '\[([A-Za-z_][A-Za-z0-9_]*)\]\(command-reference/([^)]+\.md)\)')
$documentedProcs = @($linkMatches | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)

$supportedProcs = @($publicApi.supported_procedures)
$compatProcs = @($publicApi.compatibility_procedures | ForEach-Object { [string]$_.name })

$allInventoryProcs = @($supportedProcs + $compatProcs | Sort-Object -Unique)

$missingFromSource = $allInventoryProcs | Where-Object { $allKnownProcs -notcontains $_ }
if ($missingFromSource.Count -gt 0) {
    throw "docs/public-api.json lists procedures not found in src/: $($missingFromSource -join ', ')"
}

$missingFromDocs = $allInventoryProcs | Where-Object { $documentedProcs -notcontains $_ }
if ($missingFromDocs.Count -gt 0) {
    throw "docs/public-api.json lists procedures with no docs/COMMAND_REFERENCE.md entry: $($missingFromDocs -join ', ')"
}

foreach ($module in @($publicApi.optional_modules)) {
    $moduleProcs = $moduleProcsByName[[string]$module.name]
    $notInModuleSource = @($module.procedures) | Where-Object { $moduleProcs -notcontains $_ }
    if ($notInModuleSource.Count -gt 0) {
        $sourceDesc = (@($module.source) -join ', ')
        throw "docs/public-api.json optional_modules['$($module.name)'] lists procedures not found in $sourceDesc`: $($notInModuleSource -join ', ')"
    }
}

Write-Host "verify_public_api.ps1: $($allInventoryProcs.Count) inventoried procedures all exist in src/ and are documented"

# --- 3. Public struct inventory ---

$sdfPath = Join-Path $srcDir "quaids.sdf"
if (-not (Test-Path -LiteralPath $sdfPath)) {
    throw "src/quaids.sdf not found"
}
$sdfText = Get-Content -LiteralPath $sdfPath -Raw
$structNames = @($publicApi.public_structs.names)
$missingStructs = @()
foreach ($structName in $structNames) {
    if ($sdfText -notmatch "(?m)^\s*struct\s+$([regex]::Escape($structName))\s*\{") {
        $missingStructs += $structName
    }
}
if ($missingStructs.Count -gt 0) {
    throw "docs/public-api.json lists structs not found in src/quaids.sdf: $($missingStructs -join ', ')"
}

Write-Host "verify_public_api.ps1: $($structNames.Count) inventoried structs all exist in src/quaids.sdf"
Write-Host "verify_public_api.ps1: PASS"
