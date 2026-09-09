# verify_release_artifact.ps1
#
# Milestone 7: sanity-checks a built release .zip -- name includes the
# package version, CHANGELOG.md has a matching entry, no stale/pre-version
# artifacts sit in the repo root, and the archive actually contains the
# files package.json's src array promises plus the small set of root
# files/dirs this repo currently ships. Adapted from gauss-qardl's
# scripts/verify_release_artifact.ps1.
#
# Milestone 8: README.md and docs/COMMAND_REFERENCE.md now exist and are
# required entries below, closing the gap this file's Milestone 7 comment
# flagged.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$ArtifactPath = ""
)

$packagePath = Join-Path $RepoRoot "package.json"
if (-not (Test-Path -LiteralPath $packagePath)) {
    throw "package.json not found at $packagePath"
}

$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$version = [string]$pkg.version
$packageName = [string]$pkg.name

if ([string]::IsNullOrWhiteSpace($version)) {
    throw "package.json version is empty"
}

if ([string]::IsNullOrWhiteSpace($packageName)) {
    throw "package.json name is empty"
}

if ([string]::IsNullOrWhiteSpace($ArtifactPath)) {
    $ArtifactPath = Join-Path $RepoRoot "$packageName $version.zip"
}

if (-not (Test-Path -LiteralPath $ArtifactPath)) {
    throw "release artifact not found at $ArtifactPath"
}

$artifactName = Split-Path -Leaf $ArtifactPath
if ($artifactName -notmatch [regex]::Escape($version)) {
    throw "release artifact name '$artifactName' does not include package version '$version'"
}

$changeLogPath = Join-Path $RepoRoot "CHANGELOG.md"
if (Test-Path -LiteralPath $changeLogPath) {
    $changeLog = Get-Content -LiteralPath $changeLogPath -Raw
    if ($changeLog -notmatch "(?m)^##\s+$([regex]::Escape($version))\s") {
        throw "CHANGELOG.md does not contain a top-level entry for package version $version"
    }
}

$artifactFullPath = (Resolve-Path -LiteralPath $ArtifactPath).Path
$artifactPatterns = @(
    (Join-Path $RepoRoot "$packageName *.zip"),
    (Join-Path $RepoRoot "$packageName`_*.zip")
)
$staleArtifacts = Get-ChildItem -Path $artifactPatterns -File -ErrorAction SilentlyContinue |
    Where-Object { $_.FullName -ne $artifactFullPath -and $_.Name -notmatch [regex]::Escape($version) }
if ($staleArtifacts.Count -gt 0) {
    $names = ($staleArtifacts | ForEach-Object { $_.Name }) -join ", "
    throw "stale package artifacts found in repo root: $names"
}

Add-Type -AssemblyName System.IO.Compression.FileSystem
$zip = [System.IO.Compression.ZipFile]::OpenRead((Resolve-Path -LiteralPath $ArtifactPath).Path)
try {
    $entryNames = @($zip.Entries | ForEach-Object { $_.FullName -replace "\\", "/" })

    $requiredEntries = @(
        "package.json",
        "README.md",
        "CHANGELOG.md",
        "CITATION.cff",
        "LICENSE",
        "src/quaids.sdf",
        "examples/00_real_data_quickstart.e",
        "examples/01_basic_estimation.e",
        "tests/fixtures/published/blanciforti86_food32.csv",
        "scripts/build_package.ps1",
        "scripts/verify_release_artifact.ps1",
        "docs/COMMAND_REFERENCE.md",
        "docs/USAGE_GUIDE.md",
        "docs/METHODOLOGY_NOTES.md",
        "docs/FEATURE_SUPPORT_MATRIX.md"
    )

    $missingEntries = $requiredEntries | Where-Object { $entryNames -notcontains $_ }
    if ($missingEntries.Count -gt 0) {
        throw "release artifact is missing required entries: $($missingEntries -join ', ')"
    }

    # PR-405: this pattern is the authoritative "must never ship" list --
    # kept in sync by hand with build_package.ps1's own $generatedTestFiles
    # cleanup list (which removes these from the staging directory before
    # zipping) and tests/run_examples_smoke.ps1's own post-run cleanup
    # (which removes them from the repo working tree after a local smoke
    # run). This check is the actual gate; the other two are best-effort
    # prevention -- this is what fails a release if either one is ever
    # incomplete.
    $badTempEntries = $entryNames | Where-Object {
        $_ -match "^tests/(pubtable_test_coef\.(tex|md|csv)|schema_test_quaids_wrapper_out|print_format_probe_(elas|shares|noint)\.txt|.*\.log)$" -or
        $_ -match "^examples/(quaids_coefficients\.(tex|md|csv)|quaids_income_elasticities\.md|quaids_uncompensated_elasticities\.tex|quaids_compensated_elasticities\.csv|quaids_workflow.*|blanciforti_results\.txt)$"
    }
    if ($badTempEntries.Count -gt 0) {
        throw "release artifact includes generated test-run artifacts: $($badTempEntries -join ', ')"
    }

    $pkgEntry = $zip.GetEntry("package.json")
    if ($null -eq $pkgEntry) {
        throw "release artifact does not include package.json"
    }

    $reader = [System.IO.StreamReader]::new($pkgEntry.Open())
    try {
        $artifactPkg = $reader.ReadToEnd() | ConvertFrom-Json
    } finally {
        $reader.Dispose()
    }

    if ([string]$artifactPkg.name -ne $packageName) {
        throw "artifact package name '$($artifactPkg.name)' does not match source package name '$packageName'"
    }

    if ([string]$artifactPkg.version -ne $version) {
        throw "artifact package version '$($artifactPkg.version)' does not match source package version '$version'"
    }

    foreach ($srcEntry in @($pkg.src)) {
        $entryPath = "src/$srcEntry"
        if ($entryNames -notcontains $entryPath) {
            throw "release artifact is missing source file listed in package.json: $entryPath"
        }
    }
} finally {
    $zip.Dispose()
}

Write-Host "verify_release_artifact.ps1: PASS"
