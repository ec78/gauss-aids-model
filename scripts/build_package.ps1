# build_package.ps1
#
# Milestone 7: stages a distributable copy of the repo (package.json plus
# whichever of README.md/CHANGELOG.md/CITATION.cff/CITATION.md/LICENSE/
# llms.txt exist) and zips it as "<name> <version>.zip" in the repo root.
# Adapted from gauss-qardl's scripts/build_package.ps1. Root files and
# directories are copied only if present, so this script does not need to
# change when README.md/docs/ land at Milestone 8 -- they'll just start
# being included automatically.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$OutputDir = "",
    [switch]$Force,
    [switch]$NoTests
)

if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = $RepoRoot
}

$packagePath = Join-Path $RepoRoot "package.json"
if (-not (Test-Path -LiteralPath $packagePath)) {
    throw "package.json not found at $packagePath"
}

$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
$packageName = [string]$pkg.name
$version = [string]$pkg.version

if ([string]::IsNullOrWhiteSpace($packageName) -or [string]::IsNullOrWhiteSpace($version)) {
    throw "package.json must define name and version"
}

if (-not (Test-Path -LiteralPath $OutputDir)) {
    New-Item -ItemType Directory -Path $OutputDir | Out-Null
}

$artifactPath = Join-Path $OutputDir "$packageName $version.zip"
if ((Test-Path -LiteralPath $artifactPath) -and -not $Force) {
    throw "release artifact already exists at $artifactPath. Re-run with -Force to replace it."
}

$stageRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_pkg_" + [System.Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $stageRoot | Out-Null

try {
    $rootFiles = @(
        "package.json",
        "README.md",
        "CHANGELOG.md",
        "CITATION.cff",
        "CITATION.md",
        "LICENSE",
        "llms.txt",
        "SUPPORT.md",
        "CONTRIBUTING.md"
    )

    foreach ($file in $rootFiles) {
        $srcPath = Join-Path $RepoRoot $file
        if (Test-Path -LiteralPath $srcPath) {
            Copy-Item -LiteralPath $srcPath -Destination (Join-Path $stageRoot $file)
        }
    }

    $dirs = @("src", "docs", "examples", "scripts")
    if (-not $NoTests) {
        $dirs += "tests"
    }

    foreach ($dir in $dirs) {
        $srcPath = Join-Path $RepoRoot $dir
        if (Test-Path -LiteralPath $srcPath) {
            Copy-Item -LiteralPath $srcPath -Destination (Join-Path $stageRoot $dir) -Recurse
        }
    }

    # Strip generated/gitignored run artifacts that may exist locally
    # (test/example table exports, tgauss run logs) so the staged package
    # only contains source-controlled content. Named explicitly and removed
    # by literal path rather than via Get-ChildItem -Include -Recurse --
    # that combination silently ignores -Include when the base path is not
    # itself a wildcard, which previously caused this block to delete every
    # file under tests/ and examples/, not just the generated ones (caught
    # by verify_release_artifact.ps1 failing with "missing required entry:
    # examples/quaids_example.e" the first time this script ran for real).
    $generatedTestFiles = @(
        "tests\pubtable_test_coef.tex",
        "tests\pubtable_test_coef.md",
        "tests\pubtable_test_coef.csv",
        "tests\schema_test_quaids_wrapper_out",
        "tests\print_format_probe_elas.txt",
        "tests\print_format_probe_shares.txt",
        "tests\print_format_probe_noint.txt",
        "examples\quaids_coefficients.tex",
        "examples\quaids_coefficients.md",
        "examples\quaids_coefficients.csv",
        "examples\quaids_income_elasticities.md",
        "examples\quaids_uncompensated_elasticities.tex",
        "examples\quaids_compensated_elasticities.csv",
        "examples\quaids_workflow",
        "examples\blanciforti_results.txt"
    )
    foreach ($relPath in $generatedTestFiles) {
        $fullPath = Join-Path $stageRoot $relPath
        if (Test-Path -LiteralPath $fullPath) {
            Remove-Item -LiteralPath $fullPath -Force
        }
    }
    Get-ChildItem -LiteralPath $stageRoot -Recurse -File -Filter "*.log" |
        Remove-Item -Force

    Get-ChildItem -LiteralPath $stageRoot -Recurse -File -Filter "*.zip" |
        Remove-Item -Force

    # Compress-Archive (and .NET's own ZipFile.CreateFromDirectory, under
    # the .NET Framework this Windows PowerShell 5.1 host runs on) both
    # write entry names using the platform path separator, i.e. backslashes
    # on Windows -- verified directly against a scratch fixture, not
    # assumed. The ZIP spec calls for forward slashes; a strict or
    # cross-platform unzip implementation (as the GAUSS package installer
    # may use) can fail to recognize backslash-separated entries as nested
    # paths at all, extracting them as flat, literally-backslashed
    # filenames instead of populating src/docs/tests/examples
    # subdirectories. Build the archive entry-by-entry instead, forcing
    # forward slashes explicitly.
    $tmpArtifact = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_artifact_" + [System.Guid]::NewGuid().ToString("N") + ".zip")
    Add-Type -AssemblyName System.IO.Compression
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $fs = [System.IO.File]::Open($tmpArtifact, [System.IO.FileMode]::Create)
    try {
        $archive = New-Object System.IO.Compression.ZipArchive($fs, [System.IO.Compression.ZipArchiveMode]::Create)
        try {
            $bslash = [char]92
            $fslash = [char]47
            Get-ChildItem -LiteralPath $stageRoot -Recurse -File | ForEach-Object {
                $relPath = $_.FullName.Substring($stageRoot.Length + 1).Replace($bslash, $fslash)
                [System.IO.Compression.ZipFileExtensions]::CreateEntryFromFile($archive, $_.FullName, $relPath) | Out-Null
            }
        } finally {
            $archive.Dispose()
        }
    } finally {
        $fs.Dispose()
    }
    [System.IO.File]::Copy($tmpArtifact, $artifactPath, $true)
    Remove-Item -LiteralPath $tmpArtifact -Force -ErrorAction SilentlyContinue
} finally {
    if (Test-Path -LiteralPath $stageRoot) {
        Remove-Item -LiteralPath $stageRoot -Recurse -Force
    }
}

& (Join-Path $PSScriptRoot "verify_release_artifact.ps1") -RepoRoot $RepoRoot -ArtifactPath $artifactPath
if (-not $?) {
    exit 1
}

Write-Host "build_package.ps1: wrote $artifactPath"
