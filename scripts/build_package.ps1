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
        "llms.txt"
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
        "examples\quaids_coefficients.tex",
        "examples\quaids_coefficients.md",
        "examples\quaids_coefficients.csv",
        "examples\quaids_income_elasticities.md",
        "examples\quaids_uncompensated_elasticities.tex",
        "examples\quaids_compensated_elasticities.csv"
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

    $tmpArtifact = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_artifact_" + [System.Guid]::NewGuid().ToString("N") + ".zip")
    Compress-Archive -Path (Join-Path $stageRoot "*") -DestinationPath $tmpArtifact -Force
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
