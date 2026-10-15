# run_sym_stage.ps1 -- run one harness stage of SYMMETRY_PROTOCOL.md v1.
# PRIVATE EXPERIMENT ONLY.
#
# Stages (seeds are hard-coded in sym_harness.src, never passed in):
#   fixture       seed 602 (pipeline smoke; implementation-check data)
#   dev           1001-1030 (stage D)          dev_det     1001, 1030
#   conf_interim  2001-2020 (stage C, then the interim analysis)
#   conf_rest     2021-2170 (appends to the same conf\ run)
#   conf_det      2001, 2085, 2170             real        Blanciforti86 (G6)
#
# Refusals (protocol section 8):
#   - an existing, non-empty output directory is never overwritten
#     (conf_rest is the only stage that appends, only after the interim
#     analysis exists, and only once);
#   - dev and dev_det need <OutRoot>\IMPLEMENTATION_CHECKS_PASSED.txt from
#     run_sym_checks.ps1 -Full, with every hash still matching;
#   - conf_*, and real need a verified freeze manifest
#     (sym_analyze.py verify-manifest); real also needs conf_out\results.json
#     (G6 runs after stage C).
#
# Usage: powershell -ExecutionPolicy Bypass -File run_sym_stage.ps1 -Stage dev
#
# -OutRoot defaults to the symmetry\ directory next to this script. The
# default is resolved in the script body, after parameter binding: under
# Windows PowerShell 5.1, $PSScriptRoot is empty inside the parameter
# defaults of an advanced script like this one, which made the original
# default fail before any work (Stage D, 2026-10-13; repaired afterwards,
# numerically irrelevant). -Preflight prints the resolved paths and exits
# before any guard, directory creation or solve (run_sym_checks.ps1 tests
# it).

param(
    [Parameter(Mandatory = $true)]
    [ValidateSet("fixture", "dev", "dev_det", "conf_interim", "conf_rest", "conf_det", "real")]
    [string]$Stage,
    [string]$OutRoot = "",
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$Preflight
)

$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrEmpty($OutRoot)) { $OutRoot = Join-Path $here "symmetry" }
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$testsDir = Join-Path $repo "tests"
$dirName = @{ fixture = "fixture"; dev = "dev"; dev_det = "dev_det"; conf_interim = "conf";
              conf_rest = "conf"; conf_det = "conf_det"; real = "real" }[$Stage]
$outDir = Join-Path $OutRoot $dirName
if ($Preflight) {
    Write-Host "PREFLIGHT OutRoot=$OutRoot"
    Write-Host "PREFLIGHT outDir=$outDir"
    exit 0
}

if ($Stage -in @("dev", "dev_det")) {
    $stamp = Join-Path $OutRoot "IMPLEMENTATION_CHECKS_PASSED.txt"
    if (-not (Test-Path $stamp)) { throw "REFUSED: run run_sym_checks.ps1 -Full first (no $stamp)" }
    foreach ($line in Get-Content $stamp) {
        $parts = $line -split "  ", 2
        $h = (Get-FileHash -Algorithm SHA256 (Join-Path $repo $parts[1])).Hash.ToLower()
        if ($h -ne $parts[0]) { throw "REFUSED: $($parts[1]) changed since the implementation checks passed" }
    }
}
if ($Stage -in @("conf_interim", "conf_rest", "conf_det", "real")) {
    Push-Location $here
    try { & python -B sym_analyze.py verify-manifest; $code = $LASTEXITCODE } finally { Pop-Location }
    if ($code -ne 0) { throw "REFUSED: the freeze manifest does not verify (protocol section 8 step 4)" }
}
if ($Stage -eq "real" -and -not (Test-Path (Join-Path $OutRoot "conf_out\results.json"))) {
    throw "REFUSED: G6 runs after the stage C analysis (no conf_out\results.json)"
}
if ($Stage -eq "conf_rest") {
    if (-not (Test-Path (Join-Path $OutRoot "conf_interim_out\results.json"))) {
        throw "REFUSED: conf_rest needs the interim analysis first (conf_interim_out\results.json)"
    }
    if (Test-Path (Join-Path $outDir "console_conf_rest.txt")) {
        throw "REFUSED: conf_rest has already run in $outDir; stage C is run once (protocol section 8 step 5)"
    }
} elseif ((Test-Path $outDir) -and @(Get-ChildItem $outDir).Count -gt 0) {
    throw "REFUSED: $outDir exists and is not empty; every run's outputs are kept, never overwritten"
}

New-Item -ItemType Directory -Force $outDir | Out-Null
$outFwd = ((Resolve-Path $outDir).Path) -replace "\\", "/"
$driver = Join-Path $outDir "run_$Stage.e"
@"
new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirectsym.src;
#include quaidsfixtures.src;
#include ../dev/prototypes/aids_direct_lm/sym_harness.src;

_symRunStage("$Stage", "$outFwd");
"@ | Out-File -Encoding ascii $driver

Push-Location $testsDir
try {
    $ErrorActionPreference = "Continue"
    $out = & $GaussExe -b -x $driver 2>&1 | Out-String
} finally {
    Pop-Location
}
$out | Out-File -Encoding utf8 (Join-Path $outDir "console_$Stage.txt")
if ($out -notmatch "sym_harness: stage $Stage done") {
    Write-Host $out
    throw "FAILED: stage $Stage did not complete (see console_$Stage.txt)"
}
Write-Host "stage $Stage complete: $outDir"
