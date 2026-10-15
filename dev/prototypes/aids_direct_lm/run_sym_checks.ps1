# run_sym_checks.ps1 -- implementation checks I1-I9 of SYMMETRY_PROTOCOL.md v1
# (this directory). PRIVATE EXPERIMENT ONLY; solves fixture data only.
#
#   I5 (byte identity)  quaidsdirect.src unchanged since d8b4318
#   I1-I5, I7, I8       quaidsdirectsym_test.e (tgauss, from tests/)
#   I8 (guards)         the five guard_error_cases/directlm_sym_*.e scripts
#   I6 + analyzer       python -m unittest test_sym_analyze
#   existing prototype  quaidsdirect_test.e and the five directlm_*.e guards
#                       (-WithEvidence adds quaidsdirect_evidence_test.e,
#                       which needs optmt; -Full implies -WithEvidence,
#                       since protocol I9 names the evidence test)
#   I9 (-Full only)     tests\run_source_tests.ps1 -SkipTVPKalman
#                       -SkipBootstrap, which also re-runs quaidsFit()'s own
#                       tests (quaidsFit() output unchanged). -SkipTVPKalman
#                       drops every SSLib-dependent TVP test, including
#                       tests\quaidstvpfit_test.e.
#
# Only a -Full pass with every check green writes the stamp
# <OutRoot>\IMPLEMENTATION_CHECKS_PASSED.txt (SHA-256 of every freeze file),
# which run_sym_stage.ps1 requires before stage D (protocol section 8 step 2).
# Every -Full run first deletes any existing stamp, so a failed or
# interrupted -Full run never leaves an earlier stamp in force.
#
#   wrapper preflight   run_sym_stage.ps1 -Preflight with the default and an
#                       explicit -OutRoot (path resolution only; no solve)
#
# Usage: powershell -ExecutionPolicy Bypass -File run_sym_checks.ps1 [-Full] [-WithEvidence]

param(
    [switch]$Full,
    [switch]$WithEvidence,
    [string]$OutRoot = "",
    [string]$GaussExe = "C:\gauss26\tgauss.exe"
)

$ErrorActionPreference = "Stop"
# defaults resolved after parameter binding (see run_sym_stage.ps1's header)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrEmpty($OutRoot)) { $OutRoot = Join-Path $here "symmetry" }
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$testsDir = Join-Path $repo "tests"
$failed = @()
if ($Full) {
    $WithEvidence = $true
    $stampPath = Join-Path $OutRoot "IMPLEMENTATION_CHECKS_PASSED.txt"
    if (Test-Path $stampPath) {
        Remove-Item -LiteralPath $stampPath -Force
        Write-Host "removed the existing stamp; only this run's full pass can rewrite it"
    }
}

function Run-Gauss([string]$script) {
    Push-Location $testsDir
    try {
        $ErrorActionPreference = "Continue"
        $out = & $GaussExe -b -x $script 2>&1 | Out-String
    } finally {
        Pop-Location
    }
    return $out
}

# Wrapper preflight (no solve, no directory creation): run_sym_stage.ps1's
# default -OutRoot resolves, from any working directory, to the symmetry\
# directory next to it, and an explicit -OutRoot is used as given. The
# default symmetry\ tree must be untouched by both calls.
$stageScript = Join-Path $here "run_sym_stage.ps1"
$symDefault = Join-Path $here "symmetry"
$probeRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("sym_preflight_probe_" + [guid]::NewGuid().ToString("N"))
function Tree-State([string]$dir) {
    if (-not (Test-Path $dir)) { return "<absent>" }
    return (Get-ChildItem -LiteralPath $dir -Recurse -Force |
        ForEach-Object { "$($_.FullName)|$($_.Length)|$($_.LastWriteTimeUtc.Ticks)" }) -join "`n"
}
$before = Tree-State $symDefault
$cases = @(
    @{ Name = "default OutRoot (cwd tests)"; Cwd = $testsDir; Extra = @(); Expect = $symDefault },
    @{ Name = "default OutRoot (cwd repo)"; Cwd = $repo; Extra = @(); Expect = $symDefault },
    @{ Name = "explicit OutRoot"; Cwd = $testsDir; Extra = @("-OutRoot", $probeRoot); Expect = $probeRoot }
)
foreach ($c in $cases) {
    Push-Location $c.Cwd
    try {
        $ErrorActionPreference = "Continue"
        $pf = & powershell -NoProfile -ExecutionPolicy Bypass -File $stageScript -Stage dev -Preflight @($c.Extra) 2>&1 | Out-String
        $pfCode = $LASTEXITCODE
    } finally {
        Pop-Location
        $ErrorActionPreference = "Stop"
    }
    $want = "PREFLIGHT OutRoot=$($c.Expect)"
    $wantDir = "PREFLIGHT outDir=$(Join-Path $c.Expect 'dev')"
    if ($pfCode -eq 0 -and $pf.Contains($want) -and $pf.Contains($wantDir)) {
        Write-Host "PASS  preflight: $($c.Name) -> $($c.Expect)"
    } else {
        Write-Host "FAIL  preflight: $($c.Name)"; Write-Host $pf; $failed += "preflight $($c.Name)"
    }
}
if ((Tree-State $symDefault) -ne $before) {
    Write-Host "FAIL  preflight: the default symmetry tree changed"; $failed += "preflight changed outputs"
} elseif (Test-Path $probeRoot) {
    Write-Host "FAIL  preflight: created the explicit OutRoot"; $failed += "preflight created a directory"
} else {
    Write-Host "PASS  preflight: no directory created and the symmetry tree is unchanged"
}

# I5: byte identity of the H prototype
& git -C $repo diff --quiet d8b4318 -- dev/prototypes/aids_direct_lm/quaidsdirect.src
if ($LASTEXITCODE -eq 0) { Write-Host "PASS  I5: quaidsdirect.src byte-identical to d8b4318" }
else { Write-Host "FAIL  I5: quaidsdirect.src differs from d8b4318"; $failed += "I5 byte identity" }

$tests = @("..\dev\prototypes\aids_direct_lm\quaidsdirectsym_test.e",
           "..\dev\prototypes\aids_direct_lm\quaidsdirect_test.e")
if ($WithEvidence) { $tests += "..\dev\prototypes\aids_direct_lm\quaidsdirect_evidence_test.e" }
foreach ($t in $tests) {
    $out = Run-Gauss $t
    $m = [regex]::Match($out, "ALL (\d+) CHECKS PASSED")
    if ($m.Success) { Write-Host "PASS  $t ($($m.Groups[1].Value) checks)" }
    else { Write-Host "FAIL  $t"; Write-Host $out; $failed += $t }
}

$guards = @(
    @{ S = "directlm_quaids_rejected.e"; E = "_quaidsDirectFit: aCtl.linear must be 1 (AIDS only; QUAIDS is not supported by this prototype)." },
    @{ S = "directlm_homogenous_zero.e"; E = "_quaidsDirectFit: aCtl.homogenous must be 1 (homogeneity-constrained estimation only)." },
    @{ S = "directlm_user_b0.e"; E = "_quaidsDirectFit: aCtl.b0 must be scalar 0 (user starting values are not supported by this prototype)." },
    @{ S = "directlm_bad_shape.e"; E = "_quaidsDirectFit: inputs must be conformable" },
    @{ S = "directlm_missing_data.e"; E = "_quaidsDirectFit: input data must not contain missing values." }
)
foreach ($g in $guards) {
    foreach ($prefix in @("directlm_", "directlm_sym_")) {
        $name = $g.S -replace "^directlm_", $prefix
        $out = Run-Gauss "..\dev\prototypes\aids_direct_lm\guard_error_cases\$name"
        if ($out.Contains($g.E)) { Write-Host "PASS  guard $name" }
        else { Write-Host "FAIL  guard $name"; Write-Host $out; $failed += "guard $name" }
    }
}

Push-Location $here
try {
    $ErrorActionPreference = "Continue"
    $py = & python -B -m unittest test_sym_analyze 2>&1 | Out-String
    $pyCode = $LASTEXITCODE
} finally {
    Pop-Location
}
if ($pyCode -eq 0 -and $py -match "OK") { Write-Host "PASS  python test_sym_analyze (I6 and analyzer)" }
else { Write-Host "FAIL  python test_sym_analyze"; Write-Host $py; $failed += "python" }

if ($Full) {
    & powershell -ExecutionPolicy Bypass -File (Join-Path $testsDir "run_source_tests.ps1") -SkipTVPKalman -SkipBootstrap
    if ($LASTEXITCODE -eq 0) { Write-Host "PASS  I9: run_source_tests.ps1 -SkipTVPKalman -SkipBootstrap" }
    else { Write-Host "FAIL  I9: run_source_tests.ps1"; $failed += "I9 runner" }
}

if ($failed.Count -gt 0) {
    Write-Host ""
    Write-Host "SYMMETRY IMPLEMENTATION CHECKS FAILED: $($failed -join '; ')"
    exit 1
}
Write-Host ""
if ($Full) {
    $cfg = Get-Content -Raw (Join-Path $here "sym_config.json") | ConvertFrom-Json
    New-Item -ItemType Directory -Force $OutRoot | Out-Null
    $lines = foreach ($rel in $cfg.freeze_files) {
        $h = (Get-FileHash -Algorithm SHA256 (Join-Path $repo $rel)).Hash.ToLower()
        "$h  $rel"
    }
    $lines | Out-File -Encoding ascii (Join-Path $OutRoot "IMPLEMENTATION_CHECKS_PASSED.txt")
    Write-Host "ALL SYMMETRY IMPLEMENTATION CHECKS PASSED (I1-I9); stamp written to $OutRoot"
} else {
    Write-Host "ALL TARGETED SYMMETRY CHECKS PASSED (I9 not run: pass -Full before stage D)"
}
exit 0
