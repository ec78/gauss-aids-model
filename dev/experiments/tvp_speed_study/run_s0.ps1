<#
.SYNOPSIS
  TVP speed study Phase S0 driver (analytic-gradient feasibility and speed
  pilot). See PROTOCOL.md.
.DESCRIPTION
  Runs the frozen pilot stages strictly sequentially, one GAUSS process at
  a time: A (n1=2), B (n1=3), then C (n1=4) only if the frozen inclusion
  rule (decided after C's preflight and derivative job, before any C
  fit) projects it safely. Before every stage: frozen file hashes (this
  study and the Phase C harness it reuses), the unchanged Phase C
  evidence file, clean recorded HEAD (only results\ may differ), no other
  GAUSS computation running (Test-OtherGauss), strict sslib pin, TVP Stage 6 30/30, P5 12/12 and
  the S0 harness-equivalence gate (s0_equiv.e). Each stage then runs its
  derivative-gate jobs, then its fit cells; every gate is checked as soon
  as its job finishes and any failure halts the study (no settings are
  changed and rerun). -DryRun prints the plan. -Freeze writes
  FREEZE.sha256 and prints the freeze hash. -Smoke (requires -ResultsDir
  outside this directory, -Stage A) replaces Stage A with one Stage 6
  fixture derivative job and one fixture fit cell, skips the HEAD/clean
  check and E0 (no Phase C row exists for the fixture); smoke output is
  never study evidence.
#>
param(
    [ValidateSet("ALL", "A", "B", "C")]
    [string]$Stage = "ALL",
    [switch]$DryRun,
    [switch]$Freeze,
    [switch]$Smoke,
    [string]$ResultsDir = "",
    [string]$GaussExe = "C:\gauss26\tgauss.exe"
)

$ErrorActionPreference = "Stop"
$inv = [System.Globalization.CultureInfo]::InvariantCulture
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $here "..\..\..")).Path
$convDir = Join-Path $here "..\tvp_conv_study"
$cfgDir = Join-Path $repoRoot "tests\gauss26_cfg_override"
$resultsDir = if ($ResultsDir) { $ResultsDir } else { Join-Path $here "results" }
$logDir = Join-Path $resultsDir "logs"
$genDir = Join-Path $resultsDir "gen"
$derivCsv = Join-Path $resultsDir "deriv.csv"
$profCsv = Join-Path $resultsDir "profile.csv"
$runsCsv = Join-Path $resultsDir "runs.csv"
$cellsCsv = Join-Path $resultsDir "cells.csv"
$stateFile = Join-Path $resultsDir "run_state.json"
$runLog = Join-Path $resultsDir "run_log.txt"

# ---------------------------------------------------------------- frozen design
$QSCALE = 0.01
$BUDGET_MIN = 150.0
$SEED = 1001
$MAXITERS = 500
$MAXTIME = 600
$TIMEOUT_SEC = @{ 2 = 900; 3 = 1800; 4 = 2400 }
$DERIV_TIMEOUT_SEC = 900
$C_MIN_FREE_GB = 4.0
$STAGE_ORDER = @("A", "B", "C")
$FROZEN_FILES = @("PROTOCOL.md", "s0_lib.src", "s0_includes.inc", "s0_equiv.e", "s0_deriv.template.e",
    "s0_fit.template.e", "run_s0.ps1", "analyze_s0.ps1")
# Phase C evidence this study reads; its SHA-256 (CR removed) at commit cbb3130.
$PHASEC_SUMMARY = Join-Path $convDir "results\summary.csv"
$PHASEC_SUMMARY_SHA = "6e1a5647c14f514a305c3bbba3d2d5d13c475929da9ad08d7a3071e108a51f92"
$GATE_DERIV_ERR = 1e-5      # max |gA - gFD| / max(||gFD||_inf, 1), step 1e-4*p
$GATE_DERIV_DLL = 1e-12     # |LL(gradient filter) - LL(plain filter)| / max(|LL|, 1)
$FLAG_DQ_REL = 1e-3         # Q-AGREE flag: max|QA - QN| / max|QN|
$FLAG_DB_ABS = 1e-5         # B-AGREE flag: max|bFinalA - bFinalN|
$SPEED_BAND = 1.10          # per-n1 classification band on median MLE-time ratio N/A

function Get-Datasets([string]$stg) {
    if ($Smoke) { return @([pscustomobject]@{ n1 = 2; tobs = 120; harm = 1 }) }
    switch ($stg) {
        "A" { return @([pscustomobject]@{ n1 = 2; tobs = 200; harm = 1 }, [pscustomobject]@{ n1 = 2; tobs = 200; harm = 2 }, [pscustomobject]@{ n1 = 2; tobs = 500; harm = 1 }) }
        "B" { return @([pscustomobject]@{ n1 = 3; tobs = 200; harm = 1 }, [pscustomobject]@{ n1 = 3; tobs = 200; harm = 2 }) }
        "C" { return @([pscustomobject]@{ n1 = 4; tobs = 200; harm = 1 }) }
    }
}
function Get-Cells([string]$stg) {
    $c = @()
    if ($Smoke) { return @([pscustomobject]@{ n1 = 2; tobs = 120; harm = 1; qmult = 1 }) }
    switch ($stg) {
        "A" {
            foreach ($q in @(1, 10, 0.1)) { $c += [pscustomobject]@{ n1 = 2; tobs = 200; harm = 1; qmult = $q } }
            $c += [pscustomobject]@{ n1 = 2; tobs = 200; harm = 2; qmult = 1 }
            $c += [pscustomobject]@{ n1 = 2; tobs = 500; harm = 1; qmult = 1 }
        }
        "B" {
            foreach ($q in @(1, 10, 0.1)) { $c += [pscustomobject]@{ n1 = 3; tobs = 200; harm = 1; qmult = $q } }
            $c += [pscustomobject]@{ n1 = 3; tobs = 200; harm = 2; qmult = 1 }
        }
        "C" { foreach ($q in @(1, 10, 0.1)) { $c += [pscustomobject]@{ n1 = 4; tobs = 200; harm = 1; qmult = $q } } }
    }
    return $c
}

function Get-NormHash([string]$path) {
    $txt = [System.IO.File]::ReadAllText($path) -replace "`r", ""
    $sha = [System.Security.Cryptography.SHA256]::Create()
    return (($sha.ComputeHash([System.Text.Encoding]::UTF8.GetBytes($txt)) | ForEach-Object { $_.ToString("x2") }) -join "")
}
function D([string]$s) { return [double]::Parse($s, $inv) }
function Write-RunLog([string]$msg) {
    $line = "$((Get-Date).ToString('s')) $msg"
    Write-Host $line
    if (-not $DryRun) { $line | Out-File -Encoding utf8 -Append $runLog }
}
function Get-State {
    if (Test-Path $stateFile) {
        $s = Get-Content -Raw $stateFile | ConvertFrom-Json
        $h = @{}
        foreach ($p in $s.PSObject.Properties) { $h[$p.Name] = $p.Value }
        return $h
    }
    return @{ stages_done = @(); stopped = $false }
}
function Save-State($st) { if (-not $DryRun) { ($st | ConvertTo-Json -Depth 5) | Out-File -Encoding utf8 $stateFile } }
function Get-ElapsedMin($st) {
    if (-not $st.budget_start) { return 0.0 }
    return ((Get-Date) - [datetime]::Parse($st.budget_start, $inv)).TotalMinutes
}
function Stop-Study($st, [string]$reason) {
    $st.stopped = $true
    $st.stop_reason = $reason
    Save-State $st
    Write-RunLog "STUDY HALTED: $reason"
    $reason | Out-File -Encoding utf8 (Join-Path $resultsDir "STOP.txt")
    exit 2
}

# One tgauss job under an outer timeout (same mechanism as run_conv.ps1).
function Invoke-GaussJob([string]$scriptRel, [string]$workDir, [int]$timeoutSec) {
    $tmp = Join-Path $resultsDir ("_job_" + [guid]::NewGuid().ToString("N"))
    $outFile = "$tmp.out"; $errFile = "$tmp.err"
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = "cmd.exe"
    $psi.Arguments = "/c `"`"$GaussExe`" -b -x `"$scriptRel`" > `"$outFile`" 2> `"$errFile`"`""
    $psi.UseShellExecute = $false
    $psi.WorkingDirectory = $workDir
    $psi.EnvironmentVariables["GAUSS26_CFG"] = $cfgDir
    $started = Get-Date
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $proc = [System.Diagnostics.Process]::Start($psi)
    $timedOut = $false
    if (-not $proc.WaitForExit($timeoutSec * 1000)) {
        $timedOut = $true
        & taskkill.exe /PID $proc.Id /T /F 2>&1 | Out-Null
        [void]$proc.WaitForExit(30000)
    }
    $sw.Stop()
    $text = ""
    if (Test-Path $outFile) { $text += (Get-Content -Raw $outFile) }
    if (Test-Path $errFile) { $e = Get-Content -Raw $errFile; if ($e) { $text += "`n--- stderr ---`n" + $e } }
    Remove-Item -ErrorAction SilentlyContinue $outFile, $errFile
    return [pscustomobject]@{ text = $text; timedOut = $timedOut; wall = [math]::Round($sw.Elapsed.TotalSeconds, 2); started = $started.ToString('s') }
}
function Get-Tag([string]$text, [string]$tag) {
    $m = [regex]::Match($text, "(?m)^$tag[ \t]+(\S+)")
    if ($m.Success) { return $m.Groups[1].Value } else { return "" }
}
function Test-GaussError([string]$text) { return [regex]::IsMatch($text, '(?i)\bG\d{4}\b|error G|Undefined symbol') }
# Another GAUSS computation: any tgauss.exe, or an interactive gauss.exe
# that used more than 0.3 s of CPU over a 3 s window (an idle IDE is allowed).
function Test-OtherGauss {
    if (@(Get-Process -ErrorAction SilentlyContinue -Name tgauss).Count -gt 0) { return $true }
    $g = @(Get-Process -ErrorAction SilentlyContinue -Name gauss)
    if ($g.Count -eq 0) { return $false }
    $c0 = ($g | Measure-Object -Property CPU -Sum).Sum
    Start-Sleep -Seconds 3
    $c1 = (@(Get-Process -ErrorAction SilentlyContinue -Name gauss) | Measure-Object -Property CPU -Sum).Sum
    return ($c1 - $c0) -gt 0.3
}

function Get-DatasetId($d) { $h = if ($d.harm -eq 1) { "Horacle" } else { "Hstatic" }; return "n1$($d.n1)-t$($d.tobs)-s$SEED-$h" }
function Get-CellId($c) { return "$(Get-DatasetId $c)-q$($c.qmult.ToString($inv))" }

# ---------------------------------------------------------------- Phase C reference
$phaseC = @()
if (Test-Path $PHASEC_SUMMARY) { $phaseC = @(Import-Csv $PHASEC_SUMMARY | Where-Object { $_.replicate -eq "r1" -and $_.mode -eq "1" }) }
function Get-PhaseCRow($c) {
    return $phaseC | Where-Object { [int]$_.n1 -eq $c.n1 -and [int]$_.tobs -eq $c.tobs -and [int]$_.seed_study -eq $SEED -and [int]$_.harm -eq $c.harm -and (D $_.qmult) -eq $c.qmult } | Select-Object -First 1
}
function Get-PhaseCBest($c) {
    $rows = $phaseC | Where-Object { [int]$_.n1 -eq $c.n1 -and [int]$_.tobs -eq $c.tobs -and [int]$_.seed_study -eq $SEED -and [int]$_.harm -eq $c.harm }
    return ($rows | ForEach-Object { D $_.llfit } | Measure-Object -Maximum).Maximum
}
# Phase C category rules (analyze_conv.ps1), against the Phase C best-found LL.
function Get-Category([double]$rc, [double]$ll, [double]$llTrue, [double]$best, [double]$tol) {
    if ($rc -eq 0) {
        if ($ll -lt $llTrue - $tol) { return "CONV-FAILC3" }
        if ($ll -lt $best - $tol) { return "CONV-LOCAL" }
        return "CONV-GLOBAL"
    }
    if ($rc -eq 2) { return "CAPPED" }
    return "OTHER-RC($rc)"
}

# ---------------------------------------------------------------- projection
# Projected seconds of one fit-cell job: 3 numerical fits at 1.5x the Phase C
# wall time of that cell, plus 3 analytic fits at 1.5 x (Phase C iterations)
# x (k_states + 5) x (seconds per gradient evaluation, from this stage's
# derivative job for that dataset).
function Get-JobProjectionSec($c, $st) {
    if ($Smoke) { return 0.0 }
    $pc = Get-PhaseCRow $c
    $k = 2 * $c.n1 + $c.n1 * ($c.n1 + 1) / 2
    $sg = $st.("secgrad_" + (Get-DatasetId $c))
    if (-not $sg) { $sg = 0.0 }
    return 3 * 1.5 * (D $pc.wall_s) + 3 * 1.5 * (D $pc.iters) * ($k + 5) * [double]$sg
}

# ---------------------------------------------------------------- preflight
function Invoke-Preflight([string]$stg, $st) {
    Write-RunLog "[$stg] preflight start"
    $manifest = Get-Content (Join-Path $here "FREEZE.sha256")
    foreach ($line in $manifest) {
        if ($line -match '^([0-9a-f]{64})\s+(\S+)$') {
            if ((Get-NormHash (Join-Path $here $Matches[2])) -ne $Matches[1]) { Stop-Study $st "[$stg] preflight: frozen file $($Matches[2]) changed" }
        }
    }
    foreach ($f in $FROZEN_FILES) { if (-not ($manifest -match [regex]::Escape($f))) { Stop-Study $st "[$stg] preflight: $f missing from FREEZE.sha256" } }
    foreach ($line in (Get-Content (Join-Path $convDir "FREEZE.sha256"))) {
        if ($line -match '^([0-9a-f]{64})\s+(\S+)$') {
            if ((Get-NormHash (Join-Path $convDir $Matches[2])) -ne $Matches[1]) { Stop-Study $st "[$stg] preflight: Phase C frozen file $($Matches[2]) changed" }
        }
    }
    if ((Get-NormHash $PHASEC_SUMMARY) -ne $PHASEC_SUMMARY_SHA) { Stop-Study $st "[$stg] preflight: Phase C summary.csv changed" }
    Push-Location $repoRoot
    try {
        $head = (& git rev-parse HEAD).Trim()
        $dirty = @(& git status --porcelain --untracked-files=all | Where-Object { $_ -notmatch 'dev/experiments/tvp_speed_study/results/' })
    } finally { Pop-Location }
    if (-not $st.freeze_head) { $st.freeze_head = $head; Save-State $st }
    if ($head -ne $st.freeze_head) { Stop-Study $st "[$stg] preflight: HEAD $head differs from recorded $($st.freeze_head)" }
    if ($dirty.Count -gt 0 -and -not $Smoke) { Stop-Study $st "[$stg] preflight: working tree not clean outside results/: $($dirty -join '; ')" }
    if (Test-OtherGauss) { Stop-Study $st "[$stg] preflight: another GAUSS process is running" }
    & powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repoRoot "scripts\verify_sslib_pin.ps1") -Strict *> (Join-Path $logDir "preflight_$($stg)_sslib_pin.log")
    if ($LASTEXITCODE -ne 0) { Stop-Study $st "[$stg] preflight: verify_sslib_pin.ps1 -Strict failed" }
    $tests = Join-Path $repoRoot "tests"
    foreach ($t in @(@("quaidstvpfit_test.e", 30), @("quaidstvp_loglik_xcheck_test.e", 12))) {
        $r = Invoke-GaussJob $t[0] $tests 900
        $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_$($t[0]).log")
        if ($r.timedOut -or $r.text -notmatch "ALL $($t[1]) CHECKS PASSED") { Stop-Study $st "[$stg] preflight: $($t[0]) did not report ALL $($t[1]) CHECKS PASSED" }
    }
    $r = Invoke-GaussJob "s0_equiv.e" $here 900
    $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_s0_equiv.log")
    $ok = (Get-Tag $r.text "S0_EQ_PASS") -ne "" -and (D (Get-Tag $r.text "S0_EQ_PASS")) -eq 1 -and $r.text -match '(?m)^S0_EQ_DONE'
    if ($r.timedOut -or -not $ok -or (Test-GaussError $r.text)) { Stop-Study $st "[$stg] preflight: S0 harness equivalence gate failed" }
    Write-RunLog "[$stg] preflight PASS (HEAD $head, frozen hashes, Phase C harness+evidence unchanged, no other GAUSS, pin strict, Stage6 30/30, P5 12/12, numerical arm == quaidsTVPFit exactly)"
}

# ---------------------------------------------------------------- derivative gate
function Invoke-Deriv($d, [string]$stg, $st) {
    $id = "deriv-" + (Get-DatasetId $d)
    $tmpl = Get-Content -Raw (Join-Path $here "s0_deriv.template.e")
    $src = $tmpl -replace '__FIXTURE__', $(if ($Smoke) { '1' } else { '0' }) -replace '__N1__', $d.n1 -replace '__TOBS__', $d.tobs -replace '__SEED__', $SEED `
        -replace '__HARM__', $d.harm -replace '__QSCALE__', $QSCALE.ToString($inv)
    [System.IO.File]::WriteAllText((Join-Path $genDir "$id.e"), $src)
    if (Test-OtherGauss) { Stop-Study $st "[$stg] another GAUSS process is running before $id" }
    $r = Invoke-GaussJob (Join-Path $genDir "$id.e") $here $DERIV_TIMEOUT_SEC
    "# job=$id head=$($st.freeze_head) wall_s=$($r.wall) timed_out=$($r.timedOut) started=$($r.started)`n$($r.text)" | Out-File -Encoding utf8 (Join-Path $logDir "$id.log")
    $t = $r.text
    if ($r.timedOut) { Stop-Study $st "[$stg] outer timeout in $id" }
    if ((Test-GaussError $t) -or $t -notmatch '(?m)^S0D_DONE') { Stop-Study $st "[$stg] GAUSS error in $id" }
    $blocks = [regex]::Split($t, '(?m)^S0D_POINT ') | Select-Object -Skip 1
    $allOk = $true
    foreach ($b in $blocks) {
        $pt = ($b -split "`n")[0].Trim()
        $err4 = D (Get-Tag $b "S0D_ERR4"); $dll = D (Get-Tag $b "S0D_DLL_REL")
        $ok = ($err4 -le $GATE_DERIV_ERR) -and ($dll -le $GATE_DERIV_DLL)
        if (-not $ok) { $allOk = $false }
        $vals = @((Get-DatasetId $d), $stg, $pt, (Get-Tag $b "S0D_LLPLAIN"), (Get-Tag $b "S0D_LLGRAD"), (Get-Tag $b "S0D_DLL_REL"), (Get-Tag $b "S0D_GNORM"),
            (Get-Tag $b "S0D_ERR4"), (Get-Tag $b "S0D_ERR3"), $ok)
        ($vals -join ",") | Out-File -Encoding ascii -Append $derivCsv
    }
    if ($blocks.Count -ne 5) { Stop-Study $st "[$stg] $id reported $($blocks.Count) points, expected 5" }
    $sp = Get-Tag $t "S0D_SEC_PLAIN_EVAL"; $sg = Get-Tag $t "S0D_SEC_GRAD_EVAL"
    (@((Get-DatasetId $d), $stg, $d.n1, $d.tobs, $d.harm, (Get-Tag $t "S0D_MINW"), $sp, $sg, (Get-Tag $t "S0D_SEC_QUPDATE"), (Get-Tag $t "S0D_GRAD_OVER_PLAIN")) -join ",") | Out-File -Encoding ascii -Append $profCsv
    $st.("secgrad_" + (Get-DatasetId $d)) = (D $sg); Save-State $st
    Write-RunLog "[$stg] $id wall=$($r.wall)s derivative gate=$allOk sec/plain=$sp sec/grad=$sg"
    if (-not $allOk) { Stop-Study $st "[$stg] derivative gate FAILED in $id" }
    if ((D (Get-Tag $t "S0D_MINW")) -le 0) { Stop-Study $st "[$stg] $id dataset not admissible (share <= 0)" }
}

# ---------------------------------------------------------------- fit cell
function Invoke-Cell($c, [string]$stg, $st) {
    $id = Get-CellId $c
    $pc = if ($Smoke) { [pscustomobject]@{ rc = "0"; iters = "0"; llfit = "0"; lltrue = "0"; wall_s = "0"; qfit_sum = "0" } } else { Get-PhaseCRow $c }
    if (-not $pc) { Stop-Study $st "[$stg] no Phase C reference row for $id" }
    $timeout = $TIMEOUT_SEC[[int]$c.n1]
    $proj = Get-JobProjectionSec $c $st
    if ((Get-ElapsedMin $st) + $proj / 60.0 -gt $BUDGET_MIN) { Stop-Study $st "[$stg] budget: $id projected $([math]::Round($proj/60,1)) min would exceed $BUDGET_MIN min" }
    $tmpl = Get-Content -Raw (Join-Path $here "s0_fit.template.e")
    $src = $tmpl -replace '__N1__', $c.n1 -replace '__TOBS__', $c.tobs -replace '__SEED__', $SEED -replace '__HARM__', $c.harm `
        -replace '__QMULT__', $c.qmult.ToString($inv) -replace '__MAXITERS__', $MAXITERS -replace '__MAXTIME__', $MAXTIME -replace '__QSCALE__', $QSCALE.ToString($inv)
    if ($Smoke) { $src = $src -replace '(?m)^seed = 1001;', 'seed = 7;' -replace '(?m)^qscale = 0.01;', 'qscale = 1;' }
    [System.IO.File]::WriteAllText((Join-Path $genDir "fit-$id.e"), $src)
    if (Test-OtherGauss) { Stop-Study $st "[$stg] another GAUSS process is running before $id" }
    $r = Invoke-GaussJob (Join-Path $genDir "fit-$id.e") $here $timeout
    $pin = (Get-Content -Raw (Join-Path $repoRoot 'sslib.pin.json') | ConvertFrom-Json).pinnedCommit
    "# cell=$id stage=$stg head=$($st.freeze_head) outer_timeout_s=$timeout timed_out=$($r.timedOut) wall_s=$($r.wall) sslib_pin=$pin started=$($r.started) projected_s=$([math]::Round($proj,1))`n$($r.text)" |
        Out-File -Encoding utf8 (Join-Path $logDir "fit-$id.log")
    $t = $r.text
    if ($r.timedOut) { Stop-Study $st "[$stg] outer timeout in $id" }
    if ((Test-GaussError $t) -or $t -notmatch '(?m)^S0F_DONE') { Stop-Study $st "[$stg] GAUSS error in $id" }

    # Run rows: arm rep rc iters nUpd nGrad mleSec chainSec fct llFit maxgrad qfitSum
    $runs = @()
    foreach ($m in [regex]::Matches($t, '(?m)^S0F_RUN\s+(.+)$')) {
        $v = @($m.Groups[1].Value.Trim() -split '\s+' | ForEach-Object { D $_ })
        $runs += [pscustomobject]@{ arm = $(if ($v[0] -eq 0) { "N" } else { "A" }); rep = [int]$v[1]; rc = $v[2]; iters = $v[3]; nupd = $v[4]; ngrad = $v[5]
            mle_s = $v[6]; chain_s = $v[7]; fct = $v[8]; llfit = $v[9]; maxgrad = $v[10]; qsum = $v[11] }
    }
    if ($runs.Count -ne 6) { Stop-Study $st "[$stg] $id reported $($runs.Count) runs, expected 6" }
    foreach ($x in $runs) {
        (@($id, $stg, $x.arm, $x.rep, $x.rc, $x.iters, $x.nupd, $x.ngrad, $x.mle_s.ToString("R", $inv), $x.chain_s.ToString("R", $inv),
            $x.fct.ToString("R", $inv), $x.llfit.ToString("R", $inv), $x.maxgrad.ToString("R", $inv), $x.qsum.ToString("R", $inv)) -join ",") | Out-File -Encoding ascii -Append $runsCsv
    }
    $n1r = $runs | Where-Object { $_.arm -eq "N" -and $_.rep -eq 1 }
    $a1r = $runs | Where-Object { $_.arm -eq "A" -and $_.rep -eq 1 }
    $tol = D (Get-Tag $t "S0F_TOL"); $llTrue = D (Get-Tag $t "S0F_LLTRUE")
    $llN = D (Get-Tag $t "S0F_LLN"); $llA = D (Get-Tag $t "S0F_LLA")
    $best = if ($Smoke) { [double]::NegativeInfinity } else { Get-PhaseCBest $c }
    if ($Smoke) { $pc = [pscustomobject]@{ rc = $n1r.rc.ToString("R", $inv); iters = $n1r.iters.ToString("R", $inv); llfit = $llN.ToString("R", $inv); lltrue = $llTrue.ToString("R", $inv); qfit_sum = $n1r.qsum.ToString("R", $inv) } }
    $catC = Get-Category (D $pc.rc) (D $pc.llfit) (D $pc.lltrue) $best $tol
    $catN = Get-Category $n1r.rc $llN $llTrue $best $tol
    $catA = Get-Category $a1r.rc $llA $llTrue $best $tol
    # E0: numerical arm reproduces the Phase C cell.
    $e0 = ($n1r.rc -eq (D $pc.rc)) -and ($n1r.iters -eq (D $pc.iters)) -and ([math]::Abs($llN - (D $pc.llfit)) -le $tol) -and ($catN -eq $catC)
    $e0exact = ($llN -eq (D $pc.llfit)) -and ($n1r.qsum -eq (D $pc.qfit_sum))
    # E1: same retcode and category; E2: LL within tolerance; E3: determinism.
    $e1 = ($a1r.rc -eq $n1r.rc) -and ($catA -eq $catN)
    $e2 = [math]::Abs($llA - $llN) -le $tol
    $e3 = ((D (Get-Tag $t "S0F_DET_N")) -eq 1) -and ((D (Get-Tag $t "S0F_DET_A")) -eq 1)
    $dqRel = D (Get-Tag $t "S0F_DQ_REL"); $dbAbs = D (Get-Tag $t "S0F_DB_ABS")
    $qAgree = $dqRel -le $FLAG_DQ_REL; $bAgree = $dbAbs -le $FLAG_DB_ABS
    $medN = (@($runs | Where-Object { $_.arm -eq "N" } | ForEach-Object { $_.mle_s } | Sort-Object))[1]
    $medA = (@($runs | Where-Object { $_.arm -eq "A" } | ForEach-Object { $_.mle_s } | Sort-Object))[1]
    $minN = (@($runs | Where-Object { $_.arm -eq "N" } | ForEach-Object { $_.mle_s } | Sort-Object))[0]
    $minA = (@($runs | Where-Object { $_.arm -eq "A" } | ForEach-Object { $_.mle_s } | Sort-Object))[0]
    $ratio = if ($medA -gt 0) { $medN / $medA } else { [double]::NaN }
    $vals = @($id, $stg, $c.n1, $c.tobs, $c.harm, $c.qmult.ToString($inv), $r.wall.ToString($inv), $catC, $catN, $catA, $n1r.rc, $a1r.rc, $n1r.iters, $a1r.iters,
        $n1r.nupd, $n1r.ngrad, $a1r.nupd, $a1r.ngrad, $llN.ToString("R", $inv), $llA.ToString("R", $inv), ($llA - $llN).ToString("R", $inv), $tol.ToString("R", $inv),
        $best.ToString("R", $inv), $llTrue.ToString("R", $inv), $dqRel.ToString("R", $inv), $dbAbs.ToString("R", $inv), $e0, $e0exact, $e1, $e2, $e3, $qAgree, $bAgree,
        $medN.ToString("R", $inv), $medA.ToString("R", $inv), $minN.ToString("R", $inv), $minA.ToString("R", $inv), $ratio.ToString("R", $inv), "logs/fit-$id.log")
    ($vals -join ",") | Out-File -Encoding ascii -Append $cellsCsv
    Write-RunLog "[$stg] $id wall=$($r.wall)s cat C/N/A=$catC/$catN/$catA iters N/A=$($n1r.iters)/$($a1r.iters) dLL=$($llA - $llN) E0=$e0 (exact=$e0exact) E1=$e1 E2=$e2 E3=$e3 Qagree=$qAgree Bagree=$bAgree medMLE N/A=$medN/$medA ratio=$([math]::Round($ratio,3))"
    if (-not $e0) { Stop-Study $st "[$stg] E0 (numerical arm reproduces Phase C) FAILED in $id" }
    if (-not ($e1 -and $e2 -and $e3)) { Stop-Study $st "[$stg] endpoint gate FAILED in $id (E1=$e1 E2=$e2 E3=$e3)" }
}

# ---------------------------------------------------------------- stage C inclusion rule
# Decided after Stage C's preflight and derivative job (so the measured
# seconds per gradient evaluation at n1 = 4 is known) and before any C fit.
function Test-IncludeC($st) {
    $el = Get-ElapsedMin $st
    $freeGB = [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB, 2)
    $proj = 0.0
    foreach ($c in (Get-Cells "C")) {
        $p = Get-JobProjectionSec $c $st
        $proj += $p
        if ($p -gt 0.5 * $TIMEOUT_SEC[4]) { Write-RunLog "[C inclusion] EXCLUDED: $(Get-CellId $c) projected $([math]::Round($p,0)) s > half its outer timeout"; return $false }
    }
    Write-RunLog "[C inclusion] elapsed=$([math]::Round($el,1)) min, projected C fits=$([math]::Round($proj/60,1)) min, free memory=$freeGB GB"
    if ($el + $proj / 60.0 -gt $BUDGET_MIN) { Write-RunLog "[C inclusion] EXCLUDED: budget"; return $false }
    if ($freeGB -lt $C_MIN_FREE_GB) { Write-RunLog "[C inclusion] EXCLUDED: free memory $freeGB GB < $C_MIN_FREE_GB GB"; return $false }
    Write-RunLog "[C inclusion] INCLUDED"
    return $true
}

# ---------------------------------------------------------------- stage runner
function Invoke-Stage([string]$stg, $st) {
    if ($st.stopped) { Write-RunLog "study already halted: $($st.stop_reason)"; exit 2 }
    if (@($st.stages_done) -contains $stg) { Write-RunLog "[$stg] already done; skipping"; return }
    if (-not $st.budget_start) { $st.budget_start = (Get-Date).ToString('o'); Save-State $st }
    Invoke-Preflight $stg $st
    foreach ($d in (Get-Datasets $stg)) { Invoke-Deriv $d $stg $st }
    if ($stg -eq "C") {
        if (-not (Test-IncludeC $st)) { $st.c_excluded = $true; $st.stages_done = @(@($st.stages_done) + $stg); Save-State $st; return }
    }
    foreach ($c in (Get-Cells $stg)) { Invoke-Cell $c $stg $st }
    $st.stages_done = @(@($st.stages_done) + $stg)
    Save-State $st
    Write-RunLog "[$stg] stage complete; elapsed $([math]::Round((Get-ElapsedMin $st),1)) min"
}

# ---------------------------------------------------------------- main
if ($Freeze) {
    $lines = foreach ($f in $FROZEN_FILES) { "$(Get-NormHash (Join-Path $here $f))  $f" }
    ($lines -join "`n") + "`n" | Out-File -Encoding ascii -NoNewline (Join-Path $here "FREEZE.sha256")
    Write-Host ($lines -join "`n")
    Write-Host "FREEZE HASH $(Get-NormHash (Join-Path $here 'FREEZE.sha256'))"
    exit 0
}
if ($DryRun) {
    foreach ($s in $STAGE_ORDER) {
        Write-Host "Stage $s : datasets $((Get-Datasets $s | ForEach-Object { Get-DatasetId $_ }) -join ', ')"
        foreach ($c in (Get-Cells $s)) { $pc = Get-PhaseCRow $c; Write-Host "   $(Get-CellId $c)  PhaseC rc=$($pc.rc) iters=$($pc.iters) wall=$($pc.wall_s)s llfit=$($pc.llfit)" }
    }
    Write-Host "Phase C summary sha (CR removed): $(Get-NormHash $PHASEC_SUMMARY)"
    exit 0
}
if ($Smoke -and (-not $ResultsDir -or $Stage -ne "A")) { throw "-Smoke requires -ResultsDir (outside this directory) and -Stage A" }
New-Item -ItemType Directory -Force $resultsDir, $logDir, $genDir | Out-Null
$st = Get-State
if (-not (Test-Path $derivCsv)) { "dataset,stage,point,ll_plain,ll_grad,dll_rel,gnorm_inf,err_h1e4,err_h1e3,gate_ok" | Out-File -Encoding ascii $derivCsv }
if (-not (Test-Path $profCsv)) { "dataset,stage,n1,tobs,harm,minw,sec_plain_eval,sec_grad_eval,sec_qupdate,grad_over_plain" | Out-File -Encoding ascii $profCsv }
if (-not (Test-Path $runsCsv)) { "cell_id,stage,arm,rep,rc,iters,n_obj_eval,n_grad_eval,mle_s,chain_s,fct,llfit,maxgrad,qfit_sum" | Out-File -Encoding ascii $runsCsv }
if (-not (Test-Path $cellsCsv)) { "cell_id,stage,n1,tobs,harm,qmult,job_wall_s,cat_phaseC,cat_N,cat_A,rc_N,rc_A,iters_N,iters_A,nobj_N,ngrad_N,nobj_A,ngrad_A,ll_N,ll_A,dll_A_minus_N,tol,best_phaseC,ll_true,dq_rel,db_abs,E0,E0_exact,E1,E2,E3,Q_agree,B_agree,med_mle_N,med_mle_A,min_mle_N,min_mle_A,ratio_med_N_over_A,log" | Out-File -Encoding ascii $cellsCsv }
$todo = if ($Stage -eq "ALL") { $STAGE_ORDER } else { @($Stage) }
foreach ($s in $todo) { Invoke-Stage $s $st }
Write-RunLog "run finished: stages done = $(@($st.stages_done) -join ',') C excluded=$([bool]$st.c_excluded) elapsed=$([math]::Round((Get-ElapsedMin $st),1)) min"
