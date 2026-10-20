<#
.SYNOPSIS
  TVP-AIDS Q-identification / boundary-geometry study (Phase Q0) driver.
  See PROTOCOL.md.
.DESCRIPTION
  Stages, strictly sequential, one GAUSS process at a time: A (n1=2),
  B (n1=3), ANCHOR (P6 n1=3 anchor), C (n1=4, only if the frozen
  memory/runtime inclusion gate passes). Each endpoint runs: three start
  fits (Phase C's exact cell; each must reproduce Phase C's rc,
  iterations, LL(Qfit) and sum(Qfit) text exactly), one geometry job
  (points, fixed-path slices, log-Q Hessians) and one profile job per
  variance group (true reoptimized profiles). Before every stage:
  frozen hashes, clean recorded HEAD (only results\ may differ), Phase C
  harness files unchanged, strict sslib pin, TVP Stage 6 30/30, P5 12/12,
  q0_equiv.e PASS. Before every job: no other GAUSS computation and at
  least 1 GB free memory (waits up to 10 min), and the job's outer
  timeout must fit in the remaining budget. During a job, the tgauss
  working set is polled and the job is killed above 6 GB. Stop rules
  (halt the study): preflight failure, GAUSS error, outer timeout, memory
  cap, Phase C reproduction mismatch, budget, memory wait exhausted.
  Tagged GAUSS rows are extracted into results\*.csv.
  -DryRun prints the plan. -Freeze writes FREEZE.sha256 and prints the
  freeze hash. -Smoke runs one reserve-seed (1101) n1=2 endpoint end to
  end with -ResultsDir outside this directory (not study evidence).
#>
param(
    [ValidateSet("ALL", "A", "B", "ANCHOR", "C")]
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
$cfgDir = Join-Path $repoRoot "tests\gauss26_cfg_override"
$resultsDir = if ($ResultsDir) { $ResultsDir } else { Join-Path $here "results" }
$logDir = Join-Path $resultsDir "logs"
$genDir = Join-Path $resultsDir "gen"
$workDir = Join-Path $resultsDir "work"
$stateFile = Join-Path $resultsDir "run_state.json"
$runLog = Join-Path $resultsDir "run_log.txt"
$PHASEC_DIR = Join-Path $here "..\tvp_conv_study"
$PHASEC_SUMMARY = Join-Path $PHASEC_DIR "results\summary.csv"

# ---------------------------------------------------------------- frozen design
$QSCALE = 0.01
$PHI = -12                       # log10 floor multiplier on qRef
$HSTEP1 = 0.01; $HSTEP2 = 0.03   # log-Q Hessian steps
$BUDGET_MIN = 300.0
$MIN_FREE_GB = 1.0               # every job; wait up to 10 min
$MEM_WAIT_SEC = 600
$C_MIN_FREE_GB = 4.0             # stage C inclusion
# Stage C worst-case projection (minutes): 3 fits x (Phase C max wall 81 s x 1.5)
# + geometry 5 min + 3 profile jobs x 2 levels x (300+1) iterations x 2.8 s
# (Phase C n1=4 max s/iter) + preflight 5 min, rounded up.
$C_PROJ_MIN = 100.0
$PEAK_WS_CAP_MB = 6144           # kill a tgauss job above this working set
$STARTS = @(1, 10, 0.1)
$STAGE_ORDER = @("A", "B", "ANCHOR", "C")
$FROZEN_FILES = @("PROTOCOL.md", "q0_lib.src", "q0_includes.inc", "q0_data.inc", "q0_design.inc",
    "q0_fit.template.e", "q0_geom.template.e", "q0_prof.template.e", "q0_equiv.e", "run_q0.ps1", "analyze_q0.ps1")
# Phase C files this study reuses unchanged (hashes must equal Phase C's FREEZE.sha256).
$PHASEC_REUSED = @("conv_lib.src", "conv_chain.inc")

# Endpoint and per-job settings by stage.
$LV_DGP = "$PHI|-2|-1|0|1"
$SPEC = @{
    A      = @{ mode = 1; n1 = 2; tobs = 200; harm = 1; seeds = @(1001, 1002); fitIters = 500; levels = $LV_DGP; levelsB = $LV_DGP; profIters = 500; toFit = 600; toGeom = 600; toProf = 1800 }
    B      = @{ mode = 1; n1 = 3; tobs = 200; harm = 1; seeds = @(1001, 1002); fitIters = 500; levels = $LV_DGP; levelsB = $LV_DGP; profIters = 500; toFit = 900; toGeom = 900; toProf = 2700 }
    ANCHOR = @{ mode = 2; n1 = 3; tobs = 500; harm = 2; seeds = @(204); fitIters = 2000; levels = "$PHI|-4|-3|-2|-1"; levelsB = "$PHI"; profIters = 500; toFit = 2400; toGeom = 1800; toProf = 3600 }
    C      = @{ mode = 1; n1 = 4; tobs = 200; harm = 1; seeds = @(1001); fitIters = 500; levels = "$PHI|0"; levelsB = "99"; profIters = 300; toFit = 1800; toGeom = 1800; toProf = 2400 }
}
if ($Smoke) { $SPEC.A.seeds = @(1101) }

function Get-NormHash([string]$path) {
    $txt = [System.IO.File]::ReadAllText($path) -replace "`r", ""
    $sha = [System.Security.Cryptography.SHA256]::Create()
    return (($sha.ComputeHash([System.Text.Encoding]::UTF8.GetBytes($txt)) | ForEach-Object { $_.ToString("x2") }) -join "")
}
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
    if (-not $DryRun) { $reason | Out-File -Encoding utf8 (Join-Path $resultsDir "STOP.txt") }
    exit 2
}
function Get-FreeGB { return [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB, 2) }
function Test-OtherGauss {
    if (@(Get-Process -ErrorAction SilentlyContinue -Name tgauss).Count -gt 0) { return $true }
    $g = @(Get-Process -ErrorAction SilentlyContinue -Name gauss)
    if ($g.Count -eq 0) { return $false }
    $c0 = ($g | Measure-Object -Property CPU -Sum).Sum
    Start-Sleep -Seconds 3
    $c1 = (@(Get-Process -ErrorAction SilentlyContinue -Name gauss) | Measure-Object -Property CPU -Sum).Sum
    return ($c1 - $c0) -gt 0.3
}
# Wait (up to $MEM_WAIT_SEC) until no other GAUSS computation runs and
# free memory is at least $MIN_FREE_GB; returns the free GB, or stops.
function Wait-JobGate($st, [string]$id) {
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    while ($true) {
        $other = Test-OtherGauss
        $free = Get-FreeGB
        if (-not $other -and $free -ge $MIN_FREE_GB) { return $free }
        if ($sw.Elapsed.TotalSeconds -gt $MEM_WAIT_SEC) { Stop-Study $st "job gate for ${id}: other GAUSS=$other, free memory $free GB after $MEM_WAIT_SEC s wait" }
        Start-Sleep -Seconds 15
    }
}

# Run one tgauss job under an outer timeout and a working-set cap.
function Invoke-GaussJob([string]$scriptRel, [string]$wd, [int]$timeoutSec) {
    $tmp = Join-Path $resultsDir ("_job_" + [guid]::NewGuid().ToString("N"))
    $outFile = "$tmp.out"; $errFile = "$tmp.err"
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = "cmd.exe"
    $psi.Arguments = "/c `"`"$GaussExe`" -b -x `"$scriptRel`" > `"$outFile`" 2> `"$errFile`"`""
    $psi.UseShellExecute = $false
    $psi.WorkingDirectory = $wd
    $psi.EnvironmentVariables["GAUSS26_CFG"] = $cfgDir
    $started = Get-Date
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $proc = [System.Diagnostics.Process]::Start($psi)
    $timedOut = $false; $memKill = $false; $peak = 0
    while (-not $proc.HasExited) {
        foreach ($x in @(Get-CimInstance Win32_Process -Filter "ParentProcessId=$($proc.Id)" -ErrorAction SilentlyContinue | Where-Object Name -eq 'tgauss.exe')) {
            $pp = Get-Process -Id $x.ProcessId -ErrorAction SilentlyContinue
            if ($pp -and $pp.PeakWorkingSet64 -gt $peak) { $peak = $pp.PeakWorkingSet64 }
        }
        if ($peak / 1MB -gt $PEAK_WS_CAP_MB) { $memKill = $true }
        if ($sw.Elapsed.TotalSeconds -gt $timeoutSec) { $timedOut = $true }
        if ($memKill -or $timedOut) { & taskkill.exe /PID $proc.Id /T /F 2>&1 | Out-Null; [void]$proc.WaitForExit(30000); break }
        Start-Sleep -Milliseconds 1000
    }
    $sw.Stop()
    $text = ""
    if (Test-Path $outFile) { $text += (Get-Content -Raw $outFile) }
    if (Test-Path $errFile) { $e = Get-Content -Raw $errFile; if ($e) { $text += "`n--- stderr ---`n" + $e } }
    Remove-Item -ErrorAction SilentlyContinue $outFile, $errFile
    return [pscustomobject]@{ text = $text; timedOut = $timedOut; memKill = $memKill; peakMB = [math]::Round($peak / 1MB, 0); wall = [math]::Round($sw.Elapsed.TotalSeconds, 2); started = $started.ToString('s') }
}
function Get-Tag([string]$text, [string]$tag) {
    $m = [regex]::Match($text, "(?m)^$tag[ \t]+(\S+)")
    if ($m.Success) { return $m.Groups[1].Value } else { return "" }
}
function Test-GaussError([string]$text) { return [regex]::IsMatch($text, '(?i)\bG\d{4}\b|error G|Undefined symbol') }

# Append every "TAG v1 v2 ..." row of $text to $csv as "endpoint,v1,v2,...".
function Export-Rows([string]$text, [string]$tag, [string]$csv, [string]$header, [string]$ep) {
    $p = Join-Path $resultsDir $csv
    if (-not (Test-Path $p)) { $header | Out-File -Encoding ascii $p }
    $rows = foreach ($m in [regex]::Matches($text, "(?m)^$tag[ \t]+(.+?)\s*$")) { "$ep," + (($m.Groups[1].Value.Trim() -split '\s+') -join ",") }
    if ($rows) { $rows | Out-File -Encoding ascii -Append $p }
}

$jobsHeader = "job_id,stage,kind,endpoint,started,wall_s,outer_timeout_s,timed_out,mem_killed,peak_ws_mb,free_gb_at_start,gauss_error,done,log"
function Invoke-Job([string]$id, [string]$kind, [string]$ep, [string]$stg, [string]$src, [int]$timeout, [string]$doneTag, $st) {
    if ((Get-ElapsedMin $st) + $timeout / 60.0 -gt $BUDGET_MIN) { Stop-Study $st "budget: job $id outer timeout $timeout s would pass $BUDGET_MIN min (elapsed $([math]::Round((Get-ElapsedMin $st),1)))" }
    $free = Wait-JobGate $st $id
    $src | Out-File -Encoding ascii -NoNewline (Join-Path $genDir "$id.e")
    $r = Invoke-GaussJob (Join-Path $genDir "$id.e") $here $timeout
    $pin = (Get-Content -Raw (Join-Path $repoRoot 'sslib.pin.json') | ConvertFrom-Json).pinnedCommit
    $logRel = "logs/$id.log"
    ("# job=$id stage=$stg head=$($st.freeze_head) outer_timeout_s=$timeout timed_out=$($r.timedOut) mem_killed=$($r.memKill) peak_ws_mb=$($r.peakMB) wall_s=$($r.wall) sslib_pin=$pin started=$($r.started)`n" + $r.text) |
        Out-File -Encoding utf8 (Join-Path $resultsDir $logRel)
    $done = $r.text -match "(?m)^$doneTag"
    $gerr = (Test-GaussError $r.text) -or (-not $done -and -not $r.timedOut -and -not $r.memKill)
    $p = Join-Path $resultsDir "jobs.csv"
    if (-not (Test-Path $p)) { $jobsHeader | Out-File -Encoding ascii $p }
    (@($id, $stg, $kind, $ep, $r.started, $r.wall.ToString($inv), $timeout, $r.timedOut, $r.memKill, $r.peakMB, $free.ToString($inv), $gerr, $done, $logRel) -join ",") | Out-File -Encoding ascii -Append $p
    Write-RunLog "[$stg] $id wall=$($r.wall)s peak=$($r.peakMB)MB free@start=${free}GB done=$done gerr=$gerr timeout=$($r.timedOut) memkill=$($r.memKill)"
    if ($r.timedOut) { Stop-Study $st "outer timeout in $id" }
    if ($r.memKill) { Stop-Study $st "memory cap ($PEAK_WS_CAP_MB MB) exceeded in $id" }
    if ($gerr) { Stop-Study $st "GAUSS error in $id" }
    return $r.text
}

function Get-EpId($sp, [int]$seed) {
    if ($sp.mode -eq 2) { return "ANCHOR-n1$($sp.n1)-t$($sp.tobs)-s$seed" }
    $h = if ($sp.harm -eq 1) { "Horacle" } else { "Hstatic" }
    return "n1$($sp.n1)-t$($sp.tobs)-s$seed-$h"
}
function Expand-Data([string]$tmpl, $sp, [int]$seed) {
    return $tmpl -replace '__MODE__', $sp.mode -replace '__N1__', $sp.n1 -replace '__TOBS__', $sp.tobs -replace '__SEED__', $seed `
        -replace '__HARM__', $sp.harm -replace '__QSCALE__', $QSCALE.ToString($inv)
}
function Get-WorkPath([string]$name) { return ((Join-Path $workDir $name) -replace '\\', '/') }

# ---------------------------------------------------------------- Phase C reference
$phaseC = @()
if (Test-Path $PHASEC_SUMMARY) { $phaseC = @(Import-Csv $PHASEC_SUMMARY | Where-Object { $_.replicate -eq "r1" }) }
function Get-PhaseCRow($sp, [int]$seed, [double]$qm) {
    return $phaseC | Where-Object { [int]$_.mode -eq $sp.mode -and [int]$_.n1 -eq $sp.n1 -and [int]$_.tobs -eq $sp.tobs -and [int]$_.seed_used -eq $seed -and [int]$_.harm -eq $sp.harm -and [double]::Parse($_.qmult, $inv) -eq $qm } | Select-Object -First 1
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
    $pcManifest = Get-Content (Join-Path $PHASEC_DIR "FREEZE.sha256")
    foreach ($f in $PHASEC_REUSED) {
        $line = $pcManifest | Where-Object { $_ -match "\s$([regex]::Escape($f))$" }
        if (-not $line -or ($line -split '\s+')[0] -ne (Get-NormHash (Join-Path $PHASEC_DIR $f))) { Stop-Study $st "[$stg] preflight: Phase C file $f differs from Phase C's frozen hash" }
    }
    Push-Location $repoRoot
    try {
        $head = (& git rev-parse HEAD).Trim()
        $rel = if ($resultsDir.StartsWith($repoRoot)) { ($resultsDir.Substring($repoRoot.Length + 1) -replace '\\', '/') + '/' } else { "::outside::" }
        $dirty = @(& git status --porcelain --untracked-files=all | Where-Object { -not $_.Substring(3).StartsWith($rel) })
    } finally { Pop-Location }
    if (-not $st.freeze_head) { $st.freeze_head = $head; Save-State $st }
    if ($head -ne $st.freeze_head) { Stop-Study $st "[$stg] preflight: HEAD $head differs from recorded $($st.freeze_head)" }
    if (-not $Smoke -and $dirty.Count -gt 0) { Stop-Study $st "[$stg] preflight: working tree not clean outside results/: $($dirty -join '; ')" }
    [void](Wait-JobGate $st "preflight-$stg")
    & powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repoRoot "scripts\verify_sslib_pin.ps1") -Strict *> (Join-Path $logDir "preflight_$($stg)_sslib_pin.log")
    if ($LASTEXITCODE -ne 0) { Stop-Study $st "[$stg] preflight: verify_sslib_pin.ps1 -Strict failed" }
    $tests = Join-Path $repoRoot "tests"
    foreach ($t in @(@("quaidstvpfit_test.e", 30), @("quaidstvp_loglik_xcheck_test.e", 12))) {
        $r = Invoke-GaussJob $t[0] $tests 900
        $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_$($t[0]).log")
        if ($r.timedOut -or $r.text -notmatch "ALL $($t[1]) CHECKS PASSED") { Stop-Study $st "[$stg] preflight: $($t[0]) did not report ALL $($t[1]) CHECKS PASSED" }
    }
    $r = Invoke-GaussJob "q0_equiv.e" $here 900
    $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_q0_equiv.log")
    $ok = (Get-Tag $r.text "Q0_EQ_PASS") -ne "" -and [double]::Parse((Get-Tag $r.text "Q0_EQ_PASS"), $inv) -eq 1
    if ($r.timedOut -or -not $ok -or (Test-GaussError $r.text)) { Stop-Study $st "[$stg] preflight: q0_equiv.e gate failed" }
    Write-RunLog "[$stg] preflight PASS (HEAD $head, frozen hashes, Phase C reused files unchanged, no other GAUSS, pin strict, Stage6 30/30, P5 12/12, reduced fit == quaidsTVPFit exactly, fixed element exact, Fct = LL)"
}

# ---------------------------------------------------------------- one endpoint
function Invoke-Endpoint([string]$stg, $sp, [int]$seed, $st) {
    $ep = Get-EpId $sp $seed
    $qfiles = @()
    $fitT = Get-Content -Raw (Join-Path $here "q0_fit.template.e")
    foreach ($qm in $STARTS) {
        $qs = $qm.ToString($inv)
        $qf = Get-WorkPath "$ep-q$($qs -replace '\.', 'p').fmt"   # GAUSS save truncates a name at its first dot
        $qfiles += $qf
        $src = (Expand-Data $fitT $sp $seed) -replace '__QMULT__', $qs -replace '__MAXITERS__', $sp.fitIters -replace '__MAXTIME__', 600 -replace '__QFILE__', $qf
        $t = Invoke-Job "$stg-$ep-fit-q$qs" "fit" $ep $stg $src $sp.toFit "Q0F_DONE" $st
        $vals = @($ep, $qs, (Get-Tag $t "Q0F_RC"), (Get-Tag $t "Q0F_ITERS"), (Get-Tag $t "Q0F_CMLMT_MIN"), (Get-Tag $t "Q0F_FCT"), (Get-Tag $t "Q0F_LLFIT"),
            (Get-Tag $t "Q0F_MAXGRAD"), (Get-Tag $t "Q0F_NBOUNDARY"), (Get-Tag $t "Q0F_QFIT_SUM"), (Get-Tag $t "Q0F_SAVE_DIFF"))
        $match = ""
        if (-not $Smoke) {
            $pc = Get-PhaseCRow $sp $seed $qm
            if (-not $pc) { Stop-Study $st "no Phase C row for $ep q$qs" }
            $match = ($pc.rc -ceq $vals[2]) -and ($pc.iters -ceq $vals[3]) -and ($pc.llfit -ceq $vals[6]) -and ($pc.qfit_sum -ceq $vals[9])
        }
        $p = Join-Path $resultsDir "fits.csv"
        if (-not (Test-Path $p)) { "endpoint,qmult,rc,iters,cmlmt_min,fct,llfit,maxgrad,nboundary,qfit_sum,save_diff,phasec_exact" | Out-File -Encoding ascii $p }
        ((@($vals) + $match) -join ",") | Out-File -Encoding ascii -Append $p
        if ($vals[10] -ne "" -and [double]::Parse($vals[10], $inv) -ne 0) { Stop-Study $st "Qfit save/load round trip not exact in $ep q$qs" }
        if (-not $Smoke -and -not $match) { Stop-Study $st "Phase C reproduction mismatch in $ep q$qs" }
    }
    $qhatFile = Get-WorkPath "$ep-qhat.fmt"
    $gT = Get-Content -Raw (Join-Path $here "q0_geom.template.e")
    $src = (Expand-Data $gT $sp $seed) -replace '__QFILE1__', $qfiles[0] -replace '__QFILE2__', $qfiles[1] -replace '__QFILE3__', $qfiles[2] `
        -replace '__QHATFILE__', $qhatFile -replace '__PHI__', $PHI -replace '__HSTEP1__', $HSTEP1.ToString($inv) -replace '__HSTEP2__', $HSTEP2.ToString($inv)
    $t = Invoke-Job "$stg-$ep-geom" "geom" $ep $stg $src $sp.toGeom "Q0G_DONE" $st
    Export-Rows $t "Q0G_BEST" "best.csv" "endpoint,best_start,llhat" $ep
    Export-Rows $t "Q0G_START" "starts.csv" "endpoint,start,qmult,ll,ll_minus_best,nboundary,max_dlog10_nonbnd" $ep
    Export-Rows $t "Q0G_ELEM" "elements.csv" "endpoint,j,group,qhat,qref,qtrue,qhat_over_qref,boundary,below_floor" $ep
    Export-Rows $t "Q0G_POINT" "points.csv" "endpoint,code,ll,loss" $ep
    Export-Rows $t "Q0G_SLICE" "slices.csv" "endpoint,kind,idx,log10c,ll,loss" $ep
    Export-Rows $t "Q0G_HEIG" "heig.csv" "endpoint,step,i,eig" $ep
    Export-Rows $t "Q0G_HVEC" "hvec.csv" "endpoint,step,i,v..." $ep
    Export-Rows $t "Q0G_HDIAG" "hdiag.csv" "endpoint,step,j,neg_hjj,grad_logq" $ep
    Export-Rows $t "Q0G_HINT" "hint.csv" "endpoint,step,nint,mineig,maxeig,cond" $ep
    Export-Rows $t "Q0G_CORR" "corr.csv" "endpoint,step,i,j,c" $ep
    $pT = Get-Content -Raw (Join-Path $here "q0_prof.template.e")
    foreach ($g in 1..3) {
        $src = (Expand-Data $pT $sp $seed) -replace '__QHATFILE__', $qhatFile -replace '__PHI__', $PHI -replace '__GROUP__', $g `
            -replace '__LEVELS__', $sp.levels -replace '__LEVELSB__', $sp.levelsB -replace '__MAXITERS__', $sp.profIters -replace '__MAXTIME__', 600
        $t = Invoke-Job "$stg-$ep-prof-g$g" "prof" $ep $stg $src $sp.toProf "Q0P_DONE" $st
        Export-Rows $t "Q0P_ROW" "profiles.csv" "endpoint,group,log10c,start,rc,iters,llslice,llprof,prof_minus_slice,fixdiff,nfloored,sec" $ep
        Export-Rows $t "Q0P_Q" "profile_q.csv" "endpoint,group,log10c,start,q..." $ep
    }
}

function Test-IncludeC($st) {
    $free = Get-FreeGB
    $el = Get-ElapsedMin $st
    $sp = $SPEC.C
    $worst = $C_PROJ_MIN
    Write-RunLog "[C inclusion] elapsed=$([math]::Round($el,1)) min, worst-case C=$([math]::Round($worst,1)) min, free memory=$free GB"
    if ($free -lt $C_MIN_FREE_GB) { Write-RunLog "[C inclusion] EXCLUDED: free memory $free GB < $C_MIN_FREE_GB GB"; return $false }
    if ($el + $worst -gt $BUDGET_MIN) { Write-RunLog "[C inclusion] EXCLUDED: worst case would pass $BUDGET_MIN min"; return $false }
    return $true
}

function Invoke-Stage([string]$stg, $st) {
    if ($st.stopped) { Write-RunLog "study already halted: $($st.stop_reason)"; exit 2 }
    if (@($st.stages_done) -contains $stg) { Write-RunLog "[$stg] already done; skipping"; return }
    if (-not $st.budget_start) { $st.budget_start = (Get-Date).ToString('o'); Save-State $st }
    if ($stg -eq "C" -and -not (Test-IncludeC $st)) {
        $st.stages_done = @(@($st.stages_done) + $stg); $st.c_excluded = $true; Save-State $st; return
    }
    Invoke-Preflight $stg $st
    $sp = $SPEC[$stg]
    foreach ($sd in $sp.seeds) { Invoke-Endpoint $stg $sp $sd $st }
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
if ($Smoke -and (-not $ResultsDir -or $Stage -ne "A")) { throw "-Smoke requires -ResultsDir (outside this directory) and -Stage A" }
if ($DryRun) {
    foreach ($s in $STAGE_ORDER) {
        $sp = $SPEC[$s]
        $n = $sp.seeds.Count
        Write-Host ("Stage {0}: {1} endpoint(s) {2}; per endpoint 3 fits (cap {3}) + 1 geometry + 3 profile jobs (levels {4}; start b at {5}; cap {6}); worst-case {7:N0} min" -f `
            $s, $n, (($sp.seeds | ForEach-Object { Get-EpId $sp $_ }) -join ' '), $sp.fitIters, $sp.levels, $sp.levelsB, $sp.profIters, ($n * (3 * $sp.toFit + $sp.toGeom + 3 * $sp.toProf) / 60.0))
    }
    Write-Host "free memory now: $(Get-FreeGB) GB; budget $BUDGET_MIN min"
    exit 0
}
New-Item -ItemType Directory -Force $resultsDir, $logDir, $genDir, $workDir | Out-Null
$st = Get-State
$todo = if ($Stage -eq "ALL") { $STAGE_ORDER } else { @($Stage) }
foreach ($s in $todo) { Invoke-Stage $s $st }
Write-RunLog "run finished: stages done = $(@($st.stages_done) -join ',') elapsed=$([math]::Round((Get-ElapsedMin $st),1)) min"
