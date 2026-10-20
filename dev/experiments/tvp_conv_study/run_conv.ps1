<#
.SYNOPSIS
  TVP convergence-correctness study (Phase C) driver. See PROTOCOL.md.
.DESCRIPTION
  Runs the frozen study stages strictly sequentially: A (n1=2), B (n1=3),
  C (n1=4), ANCHOR (P6 example-data anchor). Before every stage: frozen
  file hashes, clean recorded HEAD (only results\ may differ), strict
  sslib pin, TVP Stage 6 30/30, P5 cross-check 12/12, and the
  harness-equivalence/Fct gate (conv_equiv.e). DGP stages then screen
  share admissibility and map any inadmissible study seed to the next
  admissible reserve seed. Each fit is one tgauss job under an outer
  timeout; one CSV row per attempted job (results\summary.csv), full log
  per job (results\logs). Any GAUSS error, outer timeout or CONV-FAILC3
  halts the whole study; so does a job whose worst-case time would
  exceed the six-hour budget. Seed reduction is applied only by the
  projection rule. -DryRun prints the plan and projection only.
  -Freeze writes FREEZE.sha256 (SHA-256 of each frozen file with CR
  removed, so line-ending conversion cannot change it) and prints the
  freeze hash (SHA-256 of FREEZE.sha256 itself). -SmokeCells replaces
  Stage A's grid with two Stage 6 fixture-sized harness smoke cells
  (n1=2, tobs=120, seed 7) and must be combined with -ResultsDir outside
  this directory; smoke output is never study evidence.
#>
param(
    [ValidateSet("ALL", "A", "B", "C", "ANCHOR")]
    [string]$Stage = "ALL",
    [switch]$DryRun,
    [switch]$Freeze,
    [switch]$SmokeCells,
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
$summary = Join-Path $resultsDir "summary.csv"
$stateFile = Join-Path $resultsDir "run_state.json"
$runLog = Join-Path $resultsDir "run_log.txt"

# ---------------------------------------------------------------- frozen design
$QSCALE = 0.01
$BUDGET_MIN = 360.0
$PREFLIGHT_SEC = 120.0          # projected cost of one stage preflight
$OVERHEAD_PRIOR_SEC = 30.0      # projected non-CMLMT cost of one fit job
# Prior CMLMT seconds per iteration at tobs = 500, by n1 (P6 maxima,
# rounded up; n1 = 2 deliberately conservative). Scaled by tobs/500.
$RATE_PRIOR = @{ 2 = 0.30; 3 = 0.75; 4 = 2.60 }
$STARTS = @(1, 10, 0.1)         # multipliers on the shipped q0 rule
$RESERVES = 1101..1110
$TIMEOUT_SEC = @{ 2 = 900; 3 = 1800; 4 = 3600 }
$ANCHOR_TIMEOUT_SEC = 3600
$STAGE_ORDER = @("A", "B", "C", "ANCHOR")
$FROZEN_FILES = @("PROTOCOL.md", "conv_lib.src", "conv_chain.inc", "conv_includes.inc",
    "conv_fit.template.e", "conv_screen.template.e", "conv_equiv.e", "run_conv.ps1", "analyze_conv.ps1")

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
    return @{ n1_3_seeds = 6; n1_4_seeds = 4; stages_done = @(); stopped = $false }
}
function Save-State($st) { if (-not $DryRun) { ($st | ConvertTo-Json -Depth 5) | Out-File -Encoding utf8 $stateFile } }

# Grid cells for a stage (study seeds before reserve mapping).
function Get-Cells([string]$stg, $st) {
    $cells = @()
    switch ($stg) {
        "A" {
            if ($SmokeCells) {
                foreach ($q in @(1, 10)) { $cells += [pscustomobject]@{ mode = 1; n1 = 2; tobs = 120; seed = 7; harm = 1; qmult = $q; maxIters = 500; maxTime = 600 } }
                return $cells
            }
            foreach ($t in @(200, 500)) {
                foreach ($sd in 1001..1008) {
                    $harms = if ($t -eq 200) { @(1, 2) } else { @(1) }
                    foreach ($hm in $harms) { foreach ($q in $STARTS) {
                        $cells += [pscustomobject]@{ mode = 1; n1 = 2; tobs = $t; seed = $sd; harm = $hm; qmult = $q; maxIters = 500; maxTime = 600 } } }
                }
            }
        }
        "B" {
            foreach ($sd in (1001..(1000 + [int]$st.n1_3_seeds))) { foreach ($hm in @(1, 2)) { foreach ($q in $STARTS) {
                $cells += [pscustomobject]@{ mode = 1; n1 = 3; tobs = 200; seed = $sd; harm = $hm; qmult = $q; maxIters = 500; maxTime = 600 } } } }
        }
        "C" {
            foreach ($sd in (1001..(1000 + [int]$st.n1_4_seeds))) { foreach ($q in $STARTS) {
                $cells += [pscustomobject]@{ mode = 1; n1 = 4; tobs = 200; seed = $sd; harm = 1; qmult = $q; maxIters = 500; maxTime = 600 } } }
        }
        "ANCHOR" {
            foreach ($q in $STARTS) { $cells += [pscustomobject]@{ mode = 2; n1 = 3; tobs = 500; seed = 204; harm = 2; qmult = $q; maxIters = 2000; maxTime = 600 } }
            foreach ($q in $STARTS) { $cells += [pscustomobject]@{ mode = 2; n1 = 4; tobs = 500; seed = 204; harm = 2; qmult = $q; maxIters = 500; maxTime = 600 } }
        }
    }
    return $cells
}

function Get-Rate([int]$n1) {
    # max(prior, largest measured CMLMT s/iter in this study normalised to tobs 500)
    $r = $RATE_PRIOR[$n1]
    if (Test-Path $summary) {
        foreach ($row in (Import-Csv $summary)) {
            if ([int]$row.n1 -eq $n1 -and $row.s_per_iter -ne "") {
                $m = [double]::Parse($row.s_per_iter, $inv) * 500.0 / [double]$row.tobs
                if ($m -gt $r) { $r = $m }
            }
        }
    }
    return $r
}
function Get-Overhead {
    $o = $OVERHEAD_PRIOR_SEC
    if (Test-Path $summary) {
        foreach ($row in (Import-Csv $summary)) {
            if ($row.wall_s -ne "" -and $row.cmlmt_min -ne "") {
                $v = [double]::Parse($row.wall_s, $inv) - 60.0 * [double]::Parse($row.cmlmt_min, $inv)
                if ($v -gt $o) { $o = $v }
            }
        }
    }
    return $o
}
function Get-WorstSec($c) { return (Get-Overhead) + ($c.maxIters + 1) * (Get-Rate $c.n1) * $c.tobs / 500.0 }

# Projected worst-case minutes for stages: one preflight per stage plus
# every cell at its iteration cap. Determinism replicates are excluded;
# they run only when the per-job budget guard allows.
function Get-ProjectionMin([string[]]$stages, $st) {
    $sec = 0.0
    foreach ($s in $stages) {
        foreach ($c in (Get-Cells $s $st)) { $sec += Get-WorstSec $c }
        $sec += $PREFLIGHT_SEC
    }
    return $sec / 60.0
}

function Get-ElapsedMin($st) {
    if (-not $st.budget_start) { return 0.0 }
    return ((Get-Date) - [datetime]::Parse($st.budget_start, $inv)).TotalMinutes
}

# Run one tgauss job from this directory under an outer timeout.
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

function Stop-Study($st, [string]$reason) {
    $st.stopped = $true
    $st.stop_reason = $reason
    Save-State $st
    Write-RunLog "STUDY HALTED: $reason"
    if (-not $DryRun) { $reason | Out-File -Encoding utf8 (Join-Path $resultsDir "STOP.txt") }
    exit 2
}

# ---------------------------------------------------------------- preflight
function Invoke-Preflight([string]$stg, $st) {
    Write-RunLog "[$stg] preflight start"
    $manifest = Get-Content (Join-Path $here "FREEZE.sha256")
    foreach ($line in $manifest) {
        if ($line -match '^([0-9a-f]{64})\s+(\S+)$') {
            $h = Get-NormHash (Join-Path $here $Matches[2])
            if ($h -ne $Matches[1]) { Stop-Study $st "[$stg] preflight: frozen file $($Matches[2]) changed (hash mismatch)" }
        }
    }
    foreach ($f in $FROZEN_FILES) { if (-not ($manifest -match [regex]::Escape($f))) { Stop-Study $st "[$stg] preflight: $f missing from FREEZE.sha256" } }
    Push-Location $repoRoot
    try {
        $head = (& git rev-parse HEAD).Trim()
        $dirty = @(& git status --porcelain --untracked-files=all | Where-Object { $_ -notmatch 'dev/experiments/tvp_conv_study/results/' })
    } finally { Pop-Location }
    if (-not $st.freeze_head) { $st.freeze_head = $head; Save-State $st }
    if ($head -ne $st.freeze_head) { Stop-Study $st "[$stg] preflight: HEAD $head differs from recorded $($st.freeze_head)" }
    if ($dirty.Count -gt 0) { Stop-Study $st "[$stg] preflight: working tree not clean outside results/: $($dirty -join '; ')" }
    & powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $repoRoot "scripts\verify_sslib_pin.ps1") -Strict *> (Join-Path $logDir "preflight_$($stg)_sslib_pin.log")
    if ($LASTEXITCODE -ne 0) { Stop-Study $st "[$stg] preflight: verify_sslib_pin.ps1 -Strict failed" }
    $tests = Join-Path $repoRoot "tests"
    foreach ($t in @(@("quaidstvpfit_test.e", 30), @("quaidstvp_loglik_xcheck_test.e", 12))) {
        $r = Invoke-GaussJob $t[0] $tests 900
        $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_$($t[0]).log")
        if ($r.timedOut -or $r.text -notmatch "ALL $($t[1]) CHECKS PASSED") { Stop-Study $st "[$stg] preflight: $($t[0]) did not report ALL $($t[1]) CHECKS PASSED" }
    }
    $r = Invoke-GaussJob "conv_equiv.e" $here 900
    $r.text | Out-File -Encoding utf8 (Join-Path $logDir "preflight_$($stg)_conv_equiv.log")
    $okEq = (Get-Tag $r.text "CONV_EQ_PASS") -ne "" -and [double]::Parse((Get-Tag $r.text "CONV_EQ_PASS"), $inv) -eq 1
    $fd = Get-Tag $r.text "CONV_EQ_FCT_MINUS_SUM"
    $okFct = $fd -ne "" -and [double]::Parse($fd, $inv) -eq 0
    if ($r.timedOut -or -not $okEq -or -not $okFct -or (Test-GaussError $r.text)) { Stop-Study $st "[$stg] preflight: harness equivalence / Fct=sum gate failed" }
    Write-RunLog "[$stg] preflight PASS (HEAD $head, pin strict, Stage6 30/30, P5 12/12, equivalence exact, Fct = total LL)"
}

# Admissibility screen -> study seed to used seed map for one (n1, tobs).
function Get-SeedMap([int]$n1, [int]$tobs, [int[]]$seeds, [string]$stg, $st) {
    $all = @($seeds) + @($RESERVES)
    $name = "screen_$($stg)_n1$($n1)_t$tobs"
    $tmpl = Get-Content -Raw (Join-Path $here "conv_screen.template.e")
    $src = $tmpl -replace '__N1__', $n1 -replace '__TOBS__', $tobs -replace '__SEEDS__', ("{" + ($all -join ",") + "}") -replace '__QSCALE__', $QSCALE.ToString($inv)
    $src | Out-File -Encoding ascii (Join-Path $genDir "$name.e")
    $r = Invoke-GaussJob (Join-Path $genDir "$name.e") $here 600
    $r.text | Out-File -Encoding utf8 (Join-Path $logDir "$name.log")
    if ($r.timedOut -or (Test-GaussError $r.text) -or $r.text -notmatch 'CONV_SCREEN_DONE') { Stop-Study $st "[$stg] admissibility screen failed for n1=$n1 tobs=$tobs" }
    $ok = @{}
    foreach ($m in [regex]::Matches($r.text, 'seed=\s*(\S+)\s+minw=\s*(\S+)\s+maxw=\s*(\S+)\s+admissible=\s*(\S+)')) {
        $ok[[int][double]::Parse($m.Groups[1].Value, $inv)] = ([double]::Parse($m.Groups[4].Value, $inv) -eq 1)
    }
    $map = @{}
    $nextRes = 0
    foreach ($sd in $seeds) {
        if ($ok[$sd]) { $map[$sd] = $sd; continue }
        while ($nextRes -lt $RESERVES.Count -and -not $ok[$RESERVES[$nextRes]]) { $nextRes++ }
        if ($nextRes -ge $RESERVES.Count) { Stop-Study $st "[$stg] reserve seeds exhausted for n1=$n1 tobs=$tobs" }
        $map[$sd] = $RESERVES[$nextRes]; $nextRes++
        Write-RunLog "[$stg] seed $sd inadmissible for n1=$n1 tobs=$tobs -> reserve $($map[$sd])"
    }
    return $map
}

function Get-CellId($c, [string]$stg, [string]$rep) {
    $h = if ($c.harm -eq 1) { "Horacle" } else { "Hstatic" }
    return "$stg-n1$($c.n1)-t$($c.tobs)-s$($c.seed)-$h-q$($c.qmult.ToString($inv))-$rep"
}

$csvHeader = "cell_id,stage,replicate,mode,n1,tobs,seed_study,seed_used,harm,qmult,max_iters,max_time_min,qscale,outer_timeout_s,started,wall_s,timed_out,gauss_error,done,rc,iters,cmlmt_min,s_per_iter,fit_elapsed_s,fct,llkf,llfit,lltrue,llstart,tol,maxgrad,nboundary,qfit_sum,minw,maxw,failc3,warnings,log"

function Invoke-Cell($c, [string]$stg, [int]$seedUsed, [string]$rep, $st) {
    $id = Get-CellId $c $stg $rep
    $timeout = if ($c.mode -eq 2) { $ANCHOR_TIMEOUT_SEC } else { $TIMEOUT_SEC[[int]$c.n1] }
    $tmpl = Get-Content -Raw (Join-Path $here "conv_fit.template.e")
    $src = $tmpl -replace '__MODE__', $c.mode -replace '__N1__', $c.n1 -replace '__TOBS__', $c.tobs -replace '__SEED__', $seedUsed `
        -replace '__HARM__', $c.harm -replace '__QMULT__', $c.qmult.ToString($inv) -replace '__MAXITERS__', $c.maxIters `
        -replace '__MAXTIME__', $c.maxTime -replace '__QSCALE__', $QSCALE.ToString($inv)
    $src | Out-File -Encoding ascii (Join-Path $genDir "$id.e")
    $r = Invoke-GaussJob (Join-Path $genDir "$id.e") $here $timeout
    $logRel = "logs/$id.log"
    $pin = (Get-Content -Raw (Join-Path $repoRoot 'sslib.pin.json') | ConvertFrom-Json).pinnedCommit
    ("# cell=$id stage=$stg head=$($st.freeze_head) outer_timeout_s=$timeout timed_out=$($r.timedOut) wall_s=$($r.wall) sslib_pin=$pin started=$($r.started)`n" + $r.text) |
        Out-File -Encoding utf8 (Join-Path $resultsDir $logRel)
    $t = $r.text
    $done = $t -match '(?m)^CONV_DONE'
    $gerr = (Test-GaussError $t) -or (-not $done -and -not $r.timedOut)
    $iters = Get-Tag $t "CONV_ITERS"; $cm = Get-Tag $t "CONV_CMLMT_MIN"
    $spi = ""
    if ($iters -ne "" -and $cm -ne "" -and [double]::Parse($iters, $inv) -gt 0) {
        $spi = (60.0 * [double]::Parse($cm, $inv) / [double]::Parse($iters, $inv)).ToString("R", $inv)
    }
    $fc3 = if ($t -match '(?m)^CONV_FAILC3 1') { 1 } elseif ($t -match '(?m)^CONV_FAILC3 0') { 0 } else { "" }
    $warn = ([regex]::Matches($t, '(?im)^.*(complex|failed to invert|warning|not computed).*$') | ForEach-Object { $_.Value.Trim() } | Select-Object -Unique) -join " | "
    $warn = $warn -replace ',', ';' -replace '"', "'"
    $vals = @($id, $stg, $rep, $c.mode, $c.n1, $c.tobs, $c.seed, $seedUsed, $c.harm, $c.qmult.ToString($inv), $c.maxIters, $c.maxTime, $QSCALE.ToString($inv),
        $timeout, $r.started, $r.wall.ToString($inv), $r.timedOut, $gerr, $done, (Get-Tag $t "CONV_RC"), $iters, $cm, $spi, (Get-Tag $t "CONV_ELAPSED_SEC"),
        (Get-Tag $t "CONV_FCT"), (Get-Tag $t "CONV_LLKF"), (Get-Tag $t "CONV_LLFIT"), (Get-Tag $t "CONV_LLTRUE"), (Get-Tag $t "CONV_LLSTART"),
        (Get-Tag $t "CONV_TOL"), (Get-Tag $t "CONV_MAXGRAD"), (Get-Tag $t "CONV_NBOUNDARY"), (Get-Tag $t "CONV_QFIT_SUM"),
        (Get-Tag $t "CONV_MINW"), (Get-Tag $t "CONV_MAXW"), $fc3, "`"$warn`"", $logRel)
    ($vals -join ",") | Out-File -Encoding ascii -Append $summary
    Write-RunLog "[$stg] $id wall=$($r.wall)s rc=$(Get-Tag $t 'CONV_RC') iters=$iters llfit=$(Get-Tag $t 'CONV_LLFIT') lltrue=$(Get-Tag $t 'CONV_LLTRUE') failc3=$fc3 timeout=$($r.timedOut) gerr=$gerr"
    if ($r.timedOut) { Stop-Study $st "[$stg] outer timeout in $id" }
    if ($gerr) { Stop-Study $st "[$stg] GAUSS error in $id" }
    if ($fc3 -eq 1) { Stop-Study $st "[$stg] CONV-FAILC3 in $id" }
    return [pscustomobject]@{ id = $id; rc = (Get-Tag $t "CONV_RC"); text = $t }
}

# Endpoint text compared by the determinism replicate: retcode, iterations
# and everything from CONV_FCT to CONV_DONE (timing lines excluded).
function Get-Endpoint([string]$t) {
    return "rc=$(Get-Tag $t 'CONV_RC') iters=$(Get-Tag $t 'CONV_ITERS')`n" + [regex]::Match($t, '(?s)CONV_FCT.*CONV_DONE').Value
}

# ---------------------------------------------------------------- reduction rule
function Invoke-Projection([string]$after, $st) {
    $idx = [array]::IndexOf($STAGE_ORDER, $after)
    if ($idx -ge $STAGE_ORDER.Count - 1) { return }
    $remaining = @($STAGE_ORDER[($idx + 1)..($STAGE_ORDER.Count - 1)])
    $el = Get-ElapsedMin $st
    $p = Get-ProjectionMin $remaining $st
    Write-RunLog "[projection after $after] elapsed=$([math]::Round($el,1)) min, remaining worst-case=$([math]::Round($p,1)) min ($($remaining -join ',')), rates s/iter@500: n1=2 $(Get-Rate 2) n1=3 $(Get-Rate 3) n1=4 $(Get-Rate 4), overhead $(Get-Overhead) s"
    if ($el + $p -gt $BUDGET_MIN -and $remaining -contains "C" -and [int]$st.n1_4_seeds -gt 2) {
        $st.n1_4_seeds = 2; Save-State $st
        $p = Get-ProjectionMin $remaining $st
        Write-RunLog "[projection after $after] REDUCTION: n1=4 -> 2 seeds (1001-1002); remaining worst-case now $([math]::Round($p,1)) min"
    }
    if ($el + $p -gt $BUDGET_MIN -and $remaining -contains "B" -and [int]$st.n1_3_seeds -gt 4) {
        $st.n1_3_seeds = 4; Save-State $st
        $p = Get-ProjectionMin $remaining $st
        Write-RunLog "[projection after $after] REDUCTION: n1=3 -> 4 seeds (1001-1004); remaining worst-case now $([math]::Round($p,1)) min"
    }
    if ($el + $p -gt $BUDGET_MIN) { Write-RunLog "[projection after $after] still above budget after permitted reductions; per-job budget guard applies" }
}

# ---------------------------------------------------------------- stage runner
function Invoke-Stage([string]$stg, $st) {
    if ($st.stopped) { Write-RunLog "study already halted: $($st.stop_reason)"; exit 2 }
    if (@($st.stages_done) -contains $stg) { Write-RunLog "[$stg] already done; skipping"; return }
    if (-not $st.budget_start) { $st.budget_start = (Get-Date).ToString('o'); Save-State $st }
    Invoke-Preflight $stg $st
    $cells = Get-Cells $stg $st
    $maps = @{}
    if ($stg -ne "ANCHOR") {
        foreach ($t in ($cells | ForEach-Object { $_.tobs } | Sort-Object -Unique)) {
            $sds = @($cells | Where-Object { $_.tobs -eq $t } | ForEach-Object { $_.seed } | Sort-Object -Unique)
            $maps[$t] = Get-SeedMap $cells[0].n1 $t $sds $stg $st
        }
    }
    $replicated = $false
    foreach ($c in $cells) {
        $worst = Get-WorstSec $c
        if ((Get-ElapsedMin $st) + $worst / 60.0 -gt $BUDGET_MIN) { Stop-Study $st "[$stg] budget: next cell $(Get-CellId $c $stg 'r1') worst-case $([math]::Round($worst/60,1)) min would exceed $BUDGET_MIN min" }
        $su = if ($stg -eq "ANCHOR") { $c.seed } else { $maps[$c.tobs][$c.seed] }
        $res = Invoke-Cell $c $stg $su "r1" $st
        if (-not $replicated -and $res.rc -ne "" -and [double]::Parse($res.rc, $inv) -eq 2) {
            $replicated = $true
            if ((Get-ElapsedMin $st) + $worst / 60.0 -gt $BUDGET_MIN) { Write-RunLog "[$stg] determinism replicate of $($res.id) SKIPPED (budget)"; continue }
            $rep = Invoke-Cell $c $stg $su "r2" $st
            $same = (Get-Endpoint $res.text) -ceq (Get-Endpoint $rep.text)
            Write-RunLog "[$stg] determinism replicate of $($res.id): endpoint byte-identical=$same"
            "$($res.id),$same" | Out-File -Encoding ascii -Append (Join-Path $resultsDir "determinism.csv")
        }
    }
    $st.stages_done = @(@($st.stages_done) + $stg)
    Save-State $st
    Write-RunLog "[$stg] stage complete ($($cells.Count) cells)"
    Invoke-Projection $stg $st
}

# ---------------------------------------------------------------- main
if ($Freeze) {
    $lines = foreach ($f in $FROZEN_FILES) { "$(Get-NormHash (Join-Path $here $f))  $f" }
    ($lines -join "`n") + "`n" | Out-File -Encoding ascii -NoNewline (Join-Path $here "FREEZE.sha256")
    Write-Host ($lines -join "`n")
    Write-Host "FREEZE HASH $(Get-NormHash (Join-Path $here 'FREEZE.sha256'))"
    exit 0
}
if ($SmokeCells -and (-not $ResultsDir -or $Stage -ne "A")) { throw "-SmokeCells requires -ResultsDir (outside this directory) and -Stage A" }
New-Item -ItemType Directory -Force $resultsDir, $logDir, $genDir | Out-Null
$st = Get-State
if ($DryRun) {
    foreach ($s in $STAGE_ORDER) {
        $cells = Get-Cells $s $st
        Write-Host "Stage $s : $($cells.Count) cells; worst-case $([math]::Round((Get-ProjectionMin @($s) $st),1)) min"
    }
    Write-Host "Total worst-case (all stages): $([math]::Round((Get-ProjectionMin $STAGE_ORDER $st),1)) min"
    exit 0
}
if (-not (Test-Path $summary)) { $csvHeader | Out-File -Encoding ascii $summary }
if (-not (Test-Path (Join-Path $resultsDir "determinism.csv"))) { "cell_id,endpoint_identical" | Out-File -Encoding ascii (Join-Path $resultsDir "determinism.csv") }
$todo = if ($Stage -eq "ALL") { $STAGE_ORDER } else { @($Stage) }
foreach ($s in $todo) { Invoke-Stage $s $st }
Write-RunLog "run finished: stages done = $(@($st.stages_done) -join ',') elapsed=$([math]::Round((Get-ElapsedMin $st),1)) min"
