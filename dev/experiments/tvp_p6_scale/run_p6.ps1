<#
.SYNOPSIS
  TVP-AIDS P6 bounded scale re-measurement driver (dev/experiments/tvp_p6_scale).
.DESCRIPTION
  For each requested n1 value, writes gen_n1<k>_<maxIters>x<maxTime>.e from
  p6_scale_fit.template.e (nSub = n1 + 1) and runs it -Repeats times with
  tgauss: one quaidsTVPFit() call per job, under the GAUSS26_CFG tsmt
  override and an outer hard wall-clock timeout (process tree killed via
  taskkill /T on expiry). Each job's full stdout+stderr goes to
  results\<label>_r<i>.log; one summary row per job is appended to
  results\summary.csv. Jobs run strictly sequentially (no timing overlap).
#>
param(
    [int[]]$N1 = @(3, 4),
    [int]$MaxIters = 500,
    [double]$MaxTime = 5,
    [int]$Repeats = 3,
    [int]$HardTimeoutSeconds = 900,
    [string]$GaussExe = "C:\gauss26\tgauss.exe"
)

$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $here "..\..\..")
$cfgDir = Join-Path $repoRoot "tests\gauss26_cfg_override"
$resultsDir = Join-Path $here "results"
New-Item -ItemType Directory -Force $resultsDir | Out-Null
$summary = Join-Path $resultsDir "summary.csv"
if (-not (Test-Path $summary)) {
    "label,n1,max_iters,max_time_min,repeat,hard_timeout_sec,wall_sec,timed_out,fit_elapsed_sec,mle_retcode,cmlmt_iterations,cmlmt_minutes,qfit_sum,warnings,done" | Out-File -Encoding ascii $summary
}
$pin = (Get-Content -Raw (Join-Path $repoRoot 'sslib.pin.json') | ConvertFrom-Json).pinnedCommit
$template = Get-Content -Raw (Join-Path $here "p6_scale_fit.template.e")

foreach ($k in $N1) {
    $label = "n1$($k)_$($MaxIters)x$($MaxTime)"
    $script = "gen_$label.e"
    ($template -replace '__NSUB__', ($k + 1) -replace '__MAXITERS__', $MaxIters -replace '__MAXTIME__', $MaxTime) |
        Out-File -Encoding ascii (Join-Path $here $script)
    for ($i = 1; $i -le $Repeats; $i++) {
        $log = Join-Path $resultsDir "$($label)_r$i.log"
        $outFile = "$log.out"; $errFile = "$log.err"
        $psi = [System.Diagnostics.ProcessStartInfo]::new()
        $psi.FileName = "cmd.exe"
        $psi.Arguments = "/c `"`"$GaussExe`" -b -x $script > `"$outFile`" 2> `"$errFile`"`""
        $psi.UseShellExecute = $false
        $psi.WorkingDirectory = $here
        $psi.EnvironmentVariables["GAUSS26_CFG"] = $cfgDir
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        $proc = [System.Diagnostics.Process]::Start($psi)
        $timedOut = $false
        if (-not $proc.WaitForExit($HardTimeoutSeconds * 1000)) {
            $timedOut = $true
            & taskkill.exe /PID $proc.Id /T /F 2>&1 | Out-Null
            [void]$proc.WaitForExit(30000)
        }
        $sw.Stop()
        $wall = [math]::Round($sw.Elapsed.TotalSeconds, 2)
        $text = ""
        if (Test-Path $outFile) { $text += (Get-Content -Raw $outFile) }
        if (Test-Path $errFile) { $err = Get-Content -Raw $errFile; if ($err) { $text += "`n--- stderr ---`n" + $err } }
        Remove-Item -ErrorAction SilentlyContinue $outFile, $errFile
        $text = "# label=$label script=$script repeat=$i hard_timeout_sec=$HardTimeoutSeconds timed_out=$timedOut wall_sec=$wall sslib_pin=$pin started=$((Get-Date).AddSeconds(-$wall).ToString('s'))`n" + $text
        $text | Out-File -Encoding utf8 $log
        $el = if ($text -match 'P6_ELAPSED_SEC\s+(\S+)') { $Matches[1] } else { "" }
        $rc = if ($text -match 'P6_MLE_RETCODE\s+(\S+)') { $Matches[1] } else { "" }
        $qs = if ($text -match 'P6_QFIT_SUM\s+(\S+)') { $Matches[1] } else { "" }
        $it = if ($text -match 'Number of iterations\s+(\S+)') { $Matches[1] } else { "" }
        $mn = if ($text -match 'Minutes to convergence\s+(\S+)') { $Matches[1] } else { "" }
        $warn = ([regex]::Matches($text, '(?im)^.*(complex|failed to invert|warning|error).*$') | ForEach-Object { $_.Value.Trim() } | Select-Object -Unique) -join " | "
        $warn = $warn -replace ',', ';' -replace '"', "'"
        $done = $text -match 'P6_DONE'
        "$label,$k,$MaxIters,$MaxTime,$i,$HardTimeoutSeconds,$wall,$timedOut,$el,$rc,$it,$mn,$qs,`"$warn`",$done" | Out-File -Encoding ascii -Append $summary
        Write-Host "$label r$i : wall=$($wall)s timedOut=$timedOut fit=$el rc=$rc iters=$it min=$mn qsum=$qs done=$done warn=[$warn]"
    }
}
