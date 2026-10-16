# verify_runner_timeouts.ps1
#
# TVP-AIDS promotion P7 follow-up: negative control for the per-GAUSS-job
# timeout in every test/release runner that launches tgauss. It proves
# that a hung job is killed together with its whole process tree and
# reported as a failure, instead of stalling the pipeline:
#
#   1. Unit level -- the COMMITTED Invoke-GaussBatch function is parsed
#      out of each of tests/run_source_tests.ps1, tests/
#      run_examples_smoke.ps1 and scripts/run_release_verification.ps1
#      (PowerShell AST, so no copy can drift from the real code) and run
#      against a scratch GAUSS script that blocks on a grandchild
#      powershell.exe. Each must return TimedOut = $true near the limit,
#      append the TIMEOUT line, and leave neither tgauss nor the grandchild
#      running; the grandchild first writes a marker file, so "it was
#      killed" cannot pass vacuously because it never started. A quick
#      script must still pass untouched with the timeout on, and with
#      0 (disabled). The script-scope $TestTimeoutSeconds default must
#      propagate when -TimeoutSeconds is not passed. Finally, an escaped
#      pipe holder (an orphaned grandchild outside the tree that keeps
#      stdout/stderr open) must not reintroduce an unbounded wait: both a
#      timed-out and a normally exiting job must return within the
#      runner's 30 s bounded drain, with a NOTE line, the job's own tree
#      killed and its output kept. The orphan is outside the runner's
#      reach by design, so this script kills it itself.
#   2. Integration -- run_examples_smoke.ps1 itself, pointed at the scratch
#      directory (-ExamplesDir/-ExampleList), must exit 1 with the hung
#      script reported FAIL (timed out) and the quick one PASS.
#   3. Wiring -- run_release_gate.ps1 must pass -TestTimeoutSeconds to
#      run_release_verification.ps1, which must pass it to both
#      run_source_tests.ps1 and run_examples_smoke.ps1.
#
# Safe by construction: everything runs from a fresh temp directory, only
# processes whose command line carries this run's unique marker are ever
# inspected, and the runners' own taskkill targets only the PID they
# started. Takes about 4 minutes. Not part of run_source_tests.ps1 (it
# deliberately hangs jobs); run it after changing any runner's
# Invoke-GaussBatch:
#   powershell -ExecutionPolicy Bypass -File tests\verify_runner_timeouts.ps1

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [int]$TimeoutSeconds = 8
)

$ErrorActionPreference = "Stop"

$failures = 0
function Assert-That {
    param([bool]$Condition, [string]$Message)
    if ($Condition) {
        Write-Host "PASS  $Message"
    } else {
        Write-Host "FAIL  $Message"
        $script:failures++
    }
}

$token = "QRT" + [System.Guid]::NewGuid().ToString("N").Substring(0, 12)
$work = Join-Path ([System.IO.Path]::GetTempPath()) ("quaids_timeout_ctl_" + $token)
New-Item -ItemType Directory -Path $work | Out-Null

function Get-MarkedProcesses {
    @(Get-CimInstance Win32_Process |
        Where-Object { $_.CommandLine -and $_.CommandLine.Contains($token) })
}

function New-HangScript {
    param([string]$Name)
    $marker = (Join-Path $work "$Name.started") -replace "\\", "/"
    # exec() blocks until its child exits, so tgauss hangs with a live
    # grandchild -- the case taskkill /T exists for. The child writes the
    # marker first, proving the grandchild really started.
    $gauss = 'ret = exec("powershell.exe", "-NoProfile -Command Set-Content -Path ''' + $marker +
        ''' -Value x; Start-Sleep -Seconds 900 # ' + $token + '");' + "`r`n" + 'print "unreachable";'
    Set-Content -LiteralPath (Join-Path $work "$Name.e") -Value $gauss -Encoding ascii
    return (Join-Path $work "$Name.e")
}

function New-EscapeScript {
    param([string]$Name, [bool]$Hang)
    $marker = (Join-Path $work "$Name.started") -replace "\\", "/"
    # tgauss exec()s powershell A, which starts B with the inherited std
    # handles and exits at once. B is then an orphan OUTSIDE tgauss's
    # process tree (taskkill /T cannot reach it) that holds the runner's
    # stdout/stderr pipes open -- the case the runners' bounded drain
    # exists for. With $Hang, tgauss then blocks on an in-tree child too.
    $inner = "Set-Content -Path '$marker' -Value x; Start-Sleep -Seconds 900 # $token-esc"
    $a = "-NoProfile -Command Start-Process -NoNewWindow powershell.exe -ArgumentList '-NoProfile','-Command','" +
        ($inner -replace "'", "''") + "'"
    $lines = @('ret = exec("powershell.exe", "' + ($a -replace '"', '\"') + '");')
    if ($Hang) {
        $lines += 'ret = exec("powershell.exe", "-NoProfile -Command Start-Sleep -Seconds 900 # ' + $token + '-tree");'
    }
    $lines += 'print "after escape";'
    Set-Content -LiteralPath (Join-Path $work "$Name.e") -Value $lines -Encoding ascii
    return (Join-Path $work "$Name.e")
}

# Kills this run's escaped pipe holders (unreachable by the runners by
# design) and returns how many there were.
function Stop-EscapedProcesses {
    $esc = @(Get-MarkedProcesses | Where-Object { $_.CommandLine.Contains("$token-esc") })
    foreach ($p in $esc) { & taskkill.exe /PID $p.ProcessId /T /F 2>&1 | Out-Null }
    return $esc.Count
}

Set-Content -LiteralPath (Join-Path $work "quick.e") -Value 'print "quick done";' -Encoding ascii
$quick = Join-Path $work "quick.e"

try {
    # ---- 1. Unit level: each runner's committed Invoke-GaussBatch ----
    $runners = @(
        "tests\run_source_tests.ps1",
        "tests\run_examples_smoke.ps1",
        "scripts\run_release_verification.ps1"
    )

    foreach ($runner in $runners) {
        $path = Join-Path $RepoRoot $runner
        $parseErrors = $null
        $ast = [System.Management.Automation.Language.Parser]::ParseFile($path, [ref]$null, [ref]$parseErrors)
        Assert-That ($parseErrors.Count -eq 0) "${runner}: parses cleanly"
        $fn = $ast.Find({
            param($n)
            $n -is [System.Management.Automation.Language.FunctionDefinitionAst] -and $n.Name -eq "Invoke-GaussBatch"
        }, $true)
        if ($null -eq $fn) {
            Assert-That $false "${runner}: defines Invoke-GaussBatch"
            continue
        }

        # Each runner's function in its own module scope, so the copies
        # cannot shadow each other and $TestTimeoutSeconds is the
        # module-scope value the real script's param would provide.
        $mod = New-Module -ScriptBlock ([scriptblock]::Create(
            '$TestTimeoutSeconds = 0' + "`n" + $fn.Extent.Text + "`n" +
            'function Set-RunnerTimeout([int]$s) { $script:TestTimeoutSeconds = $s }' + "`n" +
            'Export-ModuleMember -Function Invoke-GaussBatch, Set-RunnerTimeout'))
        $hasWorkDir = $fn.Extent.Text -match '\$WorkingDirectory'
        $common = @{ Exe = $GaussExe }
        if ($hasWorkDir) { $common.WorkingDirectory = $work }

        $hang = New-HangScript ("hang_" + [System.IO.Path]::GetFileNameWithoutExtension($runner))
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        $r = & $mod { param($a, $t) Invoke-GaussBatch @a -TimeoutSeconds $t } ($common + @{ Arguments = @("-b", "-x", $hang) }) $TimeoutSeconds
        $secs = $sw.Elapsed.TotalSeconds
        Start-Sleep -Seconds 2
        Assert-That ($r.TimedOut -eq $true) "${runner}: hung job flagged TimedOut"
        Assert-That ($secs -ge $TimeoutSeconds -and $secs -lt $TimeoutSeconds + 40) ("${runner}: returned at the limit ({0:N1} s for a {1} s timeout)" -f $secs, $TimeoutSeconds)
        Assert-That ($r.Output -match "TIMEOUT: killed after $TimeoutSeconds s") "${runner}: TIMEOUT line appended to the output"
        Assert-That (Test-Path -LiteralPath ($hang -replace "\.e$", ".started")) "${runner}: grandchild process really started (marker written)"
        Assert-That ((Get-MarkedProcesses).Count -eq 0) "${runner}: no tgauss/grandchild process of this run left alive"

        & $mod { param($s) Set-RunnerTimeout $s } 5
        $hang2 = New-HangScript ("hangdefault_" + [System.IO.Path]::GetFileNameWithoutExtension($runner))
        $sw.Restart()
        $r = & $mod { param($a) Invoke-GaussBatch @a } ($common + @{ Arguments = @("-b", "-x", $hang2) })
        $secs = $sw.Elapsed.TotalSeconds
        Start-Sleep -Seconds 2
        Assert-That ($r.TimedOut -eq $true -and $secs -lt 45) ("${runner}: script-level -TestTimeoutSeconds default applies when -TimeoutSeconds is omitted ({0:N1} s)" -f $secs)
        Assert-That ((Get-MarkedProcesses).Count -eq 0) "${runner}: nothing left alive after the default-timeout kill"

        foreach ($t in @(60, 0)) {
            $r = & $mod { param($a, $t) Invoke-GaussBatch @a -TimeoutSeconds $t } ($common + @{ Arguments = @("-b", "-x", $quick) }) $t
            Assert-That ($r.TimedOut -eq $false -and $r.Output -match "quick done") "${runner}: a quick job completes normally with timeout $t"
        }
        # Escaped pipe holder: the runner must still return within the
        # timeout plus its 30 s bounded drain, not wait on the orphan.
        $drainLimit = 30
        $esc = New-EscapeScript ("esc_hang_" + [System.IO.Path]::GetFileNameWithoutExtension($runner)) $true
        $sw.Restart()
        $r = & $mod { param($a, $t) Invoke-GaussBatch @a -TimeoutSeconds $t } ($common + @{ Arguments = @("-b", "-x", $esc) }) $TimeoutSeconds
        $secs = $sw.Elapsed.TotalSeconds
        Start-Sleep -Seconds 2
        Assert-That ($r.TimedOut -eq $true -and $secs -ge $TimeoutSeconds -and $secs -lt $TimeoutSeconds + $drainLimit + 25) ("${runner}: hung job with an escaped pipe holder still returns, bounded ({0:N1} s)" -f $secs)
        Assert-That ($r.Output -match "TIMEOUT: killed after $TimeoutSeconds s" -and $r.Output -match "NOTE: stdout/stderr still held open") "${runner}: TIMEOUT and held-open NOTE lines appended"
        Assert-That (@(Get-MarkedProcesses | Where-Object { -not $_.CommandLine.Contains("$token-esc") }).Count -eq 0) "${runner}: tgauss and its in-tree child killed"
        Assert-That ((Stop-EscapedProcesses) -ge 1) "${runner}: the escaped pipe holder really existed outside the tree (non-vacuous)"

        $esc = New-EscapeScript ("esc_exit_" + [System.IO.Path]::GetFileNameWithoutExtension($runner)) $false
        $sw.Restart()
        $r = & $mod { param($a, $t) Invoke-GaussBatch @a -TimeoutSeconds $t } ($common + @{ Arguments = @("-b", "-x", $esc) }) 60
        $secs = $sw.Elapsed.TotalSeconds
        Assert-That ($r.TimedOut -eq $false -and $r.ExitCode -ne $null -and $secs -lt $drainLimit + 25) ("${runner}: normal exit with an escaped pipe holder returns, bounded ({0:N1} s)" -f $secs)
        Assert-That ($r.Output -match "after escape" -and $r.Output -match "NOTE: stdout/stderr still held open") "${runner}: that job's output is kept and the NOTE appended"
        Assert-That ((Stop-EscapedProcesses) -ge 1) "${runner}: escaped pipe holder existed for the normal-exit case"
        Start-Sleep -Seconds 2
        Assert-That ((Get-MarkedProcesses).Count -eq 0) "${runner}: nothing left alive after the escape cases"

        Remove-Module $mod
    }

    # ---- 2. Integration: run_examples_smoke.ps1 end to end ----
    $hangSmoke = New-HangScript "hang_smoke"
    $smokeOut = & powershell -ExecutionPolicy Bypass -File (Join-Path $RepoRoot "tests\run_examples_smoke.ps1") `
        -RepoRoot $RepoRoot -GaussExe $GaussExe -ExamplesDir $work `
        -ExampleList ([System.IO.Path]::GetFileName($hangSmoke)),"quick.e" -TestTimeoutSeconds $TimeoutSeconds 2>&1
    $smokeExit = $LASTEXITCODE
    $smokeText = ($smokeOut | Out-String)
    Start-Sleep -Seconds 2
    Assert-That ($smokeExit -eq 1) "run_examples_smoke.ps1: exits 1 when an example times out"
    Assert-That ($smokeText -match "FAIL  hang_smoke\.e \(timed out after $TimeoutSeconds s\)") "run_examples_smoke.ps1: the hung example is reported FAIL (timed out)"
    Assert-That ($smokeText -match "PASS  quick\.e") "run_examples_smoke.ps1: the next example still runs and passes"
    Assert-That ((Get-MarkedProcesses).Count -eq 0) "run_examples_smoke.ps1: no process of this run left alive"

    # ---- 3. Wiring through the release pipeline ----
    $verification = Get-Content -LiteralPath (Join-Path $RepoRoot "scripts\run_release_verification.ps1") -Raw
    $gate = Get-Content -LiteralPath (Join-Path $RepoRoot "scripts\run_release_gate.ps1") -Raw
    Assert-That ($verification -match 'run_source_tests\.ps1"\)[^\r\n]*-TestTimeoutSeconds \$TestTimeoutSeconds') "run_release_verification.ps1 passes -TestTimeoutSeconds to run_source_tests.ps1"
    Assert-That ($verification -match 'run_examples_smoke\.ps1"\)[^\r\n]*-TestTimeoutSeconds \$TestTimeoutSeconds') "run_release_verification.ps1 passes -TestTimeoutSeconds to run_examples_smoke.ps1"
    Assert-That ($gate -match '"-TestTimeoutSeconds", \[string\]\$TestTimeoutSeconds') "run_release_gate.ps1 passes -TestTimeoutSeconds to run_release_verification.ps1"
} finally {
    # Belt and braces: never leave a marked process behind, even if an
    # assertion above threw.
    foreach ($p in (Get-MarkedProcesses)) {
        & taskkill.exe /PID $p.ProcessId /T /F 2>&1 | Out-Null
    }
    Remove-Item -LiteralPath $work -Recurse -Force -ErrorAction SilentlyContinue
}

if ($failures -gt 0) {
    Write-Host "verify_runner_timeouts.ps1: FAIL ($failures checks failed)"
    exit 1
}
Write-Host "verify_runner_timeouts.ps1: PASS"
