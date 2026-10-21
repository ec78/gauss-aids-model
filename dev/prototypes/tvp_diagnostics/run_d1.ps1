param(
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$Equivalence,
    [switch]$ProbeN12,
    [int]$TimeoutSeconds = 900
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$cfg = Join-Path $repo "tests\gauss26_cfg_override"
$oldCfg = $env:GAUSS26_CFG
function Get-NormHash([string]$path) {
    $text = [IO.File]::ReadAllText($path) -replace "`r", ""
    $bytes = [Text.Encoding]::UTF8.GetBytes($text)
    $sha = [Security.Cryptography.SHA256]::Create()
    return (($sha.ComputeHash($bytes) | ForEach-Object { $_.ToString("x2") }) -join "")
}
function Test-D1Freeze {
    foreach ($line in Get-Content (Join-Path $here "FREEZE.sha256")) {
        $parts = $line -split "  ", 2
        if ((Get-NormHash (Join-Path $here $parts[1])) -ne $parts[0]) {
            throw "D1 frozen-file hash mismatch: $($parts[1])"
        }
    }
}
function Invoke-D1Gauss([string]$script) {
    $token = [Guid]::NewGuid().ToString("N")
    $stdout = Join-Path $env:TEMP ("tvpdiag-" + $token + ".out")
    $stderr = Join-Path $env:TEMP ("tvpdiag-" + $token + ".err")
    try {
        $p = Start-Process -FilePath $GaussExe -ArgumentList @("-b", $script) `
            -NoNewWindow -PassThru -RedirectStandardOutput $stdout `
            -RedirectStandardError $stderr
        # Windows PowerShell 5.1 does not reliably populate ExitCode after
        # WaitForExit unless the process handle is materialized first.
        $processHandle = $p.Handle
        if (-not $p.WaitForExit($TimeoutSeconds * 1000)) {
            & taskkill.exe /PID $p.Id /T /F | Out-Null
            throw "GAUSS job timed out after $TimeoutSeconds seconds: $script"
        }
        $p.WaitForExit()
        $text = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
        return [pscustomobject]@{ ExitCode = $p.ExitCode; Text = $text }
    } finally {
        if (Test-Path $stdout) { Remove-Item -LiteralPath $stdout -Force }
        if (Test-Path $stderr) { Remove-Item -LiteralPath $stderr -Force }
    }
}
try {
    $env:GAUSS26_CFG = $cfg
    Test-D1Freeze
    & powershell -ExecutionPolicy Bypass -File (Join-Path $repo "scripts\verify_sslib_pin.ps1") -Strict
    if ($LASTEXITCODE -ne 0) { throw "strict SSLib pin verification failed" }
    Push-Location $here
    try {
        $unit = Invoke-D1Gauss "tvpdiag_unit_test.e"
        Write-Host $unit.Text
        if ($unit.ExitCode -ne 0 -or $unit.Text -notmatch "TVPDIAG_FAILURES\s+0" -or
            $unit.Text -notmatch "TVPDIAG_UNIT_PASS") {
            throw "tvpdiag_unit_test.e failed"
        }
        if ($Equivalence) {
            $equiv = Invoke-D1Gauss "tvpdiag_equiv.e"
            Write-Host $equiv.Text
            if ($equiv.ExitCode -ne 0 -or $equiv.Text -notmatch "TVPDIAG_EQ_PASS\s+1") {
                throw "tvpdiag_equiv.e failed"
            }
        }
        if ($ProbeN12) {
            $probe1 = Invoke-D1Gauss "results\gen\probe_n12.e"
            $probe2 = Invoke-D1Gauss "results\gen\probe_n12.e"
            $tag1 = (($probe1.Text -split "`r?`n") |
                Where-Object { $_ -match "^TVPDIAG_PROBE" } |
                ForEach-Object { $_.TrimEnd() }) -join "`n"
            $tag2 = (($probe2.Text -split "`r?`n") |
                Where-Object { $_ -match "^TVPDIAG_PROBE" } |
                ForEach-Object { $_.TrimEnd() }) -join "`n"
            if ($probe1.ExitCode -ne 0 -or $probe2.ExitCode -ne 0 -or
                $tag1 -notmatch "TVPDIAG_PROBE_DONE" -or $tag1 -ne $tag2) {
                throw "n1=2 probe failed deterministic tagged-output check"
            }
            $resultDir = Join-Path $here "results"
            $utf8NoBom = New-Object Text.UTF8Encoding($false)
            [IO.File]::WriteAllText((Join-Path $resultDir "probe_n12_run1.txt"), $tag1 + "`n", $utf8NoBom)
            [IO.File]::WriteAllText((Join-Path $resultDir "probe_n12_run2.txt"), $tag2 + "`n", $utf8NoBom)
            $h1 = (Get-FileHash (Join-Path $resultDir "probe_n12_run1.txt") -Algorithm SHA256).Hash.ToLowerInvariant()
            $h2 = (Get-FileHash (Join-Path $resultDir "probe_n12_run2.txt") -Algorithm SHA256).Hash.ToLowerInvariant()
            [IO.File]::WriteAllText((Join-Path $resultDir "probe_n12.sha256"),
                "$h1  probe_n12_run1.txt`n$h2  probe_n12_run2.txt`n", $utf8NoBom)
            Write-Host "TVPDIAG_PROBE_DETERMINISTIC PASS $h1"
            Write-Host $tag1
        }
    } finally { Pop-Location }
} finally {
    $env:GAUSS26_CFG = $oldCfg
}
