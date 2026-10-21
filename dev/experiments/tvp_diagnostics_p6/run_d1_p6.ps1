param(
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [int]$TimeoutSeconds = 3600
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$d1 = Join-Path $repo "dev\prototypes\tvp_diagnostics"
$cfg = Join-Path $repo "tests\gauss26_cfg_override"
$oldCfg = $env:GAUSS26_CFG

function Get-NormHash([string]$path) {
    $content = [IO.File]::ReadAllText($path) -replace "`r", ""
    $bytes = [Text.Encoding]::UTF8.GetBytes($content)
    $sha = [Security.Cryptography.SHA256]::Create()
    return (($sha.ComputeHash($bytes) | ForEach-Object { $_.ToString("x2") }) -join "")
}

function Test-NormManifest([string]$manifest, [string]$root) {
    foreach ($line in Get-Content $manifest) {
        $parts = $line -split "  ", 2
        if ((Get-NormHash (Join-Path $root $parts[1])) -ne $parts[0]) {
            throw "normalized hash mismatch: $($parts[1])"
        }
    }
}

function Test-D1Baseline {
    foreach ($line in Get-Content (Join-Path $here "D1_BASELINE.sha256")) {
        $parts = $line -split "  ", 2
        $actual = (Get-FileHash (Join-Path $repo $parts[1]) -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($actual -ne $parts[0]) { throw "D1 byte baseline mismatch: $($parts[1])" }
    }
}

function Test-Preflight {
    Test-NormManifest (Join-Path $here "FREEZE.sha256") $here
    Test-NormManifest (Join-Path $d1 "FREEZE.sha256") $d1
    Test-D1Baseline
    & powershell -ExecutionPolicy Bypass -File (Join-Path $repo "scripts\verify_sslib_pin.ps1") -Strict
    if ($LASTEXITCODE -ne 0) { throw "strict SSLib pin verification failed" }
}

function Invoke-Gauss([string]$script) {
    $token = [Guid]::NewGuid().ToString("N")
    $stdout = Join-Path $env:TEMP ("tvpdiag-p6-" + $token + ".out")
    $stderr = Join-Path $env:TEMP ("tvpdiag-p6-" + $token + ".err")
    try {
        $p = Start-Process -FilePath $GaussExe -ArgumentList @("-b", $script) `
            -NoNewWindow -PassThru -RedirectStandardOutput $stdout `
            -RedirectStandardError $stderr
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
    Test-Preflight
    Push-Location $d1
    try {
        $unit = Invoke-Gauss "tvpdiag_unit_test.e"
        Write-Host $unit.Text
        if ($unit.ExitCode -ne 0 -or $unit.Text -notmatch "TVPDIAG_FAILURES\s+0" -or
            $unit.Text -notmatch "TVPDIAG_UNIT_PASS") {
            throw "accepted D1 unit test failed"
        }

        $probe = "..\..\experiments\tvp_diagnostics_p6\probe_p6.e"
        $run1 = Invoke-Gauss $probe
        $run2 = Invoke-Gauss $probe
        $tag1 = (($run1.Text -split "`r?`n") |
            Where-Object { $_ -match "^TVPDIAG_P6" } |
            ForEach-Object { $_.TrimEnd() }) -join "`n"
        $tag2 = (($run2.Text -split "`r?`n") |
            Where-Object { $_ -match "^TVPDIAG_P6" } |
            ForEach-Object { $_.TrimEnd() }) -join "`n"
        if ($run1.ExitCode -ne 0 -or $run2.ExitCode -ne 0 -or
            $tag1 -notmatch "TVPDIAG_P6_DONE" -or
            $tag1 -notmatch "TVPDIAG_P6_ACCEPT\s+1" -or
            $tag1 -ne $tag2) {
            throw "P6 probe failed acceptance or deterministic payload check"
        }
    } finally {
        Pop-Location
    }

    $resultDir = Join-Path $here "results"
    if (-not (Test-Path $resultDir)) { New-Item -ItemType Directory -Path $resultDir | Out-Null }
    $utf8NoBom = New-Object Text.UTF8Encoding($false)
    $path1 = Join-Path $resultDir "probe_p6_run1.txt"
    $path2 = Join-Path $resultDir "probe_p6_run2.txt"
    [IO.File]::WriteAllText($path1, $tag1 + "`n", $utf8NoBom)
    [IO.File]::WriteAllText($path2, $tag2 + "`n", $utf8NoBom)
    $h1 = (Get-FileHash $path1 -Algorithm SHA256).Hash.ToLowerInvariant()
    $h2 = (Get-FileHash $path2 -Algorithm SHA256).Hash.ToLowerInvariant()
    [IO.File]::WriteAllText((Join-Path $resultDir "probe_p6.sha256"),
        "$h1  probe_p6_run1.txt`n$h2  probe_p6_run2.txt`n", $utf8NoBom)

    Test-Preflight
    Write-Host "TVPDIAG_P6_DETERMINISTIC PASS $h1"
    Write-Host $tag1
} finally {
    $env:GAUSS26_CFG = $oldCfg
}
