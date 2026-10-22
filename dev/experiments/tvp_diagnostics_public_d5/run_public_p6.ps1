param(
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [int]$TimeoutSeconds = 2400
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$private = Join-Path $repo "dev\experiments\tvp_diagnostics_p6"
$privatePayload = Join-Path $private "results\probe_p6_run1.txt"
$expectedPrivateHash = "b33836572df51a196449208059026b6f26eaa295cdd74bc2742af9a8761090d4"
$cfg = Join-Path $repo "tests\gauss26_cfg_override"
$oldCfg = $env:GAUSS26_CFG

function Get-NormHash([string]$path) {
    $content = [IO.File]::ReadAllText($path) -replace "`r", ""
    $bytes = [Text.Encoding]::UTF8.GetBytes($content)
    $sha = [Security.Cryptography.SHA256]::Create()
    return (($sha.ComputeHash($bytes) | ForEach-Object { $_.ToString("x2") }) -join "")
}

function Test-Freeze {
    foreach ($line in Get-Content (Join-Path $here "FREEZE.sha256")) {
        $parts = $line -split "  ", 2
        if ((Get-NormHash (Join-Path $here $parts[1])) -ne $parts[0]) {
            throw "D5 normalized hash mismatch: $($parts[1])"
        }
    }
    $actual = Get-NormHash $privatePayload
    if ($actual -ne $expectedPrivateHash) {
        throw "frozen private P6 payload mismatch: $actual"
    }
    & git -C $repo diff --exit-code -- dev/experiments/tvp_diagnostics_p6
    if ($LASTEXITCODE -ne 0) { throw "frozen private P6 evidence has a tracked diff" }
}

function Test-Preflight {
    Test-Freeze
    & powershell -ExecutionPolicy Bypass -File (Join-Path $repo "scripts\verify_sslib_pin.ps1") -Strict
    if ($LASTEXITCODE -ne 0) { throw "strict SSLib pin verification failed" }
}

function Invoke-Gauss([string]$script) {
    $token = [Guid]::NewGuid().ToString("N")
    $stdout = Join-Path $env:TEMP ("tvpdiag-d5-" + $token + ".out")
    $stderr = Join-Path $env:TEMP ("tvpdiag-d5-" + $token + ".err")
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
    Push-Location $here
    try {
        $run = Invoke-Gauss "probe_public_p6.e"
        Write-Host $run.Text
        $tagged = (($run.Text -split "`r?`n") |
            Where-Object { $_ -match "^TVPDIAG_D5" } |
            ForEach-Object { $_.TrimEnd() }) -join "`n"
        if ($run.ExitCode -ne 0 -or $tagged -notmatch "TVPDIAG_D5_DONE" -or
            $tagged -notmatch "TVPDIAG_D5_ACCEPT\s+1") {
            throw "D5 public P6 probe failed acceptance"
        }
    } finally {
        Pop-Location
    }

    $resultDir = Join-Path $here "results"
    if (-not (Test-Path $resultDir)) { New-Item -ItemType Directory -Path $resultDir | Out-Null }
    $utf8NoBom = New-Object Text.UTF8Encoding($false)
    $result = Join-Path $resultDir "public_p6.txt"
    [IO.File]::WriteAllText($result, $tagged + "`n", $utf8NoBom)
    $hash = (Get-FileHash $result -Algorithm SHA256).Hash.ToLowerInvariant()
    [IO.File]::WriteAllText((Join-Path $resultDir "public_p6.sha256"),
        "$hash  public_p6.txt`n", $utf8NoBom)

    Test-Preflight
    Write-Host "TVPDIAG_D5_PUBLIC_PASS $hash"
    Write-Host $tagged
} finally {
    $env:GAUSS26_CFG = $oldCfg
}
