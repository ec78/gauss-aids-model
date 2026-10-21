param(
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$Equivalence
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = (Resolve-Path (Join-Path $here "..\..\..")).Path
$cfg = Join-Path $repo "tests\gauss26_cfg_override"
$oldCfg = $env:GAUSS26_CFG
function Invoke-D1Gauss([string]$script) {
    $token = [Guid]::NewGuid().ToString("N")
    $stdout = Join-Path $env:TEMP ("tvpdiag-" + $token + ".out")
    $stderr = Join-Path $env:TEMP ("tvpdiag-" + $token + ".err")
    try {
        $p = Start-Process -FilePath $GaussExe -ArgumentList @("-b", $script) `
            -NoNewWindow -Wait -PassThru -RedirectStandardOutput $stdout `
            -RedirectStandardError $stderr
        $text = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
        return [pscustomobject]@{ ExitCode = $p.ExitCode; Text = $text }
    } finally {
        if (Test-Path $stdout) { Remove-Item -LiteralPath $stdout -Force }
        if (Test-Path $stderr) { Remove-Item -LiteralPath $stderr -Force }
    }
}
try {
    $env:GAUSS26_CFG = $cfg
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
    } finally { Pop-Location }
} finally {
    $env:GAUSS26_CFG = $oldCfg
}
