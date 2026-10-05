# pilot_resume_test.ps1 -- V6 resume check for the pilot harness.
# Run from dev/experiments/convergence-anderson:
#   powershell -ExecutionPolicy Bypass -File pilot_resume_test.ps1
$ErrorActionPreference = 'Stop'
$dir = 'pilot_resume_test'
$csv = Join-Path $dir 'cells_m0.csv'
if (Test-Path $dir) { Remove-Item -Recurse -Force $dir }
New-Item -ItemType Directory $dir | Out-Null
$tg = 'C:\gauss26\tgauss.exe'
function RunJob($tag) {
    & $tg -b -x pilot_resume_test.e | Out-File -Encoding utf8 (Join-Path $dir "run_$tag.txt")
    $done = Select-String -Path (Join-Path $dir "run_$tag.txt") -Pattern 'PILOT CORE DONE'
    if (-not $done) { throw "run $tag did not finish" }
    return $done.Line.Trim()
}
$fail = 0
function Check($name, $cond, $detail) {
    if ($cond) { Write-Output "PASS V6 $name  $detail" } else { Write-Output "FAIL V6 $name  $detail"; $script:fail++ }
}
$r1 = RunJob 'first'
$a = Get-Content $csv
Check 'first run writes 6 cells' ($a.Count -eq 6) "rows=$($a.Count); $r1"
$r2 = RunJob 'rerun'
$b = Get-Content $csv
Check 'rerun adds no rows' ($b.Count -eq 6) "rows=$($b.Count); $r2"
Check 'rerun leaves file byte-identical' ((($a -join "`n") -eq ($b -join "`n"))) ''
# simulate an interruption: drop the last 2 finished cells
$a[0..3] | Set-Content -Encoding ascii $csv
$r3 = RunJob 'resume'
$c = Get-Content $csv
Check 'resume restores 6 rows' ($c.Count -eq 6) "rows=$($c.Count); $r3"
$keys = $c | ForEach-Object { ($_ -split ',')[0..3] -join ',' }
Check 'keys unique after resume' (($keys | Sort-Object -Unique).Count -eq 6) ''
# recomputed cells identical apart from wallTime (column 25, 1-based)
$same = $true
for ($i = 4; $i -le 5; $i++) {
    $x = $a[$i] -split ','; $y = $c[$i] -split ','
    for ($j = 0; $j -lt $x.Count; $j++) { if ($j -ne 24 -and $x[$j] -ne $y[$j]) { $same = $false } }
    if ($x.Count -ne $y.Count) { $same = $false }
}
Check 'recomputed cells identical except wallTime' $same ''
if ($fail -eq 0) { Write-Output 'RESUME TEST: ALL 6 CHECKS PASSED' } else { Write-Output "RESUME TEST: $fail CHECKS FAILED" }
