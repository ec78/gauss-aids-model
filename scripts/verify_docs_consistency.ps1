# verify_docs_consistency.ps1
#
# Public release roadmap PR-301 acceptance evidence: "A targeted
# documentation-consistency test covers defaults, dependency status, and
# feature availability" and "Search-based review finds no conflicting
# current-state claims."
#
# This is deliberately narrower than PR-303's later, broader documentation
# quality-gate work (link/heading/section checks) -- it targets exactly the
# three real, confirmed contradictions found during PR-301's own review:
#
# 1. Defaults: docs/command-reference/quaidsControlCreate.md's own defaults
#    table is cross-checked against quaidsControlCreate()'s actual coded
#    defaults in src/quaidsutil.src -- this would have caught
#    docs/command-reference/quaidsZeroFit.md's stale claim that
#    `aCtl.homogenous = 0` is the default (the real default is `1`, set at
#    Milestone 28/PR-301's own review). Any other doc page making the same
#    kind of "`aCtl.homogenous = N` (default)" claim is checked the same way,
#    not just the one file that happened to be wrong this time.
# 2. Feature availability: docs/USAGE_GUIDE.md's Zero Budget Shares section
#    must not claim homogeneity/symmetry imposition is "out of scope" for
#    quaidsZeroFit() -- Milestone 30 added it.
# 3. Dependency/feature status: docs/USAGE_GUIDE.md's survey-workflow
#    section must not claim replicate-weight (BRR/jackknife) variance
#    "remains a roadmap item" -- quaidsReplicateWeightFit() (Milestone 27)
#    already ships it; only formal strata and finite-population correction
#    are still open.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
)

$utilPath = Join-Path $RepoRoot "src\quaidsutil.src"
$controlDocPath = Join-Path $RepoRoot "docs\command-reference\quaidsControlCreate.md"
$usageGuidePath = Join-Path $RepoRoot "docs\USAGE_GUIDE.md"

foreach ($p in @($utilPath, $controlDocPath, $usageGuidePath)) {
    if (-not (Test-Path -LiteralPath $p)) {
        throw "required file not found: $p"
    }
}

# --- 1. Defaults: quaidsControlCreate()'s coded defaults vs. its own documented table ---

$utilText = Get-Content -LiteralPath $utilPath -Raw
$createBodyMatch = [regex]::Match($utilText, '(?s)proc\s*\(struct quaidsControl\)\s*=\s*quaidsControlCreate\(\);.*?endp;')
if (-not $createBodyMatch.Success) {
    throw "could not locate quaidsControlCreate() body in src/quaidsutil.src"
}
$createBody = $createBodyMatch.Value

$fieldsToCheck = @("linear", "maxiter", "homogenous", "alpha0", "err", "relax")
$codedDefaults = @{}
foreach ($field in $fieldsToCheck) {
    $m = [regex]::Match($createBody, "aCtl\.$field\s*=\s*([^;]+);")
    if (-not $m.Success) {
        throw "quaidsControlCreate() in src/quaidsutil.src does not assign aCtl.$field"
    }
    $codedDefaults[$field] = $m.Groups[1].Value.Trim()
}

$controlDocText = Get-Content -LiteralPath $controlDocPath -Raw
foreach ($field in $fieldsToCheck) {
    $rowMatch = [regex]::Match($controlDocText, "\|\s*``$field``\s*\|\s*``([^``]+)``\s*\|")
    if (-not $rowMatch.Success) {
        throw "docs/command-reference/quaidsControlCreate.md has no defaults-table row for '$field'"
    }
    $documented = $rowMatch.Groups[1].Value.Trim()
    if ($documented -ne $codedDefaults[$field]) {
        throw "default mismatch for aCtl.${field}: src/quaidsutil.src sets '$($codedDefaults[$field])', docs/command-reference/quaidsControlCreate.md documents '$documented'"
    }
}
Write-Host "verify_docs_consistency.ps1: quaidsControlCreate() defaults table matches src/quaidsutil.src for $($fieldsToCheck.Count) fields"

# --- 2. No other doc page claims a stale aCtl.homogenous default ---

$docsToScan = @(Get-ChildItem -LiteralPath (Join-Path $RepoRoot "docs") -Filter "*.md" -Recurse -File)
$docsToScan += Get-Item -LiteralPath (Join-Path $RepoRoot "README.md")

$codedHomogenousDefault = $codedDefaults["homogenous"]
$staleDefaultClaims = @()
foreach ($file in $docsToScan) {
    $text = Get-Content -LiteralPath $file.FullName -Raw
    $claimMatches = [regex]::Matches($text, "``aCtl\.homogenous\s*=\s*(\d)``[^\r\n]{0,25}\(default\)")
    foreach ($cm in $claimMatches) {
        $claimedValue = $cm.Groups[1].Value
        if ($claimedValue -ne $codedHomogenousDefault) {
            $staleDefaultClaims += "$($file.Name): claims default '$claimedValue', actual coded default is '$codedHomogenousDefault'"
        }
    }
}
if ($staleDefaultClaims.Count -gt 0) {
    throw "stale aCtl.homogenous default claim(s) found:`n$($staleDefaultClaims -join "`n")"
}
Write-Host "verify_docs_consistency.ps1: no doc page claims a stale aCtl.homogenous default"

# --- 3. Feature-availability regression guards (PR-301's own confirmed fixes) ---

$usageGuideText = Get-Content -LiteralPath $usageGuidePath -Raw

if ($usageGuideText -match "(?s)out\s+of\s+scope[\s\S]{0,20}\(no\s+homogeneity") {
    throw "docs/USAGE_GUIDE.md still claims quaidsZeroFit() homogeneity/symmetry imposition is out of scope -- stale since Milestone 30 added it"
}

if ($usageGuideText -match "(?s)replicate-weight[\s\S]{0,120}remain\s+roadmap\s+items") {
    throw "docs/USAGE_GUIDE.md still claims replicate-weight variance remains a roadmap item -- stale since quaidsReplicateWeightFit() (Milestone 27) already ships it"
}

Write-Host "verify_docs_consistency.ps1: no stale zero-share-scope or replicate-weight-scope claims found"
Write-Host "verify_docs_consistency.ps1: PASS"
