# verify_docs_quality.ps1
#
# Public release roadmap PR-303 acceptance evidence: "Documentation checks
# run in CI and fail on a deliberately introduced stale default, missing
# command page, or archive-only broken link."
#
# The "stale default" and "missing command page" dimensions are already
# covered by scripts/verify_docs_consistency.ps1 (PR-301) and
# scripts/verify_public_api.ps1 (PR-001/PR-003) respectively -- this script
# does not duplicate them. It covers the two PR-303 Work-bullet dimensions
# nothing else in this repo checks:
#
# 1. **Required command-page sections / heading structure**: every
#    docs/command-reference/*.md page must have the exact
#    Purpose/Format/Parameters/Returns/Remarks/Examples/Source/See Also
#    heading sequence (this repo's established template, confirmed against
#    every existing page before writing this check -- found and fixed 4
#    real deviations, all "## Example" instead of "## Examples", during
#    that confirmation pass).
# 2. **Internal link and anchor integrity**: every relative markdown link
#    in README.md/docs/*.md/docs/command-reference/*.md must resolve to a
#    real file, and every `#anchor` fragment must match a real heading in
#    the target file (or the current file, for a same-page anchor), using
#    GitHub's own heading-to-anchor slug algorithm.
# 3. **Keyword-argument spelling in code snippets**: for every ```gauss
#    fenced code block, a call to a known public procedure using
#    `name=value` keyword-argument syntax has each keyword name checked
#    against that procedure's actual declared parameter list in src/ --
#    catches a typo'd keyword argument in a documented example, the
#    "validate code snippets... especially signatures with keyword
#    arguments" Work bullet.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
)

$docsDir = Join-Path $RepoRoot "docs"
$cmdRefDir = Join-Path $docsDir "command-reference"
$srcDir = Join-Path $RepoRoot "src"

# --- Helper: GitHub-style heading-to-anchor slug ---
function ConvertTo-Slug {
    param([string]$Heading)
    $s = $Heading.ToLowerInvariant()
    $s = $s -replace '[^a-z0-9 \-]', ''
    $s = $s -replace ' ', '-'
    return $s
}

function Get-Headings {
    param([string]$Text)
    $matches_ = [regex]::Matches($Text, '(?m)^#{1,6}\s+(.+?)\s*$')
    $slugs = @{}
    $result = New-Object System.Collections.Generic.List[string]
    foreach ($m in $matches_) {
        $slug = ConvertTo-Slug -Heading $m.Groups[1].Value
        if ($slugs.ContainsKey($slug)) {
            $slugs[$slug] += 1
            $slug = "$slug-$($slugs[$slug])"
        } else {
            $slugs[$slug] = 0
        }
        $result.Add($slug)
    }
    return $result
}

# --- 1. Command-reference required section structure ---

$requiredHeadings = @("Purpose", "Format", "Parameters", "Returns", "Remarks", "Examples", "Source", "See Also")
$cmdRefFiles = @(Get-ChildItem -LiteralPath $cmdRefDir -Filter "*.md" -File)
if ($cmdRefFiles.Count -eq 0) {
    throw "no command-reference pages found under $cmdRefDir"
}

$structureErrors = @()
foreach ($file in $cmdRefFiles) {
    $text = Get-Content -LiteralPath $file.FullName -Raw
    $h2Matches = [regex]::Matches($text, '(?m)^##\s+(.+?)\s*$')
    $actual = @($h2Matches | ForEach-Object { $_.Groups[1].Value })
    if (@(Compare-Object -ReferenceObject $requiredHeadings -DifferenceObject $actual -SyncWindow 0).Count -gt 0) {
        $structureErrors += "$($file.Name): expected [$($requiredHeadings -join ', ')], found [$($actual -join ', ')]"
    }
}
if ($structureErrors.Count -gt 0) {
    throw "command-reference heading structure violation(s):`n$($structureErrors -join "`n")"
}
Write-Host "verify_docs_quality.ps1: all $($cmdRefFiles.Count) command-reference pages follow the required heading structure"

# --- 2. Internal link and anchor integrity ---

$docFiles = @(Get-Item -LiteralPath (Join-Path $RepoRoot "README.md"))
$docFiles += @(Get-ChildItem -LiteralPath $docsDir -Filter "*.md" -File)
$docFiles += $cmdRefFiles

$headingCache = @{}
function Get-HeadingsCached {
    param([string]$Path)
    if (-not $headingCache.ContainsKey($Path)) {
        if (Test-Path -LiteralPath $Path) {
            $headingCache[$Path] = Get-Headings -Text (Get-Content -LiteralPath $Path -Raw)
        } else {
            $headingCache[$Path] = $null
        }
    }
    return $headingCache[$Path]
}

$linkErrors = @()
foreach ($file in $docFiles) {
    $text = Get-Content -LiteralPath $file.FullName -Raw
    $linkMatches = [regex]::Matches($text, '\[[^\]]+\]\(([^)]+)\)')
    foreach ($lm in $linkMatches) {
        $target = $lm.Groups[1].Value.Trim()
        if ($target -match '^(https?:|mailto:)') { continue }

        $pathPart = $target
        $anchor = $null
        $hashIdx = $target.IndexOf('#')
        if ($hashIdx -ge 0) {
            $pathPart = $target.Substring(0, $hashIdx)
            $anchor = $target.Substring($hashIdx + 1)
        }

        if ([string]::IsNullOrEmpty($pathPart)) {
            $targetFile = $file.FullName
        } else {
            $targetFile = (Join-Path $file.DirectoryName $pathPart)
            if (-not (Test-Path -LiteralPath $targetFile)) {
                $linkErrors += "$($file.Name): broken link target '$target' (resolved path does not exist: $targetFile)"
                continue
            }
            $targetFile = (Resolve-Path -LiteralPath $targetFile).Path
        }

        if ($anchor) {
            if ([System.IO.Path]::GetExtension($targetFile) -eq ".md") {
                $headings = Get-HeadingsCached -Path $targetFile
                if ($headings -and ($headings -notcontains $anchor)) {
                    $linkErrors += "$($file.Name): anchor '#$anchor' not found in $(Split-Path -Leaf $targetFile) (target has no matching heading)"
                }
            }
        }
    }
}
if ($linkErrors.Count -gt 0) {
    throw "broken internal link(s)/anchor(s) found:`n$($linkErrors -join "`n")"
}
Write-Host "verify_docs_quality.ps1: all internal doc links and anchors resolve"

# --- 3. Keyword-argument spelling in ```gauss code snippets ---

$srcEntries = @(Get-ChildItem -LiteralPath $srcDir -Filter "*.src" -File)
$paramListByProc = @{}
foreach ($entry in $srcEntries) {
    $text = Get-Content -LiteralPath $entry.FullName -Raw
    $declMatches = [regex]::Matches($text, '(?m)^\s*proc\s*(?:\([^)]*\)\s*=\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)\s*;')
    foreach ($dm in $declMatches) {
        $procName = $dm.Groups[1].Value
        $paramListByProc[$procName] = $dm.Groups[2].Value
    }
}

$keywordErrors = @()
foreach ($file in $docFiles) {
    $text = Get-Content -LiteralPath $file.FullName -Raw
    $fenceMatches = [regex]::Matches($text, '(?s)```gauss\r?\n(.*?)```')
    foreach ($fence in $fenceMatches) {
        $code = $fence.Groups[1].Value
        $callMatches = [regex]::Matches($code, '([A-Za-z_][A-Za-z0-9_]*)\(((?:[^()]|\([^()]*\))*)\)')
        foreach ($cm in $callMatches) {
            $procName = $cm.Groups[1].Value
            if (-not $paramListByProc.ContainsKey($procName)) { continue }
            $argText = $cm.Groups[2].Value
            $paramNames = @([regex]::Matches($paramListByProc[$procName], '([A-Za-z_][A-Za-z0-9_]*)\s*(?:=|,|$)') | ForEach-Object { $_.Groups[1].Value }) |
                Where-Object { $_ -ne "struct" }
            $kwMatches = [regex]::Matches($argText, '(?<![=<>!])\b([A-Za-z_][A-Za-z0-9_]*)\s*=(?!=)')
            foreach ($kw in $kwMatches) {
                $kwName = $kw.Groups[1].Value
                if ($paramNames.Count -gt 0 -and ($paramNames -notcontains $kwName)) {
                    $keywordErrors += "$($file.Name): '$procName(...)' example uses keyword '$kwName', not found in its actual parameter list ($($paramNames -join ', '))"
                }
            }
        }
    }
}
if ($keywordErrors.Count -gt 0) {
    throw "keyword-argument spelling error(s) found in documented examples:`n$($keywordErrors -join "`n")"
}
Write-Host "verify_docs_quality.ps1: all keyword arguments in gauss code snippets match real procedure signatures"

Write-Host "verify_docs_quality.ps1: PASS"
