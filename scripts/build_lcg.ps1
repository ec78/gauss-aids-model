# build_lcg.ps1
#
# Milestone 7: writes lib/<name>.lcg, the plain-text symbol-location catalog
# GAUSS's `library` mechanism reads to resolve/lazily compile procs from an
# installed package (this is the real, functional catalog format -- verified
# against the actual installed c:\gauss26\pkgs\qardl\lib\qardl.lcg and
# c:\gauss26\pkgs\pubtable\lib\pubtable.lcg, not a stub). Adapted from
# gauss-qardl's scripts/build_lcg.ps1; package-agnostic, so kept unchanged
# aside from comments.
#
# Run against an INSTALLED copy of the package (i.e. PackageRoot should be
# a directory containing package.json and src/, such as the staging area
# produced by build_package.ps1 + Expand-Archive, or an existing
# c:\gauss26\pkgs\<name> install) -- not against this git repo directly,
# since the repo's own package.json src array intentionally omits optional
# files like pubtable_quaids.src that should not gate `library quaids;`.

param(
    [string]$PackageRoot,
    [string]$PackageName = ""
)

if ([string]::IsNullOrWhiteSpace($PackageRoot)) {
    throw "PackageRoot is required"
}

$packagePath = Join-Path $PackageRoot "package.json"
if (-not (Test-Path -LiteralPath $packagePath)) {
    throw "package.json not found at $packagePath"
}

$pkg = Get-Content -LiteralPath $packagePath -Raw | ConvertFrom-Json
if ([string]::IsNullOrWhiteSpace($PackageName)) {
    $PackageName = [string]$pkg.name
}

if ([string]::IsNullOrWhiteSpace($PackageName)) {
    throw "package name is empty"
}

$libDir = Join-Path $PackageRoot "lib"
if (-not (Test-Path -LiteralPath $libDir)) {
    New-Item -ItemType Directory -Path $libDir | Out-Null
}

$catalogPath = Join-Path $libDir "$PackageName.lcg"
$catalog = New-Object System.Collections.Generic.List[string]

$catalog.Add("/*")
$catalog.Add("** Package: $PackageName")
$catalog.Add("** Version: $($pkg.version)")
$catalog.Add("** Author: $($pkg.author)")
$catalog.Add("** Description: $($pkg.description)")
$catalog.Add("*/")
$catalog.Add("")

foreach ($entry in @($pkg.src)) {
    $srcPath = Join-Path (Join-Path $PackageRoot "src") $entry
    if (-not (Test-Path -LiteralPath $srcPath)) {
        throw "package source file listed in package.json was not found: $srcPath"
    }

    $gaussPath = ((Resolve-Path -LiteralPath $srcPath).Path -replace "\\", "/").ToLowerInvariant()
    $catalog.Add($gaussPath)

    $lines = Get-Content -LiteralPath $srcPath
    for ($ii = 0; $ii -lt $lines.Count; $ii++) {
        $lineNo = $ii + 1
        $line = $lines[$ii]

        if ($line -match '^\s*struct\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{') {
            $catalog.Add(("    {0,-34} : definition : {1}" -f ("struct " + $matches[1]), $lineNo))
            continue
        }

        # Matches all three GAUSS proc-declaration forms seen in this
        # codebase: "proc (struct X) = name(...)" / "proc (N) = name(...)"
        # (paren-enclosed, captured to $matches[1] for the typed_returns
        # check below), "proc N = name(...)" (bare digit return count,
        # e.g. "proc 3 = quaidsElas_(...)"), and "proc name(...)" (no
        # return spec at all, e.g. "proc quantile(x, s);"). The original
        # version of this regex only matched the first form -- it silently
        # dropped every bare-digit and bare-name proc from the catalog
        # (quaidsSlutzky, _quaidsIVFirstStage, quaids(), printQuaids,
        # quaidsElas_, printQuaidsElas, quaidsElas, and a private
        # quantile() helper), which surfaced as "Undefined symbol" errors
        # from `library quaids;` despite every source-tree #include-based
        # test passing -- caught by actually installing the package and
        # running tests/package_public_api.e, not by re-reading the script.
        if ($line -match '^\s*proc\s*(?:\(([^)]*)\)|\d+)?\s*(?:=\s*)?([A-Za-z_][A-Za-z0-9_]*)\s*\(') {
            $procName = $matches[2]
            $typedSuffix = ""
            if ($matches[1] -match '^\s*struct\b') {
                $typedSuffix = " : typed_returns"
            }
            $keywordSuffix = ""
            $procCallStart = $line.IndexOf($procName + "(")
            if ($procCallStart -ge 0) {
                $argPart = $line.Substring($procCallStart + $procName.Length)
                if ($argPart -match '=') {
                    $keywordSuffix = " : keywords"
                }
            }
            $catalog.Add(("    {0,-34} : proc : {1}{2}{3}" -f $procName, $lineNo, $typedSuffix, $keywordSuffix))
        }
    }

    $catalog.Add("")
}

Set-Content -LiteralPath $catalogPath -Value $catalog -Encoding ASCII
Write-Host "build_lcg.ps1: wrote $catalogPath"
