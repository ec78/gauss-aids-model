# run_source_tests.ps1
#
# Milestone 7: runs the package-manifest consistency check, then every
# source-tree test (#include-based, not library-based) in tests/, in one
# shot. Adapted from gauss-qardl's tests/run_source_tests.ps1.
#
# Public release roadmap PR-001/PR-003: also runs
# scripts/verify_public_api.ps1, which checks release-metadata version
# consistency (package.json/CITATION.cff/docs/public-api.json/CHANGELOG.md)
# and reconciles docs/public-api.json's procedure/struct inventory against
# src/ and docs/COMMAND_REFERENCE.md.
#
# Public release roadmap PR-301: also runs scripts/verify_docs_consistency.ps1,
# a targeted documentation-consistency check (defaults table vs. actual coded
# defaults, plus regression guards for two real, confirmed doc contradictions
# found during that milestone's review).
#
# Public release roadmap PR-303: also runs scripts/verify_docs_quality.ps1,
# broader documentation quality gates -- command-reference heading structure,
# internal link/anchor integrity, and keyword-argument spelling in ```gauss
# code snippets.
#
# This repo's tests print their own "PASS"/"FAIL" line per check and a
# final "...: ALL N CHECKS PASSED" (or "N CHECKS FAILED") summary line --
# CLAUDE.md documents that tgauss's process exit code is NOT a reliable
# pass/fail signal for this harness, so this runner checks the printed
# summary line (and any GAUSS-level compile/execute error) rather than
# relying on exit code alone.
#
# quaids_pubtable_test.e requires the pubtable package to be installed
# (this machine has it at c:\gauss26\pkgs\pubtable) -- pass -SkipPubtable
# to skip it on a machine without pubtable.
#
# quaids_curvature_test.e requires the optmt package to be installed
# (this machine has it at c:\gauss26\pkgs\optmt, and it is now a real
# package.json dependency -- see "Milestone 10" in CLAUDE.md) -- pass
# -SkipCurvature to skip it on a machine without optmt.
#
# quaids_curvature_bootstrap_test.e (Milestone 15) refits the whole
# quaidsFit()+quaidsCurvatureFit() pipeline B times per model, which adds
# roughly 45-50s to this script's runtime even at the small B values used
# there -- close to doubling this file's own ~30s baseline. Opt-in via
# -SkipBootstrap's ABSENCE is the default (i.e. it runs unless skipped),
# but the automatic push-triggered CI workflow (.github/workflows/tests.yml)
# passes -SkipBootstrap so routine pushes stay fast; run without the flag
# locally, or as part of release verification, to exercise it.
#
# quaids_robust_bootstrap_test.e (Milestone 20) has the same "refits the
# whole pipeline B times" cost, so it shares the same -SkipBootstrap gate
# rather than introducing a second flag.
#
# TVP-AIDS initiative, Stage 2: quaidstvp_kalman_test.e and the two
# tvp_bad_*_shape.e guard cases require `library cmlmt, tsmt, sslib;` --
# sslib is NOT a package.json dependency (deliberately; see
# src/quaidstvp.src's own Stage 2 header) and was found missing from this
# very machine's C:\gauss26\pkgs once already this initiative (reinstalled
# from gauss-state-space's pinned commit via `git archive`, not assumed
# present) -- so its presence here is not guaranteed the way optmt/pubtable
# are. Pass -SkipTVPKalman to skip these on a machine without sslib
# installed; the automatic push-triggered CI workflow
# (.github/workflows/tests.yml) passes it for exactly that reason. These
# three scripts also need the GAUSS26_CFG override documented in
# CLAUDE.md (tsmt package-shadowing on this machine) -- this script points
# GAUSS26_CFG at tests/gauss26_cfg_override for just these invocations, not
# for any of the others.
#
# TVP-AIDS initiative, Stage 3: quaidstvp_mle_test.e and its four
# tvp_mle_*.e guard cases share this same sslib dependency (and the same
# -SkipTVPKalman flag -- not a new flag, since the underlying reason to
# skip is identical: sslib's presence isn't guaranteed) and the same
# GAUSS26_CFG override.
#
# TVP-AIDS initiative, Stage 4: quaidstvp_smooth_test.e and its two
# tvp_smooth_bad_state_*.e guard cases share the same sslib dependency,
# the same -SkipTVPKalman flag, and the same GAUSS26_CFG override.
#
# TVP-AIDS initiative, Stage 5: quaidstvp_elas_test.e and its three
# tvp_elas_bad_*.e guard cases have NO sslib dependency (a state here is a
# plain vector, not any sslib struct -- see src/quaidstvpelas.src's own
# header) and so are NOT gated behind -SkipTVPKalman/GAUSS26_CFG, unlike
# Stages 2-4.
#
# TVP-AIDS initiative, Stage 6: quaidstvpfit_test.e and its two
# tvpfit_bad_*.e guard cases share the sslib dependency/-SkipTVPKalman
# flag/GAUSS26_CFG override of Stages 2-4 (quaidsTVPFit() calls
# _quaidsTVPMLEFit() internally). Also runs scripts/verify_sslib_pin.ps1
# once, right before the sslib-dependent tests, printing a warning (not a
# hard failure -- see that script's own header) if the installed sslib
# copy no longer matches sslib.pin.json's documented commit.

param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [string]$GaussExe = "C:\gauss26\tgauss.exe",
    [switch]$SkipPubtable,
    [switch]$SkipCurvature,
    [switch]$SkipBootstrap,
    [switch]$SkipTVPKalman
)

$testsDir = Join-Path $RepoRoot "tests"
$scriptsDir = Join-Path $RepoRoot "scripts"

& powershell -ExecutionPolicy Bypass -File (Join-Path $testsDir "verify_package_manifest.ps1") -RepoRoot $RepoRoot
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "verify_public_api.ps1") -RepoRoot $RepoRoot
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "verify_docs_consistency.ps1") -RepoRoot $RepoRoot
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "verify_docs_quality.ps1") -RepoRoot $RepoRoot
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

$gaussTests = @(
    "quaids_schema_test.e",
    "quaids_formula_parity_test.e",
    "quaids_synthetic_validation_test.e",
    "quaids_published_validation_test.e",
    "quaids_hypothesis_tests_test.e",
    "quaids_elasticities_test.e",
    "quaids_shares_test.e",
    "quaids_welfare_test.e",
    "quaids_reliability_regression_test.e",
    "quaids_zero_test.e",
    "quaids_robust_test.e",
    "quaids_preflight_test.e",
    "quaids_workflow_test.e",
    "quaids_survey_workflow_test.e",
    "quaids_survey_test.e",
    "quaids_replicate_test.e",
    "quaids_compatibility_test.e",
    "quaids_print_format_test.e",
    "quaidstrend_test.e",
    "quaidstvp_test.e",
    "quaidstvp_elas_test.e"
)

if (-not $SkipPubtable) {
    $gaussTests += "quaids_pubtable_test.e"
}

if (-not $SkipCurvature) {
    $gaussTests += "quaids_curvature_test.e"
}

if (-not $SkipBootstrap) {
    $gaussTests += "quaids_curvature_bootstrap_test.e"
    $gaussTests += "quaids_robust_bootstrap_test.e"
}

$sslibTests = @()
if (-not $SkipTVPKalman) {
    & powershell -ExecutionPolicy Bypass -File (Join-Path $scriptsDir "verify_sslib_pin.ps1") -RepoRoot $RepoRoot
    # Advisory only (verify_sslib_pin.ps1 exits 0 on drift unless -Strict) --
    # a warning here is a heads-up before the sslib-dependent tests below,
    # not a gate on this script's own exit code.

    $sslibTests += "quaidstvp_kalman_test.e"
    $gaussTests += "quaidstvp_kalman_test.e"
    $sslibTests += "quaidstvp_mle_test.e"
    $gaussTests += "quaidstvp_mle_test.e"
    $sslibTests += "quaidstvp_smooth_test.e"
    $gaussTests += "quaidstvp_smooth_test.e"
    $sslibTests += "quaidstvpfit_test.e"
    $gaussTests += "quaidstvpfit_test.e"
}

$gaussCfgOverride = Join-Path $testsDir "gauss26_cfg_override"

function Invoke-GaussBatch {
    param(
        [string]$Exe,
        [string[]]$Arguments,
        [string]$Gauss26Cfg = $null
    )

    # Milestone 15 finding: reading stdout fully (ReadToEnd()) before
    # touching stderr is a classic .NET Process deadlock -- if the child
    # writes enough to BOTH streams to fill their OS pipe buffers before
    # either is drained, the child blocks mid-write while this script
    # blocks reading the other stream, and neither side ever proceeds.
    # quaids_curvature_bootstrap_test.e's QUAIDS block routinely prints
    # dozens to hundreds of "Optmt: function evaluation failed" lines to
    # stderr (a normal, expected side effect of optmt hitting a bad
    # bootstrap resample -- see that test's own header), enough volume to
    # hit exactly this deadlock; a run through this function hung for
    # hours where running the same file directly (tgauss -b -x ...) never
    # did, since a direct console run has no pipe buffer to fill. Fixed by
    # draining both streams asynchronously via events instead of
    # sequential ReadToEnd() calls.
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $Exe
    $psi.Arguments = (($Arguments | ForEach-Object {
        if ($_ -match '[\s"]') {
            '"' + ($_ -replace '"', '\"') + '"'
        } else {
            $_
        }
    }) -join " ")
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.WorkingDirectory = $testsDir
    if ($Gauss26Cfg) {
        # sslib-dependent scripts only (see the TVP-AIDS Stage 2 comment
        # above param()) -- the tsmt package-shadowing fix from CLAUDE.md,
        # scoped to just this child process's environment so it cannot
        # affect any other test invoked by this script.
        $psi.EnvironmentVariables["GAUSS26_CFG"] = $Gauss26Cfg
    }

    $proc = [System.Diagnostics.Process]::new()
    $proc.StartInfo = $psi

    $outputBuilder = [System.Text.StringBuilder]::new()
    $errorBuilder = [System.Text.StringBuilder]::new()

    $outputEvent = Register-ObjectEvent -InputObject $proc -EventName OutputDataReceived -Action {
        if ($null -ne $EventArgs.Data) { [void]$Event.MessageData.AppendLine($EventArgs.Data) }
    } -MessageData $outputBuilder
    $errorEvent = Register-ObjectEvent -InputObject $proc -EventName ErrorDataReceived -Action {
        if ($null -ne $EventArgs.Data) { [void]$Event.MessageData.AppendLine($EventArgs.Data) }
    } -MessageData $errorBuilder

    try {
        [void]$proc.Start()
        $proc.BeginOutputReadLine()
        $proc.BeginErrorReadLine()
        $proc.WaitForExit()
    } finally {
        Unregister-Event -SourceIdentifier $outputEvent.Name -ErrorAction SilentlyContinue
        Unregister-Event -SourceIdentifier $errorEvent.Name -ErrorAction SilentlyContinue
    }

    [pscustomobject]@{
        ExitCode = $proc.ExitCode
        Output = ($outputBuilder.ToString() + $errorBuilder.ToString())
    }
}

$failed = @()

$guardTests = @(
    [pscustomobject]@{
        Script = "guard_error_cases\robust_nonconverged_qout.e"
        Expected = "quaidsRobustFit: qOut must come from a converged quaidsFit() result."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\robust_bad_cluster_length.e"
        Expected = "quaidsRobustFit: clusterId must be scalar 0 or a Tx1 vector matching the sample."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\robust_one_cluster.e"
        Expected = "quaidsRobustFit: cluster-robust SE require at least two clusters."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\curvature_nonconverged_qout.e"
        Expected = "quaidsCurvatureFit: qOut must come from a converged quaidsFit() result."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\curvature_invalid_sym.e"
        Expected = "quaidsCurvatureFit: qOut must have a valid homogeneity+symmetry-constrained estimate (qOut.symValid=1)."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\quaids_bad_b0_shape.e"
        Expected = "quaidsFit: aCtl.b0 must be scalar 0 or an ng x n reduced raw coefficient matrix matching qOut.homogB."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\quaids_scalar_weight.e"
        Expected = "quaidsFit: weight must be scalar 0 or a Tx1 vector matching the number of observations."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\robust_scalar_weight.e"
        Expected = "quaidsRobustFit: weight must be scalar 0 or a Tx1 vector matching the number of observations."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\robust_bootstrap_scalar_weight.e"
        Expected = "quaidsRobustBootstrapFit: weight must be scalar 0 or a Tx1 vector matching the number of observations."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\workflow_scalar_weight.e"
        Expected = "quaidsWorkflowFit: weight must be scalar 0 or a Tx1 vector matching the number of observations."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\replicate_bad_weights_shape.e"
        Expected = "quaidsReplicateWeightFit: replicateWeights must have one row per observation."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\replicate_negative_scale_factor.e"
        Expected = "quaidsReplicateWeightFit: scaleFactor must be positive."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\replicate_scalar_weight.e"
        Expected = "quaidsFit: weight must be scalar 0 or a Tx1 vector matching the number of observations."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\quaids_set_homogeneity_invalid.e"
        Expected = "quaidsSetHomogeneity: homogeneous must be scalar 0 or 1."
    },
    [pscustomobject]@{
        Script = "guard_error_cases\trend_requires_homogenous.e"
        Expected = "quaidsTrendFit: requires aCtl.homogenous == 1"
    },
    [pscustomobject]@{
        Script = "guard_error_cases\tvp_elas_bad_linear.e"
        Expected = "quaidsTVPElasFit: aCtl.linear must be 1"
    },
    [pscustomobject]@{
        Script = "guard_error_cases\tvp_elas_bad_state_length.e"
        Expected = "quaidsTVPElasFit: state must be a k_states x 1 vector"
    },
    [pscustomobject]@{
        Script = "guard_error_cases\tvp_elas_bad_prices_length.e"
        Expected = "quaidsTVPElasFit: prices must be an n x 1 vector of ABSOLUTE log prices"
    }
)

if (-not $SkipTVPKalman) {
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_bad_Q_shape.e"
        Expected = "_quaidsTVPBuildModel: Q must be k_states x k_states"
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_bad_H_shape.e"
        Expected = "_quaidsTVPBuildModel: H must be n1 x n1."
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_mle_bad_q0_shape.e"
        Expected = "_quaidsTVPMLEFit: q0 must be k_states x 1"
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_mle_nonpositive_q0.e"
        Expected = "_quaidsTVPMLEFit: every element of q0 (starting Q diagonal) must be strictly positive."
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_mle_bad_H_shape.e"
        Expected = "_quaidsTVPMLEFit: H must be n1 x n1."
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_mle_bad_y_shape.e"
        Expected = "_quaidsTVPMLEFit: y must be nobs x n1"
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_smooth_bad_state_rows.e"
        Expected = "_quaidsTVPSmoothFit: rslt.filtered_state must have tvpm.k_states rows"
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvp_smooth_bad_state_cols.e"
        Expected = "_quaidsTVPSmoothFit: rslt.filtered_state must have tvpm.nobs columns"
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvpfit_bad_H_shape.e"
        Expected = "_quaidsTVPMLEFit: H must be n1 x n1."
    }
    $guardTests += [pscustomobject]@{
        Script = "guard_error_cases\tvpfit_bad_q0_shape.e"
        Expected = "_quaidsTVPMLEFit: q0 must be k_states x 1"
    }
}

$sslibGuardScripts = @(
    "guard_error_cases\tvp_bad_Q_shape.e",
    "guard_error_cases\tvp_bad_H_shape.e",
    "guard_error_cases\tvp_mle_bad_q0_shape.e",
    "guard_error_cases\tvp_mle_nonpositive_q0.e",
    "guard_error_cases\tvp_mle_bad_H_shape.e",
    "guard_error_cases\tvp_mle_bad_y_shape.e",
    "guard_error_cases\tvp_smooth_bad_state_rows.e",
    "guard_error_cases\tvp_smooth_bad_state_cols.e",
    "guard_error_cases\tvpfit_bad_H_shape.e",
    "guard_error_cases\tvpfit_bad_q0_shape.e"
)

foreach ($guard in $guardTests) {
    Write-Host ""
    Write-Host "==> $($guard.Script) (expected guard error)"
    $cfgForThis = $null
    if ($sslibGuardScripts -contains $guard.Script) { $cfgForThis = $gaussCfgOverride }
    $result = Invoke-GaussBatch -Exe $GaussExe -Arguments @("-b", "-x", $guard.Script) -Gauss26Cfg $cfgForThis
    $output = $result.Output
    $output

    if ($output -match [regex]::Escape($guard.Expected)) {
        Write-Host "PASS  expected guard diagnostic observed"
    } else {
        Write-Host "FAIL  expected guard diagnostic not observed: $($guard.Expected)"
        $failed += $guard.Script
    }
}

foreach ($test in $gaussTests) {
    Write-Host ""
    Write-Host "==> $test"
    $cfgForThis = $null
    if ($sslibTests -contains $test) { $cfgForThis = $gaussCfgOverride }
    $result = Invoke-GaussBatch -Exe $GaussExe -Arguments @("-b", "-x", $test) -Gauss26Cfg $cfgForThis
    $output = $result.Output
    $output

    $hasGaussError = $output -match "Program execute failed|error G[0-9]+|Program compile failed"
    $hasPassSummary = $output -match "ALL \d+ CHECKS PASSED"
    $hasFailSummary = $output -match "\d+ CHECKS FAILED"

    if ($hasGaussError -or $hasFailSummary -or -not $hasPassSummary) {
        $failed += $test
    }
}

if ($failed.Count -gt 0) {
    Write-Host ""
    Write-Host "run_source_tests.ps1: FAIL -- $($failed -join ', ')"
    exit 1
}

Write-Host ""
Write-Host "run_source_tests.ps1: PASS"
