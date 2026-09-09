# Troubleshooting and Interpretation Guide

A symptom-to-action reference for the failure modes real users hit, plus
which result fields to check before trusting coefficients, elasticities,
welfare measures, or curvature results. If your data itself is the issue
(not an error message), see the [Data Preparation
Guide](DATA_PREPARATION_GUIDE.md) instead.

## Symptom-to-Action Table

| Symptom | Likely cause | Recommended action |
| --- | --- | --- |
| `library quaids;` fails, or any `quaidsFit()`/`quaidsControlCreate()` call throws an undefined-symbol/undefined-procedure error immediately | The package was never installed, or the install is stale/incomplete | Reinstall via **Tools > Install Application** (or `scripts/run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact` from this repo), then retry with a fresh `new;` session. See [README's Installation section](../README.md#installation). |
| `error G0507: Undefined structure 'quaidsOut'` (or similar, for any struct name) in a custom script that `#include`s individual `.src` files directly | `quaids.sdf` was not included before the file that references its structs -- GAUSS does not guarantee struct registration order across separately `#include`d files | Add `#include quaids.sdf` (or `#include quaidsutil.src`, which already includes it) as the *first* line of your script, before any other `quaids*.src` include. Every file in `package.json`'s own `src` array already self-includes `quaids.sdf` for exactly this reason. |
| `error G0025: Undefined symbol: 'X'` for a proc you expect to exist (e.g. `quaidsCurvatureFit`, `ptFromQuaids`) | An optional-adapter module (curvature or `pubtable` reporting) was not `#include`d after its `library` statement -- these are deliberately not part of `library quaids;`'s own lazy-load catalog | Add the module's own `#include` line: `#include quaidscurvature.src` (after `library optmt, quaids;`) or `#include pubtable_quaids.src` (after `library pubtable, quaids;` and `#include quaids.sdf`). See [Curvature Imposition](../README.md#curvature-imposition-optional-optmt) or [Reporting](../README.md#reporting-optional-pubtable) in the README. |
| `error G0025` for a core proc (`quaidsFit`, `quaidsElasFit`, ...) that should always be available under `library quaids;` | A stale or corrupted `.lcg` catalog on the installed package -- confirmed to happen in practice (a real, previously-encountered finding, not hypothetical) | Reinstall the package fresh (see the first row above); this regenerates the catalog. |
| A dataframe column selection throws `error G0472: Invalid name` (or similar) inside `quaidsFull()` | A name in `shareVars`/`priceVars`/`totexpVar`/`instrVars`/`extraVars` does not match an actual column in your loaded dataframe -- often a typo, or a transformed column (e.g. a log price) that was never added via `dfaddcol()` | Double-check every column name against `data`'s actual columns, and confirm required transforms were applied first -- see [Data Preparation Guide, Section 2](DATA_PREPARATION_GUIDE.md#2-price-and-expenditure-transformations-and-unit-consistency). |
| A raw GAUSS shape error (`Columns don't match`, `Rows don't match`) inside `quaidsFit()`/`quaidsFull()` | `w`/`intcpt`/`prices`/`totexp`/`instr` have inconsistent row counts, or `w`/`prices` have a different number of goods (columns) | Run [quaidsPreflight](command-reference/quaidsPreflight.md) first -- its `dimensionsOk` check catches this *before* fitting, with a clear diagnostic instead of a raw GAUSS shape error. |
| A specific, named "must be..." error at the top of the call (`aCtl.b0 must be scalar 0 or an ng x n reduced raw coefficient matrix...`, `replicateWeights must have one row per observation`, `scaleFactor must be positive`, `weight must be scalar 0 or a Tx1 vector...`, `clusterId must be scalar 0 or a Tx1 vector...`) | An explicit, deliberate input-validation guard rejected a malformed argument -- these are not bugs, they are fail-fast checks | The error message itself states the exact required shape/value; fix the argument accordingly. Each guard has a dedicated test in `tests/guard_error_cases/` if you want to see a minimal reproduction. |
| `quaidsRobustFit`/`quaidsCurvatureFit` errors with `"... must come from a converged quaidsFit() result"` | You passed a `qOut` where `qOut.converged == 0` into a post-estimation proc that requires a converged base fit | Check `qOut.converged` immediately after `quaidsFit()` and stop before calling any post-estimation proc if it is `0` -- see [What Establishes a Result's Validity](#what-establishes-a-results-validity) below. |
| `quaidsCurvatureFit` errors with `"... qOut.symValid=1"` | Curvature imposition requires a homogeneity+symmetry-constrained starting fit, but the supplied `qOut` was fit unconstrained (`aCtl.homogenous=0`) or its symmetry stage itself failed | Fit with `aCtl = quaidsSetHomogeneity(aCtl, 1);` and confirm `qOut.symValid == 1` before calling `quaidsCurvatureFit`. |
| `weakIV == 1` from `quaidsPreflight`, or a very large first-stage standard error | Your instrument(s) are not strongly correlated with log total expenditure | See [Data Preparation Guide, Section 5](DATA_PREPARATION_GUIDE.md#5-instrument-selection-and-weak-instrument-diagnostics) -- choose a stronger instrument if possible; treat estimates with caution otherwise. |
| `quaidsPreflight` hard-fails on `shareAddOk` or `negativeShareCount` | Real, rounded data rarely sums to floating-point-exact `1`; negative shares usually indicate a genuine data error | See [Data Preparation Guide, Section 1](DATA_PREPARATION_GUIDE.md#1-budget-shares-and-total-expenditure) (row-normalize) and [Section 4](DATA_PREPARATION_GUIDE.md#4-missing-values-invalid-observations-zeros-and-corner-solutions) (negative shares/corner solutions). |
| `qOut.converged == 0` after fitting with `aCtl.maxiter > 1` | The iterated estimator has no global-convergence guarantee -- a real, measured, and documented failure rate (58% for iterated AIDS, 76% for QUAIDS at default settings) | Try `aCtl.relax = .75` (a modest, evidence-backed mitigation), or fall back to the stable LA-AIDS baseline (`aCtl.maxiter = 1`). See [Model & Feature Support Tiers](../README.md#model--feature-support-tiers). |
| The fit reports `qOut.converged == 1` but coefficients look economically implausible (wrong signs, extreme magnitudes) | A real, distinct failure mode from non-convergence: the iteration reached its OWN tolerance at a self-consistent but wrong fixed point ("converged-but-wrong" in the measured sweep, 19-21.5% of fits at default settings) -- `qOut.converged` proves tolerance convergence only, never solution uniqueness | There is no automated multi-start/stability diagnostic in this release (tracked as a deferred roadmap item). Manually cross-check: does `aCtl.relax = .75` or a different `aCtl.b0` starting value converge to a materially different answer? Do the signs/magnitudes match economic priors (e.g. own-price elasticities generally negative)? Compare against the stable LA-AIDS fit on the same data as a sanity baseline. See [Model & Feature Support Tiers](../README.md#model--feature-support-tiers) for the full caveat. |
| `bootOut.nCompleted` is far below `nRequested` (curvature or robust bootstrap), or `rOut.nCompleted` is far below the number of replicate columns supplied (`quaidsReplicateWeightFit`) | Individual replicates/resamples are failing to converge or hitting a low-effective-sample-size guard, and are dropped rather than retried (bootstrap resamples get up to `5x` attempts to reach `B`; replicate-weight columns are fixed and never redrawn) | A low completion count is itself informative -- it suggests the base model/data combination is fragile. Check `nFailed` too; if most replicates fail, treat the resulting SE with real caution rather than as a routine result. See [quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md), [quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md), or [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md). |
| `error G0025` for `ptFromQuaids`/`ptExport`/similar, or a `pubtable`/`optmt`-not-found error | The optional `pubtable`/`optmt` package is not installed, or was not loaded/`#include`d correctly | Install the optional package separately, then follow the exact `library`/`#include` pattern shown in [Reporting](../README.md#reporting-optional-pubtable) or [Curvature Imposition](../README.md#curvature-imposition-optional-optmt) -- both require the combined `library <pkg>, quaids;` form, not two separate `library` statements (a real GAUSS quirk: a second `library` statement can break an earlier one's own cross-file symbol resolution). |

## What Establishes a Result's Validity

Every estimation/post-estimation struct in this library carries explicit
fields that tell you whether its contents are trustworthy -- check these
*before* reading coefficients, elasticities, welfare measures, or
curvature results, not after:

- **[quaidsFit](command-reference/quaidsFit.md)** (`qOut`): check
  `qOut.converged` first (see the non-convergence rows above -- and note
  it proves tolerance convergence only, not correctness). If
  `aCtl.homogenous = 1`, also check `qOut.symValid` before trusting
  `qOut.bS`/`qOut.vS`/`qOut.symStat`/`qOut.symPval` -- a failed symmetry
  restriction (e.g. a non-positive-definite intermediate matrix) leaves
  `symValid = 0` and falls back to the homogeneity-constrained estimate.
  `qOut.bestB`/`qOut.bestV` always reflect the most-constrained estimate
  that actually succeeded.
- **[quaidsPreflight](command-reference/quaidsPreflight.md)** (`pOut`):
  check `pOut.ok` before proceeding to fit at all -- `0` means a hard
  data/design problem (see the table above); a `1` with warnings is safe
  to proceed on but worth reading.
- **[quaidsWorkflowFit](command-reference/quaidsWorkflowFit.md)** /
  **[quaidsWorkflowScenarioFit](command-reference/quaidsWorkflowScenarioFit.md)**
  (`wfOut`): `wfOut.converged` gates everything else. `wfOut.postValid`
  must be `1` before reading `shares`/`incomeElas`/`priceElas`/
  `compPriceElas`; `wfOut.robustValid`/`postRobustValid` must be `1`
  before reading any `*RobustSE`/`*RobustV` field; `wfOut.welfareValid`/
  `welfareRobustValid` (only set by `quaidsWorkflowScenarioFit`) gate
  `cv`/`ev`/`seCV`/`seEV` and their robust counterparts respectively.
- **[quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)**
  (`cOut`): check `cOut.converged` -- the proc itself already refuses to
  run on a `qOut` that isn't converged with `symValid=1` (see the error
  rows above), but the curvature-imposition outer loop has its own,
  separate convergence question.
- **[quaidsZeroFit](command-reference/quaidsZeroFit.md)** (`zOut`): check
  `zOut.converged` first, then `zOut.probitConverged` (one flag per good
  -- a failed first-stage probit for any good undermines that good's
  correction specifically), then `zOut.symValid` if
  `aCtl.homogenous = 1`, mirroring `qOut`'s own three-tier check above.
- **Bootstrap/replicate-weight structs**
  ([quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md),
  [quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md),
  [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)):
  check `nCompleted` relative to `nRequested`/the number of supplied
  replicate columns (see the table above) before trusting the reported
  `seBoot`/`se`.
- **[quaidsElasFit](command-reference/quaidsElasFit.md)**,
  **[quaidsSharesFit](command-reference/quaidsSharesFit.md)**,
  **[quaidsWelfareFit](command-reference/quaidsWelfareFit.md)**: these
  are pure closed-form evaluations at a point, given an already-fitted
  `bestB`/`bestV` -- they carry no separate validity flag of their own;
  their trustworthiness is entirely inherited from whichever upstream
  `qOut`/`wfOut` you evaluated them against.

## Robust Sandwich vs. Bootstrap: Which Should I Use?

Both [quaidsRobustFit](command-reference/quaidsRobustFit.md) (closed-form
sandwich) and
[quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md)
(resampling) answer the same question -- heteroskedasticity- or
cluster-robust standard errors -- but with a real, measured tradeoff:

- **The closed-form sandwich is fast** (one extra pass over already-fitted
  data, no refitting), but its **bread is a simplified formula**
  (`inv(gg)`-based, not `quaidsFit()`'s own nonlinear-price-index-feedback-
  corrected Jacobian). This makes its `se` **dramatically more
  conservative** than `qOut`'s own classical SE -- often 10-100x larger --
  a confirmed, expected consequence of comparing a simple
  equation-by-equation sandwich against the full cross-equation-efficient
  FGLS system, not a defect in the formula.
- **The bootstrap resamples and refits the actual estimator** (`B`
  replications of `quaidsFit()` on resampled data), so its `seBoot`
  typically lands **much closer to `qOut`'s own SE** -- at the cost of
  real runtime (`B` full refits; a single QUAIDS fit is meaningfully
  slower than a single AIDS fit, so choose `B` deliberately, especially
  for QUAIDS).

**Recommendation**: if the sandwich's conservatism doesn't matter for your
use case (e.g. you only need a defensible upper bound, or speed matters
more than precision), use `quaidsRobustFit()` alone. If the gap between
`se` and `qOut`'s classical SE looks implausibly large, or you need SE
close to what the efficient estimator itself would produce, run
`quaidsRobustBootstrapFit()` and prefer `seBoot`. Both are silent,
struct-returning, and can be run in the same script.

## See Also

- [Data Preparation Guide](DATA_PREPARATION_GUIDE.md) -- for data
  problems rather than error messages or code-level failures.
- [Feature support matrix](FEATURE_SUPPORT_MATRIX.md) -- what is and is
  not supported.
- [Model & Feature Support Tiers](../README.md#model--feature-support-tiers) --
  the full convergence-reliability data referenced throughout this guide.
