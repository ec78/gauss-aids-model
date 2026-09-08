# Feature Support Matrix

This matrix summarizes the current public support surface across the three
model choices `quaidsFit()` (and its wrappers `quaids()`/`quaidsFull()`)
selects via `aCtl.linear`/`aCtl.maxiter` -- see the
[usage guide](USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids)
for the exact switch values.

## Support Tier Summary

A new user should be able to identify the stable baseline and the extra
checks experimental estimators require without reading internal roadmap
files -- this section is that summary. See the README's own [Model &
Feature Support Tiers](../README.md#model--feature-support-tiers) for the
same table in a shorter form next to the quick start.

| Component | Tier | Why |
| --- | --- | --- |
| LA-AIDS (`aCtl.maxiter = 1`) | **Stable** | One-step Stone price index -- no iteration, so no convergence-failure mode of this kind exists. Tradeoff: linear-approximation bias (see the synthetic-validation tolerance of 1.20 vs. 0.10 for the iterated models). |
| Iterated AIDS (`aCtl.linear=1`, `aCtl.maxiter>1`) | Experimental | A committed 200-seed sweep measured 58% combined failure at default settings (39% never converges, 19% converges to a self-consistent but wrong answer). Always check `qOut.converged`; `aCtl.relax=.75` measurably helps (see Notes below). |
| QUAIDS (`aCtl.linear=0`, `aCtl.maxiter>1`) -- `quaidsControlCreate()`'s actual shipped default | Experimental (highest measured risk) | Same sweep measured 76% combined failure (54.5% never converges, 21.5% wrong). This is the library's own default combination -- do not rely on the coded defaults without explicitly checking `qOut.converged`. |
| Zero-share correction (`quaidsZeroFit`) | Experimental | Inherits the base model's convergence risk above (it runs the same translog-price-index iteration), plus a simplified (non-sandwich) SE formula, approximate (not exact) adding-up in the corrected coefficients, and a known, non-trappable GAUSS `glm()` crash mode on some inputs -- see Notes below. |
| Curvature imposition (`quaidsCurvatureFit`, requires `optmt`) | Experimental | Delta-method standard errors are known-unreliable whenever the estimated Cholesky factor sits at the boundary of the negative-semidefinite cone, a common outcome in this library's own test fixtures -- prefer `quaidsCurvatureBootstrapFit`/`quaidsCurvatureBootstrapCI`. QUAIDS curvature additionally requires damping (`aCtl.relax=.25`-ish) to converge at all -- see Notes below. |
| Bootstrap / replicate-weight procedures (`quaidsCurvatureBootstrapFit`, `quaidsRobustBootstrapFit`, `quaidsReplicateWeightFit`) | Inherits the base model's tier | Each resample or caller-supplied replicate is an independent refit of the base model above -- a replicate that itself fails to converge (or falls below a minimum effective sample size) is dropped, not retried. Check the returned struct's `nCompleted`/`nFailed` (or `nRequested`/`nAttempts`) before trusting the reported SE. |

**What `qOut.converged == 1` proves, and what it does not**: it means the
iteration's relative parameter change (`err = max(abs((b-b0)/b0))`) fell
below `aCtl.err` before `aCtl.maxiter` was reached -- nothing more. It does
**not** prove the fixed point found is unique, or that it is the fixed
point you intended: the "converged-but-wrong" bucket in the sweep above is
defined as exactly this -- a fit that satisfies this tolerance test while
still landing far (>10x the normal structural tolerance) from the true
answer on a known synthetic DGP. Naive successive-substitution on this
nonlinear FGLS system can have multiple fixed points for a bad price draw;
no amount of tolerance-tightening or damping changes which basin of
attraction a given dataset's iteration falls into. For `aCtl.maxiter = 1`
(LA-AIDS), `qOut.converged` is unconditionally `1` (there is no iteration
to fail) and this caveat does not apply.

| Feature | LA-AIDS | Iterated AIDS | QUAIDS |
| --- | --- | --- | --- |
| Estimator | `quaidsFit`/`quaids`/`quaidsFull` | `quaidsFit`/`quaids`/`quaidsFull` | `quaidsFit`/`quaids`/`quaidsFull` |
| Price index | Stone (one-step) | Nonlinear translog, iterated | Nonlinear translog, iterated |
| Quadratic log-expenditure term | No | No | Yes |
| IV treatment of total expenditure | Always (control-function) | Always (control-function) | Always (control-function) |
| Overidentification test | Yes, if `ninst > nu` | Yes, if `ninst > nu` | Yes, if `ninst > nu` |
| Homogeneity imposition | Yes (`aCtl.homogenous=1`) | Yes (`aCtl.homogenous=1`) | Yes (`aCtl.homogenous=1`) |
| Symmetry imposition (given homogeneity) | Yes | Yes | Yes |
| Symmetry-given-homogeneity test | Yes (built into `quaidsFit`) | Yes (built into `quaidsFit`) | Yes (built into `quaidsFit`) |
| Standalone homogeneity test | Yes (`quaidsHomogeneityTest`, needs `aCtl.homogenous=0` fit) | Yes | Yes |
| Standalone joint homogeneity+symmetry test | Yes (`quaidsJointTest`) | Yes | Yes |
| Quadratic-term (AIDS-vs-QUAIDS) specification test | Not applicable (no quadratic term to test) | Not applicable | Yes (`quaidsQuadraticTest`) |
| Elasticities at arbitrary points | Yes (`quaidsElasFit`) | Yes | Yes |
| Delta-method elasticity standard errors | Yes | Yes | Yes |
| Predicted budget shares at arbitrary points | Yes (`quaidsSharesFit`, no extra dependency) | Yes | Yes |
| Exact algebraic identity validation (Engel/Cournot/homogeneity) | Yes | Yes | Yes |
| Slutzky negativity diagnostic | Yes (`quaidsSlutzky`) | Yes | Yes |
| Preflight data/design diagnostics | Yes (`quaidsPreflight`) | Yes | Yes |
| Welfare measures (exact CV/EV) | Yes (`quaidsWelfareFit`, no extra dependency) | Yes | Yes |
| Curvature imposition | Yes (`quaidsCurvatureFit`, sample mean, requires `optmt` -- see Notes) | Yes (same) | Yes, requires `aCtl.relax` -- see Notes |
| Curvature bootstrap standard errors | Yes (`quaidsCurvatureBootstrapFit` -- see Notes) | Yes (same) | Yes (same) |
| Curvature bootstrap percentile CIs | Yes (`quaidsCurvatureBootstrapCI`) | Yes (same) | Yes (same) |
| Dataframe/column-name entry point | Yes (`quaidsFull`) | Yes | Yes |
| Formula-string (`"y ~ x"`) API | Not applicable (multi-equation system) | Not applicable | Not applicable |
| `pubtable` export (LaTeX/Markdown/CSV/...) | Yes (`src/pubtable_quaids.src`, optional) | Yes | Yes |
| Synthetic deterministic validation | Yes (`tests/quaids_synthetic_validation_test.e`) | Yes | Yes |
| Published-data cross-validation vs. R | Yes (`Blanciforti86` vs. 3SLS, `tests/quaids_published_validation_test.e`) | Yes (`Blanciforti86` vs. `method="IL"`, wider tolerance -- see Notes) | No independent reference implementation exists (see Notes) |
| Iteration convergence guarantee | Not applicable (one-step) | No -- a 200-seed sweep measured a 58% failure rate (never-converges or converges to a wrong answer) at default settings; check `qOut.converged`. `aCtl.relax=.75` measurably reduces this -- see Notes | No -- same caveat, 76% failure rate measured |
| Zero budget share correction | Yes (`quaidsZeroFit`, unconstrained or homogeneity/symmetry -- see Notes) | Yes (same) | Yes (same) |
| Robust / cluster-robust standard errors | Yes (`quaidsRobustFit`, simplified bread -- see Notes) | Yes (same) | Yes (same) |
| Robust / cluster bootstrap | Yes (`quaidsRobustBootstrapFit`) | Yes (same) | Yes (same) |
| Robust covariance propagation to shares/elasticities/welfare | Yes (`quaidsRobustCovariance`/`quaidsRobustBootstrapCovariance`) | Yes (same) | Yes (same) |
| Applied workflow bundle | Yes (`quaidsWorkflowFit` with compact preflight summary; `quaidsWorkflowScenarioFit` for CV/EV scenarios) | Yes (same) | Yes (same) |
| Sampling-weighted point estimate + weighted/clustered SE | Yes (`quaidsFit`'s optional `weight` argument -- see Notes) | Yes (same) | Yes (same) |
| Sampling-weighted workflow (estimator + evaluation point) | Yes (`quaidsSurveyWorkflowFit` -- see Notes) | Yes (same) | Yes (same) |
| Replicate-weight (jackknife/BRR-style) standard errors | Yes (`quaidsReplicateWeightFit`, caller-supplied design only -- see Notes) | Yes (same) | Yes (same) |
| Installed-package (`library quaids;`) support | Yes | Yes | Yes |

## Notes

- **Iteration convergence guarantee**: a committed 200-seed x 2-model
  sweep (`tests/quaids_convergence_sweep.e`, default settings
  `aCtl.relax=1`, `aCtl.err=.0001`, `aCtl.maxiter=100`) measures: iterated
  AIDS never-converges 39% of the time and converges to a wrong answer (a
  self-consistent but incorrect fixed point, distinct from simply running
  out of iterations) another 19% (58% combined failure); QUAIDS
  never-converges 54.5% and converges wrong another 21.5% (76%
  combined). `aCtl.relax` under-relaxes the fixed-point update; `relax=.75`
  measurably improves the correct-convergence rate (iterated AIDS to 43%,
  QUAIDS to 26.5%), but more aggressive damping (`.5`, `.3`) does not help
  further and often makes things worse. This is a modest, evidence-backed
  mitigation, not a solved problem -- see
  [Usage guide](USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids)
  and `GOLD_STANDARD_TODO.md` for the full grid.
- "Always (control-function)" means `instr` is a required argument to
  every estimator entry point -- there is no exogenous-total-expenditure
  estimation mode in this library.
- Sampling-weighted estimation (`quaidsFit`'s optional `weight` argument)
  is a genuine weighted point estimate: every cross-product in the
  starting value, iteration loop, Jacobian-corrected variance, and
  overidentification test is pre-scaled by `sqrt(weight)`, the standard
  survey-WLS trick, an exact no-op when `weight` is uniform. `weight` is
  renormalized internally to sum to `nobs`. `quaidsRobustFit`/
  `quaidsRobustBootstrapFit` accept the same optional `weight` for a
  matching Horvitz-Thompson pweight-robust sandwich -- **a different
  scaling convention from the point estimate's own `sqrt(weight)`**: the
  bread keeps `sqrt(weight)`, but the per-observation score contribution
  uses plain `weight`. `quaidsPreflight` validates the same weight (a
  required positional argument there, mirroring `clusterId`'s convention)
  and `quaidsWorkflowFit` threads an optional `weight` through its own
  sub-calls. `quaidsSurveyWorkflowFit`'s own `weight` argument both fits
  the weighted estimator (via `quaidsWorkflowFit`'s argument) and
  recomputes the workflow's representative evaluation point as the
  weighted mean of intercept shifters, prices, and total expenditure.
  Formal strata as a concept distinct from clustering, and
  finite-population correction, remain future survey/microdata work.
- Replicate-weight (jackknife/BRR-style) standard errors
  (`quaidsReplicateWeightFit`) implement the shared linear form underlying
  every linearized replication variance estimator,
  `V = sum_r c_r * vec(b_r - b_full) * vec(b_r - b_full)'`, from a
  caller-supplied `TxR` matrix of replicate weight columns and a
  scale factor (scalar or `Rx1`) -- **no specific design (JK1, JKn, BRR,
  Fay's BRR) is implemented or auto-detected**; both inputs are always
  required. Unlike the bootstrap procs it otherwise resembles, there is no
  resampling loop, no `seed`, and no retry -- a failed replicate (fixed,
  caller-supplied, not random) is simply dropped from the sum, since the
  formal jackknife/BRR literature does not define a missing-replicate
  adjustment this library implements. `rOut.b`/`rOut.v` are already in
  `quaidsFit()`'s own full `bestB` basis, so -- unlike `quaidsRobustFit()`
  -- no expansion helper is needed before
  `quaidsSharesFit()`/`quaidsElasFit()`/`quaidsWelfareFit()`. A replicate
  weight concentrated on too few effectively-weighted rows is skipped
  before it can drive `quaidsFit()`'s iteration into a rank-deficient,
  crashing state. See
  [Methodology Notes](METHODOLOGY_NOTES.md#replicate-weight-jackknifebrr-variance-estimation).
- Welfare measures (`quaidsWelfareFit`) are exact, not approximated, for
  all three model choices -- unlike curvature imposition, computing CV/EV
  needs no new estimation, only a closed-form evaluation of the
  already-fitted expenditure function at two points, so QUAIDS needs no
  separate scoping. See [Methodology Notes](METHODOLOGY_NOTES.md#welfare-measures).
- The published-data cross-validation (`Blanciforti86` vs. R's
  `micEconAids`) covers both LA-AIDS (`aCtl.linear=1, aCtl.maxiter=1`, vs.
  `aidsEst(..., instNames=...)`, 3SLS -- max abs difference ~0.021) and
  iterated AIDS (`aCtl.linear=1, aCtl.maxiter>1`, vs.
  `aidsEst(method="IL", ...)`, the Iterated Linear Least Squares Estimator
  -- max abs difference ~0.11, tolerance `0.15`). The iterated-AIDS
  comparison has a wider gap for a real, understood reason, not
  approximation slop: `micEconAids`'s `method="IL"` does not support
  instrumental variables (combining it with `instNames` segfaults R's
  `aidsEst` rather than erroring cleanly), so that reference is
  SUR-estimated, while GAUSS's iterated fit always instruments log total
  expenditure. The comparison therefore spans both a different estimation
  algorithm *and* an IV-vs-no-IV difference. **QUAIDS has no independent
  reference implementation available**: `micEconAids` does not implement a
  quadratic log-expenditure term at all, and no other comparably-established
  QUAIDS implementation was found (see `GOLD_STANDARD_TODO.md` on the
  Python from-scratch replica, kept as supplementary evidence only).
  QUAIDS's validation is therefore the known-true synthetic-DGP recovery
  in `tests/quaids_synthetic_validation_test.e` -- a real, non-circular
  check (independently-generated data with known-true parameters, not
  just re-running the estimator on its own output), but a different,
  weaker tier of evidence than cross-implementation agreement on real
  published data. Documented here rather than silently claimed as
  equivalent.
- Curvature imposition (Diewert-Wales Cholesky reparametrization,
  `quaidsCurvatureFit`) is available for LA-AIDS/AIDS (`aCtl.linear=1`)
  and QUAIDS (`aCtl.linear=0`), imposed locally at the sample mean,
  requiring the `optmt` package -- an opt-in adapter
  (`src/quaidscurvature.src`, not in `package.json`'s `src` array, same
  treatment as the optional `pubtable` reporting adapter), not a
  `library quaids;` dependency; see `docs/public-api.json`'s
  `optional_modules` entry. QUAIDS's curvature outer loop is measurably
  less stable than AIDS's own: `aCtl.relax` is effectively required, not
  optional -- undamped runs on the validation fixture diverge to NaN
  within a handful of iterations. Standard errors from
  `quaidsCurvatureFit` are a simplified delta-method approximation, known
  to be unreliable when the estimated Cholesky factor has boundary
  (near-zero) entries (a standard complication of Cholesky-based
  negative-semidefinite-cone estimation) -- point estimates and the exact
  curvature property at the reference point are unaffected.
  `quaidsCurvatureBootstrapFit` closes this gap with a nonparametric
  i.i.d. row bootstrap (resample, refit the whole pipeline, take the
  empirical SE), reported alongside rather than replacing the delta-method
  SE; it has no default replication count, since a single AIDS curvature
  fit and a single QUAIDS curvature fit differ in runtime by roughly an
  order of magnitude, making a one-size-fits-all default misleading.
  `quaidsCurvatureBootstrapCI` adds percentile confidence intervals
  directly from the bootstrap's raw draws, no new resampling needed. There
  is no independent published/cross-implementation validation for the
  *imposed* estimator on either model: even the R `micEconAids` reference
  implementation used elsewhere in this library only diagnoses curvature,
  never imposes it. For QUAIDS specifically, `tests/quaids_curvature_test.e`
  validates convergence/exact negative-semidefiniteness/non-vacuousness/
  shape rather than "recovers a known true curvature-consistent gamma" the
  way the AIDS block does -- a deliberately weaker (but still real) tier
  of evidence, documented as such rather than silently equated with AIDS's.
  See [Methodology Notes](METHODOLOGY_NOTES.md#curvature-imposition-diewert-wales)
  and `GOLD_STANDARD_TODO.md` for the full history.

- Zero budget share correction (Shonkwiler-Yen, `quaidsZeroFit`)
  addresses real survey/microdata's corner solutions (zero-expenditure
  goods), which `quaidsFit()` does not model. A per-good first-stage
  probit's fitted probability `F_i` is divided into the second-stage share
  equation (`w_i/F_i = ...`), which avoids breaking the shared-design-
  matrix Kronecker-product identity every stage of `quaidsFit()` relies on
  -- a literal textbook implementation (rescaling every regressor by
  `F_i`) would not. `aCtl.homogenous=1` also imposes homogeneity/symmetry
  on the corrected model, in the same minimum-distance projection as the
  method's own diagonal-delta restriction. Standard errors are a
  simplified `S .*. inv(gg)` formula that does not correct for the
  nonlinear translog-price-index feedback or first-stage probit/IV
  generated-regressor uncertainty. Adding-up does not hold exactly for the
  corrected coefficients -- a real, known property of the method itself,
  not a bug. Validated on a synthetic fixture with a known latent
  (uncensored) DGP and a genuine, non-degenerate zero-share censoring
  pattern (`tests/quaids_zero_test.e`): the corrected fit recovers the
  true latent parameters measurably better than naively fitting
  `quaidsFit()` on the same censored data. GAUSS's built-in `glm()` (used
  for the first-stage probits, no new package dependency) can hard-crash
  on some degenerate inputs, a known non-trappable failure mode. See
  [Methodology Notes](METHODOLOGY_NOTES.md#zero-budget-share-correction-shonkwiler-yen).

- Robust / cluster-robust standard errors (`quaidsRobustFit`) generalize
  the pooled, homoskedastic `S.*.inv(gg)` sandwich every other covariance
  in this library uses to a per-observation (heteroskedasticity-robust) or
  per-cluster (cluster-robust, with a CR1 small-sample correction) score
  aggregation -- genuinely new math, since neither GAUSS's base runtime
  nor the `tsmt` package's single-equation `robustSE`/`clusterSE`
  generalize to this library's stacked multi-equation system. Robust is
  the literal `G=nobs` special case of cluster-robust, unified through one
  `clusterId` argument. Uses a **simplified bread** (`inv(gg)`-based, not
  `quaidsFit()`'s own nonlinear-feedback-corrected Jacobian), which makes
  its `se` dramatically more conservative than `qOut`'s own classical SE
  -- an expected consequence of comparing a simple sandwich against the
  full cross-equation-efficient FGLS system, not a bug.
  `quaidsRobustBootstrapFit` offers a cluster-aware nonparametric
  bootstrap alternative that resamples whole clusters and refits
  `quaidsFit()` itself, typically landing much closer to `qOut`'s own SE
  than the closed-form sandwich does. The reduced robust coefficient table
  covers only the `n1` independently-estimated equations;
  `quaidsRobustCovariance`/`quaidsRobustBootstrapCovariance` expand the
  robust or bootstrap covariance into `qOut.bestB`'s full basis for
  elasticities, shares, and welfare. See
  [Methodology Notes](METHODOLOGY_NOTES.md#robust-and-cluster-robust-standard-errors).

Related documentation:

- [Usage guide](USAGE_GUIDE.md)
- [Methodology notes](METHODOLOGY_NOTES.md)
- [Command reference](COMMAND_REFERENCE.md)
