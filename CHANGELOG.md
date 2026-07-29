# Changelog

All notable changes to this project are documented here. This project is
pre-alpha and does not yet follow strict semantic versioning guarantees
(see `GOLD_STANDARD_TODO.md` for the release roadmap); version numbers
below match `package.json` at the time each milestone landed.

## Unreleased

### Fixed
- `quaidsFit()` now correctly honors documented `aCtl.b0` starting values by
  initializing the shared Stone-index design state regardless of whether
  the start matrix is caller-supplied or internally estimated. Supplied
  `b0` matrices are now shape-validated with a clear diagnostic.
- Corrected first-stage IV and homogeneity-stage diagnostic degrees of
  freedom, which previously subtracted an extra constant even though the
  design matrices already included one.
- `quaidsRobustFit()` and `quaidsRobustBootstrapFit()` now fail fast with
  explicit diagnostics for non-converged base fits, malformed cluster
  vectors, one-cluster inputs, insufficient CR1 degrees of freedom, and
  singular robust-covariance bread matrices.
- `quaidsCurvatureFit()` now rejects non-converged or invalid
  symmetry-constrained fits before imposing curvature.
- `quaidsSharesFit()`, `quaidsElasFit()`, and `quaidsWelfareFit()` now build
  finite-difference perturbation directions elementwise, avoiding scalar
  vector-condition behavior when coefficient vectors mix zero and nonzero
  entries.
- `quaidsZeroFit()` now accepts and validates compatible `aCtl.b0` starting
  matrices instead of silently ignoring the shared control field.
- README, CLAUDE orientation notes, and command-reference pages now reflect
  package version `0.15.0`, installability, `optmt` dependency scope, and
  the refined `b0` contracts.

### Added
- `quaidsRobustCovariance()` expands `quaidsRobustFit()`'s reduced
  robust/cluster covariance into `qOut.bestB`'s full coefficient basis for
  shares, elasticities, and welfare delta-method standard errors.
- `quaidsRobustBootstrapCovariance()` computes the empirical covariance of
  `quaidsRobustBootstrapFit()` draws and expands it into the same full
  post-estimation basis.
- `quaidsWorkflowFit()` now propagates robust/cluster-robust covariance
  into mean-point share and elasticity standard errors, and
  `quaidsWorkflowScenarioFit()` now reports robust CV/EV welfare standard
  errors.
- `quaidsWorkflowFit()`/`quaidsWorkflowOut` as the first Milestone 21
  applied workflow bundle, composing `quaidsFit()`, `quaidsSharesFit()`,
  `quaidsElasFit()`, and `quaidsRobustFit()` into one silent,
  struct-returning call.
- `quaidsWorkflowFit()` now includes restriction/model-choice summary
  fields for symmetry, overidentification, and unconstrained-QUAIDS
  quadratic-term tests.
- `quaidsWorkflowScenarioFit()` fills the workflow struct's exact CV/EV
  welfare fields for one explicit price-change scenario while preserving
  `quaidsWorkflowFit()`'s shorter base signature.
- `ptTablesFromQuaidsWorkflow()` in the optional `pubtable_quaids.src`
  adapter exports workflow shares, elasticities, and welfare-scenario
  outputs without manual reshaping.
- `examples/workflow_example.e` demonstrates the installed-package workflow
  path.
- `tests/quaids_workflow_test.e` parity coverage proving workflow output
  matches explicit calls to the underlying public APIs, including welfare
  scenario parity.
- `tests/package_public_api.e` now exercises the workflow procs against
  the installed package release gate.
- Source-test harness coverage for expected guard failures, including
  invalid robust/cluster inputs, invalid curvature inputs, and malformed
  `aCtl.b0` matrices.

### Documentation
- README, usage guide, feature matrix, command reference, `CLAUDE.md`, and
  `GOLD_STANDARD_TODO.md` now describe the post-20 roadmap and Milestone
  21 Gold Standard to-dos.

## 0.15.0 - 2026-07-28

### Added
- `quaidsRobustFit()`/`printQuaidsRobust()` (`src/quaidsrobust.src`, new
  file): heteroskedasticity-robust and cluster-robust standard errors for
  an already-fitted `quaidsFit()` result, generalizing the pooled,
  homoskedastic `S.*.inv(gg)` sandwich every other covariance in this
  library uses to a per-observation or per-cluster score aggregation.
  Genuinely new math -- neither GAUSS's base runtime nor the `tsmt`
  package's single-equation `robustSE`/`clusterSE` generalize to this
  library's stacked multi-equation system. Robust and cluster-robust are
  unified through one `clusterId` argument (`0` = robust, a `Tx1` group
  vector = cluster-robust with a CR1 correction), confirmed to be
  literally the same formula (an exact-identity regression test).
- `quaidsRobustBootstrapFit()`/`printQuaidsRobustBootstrap()` (same
  file): a cluster-aware nonparametric bootstrap alternative, resampling
  whole clusters (or plain rows) and refitting `quaidsFit()` itself.
- New `quaidsRobustOut`/`quaidsRobustBootOut` structs (`src/quaids.sdf`).
- `tests/quaidsfixtures.src`: `_quaidsClusterSyntheticDGP()`, a 5-good
  fixture with a genuine within-cluster-correlated noise component,
  needed to make the "cluster-robust se exceeds naive se" check
  non-vacuous.
- `tests/quaids_robust_test.e` (17 checks) and
  `tests/quaids_robust_bootstrap_test.e` (13 checks, added to the
  existing `-SkipBootstrap`-gated group).
- `tests/package_public_api.e` extended to exercise the four new procs.
- New docs: `docs/command-reference/quaidsRobustFit.md`,
  `printQuaidsRobust.md`, `quaidsRobustBootstrapFit.md`,
  `printQuaidsRobustBootstrap.md`; new sections in
  `docs/USAGE_GUIDE.md` and `docs/METHODOLOGY_NOTES.md`; new rows in
  `docs/FEATURE_SUPPORT_MATRIX.md`.

### Found (documented, not a defect)
- `quaidsRobustFit()`'s simplified bread (`inv(gg)`-based, not
  `quaidsFit()`'s own nonlinear-feedback-corrected Jacobian) makes its
  `se` dramatically more conservative (often 10-100x larger) than
  `qOut.homogSE`/`symcSE` -- confirmed, via an independent hand-
  derivation using the same regressors/residuals, to be an expected
  consequence of comparing a simple equation-by-equation sandwich against
  the full cross-equation-efficient FGLS system, not a bug.
  `quaidsRobustBootstrapFit()`'s bootstrap SE does not share this gap.

### Fixed (found and fixed before shipping)
- An early version of the bootstrap tracked the point estimate in
  `qOut.bestB`'s full (adding-up-recovered) shape while
  `quaidsRobustFit()`'s own SE is in a reduced form with one fewer row --
  a genuine shape mismatch caught by running
  `printQuaidsRobustBootstrap()` against real data
  (`error G0058: Index out of range`), not by re-reading the code. Fixed
  by sharing one reduction helper (`_quaidsRobustReduceB()`) between both
  procs.
- `tests/quaidsfixtures.src`'s new `_quaidsClusterSyntheticDGP()` drew its
  `clusterId` via unseeded `rndu()`, unlike every other draw in this
  codebase's fixtures (`rndns(rows,cols,seed)`'s explicit-seed form) --
  `rndu()` has no such seeded form, so this one draw depended on GAUSS's
  ambient global random state, making
  `tests/quaids_robust_bootstrap_test.e`'s core check genuinely flaky
  (confirmed by running it three times in a row: pass, fail, fail).
  Fixed with one `rndseed seed;` call at the top of the fixture.
- `tests/package_public_api.e`'s `quaidsRobustBootstrapFit()` block reused
  a known non-converging fixture (fine for `quaidsRobustFit()`, which
  needs no convergence, but not for `quaidsRobustBootstrapFit()`, which
  requires its own base fit to converge) -- caught by the release-
  verification gate against the real installed package. Fixed by reusing
  the file's own already-converging seed=500 AIDS fixture instead.

This milestone completes the five-item full-demand-system-workflow
outline (predicted shares, a quadratic-term specification test, bootstrap
percentile CIs, zero-share correction, and now robust/cluster standard
errors) worked through in order since Milestone 16.

## 0.14.0 - 2026-07-27

### Added
- `quaidsZeroFit()`/`printQuaidsZero()` (`src/quaidszerocorrect.src`, new
  file): AIDS/QUAIDS estimation corrected for zero budget shares (corner
  solutions) via the Shonkwiler & Yen (1999) two-step procedure. A
  per-good first-stage probit's fitted probability is divided into the
  second-stage share equation, turning what would otherwise be a
  regressor-rescaling problem (incompatible with the shared-design-matrix
  Kronecker-product identity every stage of `quaidsFit()` relies on) into
  a dependent-variable transform plus one new shared regressor column per
  equation. A new minimum-distance restriction (mirroring `quaidsFit()`'s
  own symmetry-restriction machinery) forces the resulting hazard-
  coefficient block to be diagonal, as the method requires. Unconstrained
  only in this pass (errors if `aCtl.homogenous=1`); standard errors are a
  simplified `S .*. inv(gg)` formula; adding-up does not hold exactly for
  the corrected coefficients (a real property of the method itself, not a
  bug) -- all documented explicitly, not silently absent. Uses GAUSS's
  built-in `glm()` for the probit stage, no new package dependency.
- New `quaidsZeroOut` struct (`src/quaids.sdf`).
- `tests/quaidsfixtures.src`: `_quaidsZeroSyntheticDGP()`, a 5-good QUAIDS
  fixture with a known latent (uncensored) DGP, censored the economically
  correct way (`w_i = max(0,latent_i) / sum_j max(0,latent_j)`) to
  produce a genuine, non-degenerate zero-share pattern.
- `tests/quaids_zero_test.e` (17 checks): diagonal-delta restriction
  correctness, shape/finiteness, and the core validation -- corrected
  recovery of the true latent-DGP parameters beats naive `quaidsFit()` on
  the same censored data, on both a max- and mean-absolute-difference
  basis.
- `tests/package_public_api.e` extended to exercise `quaidsZeroFit()`/
  `printQuaidsZero()` against the installed package.
- New docs: `docs/command-reference/quaidsZeroFit.md`,
  `printQuaidsZero.md`; new sections in `docs/USAGE_GUIDE.md` and
  `docs/METHODOLOGY_NOTES.md`; new row in
  `docs/FEATURE_SUPPORT_MATRIX.md`.

### Known limitation
- GAUSS's built-in `glm()` (used for the first-stage probits) can
  hard-crash on some degenerate inputs (`Intel MKL ERROR ... DGELS`), a
  failure mode not intercepted by this codebase's usual `trap`/`scalmiss`
  guard -- the same class of non-trappable failure already documented for
  `eighv()` inside `quaidsCurvatureBootstrapFit()`. Not hardened against
  in this pass.

## 0.13.0 - 2026-07-27

### Added
- `quaidsCurvatureBootstrapCI()` (`src/quaidscurvature.src`): percentile
  bootstrap confidence intervals for a `quaidsCurvatureFit()` coefficient
  vector, computed directly from an already-computed
  `quaidsCurvatureBootstrapFit()` result's raw draws (`bootOut.bBoot`) --
  no new resampling or refitting. Returns a bare `(ciLower, ciUpper)`
  tuple, matching this codebase's existing bare-tuple hypothesis-test
  convention. `alpha` is required, no default.
- `tests/quaids_curvature_bootstrap_test.e` gained 11 checks (26 -> 37):
  shape/ordering/containment checks for the new CI proc, plus a direct
  `quantile()` cross-check at a specific flattened index verifying the
  reshape/index mapping.

### Fixed
- **A real, silent bug present since Milestone 10/15**:
  `quaidsCurvatureFit()`'s `se` and `quaidsCurvatureBootstrapFit()`'s
  `seBoot` had their individual cells scrambled relative to `b`/`bootOut.b`
  -- GAUSS's `reshape()` fills row-major, not column-major like `vec()`,
  and both procs' reshape calls had dropped the transpose needed to
  correctly invert `vec()`'s ordering. Invisible to existing shape/sign/
  finiteness checks (all permutation-invariant). Confirmed directly
  against real fitted data: `se`'s cells did not match
  `sqrt(diag(v))` at the correct positions until fixed; `v` (the full
  covariance) was never affected. Found while building
  `quaidsCurvatureBootstrapCI()`'s own ground-truth cross-check. Both
  `tests/quaids_curvature_test.e` and `tests/
  quaids_curvature_bootstrap_test.e` gained explicit cell-position
  regression guards (2 and part of the 11 checks above, respectively) so
  this bug class cannot silently return.

## 0.12.0 - 2026-07-27

### Added
- `quaidsQuadraticTest()` (`src/quaidstests.src`): a standalone Wald test
  of whether the QUAIDS quadratic log-expenditure term (`lambda`) is
  needed at all, or whether plain AIDS would fit equally well. Mirrors
  `quaidsHomogeneityTest()`/`quaidsJointTest()` exactly (same file, same
  `(stat, pval, df)` return, same unconstrained-fit guard), plus an
  additional guard requiring `qOut.linear == 0` -- an AIDS fit never
  estimates `lambda` at all, so there is nothing to test. `df = n-1`.
- `tests/quaids_hypothesis_tests_test.e` gained 3 checks (19 -> 22): size
  and power for the new test, reusing the file's existing `quadratic=1`
  fixture for power and a fresh `quadratic=0` fixture for size.

## 0.11.0 - 2026-07-26

### Added
- `quaidsSharesFit()` and `printQuaidsShares()` (`src/quaidsshares.src`),
  plus a new `quaidsSharesOut` struct (`src/quaids.sdf`): the model-
  implied predicted budget share vector, and its full delta-method
  covariance, at an arbitrary evaluation point -- useful for out-of-
  sample prediction and policy simulation without hand-deriving the
  share equation. Reports the full `n x n` covariance of `w` (not just
  marginal SEs), built via a numerical Jacobian propagated as
  `jacW*v*jacW'`. A deliberately independent third implementation of the
  share formula (already duplicated in `quaidsElas_()` and a test-only
  helper) rather than a refactor of shipped code -- see
  `GOLD_STANDARD_TODO.md`'s Milestone 16 section.
- `tests/quaids_shares_test.e` (21 checks): point-estimate correctness
  against an independent hand-evaluated formula (AIDS and QUAIDS
  fixtures), exact adding-up, SE/covariance shape and consistency,
  non-vacuousness.

### Changed
- `tests/quaids_elasticities_test.e`'s private `modelShareAt()` helper
  removed; its Engel/Cournot/homogeneity identity checks now use the
  public `quaidsSharesFit()` directly, removing a formula duplicated
  three times across the codebase down to two (the estimator's own
  internal share equation, and `quaidsSharesFit()`).

## 0.10.0 - 2026-07-26

### Added
- `quaidsCurvatureBootstrapFit()` and `printQuaidsCurvatureBootstrap()`
  (`src/quaidscurvature.src`), plus a new `quaidsCurvBootOut` struct
  (`src/quaids.sdf`): a nonparametric i.i.d. row (pairs) bootstrap
  standard error for `quaidsCurvatureFit()`'s coefficient vector,
  reported alongside (not replacing) the existing delta-method SE. Closes
  a gap documented since Milestone 10: the delta-method SE is known to be
  unreliable whenever the estimated Cholesky factor sits at the boundary
  of the negative-semidefinite cone. Resamples rows of `(w, intcpt,
  prices, totexp, instr)` with replacement and refits the whole pipeline
  (`quaidsFit()` then `quaidsCurvatureFit()`) on each resample. No default
  replication count (`B` is required) given a roughly order-of-magnitude
  timing gap between a single AIDS curvature fit (~0.9s) and a single
  QUAIDS curvature fit (~7.3s).
- `tests/quaids_curvature_bootstrap_test.e` (26 checks; opt-in via a new
  `-SkipBootstrap` flag in `tests/run_source_tests.ps1`, since it adds
  ~45-50s to that script's runtime even at a small `B`).
- `.github/workflows/tests.yml` (Milestone 14, no version bump of its own
  at the time): push-triggered CI on a self-hosted GitHub Actions runner,
  since this repo's tests require licensed GAUSS unavailable on GitHub-
  hosted runners. Triggers on `push` to `master` only (never
  `pull_request`), a deliberate mitigation for the fork/PR code-execution
  risk self-hosted runners carry on public repositories. Passes
  `-SkipBootstrap` to keep the routine per-push run fast.

### Fixed
- `quaidsCurvatureFit()`'s two internal `eighv()` calls (the warm-start
  and the Hessian-based standard error) could crash the entire calling
  job on a sufficiently degenerate input (`error G0528: More returns than
  targets`) -- a call-arity failure mode that this codebase's usual
  `trap`/`scalmiss` guard does not intercept. Found by stress-testing
  under bootstrap resampling. Fixed with an explicit pre-call finiteness
  check (both NaN and plain `Inf`, which `x .eq x` alone does not catch)
  ahead of both calls, falling back to a harmless identity decomposition
  on failure.
- `tests/run_source_tests.ps1`'s `Invoke-GaussBatch` helper could deadlock
  (child process hangs indefinitely, observed for over eight hours before
  being killed) when a test produces enough combined stdout+stderr output
  to fill both OS pipe buffers -- a classic .NET `Process` issue caused by
  reading stdout fully before stderr. `quaids_curvature_bootstrap_test.e`
  is the first test file to produce enough stderr volume (harmless
  "Optmt: function evaluation failed" diagnostics from bad bootstrap
  resamples) to trigger it. Fixed by draining both streams asynchronously
  via `OutputDataReceived`/`ErrorDataReceived` events instead of
  sequential `ReadToEnd()` calls.

### Notes
- Building the bootstrap's plausibility test found that the delta-method
  SE's boundary-inference unreliability does not stay confined to the
  specific gamma row/column tied to a boundary Cholesky entry -- it can
  inflate the reported SE anywhere in the coefficient vector, since the
  classical NLS covariance is for the whole `vech(A)` vector at once. See
  `GOLD_STANDARD_TODO.md`'s Milestone 15 section.

## 0.9.0 - 2026-07-25

### Added
- `quaidsCurvatureFit()` (`src/quaidscurvature.src`) now accepts QUAIDS
  fits (`aCtl.linear=0`), not just AIDS -- deferred at Milestone 10
  because QUAIDS's Slutzky matrix has an extra `lambda`-dependent cross-
  term. Resolved using the same lag-then-solve trick `quaidsFit()`'s own
  iteration already uses: `beta`/`lambda` are profiled out by OLS every
  outer round rather than joining `optmt`'s searched parameter set, which
  stays `vech(A)`-only, unchanged in size. `quaidsCurvatureFit()`'s
  signature is unchanged; no new proc, no new required struct field --
  version bumped anyway, since this is a real capability addition (QUAIDS
  curvature imposition going from unsupported to supported), not a bugfix.
- `tests/quaids_curvature_test.e` extended (17 -> 31 checks) with a QUAIDS
  block: convergence, exact negative-semidefiniteness at the reference
  point, non-vacuousness, shape/finiteness -- deliberately does NOT test
  "recovers a known true curvature-consistent gamma" the way the AIDS
  block does (see Notes).

### Fixed
- A real, pre-existing bug in `_quaidsCurvRecoverFull()`, latent since
  Milestone 10: the adding-up recovery for extra intercept-shifter rows
  (`aCtl`'s `intcpt` argument with `nint>0`) incorrectly applied the
  CONSTANT row's "sums to 1" formula to every intercept row, including
  shifter rows that must sum to 0. Never triggered before because the
  AIDS-only curvature fixture deliberately has no intercept shifters
  (`nint=0`); found only once a `nint>0` QUAIDS fixture was tried.

### Notes
- QUAIDS's curvature outer loop is measurably less stable than AIDS's own
  two-block version (found empirically) -- `aCtl.relax` (Milestone 12)
  is effectively required, not just optional, for QUAIDS; `maxOuterIter`
  (an internal constant) was bumped from 50 to 300 to accommodate the
  slower damped convergence this requires. AIDS is unaffected (still
  converges in ~10-20 iterations, well under either cap).
- A dedicated QUAIDS analog of the AIDS-only curvature fixture (a true
  gamma that is curvature-consistent at its own self-consistent sample
  mean) was attempted but not shipped: a broad screen found dozens of
  seeds where the construction is numerically self-consistent and
  genuinely negative-semidefinite, but every one implied economically
  implausible mean budget shares. The QUAIDS test instead validates
  against the existing general QUAIDS fixture, checking convergence/NSD/
  shape rather than true-parameter recovery -- a real, documented,
  weaker tier of evidence than the AIDS block's. See
  `GOLD_STANDARD_TODO.md`'s Milestone 13 section.

## 0.8.0 - 2026-07-23

### Added
- `aCtl.relax` (`src/quaids.sdf`/`quaidsutil.src`): optional under-
  relaxation (damping) of `quaidsFit()`'s iterated translog-price-index
  fixed-point update, `(0,1]`, default `1` (no damping, byte-identical to
  every prior release). Empirically the one mitigation of three tested
  that measurably reduced the iterated estimator's convergence-failure
  rate, at `relax=.75` -- see `GOLD_STANDARD_TODO.md`'s Milestone 12.
- `tests/quaids_convergence_sweep.e`/`tests/run_convergence_sweep.ps1`: a
  real, committed 200-seed x 2-model convergence-reliability diagnostic,
  replacing an informal 8-seed probe referenced since Milestone 3 that
  never survived as a repo artifact. Diagnostic only, not a pass/fail
  gate.
- `tests/quaids_reliability_regression_test.e` (8 checks): regression
  guard for the fixes below.

### Fixed
- A real crash: an unguarded `invpd()` in `quaidsFit()`'s symmetry-test
  block threw `error G0121: Matrix not positive definite` and aborted the
  entire call for the caller whenever a badly-diverged iterated fit
  produced a non-positive-definite variance block -- found by the new
  200-seed sweep, not a hypothetical. Now degrades gracefully
  (`qOut.symValid=0`, `qOut.converged` reflects the iteration's own
  state) instead of crashing.
- A near-zero-denominator guard on `quaidsFit()`'s convergence check
  (`b0 + (b0 .== 0)`, matching `quaidscurvature.src`'s own analogous
  guard) -- a real correctness gap (this codebase's own synthetic
  fixtures round true coefficients to the nearest 0.1, so exact zeros are
  a real occurrence), but measured via the sweep to have **zero effect**
  on the observed failure rate. Kept and documented as an honest
  non-result, not omitted.

### Notes
- The 200-seed sweep replaced "roughly half of seeds fail" with precise,
  reproducible numbers: iterated AIDS fails (never-converges or converges
  to a wrong answer) 58% of the time at default settings, QUAIDS 76% --
  see `docs/USAGE_GUIDE.md` and `docs/FEATURE_SUPPORT_MATRIX.md`. This
  milestone characterizes and partially mitigates the instability; it
  does not claim to solve it.

## 0.7.0 - 2026-07-23

### Added
- `quaidsWelfareFit()`/`printQuaidsWelfare()` (`src/quaidswelfare.src`):
  exact compensating variation (CV) and equivalent variation (EV) for a
  price change, holding nominal expenditure fixed, with delta-method
  standard errors. Works for LA-AIDS, iterated AIDS, and QUAIDS alike --
  a closed-form evaluation of the already-fitted expenditure function at
  two points, no new estimation and no new package dependency (unlike
  Milestone 10's curvature imposition).
- `quaidsWelfareOut` struct (`src/quaids.sdf`).
- `tests/quaids_welfare_test.e` (20 checks): exact zero-price-change
  identity, exact round-trip inverse-function identity, first-order
  (Marshallian) limiting-case consistency, and sign-agreement between CV
  and EV -- exercised on both a QUAIDS and an AIDS/LA-AIDS fit from the
  existing `_quaidsSyntheticDGP` fixture.

### Notes
- The QUAIDS expenditure-function inversion (Banks, Blundell & Lewbel
  1997) was verified by Shephard's-lemma cross-check against the already-
  validated QUAIDS share equation before implementation, after an initial
  from-memory attempt failed that same check (a sign/inversion error,
  caught before any code was written). See `docs/METHODOLOGY_NOTES.md`.

## 0.6.0 - 2026-07-22

### Added
- `quaidsCurvatureFit()`/`printQuaidsCurvature()` (`src/quaidscurvature.src`):
  local curvature (Slutzky negative semidefiniteness) imposition for
  LA-AIDS/AIDS via the Diewert-Wales (1987) Cholesky reparametrization,
  evaluated at the sample mean. Built on GAUSS's `optmt` package (now a
  real package dependency -- `package.json`'s `deps` array is no longer
  empty) for the small (`vech(A)`-only) profiled nonlinear IV step; the
  remaining coefficients are exactly identified by OLS once the
  reparametrized `gamma` is substituted in, reusing the same
  `moment`/`solpd` primitives as `quaidsFit()` itself.
- `quaidsCurvOut` struct (`src/quaids.sdf`).
- `tests/quaidsfixtures.src`: `_quaidsCurvatureSyntheticDGP()`, a
  synthetic dataset whose true gamma is curvature-consistent at its own
  sample mean by construction, found by direct seed screening (most seeds
  make the underlying nonlinear system numerically unstable -- the same
  kind of iteration sensitivity already documented for the main estimator).
- `tests/quaids_curvature_test.e` (17 checks): recovery, exact adding-up/
  homogeneity/symmetry, near-exact negative-semidefiniteness at the
  reference point, and a non-vacuousness check confirming the
  unconstrained fit genuinely violates curvature on this fixture.
- QUAIDS (the quadratic extension) curvature imposition is explicitly
  deferred, not silently absent -- see `docs/FEATURE_SUPPORT_MATRIX.md`
  and `GOLD_STANDARD_TODO.md`'s Milestone 10 section.

### Known limitation
- Standard errors from `quaidsCurvatureFit()` are a simplified,
  homoskedastic-NLS delta-method approximation, not a full SUR/GMM
  sandwich, and are unreliable specifically when the estimated Cholesky
  factor has boundary (near-zero) entries -- a known, documented
  complication of Cholesky-based negative-semidefinite-cone estimation,
  not unique to this implementation. Point estimates and the exact
  curvature property are unaffected.

## 0.5.0 - 2026-07-20

### Added
- `src/pubtable_quaids.src`: optional adapter onto the `pubtable` package
  for LaTeX/Markdown/CSV/RTF/HTML/XLSX export. `ptModelFromQuaids`/
  `ptFromQuaids` build coefficient tables (one comparison column per good)
  from `quaidsOut`; `ptModelFromQuaidsElas`/`ptFromQuaidsElas`/
  `ptTablesFromQuaidsElas` build income/uncompensated/compensated
  elasticity tables from `quaidsElasOut`; `ptFromQuaidsFamily` dispatches
  on either struct type. Guarded by `#ifDef QUAIDS_SDF_INCLUDED`, following
  the same pattern as `pubtable`'s own bundled `pubtable_qardl.src`, but
  kept inside this repo rather than the installed `pubtable` package (see
  `GOLD_STANDARD_TODO.md`'s Milestone 6 section for why).
- `examples/pubtable_export_example.e`: exports a coefficient table and all
  three elasticity tables to `.tex`/`.md`/`.csv`.
- `tests/quaids_pubtable_test.e` (30 checks): exact numeric parity between
  `pubtable`'s `ptModel.estimates`/`stdErrors` and the `quaidsOut`/
  `quaidsElasOut` values they're built from, plus an end-to-end export
  smoke test that reads the generated files back.
- `#ifndef`/`#define QUAIDS_SDF_INCLUDED`/`#endif` include guard on
  `src/quaids.sdf`, so `pubtable_quaids.src` has a symbol to detect.
- Package build/release tooling (`scripts/build_lcg.ps1`,
  `scripts/build_package.ps1`, `scripts/verify_release_artifact.ps1`,
  `scripts/run_release_verification.ps1`, `tests/verify_package_manifest.ps1`,
  `tests/run_source_tests.ps1`) and an installed-package public API gate
  (`tests/package_public_api.e`). The package is now actually installed and
  loadable via `library quaids;`.
- Full documentation set: `README.md`, `docs/COMMAND_REFERENCE.md` plus 18
  `docs/command-reference/*.md` pages, `docs/USAGE_GUIDE.md`,
  `docs/METHODOLOGY_NOTES.md`, `docs/FEATURE_SUPPORT_MATRIX.md`.
- Published-data validation extended to iterated AIDS (`aCtl.linear=1,
  aCtl.maxiter>1`) against R's `micEconAids::aidsEst(method="IL", ...)`,
  alongside the existing LA-AIDS check
  (`tests/quaids_published_validation_test.e`, now 19 checks).

### Changed
- `examples/quaids_example.e` and `examples/pubtable_export_example.e` now
  use `library quaids;` against the installed package, matching
  README.md/USAGE_GUIDE.md's documented usage pattern, instead of
  `#include`-ing the source tree directly (only meaningful once the
  package became installable).

### Fixed
- `tests/package_public_api.e` now actually calls `printQuaids()` and
  `quaidsElas()` -- both are required, `package.json`-listed public procs
  that the installed-package gate exercised only indirectly (via their
  split `quaidsFit`/`quaidsElasFit`/`printQuaidsElas` components) until the
  final integration pass found the gap.

No public API in `src/` changed across this tooling/documentation/testing
work, so the package version did not bump for it -- see
`GOLD_STANDARD_TODO.md`'s Milestones 7-9 for the full writeups.

## 0.4.0 - 2026-07-20

### Added
- `quaidsElasFit(b, v, intcpt, prices, totexp, aCtl)`: silent,
  struct-returning elasticity computation (`quaidsElasOut`: point estimates
  plus delta-method standard errors) at any evaluation point.
- `printQuaidsElas(elasOut)`: separated console printer for a
  `quaidsElasOut`.
- `tests/quaids_elasticities_test.e` (17 checks): parity between
  `quaidsElasFit()`/`printQuaidsElas()` and the pre-split `quaidsElas_()`,
  plus three exact algebraic identities (Engel aggregation, Cournot
  aggregation, elasticity homogeneity) checked at an out-of-sample
  observation and a synthetic counterfactual price scenario.

### Changed
- `quaidsElas(b, v, intcpt, prices, totexp, aCtl)` is now a thin wrapper
  (`quaidsElasFit()` then `printQuaidsElas()`); its signature and printed
  output are unchanged.

## 0.3.0 - 2026-07-20

### Added
- `src/quaidstests.src`: `quaidsHomogeneityTest(qOut)` and
  `quaidsJointTest(qOut)`, standalone Wald chi-squared tests against an
  unconstrained (`aCtl.homogenous == 0`) `quaidsOut` fit.
- `tests/quaids_hypothesis_tests_test.e` (19 checks): size and power for
  both new tests, a power check for the existing symmetry-given-homogeneity
  test, and the first exercise in this repo's history of the
  overidentification test (`ninst > nu`).

## 0.2.0 - 2026-07-20

### Added
- `src/quaids.src` split into `quaidsiv.src` (`_quaidsIVFirstStage()`),
  `quaidselas.src` (elasticities), `quaidsslutzky.src` (Slutzky negativity
  diagnostic), and `quaids.src` (core estimation, printing, legacy
  wrapper).
- `quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars,
  aCtl)`: dataframe/column-name entry point (`src/quaidsformula.src`).
- `tests/quaids_formula_parity_test.e` (17 checks): `quaidsFull()` vs.
  `quaidsFit()` parity, including the `extraVars == 0` path.
- `tests/quaidsfixtures.src`: shared synthetic 5-good AIDS/QUAIDS data
  generator used by the validation test suite below.
- `tests/quaids_synthetic_validation_test.e` (22 checks): recovers true DGP
  parameters within a documented tolerance across LA-AIDS/iterated-AIDS/
  QUAIDS x with/without-IV.
- `tests/quaids_published_validation_test.e` (11 checks): cross-checks
  `quaidsFit()` against an independent R (`micEconAids`) reference on
  published `Blanciforti86` food-demand data.

### Fixed
- Stone price index starting value (`quaidsFit()`'s "STARTING VALUE"
  block): was computed on a partially-relative, partially-absolute price
  matrix, silently producing wrong results for every `aCtl.maxiter==1`
  (LA-AIDS) call using default starting values, and a bad (though
  self-correcting via iteration) starting point for `aCtl.maxiter>1`
  calls. Found via the published-data cross-check against R (~5x
  discrepancy), confirmed by isolating the formula in Python. After the
  fix, GAUSS matches the R reference to within ~0.021 max absolute
  difference (was ~5x off before).
- Two dead-code bugs in `quaidsFit()`'s `intcpt == 0` branch (uninitialized
  `xnam`; a `G0071 Type mismatch` assigning a native string literal into a
  `matrix`-typed struct field), surfaced by the first-ever exercise of that
  code path in `quaids_formula_parity_test.e`.

## 0.1.0 - 2026-07-20

### Added
- Initial `quaids` package scaffold: renamed from `aids` (Milestone 0;
  QUAIDS is the more general model actually implemented, and `aids` risked
  colliding as a bare identifier). Dead code removed (`rankControl`/
  `latentControl` structs and their unused constructors). `package.json`,
  `LICENSE` (MIT), `CITATION.cff`, `.gitignore`, `CLAUDE.md` added.
- `quaidsFit(w, intcpt, prices, totexp, instr, aCtl)`: silent,
  struct-returning estimation core (`quaidsOut`, ~75 fields), split out of
  the original combined estimate-and-print `quaids()`/`aids()` proc.
- `printQuaids(qOut)`: separated console estimation-report printer.
- `tests/quaids_schema_test.e` (34 checks): `quaidsOut` field values/
  shapes, that `quaidsFit()` prints nothing, and that the legacy `quaids()`
  wrapper's returned matrices are byte-identical to the struct fields
  they're drawn from.

### Changed
- `quaids(w, intcpt, prices, totexp, instr, aCtl)` (formerly `aids()`) is
  now a thin wrapper: `quaidsFit()` then `printQuaids()` then the legacy
  elasticities-at-four-points/descriptive-statistics/Slutzky report.
  Signature and printed output verified byte-for-byte identical to the
  pre-split code at the time of the split.
