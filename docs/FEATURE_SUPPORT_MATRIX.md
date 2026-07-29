# Feature Support Matrix

This matrix summarizes the current public support surface across the three
model choices `quaidsFit()` (and its wrappers `quaids()`/`quaidsFull()`)
selects via `aCtl.linear`/`aCtl.maxiter` -- see the
[usage guide](USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids)
for the exact switch values.

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
| Quadratic-term (AIDS-vs-QUAIDS) specification test | Not applicable (no quadratic term to test) | Not applicable | Yes (`quaidsQuadraticTest`, since Milestone 17) |
| Elasticities at arbitrary points | Yes (`quaidsElasFit`) | Yes | Yes |
| Delta-method elasticity standard errors | Yes | Yes | Yes |
| Predicted budget shares at arbitrary points | Yes (`quaidsSharesFit`, since Milestone 16, no extra dependency) | Yes | Yes |
| Exact algebraic identity validation (Engel/Cournot/homogeneity) | Yes | Yes | Yes |
| Slutzky negativity diagnostic | Yes (`quaidsSlutzky`) | Yes | Yes |
| Preflight data/design diagnostics | Yes (`quaidsPreflight`, since Milestone 23) | Yes | Yes |
| Welfare measures (exact CV/EV) | Yes (`quaidsWelfareFit`, no extra dependency) | Yes | Yes |
| Curvature imposition | Yes (`quaidsCurvatureFit`, sample mean, requires `optmt` -- see Notes) | Yes (same) | Yes since Milestone 13, requires `aCtl.relax` -- see Notes |
| Curvature bootstrap standard errors | Yes (`quaidsCurvatureBootstrapFit`, since Milestone 15 -- see Notes) | Yes (same) | Yes (same) |
| Curvature bootstrap percentile CIs | Yes (`quaidsCurvatureBootstrapCI`, since Milestone 18) | Yes (same) | Yes (same) |
| Dataframe/column-name entry point | Yes (`quaidsFull`) | Yes | Yes |
| Formula-string (`"y ~ x"`) API | Not applicable (multi-equation system) | Not applicable | Not applicable |
| `pubtable` export (LaTeX/Markdown/CSV/...) | Yes (`src/pubtable_quaids.src`, optional) | Yes | Yes |
| Synthetic deterministic validation | Yes (`tests/quaids_synthetic_validation_test.e`) | Yes | Yes |
| Published-data cross-validation vs. R | Yes (`Blanciforti86` vs. 3SLS, `tests/quaids_published_validation_test.e`) | Yes (`Blanciforti86` vs. `method="IL"`, wider tolerance -- see Notes) | No independent reference implementation exists (see Notes) |
| Iteration convergence guarantee | Not applicable (one-step) | No -- a 200-seed sweep measured a 58% failure rate (never-converges or converges to a wrong answer) at default settings; check `qOut.converged`. `aCtl.relax=.75` measurably reduces this -- see Notes | No -- same caveat, 76% failure rate measured |
| Zero budget share correction | Yes (`quaidsZeroFit`, since Milestone 19, unconstrained only -- see Notes) | Yes (same) | Yes (same) |
| Robust / cluster-robust standard errors | Yes (`quaidsRobustFit`, since Milestone 20, simplified bread -- see Notes) | Yes (same) | Yes (same) |
| Robust / cluster bootstrap | Yes (`quaidsRobustBootstrapFit`, since Milestone 20) | Yes (same) | Yes (same) |
| Robust covariance propagation to shares/elasticities/welfare | Yes (`quaidsRobustCovariance`/`quaidsRobustBootstrapCovariance`, since Milestone 22) | Yes (same) | Yes (same) |
| Applied workflow bundle | Yes (`quaidsWorkflowFit` with compact preflight summary; `quaidsWorkflowScenarioFit` for CV/EV scenarios) | Yes (same) | Yes (same) |
| Installed-package (`library quaids;`) support | Yes | Yes | Yes |

## Notes

- **Iteration convergence guarantee** (Milestone 12): a real, committed
  200-seed x 2-model sweep (`tests/quaids_convergence_sweep.e`, replacing
  an informal 8-seed probe referenced since Milestone 3 that never
  survived as a repo artifact) measured, at default settings
  (`aCtl.relax=1`, `aCtl.err=.0001`, `aCtl.maxiter=100`): iterated AIDS
  never-converges 39% of the time and converges to a wrong answer (a
  self-consistent but incorrect fixed point, distinct from simply running
  out of iterations) another 19% (58% combined failure); QUAIDS
  never-converges 54.5% and converges wrong another 21.5% (76% combined).
  Two candidate fixes were tested empirically, not assumed: a near-zero-
  denominator guard on the convergence check had **zero measurable
  effect** (documented as an honest non-result); an opt-in damping field
  `aCtl.relax` (default `1`, i.e. off) measurably improved the correct-
  convergence rate at `relax=.75` (iterated AIDS 42%->43% correct, QUAIDS
  24%->26.5% correct) but more aggressive damping (`.5`, `.3`) did not
  help further and often made things worse. This is a modest, evidence-
  backed mitigation, not a solved problem -- see
  [Usage guide](USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids)
  and `GOLD_STANDARD_TODO.md`'s Milestone 12 section for the full grid
  and the separate, unrelated crash fix this same sweep also found (an
  unguarded `invpd()` in the symmetry-test block that used to abort the
  entire `quaidsFit()` call on a badly-diverged fit).
- "Always (control-function)" means `instr` is a required argument to
  every estimator entry point -- there is no exogenous-total-expenditure
  estimation mode in this library.
- Welfare measures (`quaidsWelfareFit`, Milestone 11) are exact, not
  approximated, for all three model choices -- unlike curvature
  imposition, computing CV/EV needs no new estimation, only a closed-form
  evaluation of the already-fitted expenditure function at two points, so
  there was no reason to scope QUAIDS out the way Milestone 10 did. See
  [Methodology Notes](METHODOLOGY_NOTES.md#welfare-measures).
- The published-data cross-validation (`Blanciforti86` vs. R's
  `micEconAids`) covers both LA-AIDS (`aCtl.linear=1, aCtl.maxiter=1`, vs.
  `aidsEst(..., instNames=...)`, 3SLS -- max abs difference ~0.021) and
  iterated AIDS (`aCtl.linear=1, aCtl.maxiter>1`, vs.
  `aidsEst(method="IL", ...)`, the Iterated Linear Least Squares Estimator
  -- max abs difference ~0.11, tolerance `0.15`). The iterated-AIDS
  comparison has a wider gap for a real, understood reason, not
  approximation slop: `micEconAids`'s `method="IL"` does not support
  instrumental variables (combining it with `instNames` segfaults R's
  `aidsEst` rather than erroring cleanly -- confirmed by direct testing),
  so that reference is SUR-estimated, while GAUSS's iterated fit always
  instruments log total expenditure. The comparison therefore spans both a
  different estimation algorithm *and* an IV-vs-no-IV difference. **QUAIDS
  has no independent reference implementation available**: `micEconAids`
  does not implement a quadratic log-expenditure term at all, and no other
  comparably-established QUAIDS implementation was found (see
  `GOLD_STANDARD_TODO.md`'s Milestone 3 section on the Python from-scratch
  replica, kept as supplementary evidence only). QUAIDS's validation is
  therefore the known-true synthetic-DGP recovery in
  `tests/quaids_synthetic_validation_test.e` -- a real, non-circular check
  (independently-generated data with known-true parameters, not just
  re-running the estimator on its own output), but a different, weaker
  tier of evidence than cross-implementation agreement on real published
  data. Documented here rather than silently claimed as equivalent.
- Curvature imposition (Diewert-Wales Cholesky reparametrization,
  `quaidsCurvatureFit`) is available for LA-AIDS/AIDS (`aCtl.linear=1`,
  Milestone 10) and QUAIDS (`aCtl.linear=0`, Milestone 13), imposed
  locally at the sample mean, requiring the `optmt` package
  (`package.json`'s `deps` array, no longer empty). QUAIDS was initially
  deferred at Milestone 10 -- its Slutzky matrix adds a `lambda`-
  dependent cross-term entangling three nonlinear parameter blocks
  instead of two -- but this resolved (Milestone 13) using the same
  lag-then-solve trick `quaidsFit()`'s own iteration already uses, with
  no growth in what `optmt` searches over. QUAIDS's curvature outer loop
  is measurably less stable than AIDS's own though: `aCtl.relax`
  (Milestone 12) is effectively required, not optional, for QUAIDS --
  undamped runs on the validation fixture diverge to NaN within a
  handful of iterations. Standard errors from `quaidsCurvatureFit` are a
  simplified delta-method approximation, known to be unreliable when the
  estimated Cholesky factor has boundary (near-zero) entries (a standard
  complication of Cholesky-based negative-semidefinite-cone estimation)
  -- point estimates and the exact curvature property at the reference
  point are unaffected. `quaidsCurvatureBootstrapFit` (Milestone 15) closes
  this gap with a nonparametric i.i.d. row bootstrap (resample, refit the
  whole pipeline, take the empirical SE), reported alongside rather than
  replacing the delta-method SE; it has no default replication count, since
  a single AIDS curvature fit and a single QUAIDS curvature fit differ in
  runtime by roughly an order of magnitude, making a one-size-fits-all
  default misleading. `quaidsCurvatureBootstrapCI` (Milestone 18) adds
  percentile confidence intervals directly from the bootstrap's raw
  draws, no new resampling needed. Building it surfaced a real, silent
  bug (present since Milestone 10/15): `quaidsCurvatureFit`'s `se` and
  `quaidsCurvatureBootstrapFit`'s `seBoot` had their individual cells
  scrambled relative to `b` (GAUSS's `reshape()` fills row-major, not
  column-major like `vec()`) -- invisible to shape/sign/finiteness
  checks, since those are permutation-invariant. Fixed in both places;
  `v` (the full covariance) was never affected. There is no independent
  published/cross-
  implementation validation for the *imposed* estimator on either model:
  even the R `micEconAids` reference implementation used elsewhere in
  this library only diagnoses curvature, never imposes it. For QUAIDS
  specifically, `tests/quaids_curvature_test.e` validates convergence/
  exact negative-semidefiniteness/non-vacuousness/shape rather than
  "recovers a known true curvature-consistent gamma" the way the AIDS
  block does -- a real, committed attempt to build a QUAIDS analog of the
  AIDS fixture found dozens of numerically self-consistent, genuinely
  NSD candidate seeds, but every one implied economically implausible
  mean budget shares, so this is a deliberately weaker (but still real)
  tier of evidence, documented as such rather than silently equated with
  AIDS's. See
  [Methodology Notes](METHODOLOGY_NOTES.md#curvature-imposition-diewert-wales)
  and `GOLD_STANDARD_TODO.md`'s Milestone 10 and 13 sections.

- Zero budget share correction (Shonkwiler-Yen, `quaidsZeroFit`, Milestone
  19) addresses real survey/microdata's corner solutions (zero-expenditure
  goods), which `quaidsFit()` does not model. A per-good first-stage
  probit's fitted probability `F_i` is divided into the second-stage share
  equation (`w_i/F_i = ...`), which avoids breaking the shared-design-
  matrix Kronecker-product identity every stage of `quaidsFit()` relies on
  -- a literal textbook implementation (rescaling every regressor by
  `F_i`) would not. **Unconstrained only in this pass** -- no homogeneity/
  symmetry imposition on the corrected model (errors if
  `aCtl.homogenous=1`); combining the method's own diagonal-delta
  restriction with a second, simultaneous homogeneity/symmetry
  restriction is left for a follow-up. Standard errors are a simplified
  `S .*. inv(gg)` formula that does not correct for the nonlinear
  translog-price-index feedback or first-stage probit/IV
  generated-regressor uncertainty. Adding-up does not hold exactly for the
  corrected coefficients -- a real, known property of the method itself,
  not a bug. Validated on a synthetic fixture with a known latent
  (uncensored) DGP and a genuine, non-degenerate zero-share censoring
  pattern (`tests/quaids_zero_test.e`): the corrected fit recovers the
  true latent parameters measurably better than naively fitting
  `quaidsFit()` on the same censored data. GAUSS's built-in `glm()` (used
  for the first-stage probits, no new package dependency) can hard-crash
  on some degenerate inputs, a known non-trappable failure mode not
  hardened against in this pass. See
  [Methodology Notes](METHODOLOGY_NOTES.md#zero-budget-share-correction-shonkwiler-yen).

- Robust / cluster-robust standard errors (`quaidsRobustFit`, Milestone
  20) generalize the pooled, homoskedastic `S.*.inv(gg)` sandwich every
  other covariance in this library uses to a per-observation
  (heteroskedasticity-robust) or per-cluster (cluster-robust, with a CR1
  small-sample correction) score aggregation -- genuinely new math, since
  neither GAUSS's base runtime nor the `tsmt` package's single-equation
  `robustSE`/`clusterSE` generalize to this library's stacked multi-
  equation system. Robust is the literal `G=nobs` special case of
  cluster-robust, unified through one `clusterId` argument, confirmed by
  an exact-identity regression test. Uses a **simplified bread**
  (`inv(gg)`-based, not `quaidsFit()`'s own nonlinear-feedback-corrected
  Jacobian), which was found, empirically, to make its `se` dramatically
  more conservative than `qOut`'s own classical SE -- confirmed to be an
  expected consequence of comparing a simple sandwich against the full
  cross-equation-efficient FGLS system (not a bug), via an independent
  hand-derivation using the same regressors/residuals.
  `quaidsRobustBootstrapFit` (also Milestone 20, shipped in the same pass
  rather than a later follow-up as Milestone 10/15's curvature/bootstrap
  split did) offers a cluster-aware nonparametric bootstrap alternative
  that resamples whole clusters and refits `quaidsFit()` itself, typically
  landing much closer to `qOut`'s own SE than the closed-form sandwich
  does. The reduced robust coefficient table covers only the `n1`
  independently-estimated equations, but Milestone 22 adds
  `quaidsRobustCovariance`/`quaidsRobustBootstrapCovariance` to expand the
  robust or bootstrap covariance into `qOut.bestB`'s full basis for
  elasticities, shares, and welfare. See
  [Methodology Notes](METHODOLOGY_NOTES.md#robust-and-cluster-robust-standard-errors).

Related documentation:

- [Usage guide](USAGE_GUIDE.md)
- [Methodology notes](METHODOLOGY_NOTES.md)
- [Command reference](COMMAND_REFERENCE.md)
