# Methodology Notes

This file summarizes the estimator implemented by this package. It is
intentionally concise; command-reference pages document argument syntax
and return structures.

## The Almost Ideal Demand System

The Almost Ideal Demand System (Deaton & Muellbauer, 1980) models budget
shares `w_i` for `N` goods as a function of log prices and log real total
expenditure:

```
w_i = alpha_i + sum_j(gamma_ij * log(p_j)) + beta_i * log(x / P)
```

where `P` is a price index (the translog price aggregator `a(p)`) and `x`
is total expenditure. Demand theory implies:

- **Adding-up**: `sum_i(alpha_i) = 1`, `sum_i(gamma_ij) = 0` for each `j`,
  `sum_i(beta_i) = 0`.
- **Homogeneity** (degree zero in prices and expenditure): `sum_j(gamma_ij)
  = 0` for each `i`.
- **Slutzky symmetry**: `gamma_ij = gamma_ji`.

This package always imposes adding-up by construction (the `N`-th
equation is recovered from the other `N-1` via the adding-up identity,
never separately estimated) and optionally imposes homogeneity and
symmetry via iterated FGLS with cross-equation restrictions, applied
through a minimum-distance reparametrization
(`design(vec(xpnd(seqa(1,1,k*(k+1)/2))))` builds the selection matrix `R`
such that `vec(G) = R*vech(G)` for a symmetric `G`).

## Price Index: Stone vs. Translog

- **Stone price index** (linear approximation): `log(P) = sum_i(w_i *
  log(p_i))`, using observed budget shares. Used for the one-step
  **LA-AIDS** special case (`aCtl.maxiter == 1`) -- the model is then
  linear in parameters and estimable by a single FGLS/IV step.
- **Translog price index** (nonlinear): `log(P) = alpha_0 + sum_i(alpha_i *
  log(p_i)) + 0.5 * sum_i sum_j(gamma_ij * log(p_i) * log(p_j))`. Used for
  **iterated AIDS**/**QUAIDS** (`aCtl.maxiter > 1`) -- the price index
  depends on the parameters being estimated, so the model is fit by
  iterating: form the price index from the current parameter estimates,
  re-estimate, repeat until the relative parameter change falls below
  `aCtl.err` or `aCtl.maxiter` is reached. `aCtl.alpha0` fixes the
  translog intercept (a normalization choice, not estimated).

LA-AIDS never advances past the Stone-index starting value regardless of
`aCtl.linear` -- iteration only applies when `aCtl.maxiter > 1`.

## QUAIDS: The Quadratic Extension

Banks, Blundell & Lewbel (1997) extend the linear-in-`log(x)` AIDS Engel
curve with a quadratic term:

```
w_i = alpha_i + sum_j(gamma_ij * log(p_j)) + beta_i * log(x/P) + (lambda_i / b(p)) * [log(x/P)]^2
```

where `b(p) = exp(sum_i(beta_i * log(p_i)))`. This allows goods to be
luxuries at low expenditure and necessities at high expenditure (or vice
versa), which a purely linear Engel curve cannot represent. Set
`aCtl.linear = 0` (with `aCtl.maxiter > 1`) to fit QUAIDS; `aCtl.linear =
1` drops the quadratic term (`lambda_i = 0` for all `i`), giving iterated
AIDS instead.

## Instrumental Variables: Control-Function Approach

Log total expenditure is treated as endogenous throughout (never
optional). The first stage (`_quaidsIVFirstStage()`, `src/quaidsiv.src`)
regresses log total expenditure on the supplied instruments; the
first-stage residual `u` is then included as an additional regressor in
every share equation (the control-function/residual-inclusion approach to
endogeneity, as opposed to two-stage least squares or three-stage least
squares). This is a different, but equally valid, IV algorithm from the
3SLS approach used by R's `micEconAids` package -- cross-validation
against `micEconAids` on published data found close but not bit-for-bit
agreement (~0.021 max absolute difference in `alpha`/`beta`/`gamma`),
consistent with two different valid algorithms for the same model rather
than a discrepancy in either implementation.

`ninst > nu` (more instruments than endogenous regressors) activates an
overidentification test comparing the instrumented and non-instrumented
coefficient estimates.

## Estimation Algorithm Summary

1. **Starting values**: linearized AIDS with the Stone price index (fixed
   for `aCtl.maxiter == 1`; a starting point for the iteration otherwise).
2. **IV first stage**: regress log total expenditure on instruments;
   residual becomes an added regressor.
3. **Homogeneity-constrained (or unconstrained) FGLS**, iterating the
   translog price index when `aCtl.maxiter > 1`, until relative parameter
   change (guarded against a near-zero denominator, `b0 + (b0 .== 0)`)
   falls below `aCtl.err` or `aCtl.maxiter` is reached. Each step is a
   plain fixed-point update (`b_new = relax*b + (1-relax)*b_old`,
   `aCtl.relax` default `1` = no damping); this has no global-convergence
   guarantee -- see `GOLD_STANDARD_TODO.md`'s Milestone 12 section for a
   measured failure-rate characterization and `aCtl.relax`'s effect.
4. **Overidentification test**, if `ninst > nu`.
5. **Symmetry test given homogeneity**, and a **symmetry-constrained
   re-estimation** via minimum distance, if `aCtl.homogenous == 1`.
6. **Recovery**: reparametrized/relative-price-form coefficients are
   converted back to absolute-price, full-`N`-equation form
   (`qOut.b`/`qOut.v`/`qOut.bS`/`qOut.vS`; `qOut.bestB`/`qOut.bestV` holds
   whichever is most-constrained).

Standard errors throughout are asymptotic (delta-method/FGLS), valid under
large-T, correctly-specified-instrument regularity conditions -- no claim
of finite-sample exactness.

## Elasticities

Income and price (Marshallian/Hicksian) elasticities are computed
analytically from the fitted budget-share function and its derivatives at
an arbitrary evaluation point, with delta-method standard errors from a
numerically differenced Jacobian. Given adding-up/homogeneity (always
imposed by construction), the resulting elasticities satisfy three exact
algebraic identities -- Engel aggregation, Cournot aggregation, and
elasticity homogeneity -- as consequences of the functional form, checked
to floating-point precision in `tests/quaids_elasticities_test.e`. See
[quaidsElasFit](command-reference/quaidsElasFit.md) for the formulas.

## Welfare Measures

[quaidsWelfareFit](command-reference/quaidsWelfareFit.md) computes exact
compensating variation (CV) and equivalent variation (EV) for a price
change, holding nominal expenditure fixed, for LA-AIDS, iterated AIDS,
and QUAIDS alike -- unlike curvature imposition, this needs no new
estimation, just a closed-form evaluation of the already-fitted
expenditure function at two points.

**The QUAIDS indirect utility function** (Banks, Blundell & Lewbel 1997):

```
ln V(x,p) = 1 / ( b(p)/L(p,x) + lambda(p) )
```

where `a(p)`/`b(p)` are the same translog price-index/scale terms used
throughout this library, `L(p,x) = ln(x) - a(p)` (the same `lx` computed
elsewhere), and `lambda(p) = lambda . p` (a scalar; zero when
`aCtl.linear=1`, since there is no quadratic term to weight).

**Inverting for the expenditure function** (the money needed to achieve
utility `u` at prices `p`):

```
ln e(u,p) = a(p) + b(p) / ( 1/u - lambda(p) )
```

**Verified before implementation, not assumed from memory**: an initial
derivation attempt produced a different (incorrect) formula. It was
checked -- and found wrong -- by confirming whether
`w_i = d ln(e)/d ln(p_i)` (Shephard's lemma, holding `u` fixed) reproduces
the *already validated* QUAIDS share equation
(`w_i = alpha_i + sum_j gamma_ij ln p_j + beta_i*L + (lambda_i/b(p))*L^2`).
The corrected formula above reproduces it exactly, term by term, and
collapses to the simpler, independently-verified AIDS expenditure
function `a(p) + u*b(p)` when `lambda=0`.

**Welfare measures**, comparing base prices `p0` (expenditure `x0`) to new
prices `p1`, with nominal expenditure held fixed at `x0` (sign
convention: both positive when the price change reduces welfare, zero
when `p0=p1`, matching Deaton & Muellbauer's own applied convention):

```
u0 = ln V(x0, p0)              // utility actually enjoyed at the base point
CV = e(p1, u0) - x0             // extra cost, at NEW prices, to keep the OLD utility
u1 = ln V(x0, p1)               // utility x0 actually buys at the NEW prices
EV = x0 - e(p0, u1)             // at OLD prices, cost difference between old and new utility
```

`tests/quaids_welfare_test.e` checks this via exact algebraic identities
rather than tolerance-based approximation: `CV`/`EV` are exactly zero at
`p0=p1`; feeding `e(p1,u0)` back into `ln V(.,p1)` returns `u0` exactly
(the defining inverse-function property); and `CV`/`EV` for a small price
change converge to the standard Marshallian first-order (share-weighted)
approximation as the change shrinks toward zero.

## Slutzky Negativity

Demand theory implies the Slutzky (compensated price-response) matrix
should be negative semi-definite. `quaidsSlutzky()` computes this matrix
observation by observation across a sample and reports descriptive
statistics on its eigenvalues -- a diagnostic, always available regardless
of whether curvature is imposed.

## Curvature Imposition (Diewert-Wales)

Concavity of a flexible functional form like AIDS/QUAIDS cannot be imposed
*globally* without over-restricting the functional form -- a standard
result in the demand-systems literature (Diewert & Wales, 1987), not a gap
in this implementation. [quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)
imposes it *locally*, at the sample mean, for LA-AIDS/AIDS
(`aCtl.linear = 1`) and, since Milestone 13, QUAIDS (`aCtl.linear = 0`)
too -- see [Feature Support Matrix](FEATURE_SUPPORT_MATRIX.md).

The reparametrization: write the upper-left `(n-1) x (n-1)` block of the
gamma matrix as `gamma = -A*A' - K0`, where `A` is lower triangular and
`K0` is the non-gamma part of `quaidsSlutzky()`'s own matrix, evaluated at
the sample mean. This makes the Slutzky matrix equal `-A*A'` at that
point by construction -- negative semidefinite for *any* `A`, regardless
of the data. Estimation is a profiled/concentrated nonlinear IV problem:
GAUSS's `optmt` package searches only over `A`'s free elements
(`vech(A)`), with the remaining coefficients exactly identified by OLS
once the reparametrized gamma is substituted in as a fixed offset --
within an outer iteration that mirrors `quaidsFit()`'s own translog-
price-index iteration. For QUAIDS, `beta` and `lambda` are profiled out
alongside the intercept/gamma the same way -- `beta` appears both in the
plain `beta'beta` term (shared with AIDS) and, via `b(p)=exp(beta'p)`,
inside a `lambda`-scaled cross-term; both are lagged by one outer round
to build that round's regressors, exactly mirroring how `quaidsFit()`'s
own loop lags `beta` to build its `lx2` regressor -- so `optmt`'s search
dimension (`vech(A)` only) never grows for QUAIDS. QUAIDS's outer loop is
measurably less stable than AIDS's own, found empirically: `aCtl.relax`
(the same field `quaidsFit()` uses, see "Estimation Algorithm Summary"
above) is effectively required for QUAIDS, not just optional.

Standard errors are a simplified, homoskedastic-NLS delta-method
approximation, not a full SUR/GMM sandwich, and are known to be
unreliable when the estimated Cholesky factor has boundary (near-zero)
entries -- a standard complication of Cholesky-based negative-
semidefinite-cone estimation, analogous to inference on non-negativity-
constrained variance components elsewhere in econometrics. Point
estimates and the exact curvature property are unaffected. See
[quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md) and the
[usage guide's Limitations section](USAGE_GUIDE.md#limitations).

### Bootstrap Standard Errors (Milestone 15)

[quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md)
closes the gap above with a nonparametric i.i.d. row (pairs) bootstrap:
resample `T` rows with replacement (one draw of row indices applied
identically to `w`/`intcpt`/`prices`/`totexp`/`instr`), refit the *whole*
pipeline (`quaidsFit()` then `quaidsCurvatureFit()`) on each resample, and
take the empirical standard deviation of the resulting coefficient draws.
Refitting the whole pipeline, not just perturbing the already-fitted
curvature struct, is deliberate -- it captures first-stage IV sampling
variability as well as the curvature reparametrization's own variability.
A replication is included only if both stages converge and the resulting
coefficient vector is finite; up to `5*B` resamples are attempted before
giving up on reaching `B` completed replications, matching the same
attempt-cap convention this author's `gauss-qardl` package already uses
for its own moving-block bootstrap.

No default replication count is given: a single AIDS curvature fit
measures well under a second, but a single QUAIDS curvature fit measures
several seconds (and typically needs `aCtl.relax < 1` to converge, as
above) -- at conventional bootstrap sizes (200-1000), QUAIDS alone can
take on the order of half an hour to a couple of hours. `B` is a required
argument so this tradeoff is chosen deliberately rather than defaulted.

Building this bootstrap surfaced a real, separate finding about
[quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md) itself: on a
sufficiently degenerate resample, its internal eigendecomposition calls
(`eighv()`, used both for the warm start and for the Hessian-based
standard error) can fail outright rather than merely return a poor answer
-- and this specific failure mode (GAUSS returning fewer values than a
multiple-assignment expects) is not intercepted by this codebase's usual
`trap`/`scalmiss` guard (confirmed directly: trapping the call alone did
not prevent it). `quaidsCurvatureFit()` was hardened with an explicit
pre-call finiteness check (catching both NaN, via `x .eq x`, and plain
Inf, which is *not* caught by that check since `Inf .eq Inf` is true under
IEEE 754 and requires a separate magnitude bound) ahead of both `eighv()`
calls, falling back to a harmless identity decomposition on failure --
this degrades the standard error for that one fit rather than crashing
the entire bootstrap run.

## Zero Budget Share Correction (Shonkwiler-Yen)

Real survey/microdata routinely has corner solutions: some households
report zero expenditure on some goods. The linear/log-linear AIDS/QUAIDS
share equation has no mechanism for a censored dependent variable, so
fitting `quaidsFit()` directly on such data is a known source of bias.
[quaidsZeroFit](command-reference/quaidsZeroFit.md) implements the
Shonkwiler & Yen (1999) two-step correction.

**The method**: for each good `i`, a first-stage probit estimates
`F_i = Pr(w_i > 0 | intcpt, prices, totexp)`. The second stage fits

```
w_i = F_i*(alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2]) + f_i*delta_i + e_i
```

where `f_i` is the probit density at the same point and `delta_i` is an
estimated coefficient on the resulting inverse-Mills-ratio-like hazard
term -- the mechanism that lets the corrected model account for the
selection into a positive share.

**The architectural constraint this library is built around**: every
stage of `quaidsFit()` (the GLS solve, the variance formula, the
homogeneity/symmetry minimum-distance restriction, the overidentification
test, the absolute-price recovery) relies on a Kronecker-product identity
(e.g. `S[1:n-1,1:n-1].*.gg`, `src/quaids.src`) that holds only because
every equation shares the *same* design matrix `X`. A literal Shonkwiler-
Yen implementation rescales *every* regressor in equation `i` by that
equation's own `F_i`, which breaks that shared-`X` assumption outright.

**The reformulation used here** avoids this by dividing the whole equation
by `F_i` (a known, first-stage-fitted quantity, held fixed during the
second stage):

```
w_i/F_i = alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2] + (f_i/F_i)*delta_i + e_i/F_i
```

This turns the problematic *regressor* rescaling into a **dependent-
variable transform** (`wTilde_i = w_i/F_i`, computed once from the
first-stage probit) plus **one new shared regressor column per equation**
(`h_i = f_i/F_i`) -- structurally the same kind of addition as the `u`
(IV-residual) column `quaidsFit()` already appends to its own shared `X`.
The shared-design-matrix machinery, and everything built on it (including
the whole translog-price-index outer iteration, reused unchanged with
`wTilde`/`h` substituted for `w`), survives intact.

**The diagonal-delta restriction**: appending `n` hazard columns to a
shared `X` means the one-shot GLS solve initially estimates a full
`n x n` cross-equation `delta` block (every equation's response to
*every* good's hazard term), when Shonkwiler-Yen only wants the diagonal
(good `i`'s own hazard term in good `i`'s own equation). This is imposed
by a minimum-distance restriction forcing the off-diagonal entries to
exactly zero -- the same `design()`-based selection-matrix construction
`quaidsFit()`'s own symmetry-restriction stage already uses, just with a
diagonal (not `gamma_ij=gamma_ji`) restriction pattern.

**Scope of this implementation, deliberately limited**:

- **Unconstrained only** -- no homogeneity/symmetry imposition on the
  corrected model in this pass (errors if `aCtl.homogenous = 1`).
  Combining the diagonal-delta restriction with a second, simultaneous
  homogeneity/symmetry restriction is real additional work, left for a
  follow-up.
- **Standard errors are a simplified formula**
  (`V(vec(b)) = S .*. inv(gg)`, the classic SUR-with-shared-regressors
  covariance) -- unlike `quaidsFit()`'s own variance, it does not correct
  for the nonlinear translog-price-index feedback, nor for first-stage
  probit/IV generated-regressor uncertainty.
- **Adding-up does not hold exactly** for the corrected coefficients --
  a real, known property of Shonkwiler-Yen itself (each equation is
  independently rescaled by its own good-specific `F_i`), not a bug. No
  post-hoc renormalization is applied.

**Validation**: `tests/quaids_zero_test.e` fits a synthetic 5-good QUAIDS
fixture (`_quaidsZeroSyntheticDGP`, `tests/quaidsfixtures.src`) whose
*latent* (uncensored) shares are known exactly, then censors them the
economically correct way -- `w_i = max(0, latent_i) / sum_j max(0,
latent_j)`, an accounting identity, not an ad hoc redistribution -- to
produce genuine, non-degenerate per-good zero-share fractions (found by
direct screening at `seed=1`: roughly 17-84% depending on the good). The
core check compares `quaidsZeroFit()`'s recovered coefficients against
the true latent-DGP parameters, and against a naive `quaidsFit()` fit to
the same censored data: the corrected fit recovers the truth measurably
better on both a max-absolute-difference and a mean-absolute-difference
basis. Also checks the diagonal-delta restriction holds exactly (off-
diagonal entries exactly zero, on-diagonal entries genuinely estimated).

**A known, unresolved limitation**: GAUSS's built-in `glm()` (used for the
first-stage probits, no new package dependency) can hard-crash on some
degenerate inputs (`error: Intel MKL ERROR ... DGELS`) -- a failure mode
that is *not* intercepted by this codebase's usual `trap 1,1;`/
`scalmiss()` guard idiom, the same class of non-trappable failure already
documented for `eighv()` inside
[quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md).
Confirmed empirically (some seeds trigger it, others don't); not hardened
against in this pass.

## Robust and Cluster-Robust Standard Errors

Every covariance elsewhere in this library (`quaidsFit()`'s
`qOut.homogV`/`symcV`, the IV first stage, `quaidsCurvatureFit()`'s
`cOut.v`, `quaidsZeroFit()`'s `zOut.v`) rests on a single pooled,
homoskedastic `S = sse/nobs` combined with the shared-design-matrix
`S.*.inv(gg)` Kronecker sandwich -- a real, if standard, simplification:
it assumes the same error covariance for every observation, and does not
allow for correlation within groups (e.g. multiple observations per
household or region). [quaidsRobustFit](command-reference/quaidsRobustFit.md)
generalizes this to a heteroskedasticity-robust (White) or cluster-robust
sandwich.

**Why this needed genuinely new math, not an existing utility**: GAUSS's
base runtime (`robustSE`/`clusterSE`/`hacSE` in `robust.src`) and the
separately-licensed `tsmt` package's near-identical procs are both
single-equation `(X'X)^-1 (...) (X'X)^-1` sandwiches for one dependent
variable and shared regressors -- neither generalizes to this library's
*stacked multi-equation* system, where `n1` share equations share the
same design matrix `X` but have their own residual columns. Reusing them
would require unpacking into exactly the same per-cluster score
aggregation this library builds directly, for zero net simplification --
the same conclusion this project reached about `gmmFitIV` at Milestone 2.

**The construction**: given an already-fitted `qOut` and the raw sample,

1. Rebuild per-observation model-implied fitted shares at `qOut.bestB`
   (the same formula `quaidsElas_()`/`quaidsSharesFit()` already
   duplicate independently, vectorized across the whole sample). Residuals
   `U = w[.,1:n1] - fittedW[.,1:n1]` (only the `n1` independently-
   estimated equations).
2. Rebuild the shared regressor block `X` at the converged point (one
   evaluation, not a re-iteration) -- the same reduced-form regressor set
   `quaidsFit()`'s own STARTING VALUE block uses (`intcptFull~pricesHybrid[.,1:n1]~endog~u`).
3. `Infl[.,(i-1)*k+1:i*k] = X.*U[.,i]` for each of the `n1` equations --
   the standard "per-observation score contribution."
4. Aggregate `Infl`'s rows by cluster label (or leave ungrouped for the
   robust case): `Meat = c*(InflG'InflG)/nobs`, with a CR1 finite-sample
   correction `c = (G/(G-1))*((nobs-1)/(nobs-K))` unifying the robust
   (`G=nobs`) and cluster (`G<nobs`) cases through the same formula.
5. `bread = eye(n1).*.inv(gg)`; `v = bread*Meat*bread`.

Robust is the literal `G=nobs` (every row its own cluster) special case
of cluster-robust -- confirmed by an exact-identity regression test
(`tests/quaids_robust_test.e`: `clusterId=0` and an explicit
`seqa(1,1,nobs)` per-row label give byte-identical output), not just
argued algebraically.

**A real, empirically-confirmed finding, not a theoretical worry**: this
sandwich's `bread` is `inv(gg)` -- it does *not* replicate
`quaidsFit()`'s own nonlinear-translog-price-index-feedback Jacobian
correction (its `Ji`/`J` construction). Building this milestone's test
surfaced that this makes `quaidsRobustFit()`'s `se` dramatically more
*conservative* (often more than an order of magnitude larger) than
`qOut.homogSE`/`symcSE` -- confirmed, via an independent hand-derivation
in `tests/quaids_robust_test.e`, to be entirely attributable to comparing
a simple equation-by-equation sandwich against the full cross-equation-
efficient FGLS system (the same hand-derived formula, applied to the
*same* regressors/residuals this proc itself uses, lands in the same
order of magnitude under genuine homoskedasticity) -- not a defect in the
sandwich formula itself.
[quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md)'s
bootstrap SE, which resamples and refits the actual efficient estimator,
does not share this gap and is typically much closer to `qOut`'s own SE.

**Scope, deliberately limited** (the same "new sibling, not a
modification of already-shipped code" choice as curvature/welfare/
shares/zero-correction): covers only the `n1` independently-estimated
equations in its printed coefficient table. Milestone 22 adds an explicit
linear recovery step,
[quaidsRobustCovariance](command-reference/quaidsRobustCovariance.md), that
expands this reduced covariance into `qOut.bestB`'s full basis: equation
`n` is recovered through adding-up, and under homogeneity the reference-
price gamma row/column are recovered through the same linear restrictions
used by the point estimates. The bootstrap analogue,
[quaidsRobustBootstrapCovariance](command-reference/quaidsRobustBootstrapCovariance.md),
forms the empirical covariance of `rbOut.bBoot` first, then applies the
same recovery. These helpers let
[quaidsSharesFit](command-reference/quaidsSharesFit.md),
[quaidsElasFit](command-reference/quaidsElasFit.md), and
[quaidsWelfareFit](command-reference/quaidsWelfareFit.md) consume robust
or bootstrap covariance directly, without changing point estimates or
mutating `qOut.symcV`. The cluster bootstrap still resamples whole clusters
and refits `quaidsFit()` only, not `quaidsRobustFit()`'s own sandwich,
mirroring `quaidsCurvatureBootstrapFit()`'s identical "refit the estimator,
not its SE stage, each replication" pattern.

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312-326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527-539.
- Blanciforti, L., Green, R., King, G. (1986). *U.S. Consumer Behavior Over
  the Postwar Period: An Almost Ideal Demand System Analysis*. Giannini
  Foundation Monograph No. 40. Used for published-data cross-validation --
  see `tests/fixtures/published/SOURCE.md`.
- Diewert, W. E., Wales, T. J. (1987). "Flexible Functional Forms and
  Global Curvature Conditions." *Econometrica*, 55(1), 43-68. The
  Cholesky reparametrization used by `quaidsCurvatureFit()`.
- Shonkwiler, J. S., Yen, S. T. (1999). "Two-Step Estimation of a
  Censored System of Equations." *American Journal of Agricultural
  Economics*, 81(4), 972-982. The two-step zero-budget-share correction
  used by `quaidsZeroFit()`.
