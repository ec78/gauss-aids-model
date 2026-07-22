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
   change falls below `aCtl.err` or `aCtl.maxiter` is reached.
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

## Slutzky Negativity

Demand theory implies the Slutzky (compensated price-response) matrix
should be negative semi-definite. `quaidsSlutzky()` computes this matrix
observation by observation across a sample and reports descriptive
statistics on its eigenvalues -- a diagnostic, not an imposed constraint.
Curvature **imposition** is out of scope for this release; see the
[usage guide's Limitations section](USAGE_GUIDE.md#limitations).

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312-326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527-539.
- Blanciforti, L., Green, R., King, G. (1986). *U.S. Consumer Behavior Over
  the Postwar Period: An Almost Ideal Demand System Analysis*. Giannini
  Foundation Monograph No. 40. Used for published-data cross-validation --
  see `tests/fixtures/published/SOURCE.md`.
