# quaidsSharesFit

## Purpose

Computes the model-implied predicted budget share vector -- and its
delta-method covariance/standard errors -- at one arbitrary evaluation
point. Useful for out-of-sample prediction and policy simulation without
hand-deriving the AIDS/QUAIDS share equation. Silent, no printing -- see
[printQuaidsShares](printQuaidsShares.md).

## Format

```gauss
struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(b, v, intcpt, prices, totexp, aCtl);
```

## Parameters

- `b` (*matrix*) - parameter matrix; each column corresponds to a
  specific share; row order is intcp|gamma|beta|[lambda]. Typically
  `qOut.bestB`. Do not include the "u" (IV-residual) part -- not used,
  same as [quaidsElasFit](quaidsElasFit.md)/[quaidsWelfareFit](quaidsWelfareFit.md).
- `v` (*matrix*) - variance of `vec(b)`. Typically `qOut.bestV`.
- `intcpt` (*vector*) - independent variables at the evaluation point,
  INCLUDING the leading constant (length `1+nint`) -- same convention as
  `quaidsElasFit`/`quaidsWelfareFit`; pass a scalar `1` if there are no
  extra intercept shifters.
- `prices` (*vector*) - log prices (absolute) at the evaluation point.
- `totexp` (*scalar*) - log total expenditure at the evaluation point.
- `aCtl` (*`quaidsControl` structure*).

## Returns

`sharesOut` is a `quaidsSharesOut` structure:

- `w` - `n x 1`, predicted budget shares at the point.
- `se` - `n x 1`, delta-method standard errors of `w` (equal to
  `sqrt(diag(v))`).
- `v` - `n x n`, the full delta-method covariance of `w` -- lets a caller
  test hypotheses spanning more than one good (e.g. whether two goods'
  shares differ significantly, or the variance of an aggregated group of
  goods) using the correct combined variance, not just marginal SEs.
- `n`, `wnam` - metadata (number of goods, good names).

## Remarks

**Budget shares only, not physical quantities.** `w_i` is the fraction of
total expenditure spent on good `i` -- the model's direct, unambiguous
output. Converting to a quantity demanded would require assuming `prices`/
`totexp` are logs of levels in mutually consistent units, which nothing
else in this library requires of the caller; deliberately out of scope.

**Adding-up holds exactly**: `sum(sharesOut.w) == 1` to floating-point
precision at any point, a direct consequence of `qOut.bestB`'s own
adding-up construction (`sum_i alpha_i = 1`, `sum_i gamma_ij = 0`,
`sum_i beta_i = 0`, `sum_i lambda_i = 0`) -- not something this proc
imposes itself.

**Standard errors**: a numerically differenced Jacobian (identical
stepsize/perturbation technique to `quaidsElasFit`/`quaidsWelfareFit`)
propagated as `jacW*v*jacW'` -- the same matrix-form delta method
`quaidsCurvatureFit` uses -- rather than the cheaper recursive-variance-
only formula those two procs use, since this proc reports the full `n x n`
covariance rather than only marginal SEs.

**Not a new formula**: the share equation itself already existed inside
`quaidsElas_()` (`src/quaidselas.src`) as an intermediate step toward
elasticities, and was independently duplicated in
`tests/quaids_elasticities_test.e`'s test-only `modelShareAt()` helper
(now replaced by a direct call to this proc). This is a third,
deliberately independent implementation of the same formula rather than a
refactor of `quaidsElas_()` -- see `GOLD_STANDARD_TODO.md`'s Milestone 16
section for why.

## Examples

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
intcptPt = meanc(qOut.intcptFull);
pricesPt = meanc(prices);
totexpPt = meanc(totexp);

struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

call printQuaidsShares(sharesOut);
print "sum of predicted shares:" sumc(sharesOut.w);   // exactly 1
```

## Source

`quaidsshares.src`

## See Also

[printQuaidsShares](printQuaidsShares.md), [quaidsFit](quaidsFit.md),
[quaidsElasFit](quaidsElasFit.md), [quaidsWelfareFit](quaidsWelfareFit.md)
