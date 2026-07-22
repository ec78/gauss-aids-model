# quaidsElasFit

## Purpose

Computes income and price elasticities, and their asymptotic (delta-method)
standard errors, at one arbitrary evaluation point, silently (no console
output), returning a `quaidsElasOut` structure.

## Format

```gauss
struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(b, v, intcpt, prices, totexp, aCtl);
```

## Parameters

- `b` (*matrix*) - parameter matrix; each column corresponds to a specific
  share; row order is intercept | prices | lx | lx2 | u. Typically
  `qOut.bestB`. Do not include the `u` (IV control-function residual)
  block -- it is not used to compute elasticities; extra trailing rows
  beyond `1+nint+n+nendog` are silently ignored.
- `v` (*matrix*) - variance of `vec(b)`. Typically `qOut.bestV`.
- `intcpt` (*vector*) - independent variables at the evaluation point,
  **including** the leading constant (length `1+nint`); for a fit with no
  extra intercept shifters, pass a scalar `1`.
- `prices` (*vector*) - log prices (absolute) at the evaluation point.
- `totexp` (*scalar*) - log total expenditure at the evaluation point.
- `aCtl` (*`quaidsControl` structure*).

The evaluation point (`intcpt`, `prices`, `totexp`) can be **any** triple
in the same units the model was fit in -- the sample mean, a quartile, a
specific observation, or a fully synthetic counterfactual scenario. It
does not need to come from the data at all.

## Returns

`elasOut` is a `quaidsElasOut` structure:

- `er` (*Nx1*) - income elasticities.
- `ep` (*NxN*) - uncompensated (Marshallian) price elasticities; `ep[i,j]`
  is good `i`'s elasticity with respect to price `j`.
- `epc` (*NxN*) - compensated (Hicksian) price elasticities.
- `ser`, `sep`, `sepc` - delta-method standard errors of `er`/`ep`/`epc`
  (raw numeric values, not pre-formatted display strings).
- `n`, `wnam` - metadata (number of goods, good names).

## Remarks

Given adding-up/homogeneity hold (which this estimator always imposes by
construction), the returned elasticities satisfy three *exact* algebraic
identities to floating-point precision:

- Engel aggregation: `sum_i(w_i * er_i) = 1`
- Cournot aggregation: `sum_i(w_i * ep_ij) + w_j = 0`, for each price `j`
- Elasticity homogeneity: `sum_j(ep_ij) + er_i = 0`, for each good `i`

where `w` is the **model-implied** share at the evaluation point, not a
noisy observed share -- see `tests/quaids_elasticities_test.e` for a
worked example computing the model-implied share and checking all three
identities.

For a printed report, use [printQuaidsElas](printQuaidsElas.md), or the
backward-compatible [quaidsElas](quaidsElas.md) wrapper for both in one
call.

## Examples

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
print elasOut.er;
```

## Source

`quaidselas.src`

## See Also

[printQuaidsElas](printQuaidsElas.md), [quaidsElas](quaidsElas.md),
[quaidsSlutzky](quaidsSlutzky.md), [ptModelFromQuaidsElas](ptModelFromQuaidsElas.md)
