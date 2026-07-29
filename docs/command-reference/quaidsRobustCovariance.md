# quaidsRobustCovariance

## Purpose

Expands [quaidsRobustFit](quaidsRobustFit.md)'s reduced robust or
cluster-robust covariance into [quaidsFit](quaidsFit.md)'s full
`qOut.bestB` coefficient basis. Use this when you want robust/cluster-
robust delta-method standard errors from [quaidsSharesFit](quaidsSharesFit.md),
[quaidsElasFit](quaidsElasFit.md), or [quaidsWelfareFit](quaidsWelfareFit.md).

## Format

```gauss
{ bRobust, vRobust } = quaidsRobustCovariance(qOut, rOut, aCtl);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - the converged fit passed to
  [quaidsRobustFit](quaidsRobustFit.md).
- `rOut` (*`quaidsRobustOut` structure*) - output from
  [quaidsRobustFit](quaidsRobustFit.md) for the same fit and sample.
- `aCtl` (*`quaidsControl` structure*) - the same control structure used
  for the original fit.

## Returns

- `bRobust` - equal to `qOut.bestB`, returned for convenient direct use
  with post-estimation procs.
- `vRobust` - the robust/cluster-robust covariance transformed into the
  full `vec(qOut.bestB)` basis.

## Remarks

[quaidsRobustFit](quaidsRobustFit.md) reports covariance in the reduced
basis actually estimated by the robust sandwich: the independently-
estimated `n1` equations and, under homogeneity, the independent
top-left gamma block. The existing post-estimation procs consume the full
`qOut.bestB` basis, including adding-up and homogeneity recoveries. This
helper applies that linear recovery to the covariance before delta-method
propagation.

The transformation changes standard errors only. Point estimates remain
`qOut.bestB`.

## Example

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, householdId);

{ bR, vR } = quaidsRobustCovariance(qOut, rOut, aCtl);

struct quaidsElasOut elasR;
elasR = quaidsElasFit(bR, vR, intcptPt, pricesPt, totexpPt, aCtl);
```

## Source

`quaidsrobust.src`

## See Also

[quaidsRobustFit](quaidsRobustFit.md),
[quaidsRobustBootstrapCovariance](quaidsRobustBootstrapCovariance.md),
[quaidsSharesFit](quaidsSharesFit.md), [quaidsElasFit](quaidsElasFit.md),
[quaidsWelfareFit](quaidsWelfareFit.md)
