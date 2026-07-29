# quaidsRobustBootstrapCovariance

## Purpose

Computes the empirical covariance of
[quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md)'s reduced
bootstrap coefficient draws and expands it into [quaidsFit](quaidsFit.md)'s
full `qOut.bestB` coefficient basis.

## Format

```gauss
{ bBoot, vBoot } = quaidsRobustBootstrapCovariance(qOut, rbOut, aCtl);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - the base fit for the same sample used
  by [quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md).
- `rbOut` (*`quaidsRobustBootOut` structure*) - output from
  [quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md).
- `aCtl` (*`quaidsControl` structure*) - the same control structure used
  for the original fit.

## Returns

- `bBoot` - equal to `qOut.bestB`, returned for convenient direct use with
  post-estimation procs.
- `vBoot` - the empirical bootstrap covariance transformed into the full
  `vec(qOut.bestB)` basis.

## Remarks

This proc requires at least two completed bootstrap replications. The
raw bootstrap draws in `rbOut.bBoot` are stored in the same reduced basis
as `rbOut.b` and `rbOut.seBoot`; this helper first forms the empirical
covariance of those draws, then applies the same full-basis linear
recovery used by [quaidsRobustCovariance](quaidsRobustCovariance.md).

Use this helper when the closed-form robust sandwich is too conservative
for the application and you want bootstrap uncertainty propagated through
shares, elasticities, or welfare measures.

## Example

```gauss
struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr,
    aCtl, householdId, 200, 42);

{ bB, vB } = quaidsRobustBootstrapCovariance(qOut, rbOut, aCtl);

struct quaidsSharesOut sharesB;
sharesB = quaidsSharesFit(bB, vB, intcptPt, pricesPt, totexpPt, aCtl);
```

## Source

`quaidsrobust.src`

## See Also

[quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md),
[quaidsRobustCovariance](quaidsRobustCovariance.md),
[quaidsSharesFit](quaidsSharesFit.md), [quaidsElasFit](quaidsElasFit.md),
[quaidsWelfareFit](quaidsWelfareFit.md)
