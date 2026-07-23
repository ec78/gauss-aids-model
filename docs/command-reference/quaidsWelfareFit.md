# quaidsWelfareFit

## Purpose

Computes exact compensating variation (CV) and equivalent variation (EV)
for a price change, holding nominal expenditure fixed, with delta-method
standard errors. Silent, no printing -- see
[printQuaidsWelfare](printQuaidsWelfare.md).

## Format

```gauss
struct quaidsWelfareOut wOut;
wOut = quaidsWelfareFit(b, v, intcpt, pricesPt0, pricesPt1, totexpPt0, aCtl);
```

## Parameters

- `b` (*matrix*) - parameter matrix; each column corresponds to a
  specific share; row order is intercept | prices | lx | lx2 | u.
  Typically `qOut.bestB`. Do not include the "u" part -- not used, same
  as [quaidsElasFit](quaidsElasFit.md).
- `v` (*matrix*) - variance of `vec(b)`. Typically `qOut.bestV`.
- `intcpt` (*vector*) - independent variables at the reference point,
  **including** the leading constant (length `1+nint`) -- same convention
  as [quaidsElasFit](quaidsElasFit.md).
- `pricesPt0` (*vector*) - log prices (absolute) at the **base** point.
- `pricesPt1` (*vector*) - log prices (absolute) at the **new** point.
- `totexpPt0` (*scalar*) - log total (base) expenditure. Nominal
  expenditure is held fixed at this level throughout -- the standard
  "what if prices change and income doesn't" applied setup.
- `aCtl` (*`quaidsControl` structure*).

## Returns

`wOut` is a `quaidsWelfareOut` structure:

- `cv`, `ev` (*scalar*) - compensating and equivalent variation, in
  expenditure units (money, not log-money). **Sign convention**: both are
  positive when the price change reduces welfare, negative when it
  improves welfare, and exactly zero when `pricesPt1 == pricesPt0`.
- `seCV`, `seEV` (*scalar*) - delta-method standard errors.
- `u0`, `u1` (*scalar*) - the implied log indirect-utility levels at the
  base point and at the new prices with unchanged nominal expenditure,
  for diagnostic purposes.
- `pricesPt0`, `pricesPt1`, `totexpPt0` - echoed inputs.

## Remarks

Works for LA-AIDS, iterated AIDS, and QUAIDS alike -- unlike
[quaidsCurvatureFit](quaidsCurvatureFit.md), this needs no new estimation
and no external package dependency, since it's a closed-form evaluation
of the already-fitted expenditure function at two points (`aCtl.linear`
controls whether the quadratic term is used, same as everywhere else).

**Formula** (Banks, Blundell & Lewbel 1997's QUAIDS indirect utility
function, inverted): `ln e(u,p) = a(p) + b(p) / (1/u - lambda(p))`, where
`a(p)`/`b(p)` are the same translog price-index/scale terms used
elsewhere in this library and `lambda(p)` is a weighted sum of the
QUAIDS quadratic coefficients against log prices (zero when
`aCtl.linear=1`). See
[Methodology Notes](../METHODOLOGY_NOTES.md#welfare-measures) for the
full derivation and how it was verified before implementation.

## Examples

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
intcptPt = meanc(qOut.intcptFull);
pricesPt0 = meanc(prices);
totexpPt0 = meanc(totexp);

// A hypothetical 5% price increase on good 1:
pricesPt1 = pricesPt0;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

struct quaidsWelfareOut wOut;
wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
call printQuaidsWelfare(wOut);
```

## Source

`quaidswelfare.src`

## See Also

[printQuaidsWelfare](printQuaidsWelfare.md), [quaidsElasFit](quaidsElasFit.md),
[quaidsFit](quaidsFit.md)
