# quaidsElas_

## Purpose

**Deprecated compatibility procedure.** Low-level income/price elasticity
computation at one evaluation point, returning raw matrices (not a
struct). Retained callable through the `0.x` series for source
compatibility with code written before [quaidsElasFit](quaidsElasFit.md)
existed; not part of the supported `0.x` public API surface for new code.
See [`docs/public-api.json`](../public-api.json)'s
`compatibility_procedures` entry.

## Format

```gauss
{ er, ep, epc } = quaidsElas_(b, intcpt, prices, totexp, aCtl);
```

## Parameters

Identical to [quaidsElasFit](quaidsElasFit.md)'s `b`/`intcpt`/`prices`/
`totexp`/`aCtl` parameters.

## Returns

- `er` (*Nx1*) - income elasticities.
- `ep` (*NxN*) - uncompensated (Marshallian) price elasticities.
- `epc` (*NxN*) - compensated (Hicksian) price elasticities.

No standard errors -- `quaidsElas_()` never computed them; that is the
main reason [quaidsElasFit](quaidsElasFit.md) exists.

## Remarks

Use [quaidsElasFit](quaidsElasFit.md) in new application code: it
returns a `quaidsElasOut` structure with the same point estimates plus
delta-method standard errors, and is what
[printQuaidsElas](printQuaidsElas.md)/[quaidsElas](quaidsElas.md)/
[quaidsSlutzky](quaidsSlutzky.md)-adjacent workflows are built on.
`quaidsElas_()` is not scheduled for removal before `1.0.0` (see the
[compatibility policy](../../README.md#compatibility-policy)), but
carries no further development.

## Examples

```gauss
{ er, ep, epc } = quaidsElas_(qOut.bestB, intcptPt, pricesPt, totexpPt, aCtl);
```

Prefer:

```gauss
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
print elasOut.er;
```

## Source

`quaidselas.src`

## See Also

[quaidsElasFit](quaidsElasFit.md), [quaidsElas](quaidsElas.md)
