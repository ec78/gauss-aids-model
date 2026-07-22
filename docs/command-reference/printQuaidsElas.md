# printQuaidsElas

## Purpose

Prints the income and price elasticities and standard errors captured in a
`quaidsElasOut` structure.

## Format

```gauss
call printQuaidsElas(elasOut);
```

## Parameters

- `elasOut` (*`quaidsElasOut` structure*) - the result of
  [quaidsElasFit](quaidsElasFit.md).

## Returns

Nothing (prints to the console): income/own-price elasticity table,
uncompensated price elasticity matrix, and compensated price elasticity
matrix, each with standard errors.

## Remarks

Separated from [quaidsElasFit](quaidsElasFit.md) so elasticities can be
computed silently and printed only when wanted -- mirrors the
`quaidsFit()`/`printQuaids()` split. Verified byte-for-byte identical
printed output against the pre-split combined
[quaidsElas](quaidsElas.md) call.

## Examples

```gauss
struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call printQuaidsElas(elasOut);
```

## Source

`quaidselas.src`

## See Also

[quaidsElasFit](quaidsElasFit.md), [quaidsElas](quaidsElas.md)
