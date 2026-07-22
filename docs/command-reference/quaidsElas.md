# quaidsElas

## Purpose

Computes income and price elasticities, and their asymptotic standard
errors, at one evaluation point, and prints them (backward-compatible
wrapper: calls `quaidsElasFit()` then `printQuaidsElas()`).

## Format

```gauss
call quaidsElas(b, v, intcpt, prices, totexp, aCtl);
```

## Parameters

Same as [quaidsElasFit](quaidsElasFit.md): `b`, `v`, `intcpt`, `prices`,
`totexp`, `aCtl`.

## Returns

Nothing (prints to the console). For a silent, struct-returning call, use
[quaidsElasFit](quaidsElasFit.md) directly.

## Remarks

Unchanged signature and printed output from before the Milestone 5
elasticities-generalization split; now a thin wrapper around
[quaidsElasFit](quaidsElasFit.md) and
[printQuaidsElas](printQuaidsElas.md). The legacy [quaids](quaids.md)
wrapper calls this at four fixed points (mean, Q1, median, Q3) as part of
its full console report.

## Examples

```gauss
call quaidsElas(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
```

## Source

`quaidselas.src`

## See Also

[quaidsElasFit](quaidsElasFit.md), [printQuaidsElas](printQuaidsElas.md),
[quaids](quaids.md)
