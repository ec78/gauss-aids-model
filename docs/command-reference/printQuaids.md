# printQuaids

## Purpose

Prints the estimation-stage console report (instrumental-variables first
stage, iteration summary, homogeneity-constrained coefficient table,
overidentification test, symmetry test, and symmetry-constrained
coefficient table) from an already-computed `quaidsOut` structure.

## Format

```gauss
call printQuaids(qOut);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - the result of
  [quaidsFit](quaidsFit.md) or [quaidsFull](quaidsFull.md).

## Returns

Nothing (prints to the console).

## Remarks

Does **not** print elasticities, descriptive statistics, or the Slutzky
diagnostic -- those remain separate, explicitly-callable reports
([quaidsElas](quaidsElas.md)/[printQuaidsElas](printQuaidsElas.md),
[quaidsSlutzky](quaidsSlutzky.md)), reproduced automatically by the legacy
[quaids](quaids.md) wrapper.

A known, deliberately preserved anomaly: in the symmetry-constrained
table's "Residuals of instrumental regressions" row, point estimates come
from the homogeneity-stage fit (`qOut.homogB`) while standard errors/t/p
come from the symmetry-stage fit (`qOut.symcSE`/`qOut.symcT`/
`qOut.symcPvt`) -- both stages leave this coefficient block numerically
identical, but this looks like it could be an original-authoring
inconsistency. Carried over unchanged from the pre-refactor code rather
than silently corrected.

## Examples

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaids(qOut);
```

## Source

`quaids.src`

## See Also

[quaidsFit](quaidsFit.md), [quaids](quaids.md)
