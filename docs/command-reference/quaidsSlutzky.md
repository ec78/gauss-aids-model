# quaidsSlutzky

## Purpose

Tests negativity of the Slutzky matrix observation by observation, and
prints descriptive statistics on the eigenvalues.

## Format

```gauss
call quaidsSlutzky(b, intcpt, prices, totexp, aCtl);
```

## Parameters

- `b` (*matrix*) - parameter matrix; each column corresponds to a specific
  share; row order is intercept | prices | lx | lx2 | u. Typically
  `qOut.bestB`.
- `intcpt` (*Txnint matrix*) - independent variables, **including** the
  leading constant column -- typically `qOut.intcptFull`.
- `prices` (*Txn matrix*) - log prices, the full sample (not a single
  point).
- `totexp` (*Tx1 vector*) - log total expenditure, the full sample.
- `aCtl` (*`quaidsControl` structure*).

## Returns

Nothing (prints to the console): a table of eigenvalue descriptive
statistics (mean, std dev, min, Q1, Q2, Q3, max) computed observation by
observation, from smallest to biggest. Demand theory implies the Slutzky
matrix should be negative semi-definite; negative eigenvalues are
consistent with theory, positive eigenvalues indicate a theory violation
at that observation.

## Remarks

Unlike [quaidsElasFit](quaidsElasFit.md), which evaluates at one point,
`quaidsSlutzky()` always accepts an arbitrary sample (any number of rows)
-- this did not need to change during the Milestone 5 elasticities
generalization, since it was already general.

Curvature **imposition** (as opposed to this diagnosis-only check) is
explicitly out of scope -- even the R `micEconAids` reference implementation
used for cross-validation in this library only offers curvature diagnosis,
not imposition as a constrained-estimation mode.

## Examples

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);
```

## Source

`quaidsslutzky.src`

## See Also

[quaidsFit](quaidsFit.md), [quaidsElasFit](quaidsElasFit.md)
