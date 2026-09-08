# quaids

## Purpose

Estimates a linearized AIDS, iterated AIDS, or QUAIDS demand system
(backward-compatible wrapper): fits via `quaidsFit()`, prints the full
console report (estimation tables, elasticities at the mean/quartiles,
descriptive statistics, and the Slutzky negativity diagnostic), and
returns the same four positional matrices this proc has always returned.

## Format

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
```

## Parameters

Same as [quaidsFit](quaidsFit.md): `w`, `intcpt`, `prices`, `totexp`,
`instr`, `aCtl`.

## Returns

- `b1`, `v1` (*matrix*) - if `aCtl.homogenous == 1`,
  homogeneity-constrained estimates/covariance; if `aCtl.homogenous == 0`,
  the unconstrained (reparametrized) estimates/covariance.
- `b2`, `v2` (*matrix*) - if `aCtl.homogenous == 1`,
  homogeneity+symmetry-constrained estimates/covariance; if
  `aCtl.homogenous == 0`, `0`.

Identical to `qOut.b`, `qOut.v`, `qOut.bS`, `qOut.vS` from
[quaidsFit](quaidsFit.md) -- verified exactly (not approximately) equal by
`tests/quaids_schema_test.e`.

## Remarks

For a silent call that returns a struct without printing, use
[quaidsFit](quaidsFit.md) directly. `quaids()` is a backward-compatible
wrapper that fits, prints the full console report, and returns the same
four legacy matrices this proc has always returned.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.homogenous = 1;

{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
```

## Source

`quaids.src`

## See Also

[quaidsFit](quaidsFit.md), [printQuaids](printQuaids.md),
[quaidsElas](quaidsElas.md), [quaidsSlutzky](quaidsSlutzky.md)
