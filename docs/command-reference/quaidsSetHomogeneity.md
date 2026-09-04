# quaidsSetHomogeneity

## Purpose

Correctly spelled public setter for the historical `aCtl.homogenous`
control field (note the missing second "e"). New application code should
configure homogeneity imposition through this procedure instead of
assigning the misspelled field directly, so it does not depend on a name
that may eventually be removed. See the
[compatibility policy](../../README.md#compatibility-policy) and
[`docs/public-api.json`](../public-api.json).

## Format

```gauss
aCtl = quaidsSetHomogeneity(aCtl, homogeneous);
```

## Parameters

- `aCtl` (*`quaidsControl` structure*) - a control struct, typically from
  `quaidsControlCreate()`.
- `homogeneous` (*scalar*) - `1` to impose homogeneity (and test/report
  symmetry); `0` to leave the fit unconstrained.

## Returns

`aCtl` - the same `quaidsControl` structure with its (internally still
misspelled) `homogenous` field set to `homogeneous`.

## Remarks

Errors clearly (`errorlog` + `end`) if `homogeneous` is not a scalar `0`
or `1` -- including non-finite/missing input -- rather than silently
storing an invalid control value. This mirrors this library's existing
fail-fast conventions elsewhere (e.g. `quaidsRobustFit`'s `clusterId`
validation).

Internally, this only ever writes `aCtl.homogenous` (the field GAUSS
struct compatibility requires this library to retain through the `0.x`
series) -- there is no separate, independently-tracked "correct" field
on `quaidsControl` itself. Use [quaidsGetHomogeneity](quaidsGetHomogeneity.md)
to read the value back without referencing the misspelled field name
either.

This setter only affects `quaidsControl` (the estimator's *input*
options). The *returned* `quaidsOut`/`quaidsZeroOut`/`quaidsWorkflowOut`
structs separately expose a correctly spelled `homogeneous` field
alongside the deprecated `homogenous` alias -- see those structs' own
documentation.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl = quaidsSetHomogeneity(aCtl, 1);   // impose homogeneity, test/report symmetry

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
print qOut.homogeneous;   // 1 -- correctly spelled
print qOut.homogenous;    // 1 -- deprecated alias, same value
```

## Source

`quaidsutil.src`

## See Also

[quaidsGetHomogeneity](quaidsGetHomogeneity.md), [quaidsControlCreate](quaidsControlCreate.md),
[quaidsFit](quaidsFit.md)
