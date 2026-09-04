# quaidsGetHomogeneity

## Purpose

Correctly spelled public getter for the historical `aCtl.homogenous`
control field, matching [quaidsSetHomogeneity](quaidsSetHomogeneity.md).
Lets application code read the current homogeneity setting back without
referencing the misspelled field name directly.

## Format

```gauss
homogeneous = quaidsGetHomogeneity(aCtl);
```

## Parameters

- `aCtl` (*`quaidsControl` structure*).

## Returns

`homogeneous` (*scalar*) - the current value of `aCtl.homogenous` (`1` =
homogeneity imposed, `0` = unconstrained).

## Remarks

A thin, read-only accessor -- `retp(aCtl.homogenous)` -- kept as its own
procedure (rather than expecting callers to read the field directly) for
symmetry with `quaidsSetHomogeneity()` and so application code has no
remaining reason to spell out the misspelled field name at all.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl = quaidsSetHomogeneity(aCtl, 0);

if quaidsGetHomogeneity(aCtl) == 0;
    qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
    { statH, pvalH, dfH } = quaidsHomogeneityTest(qOut);
endif;
```

## Source

`quaidsutil.src`

## See Also

[quaidsSetHomogeneity](quaidsSetHomogeneity.md), [quaidsControlCreate](quaidsControlCreate.md)
