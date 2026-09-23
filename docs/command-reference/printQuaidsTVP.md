# printQuaidsTVP

## Purpose

Print the report captured in a `quaidsTVPOut` structure (from
[quaidsTVPFit](quaidsTVPFit.md)): MLE convergence, the fitted `Q` diagonal
grouped by the state vector's own alpha/gamma/beta blocks, and the
final-period coefficient snapshot.

## Format

```gauss
call printQuaidsTVP(tvOut);
```

## Parameters

- `tvOut` (*`quaidsTVPOut` structure*) - from [quaidsTVPFit](quaidsTVPFit.md).

## Returns

None (prints to the screen/log).

## Remarks

The final-period coefficient table shows **point values only** -- no
standard errors are available for a TVP-AIDS fit yet (see
[quaidsTVPFit](quaidsTVPFit.md)'s own Remarks). The table's row labels
follow the state vector's own block layout (see `quaidstvp.src`'s header):
one `alpha` row, `n` `gamma` rows (the full, symmetric, homogeneity-
respecting absolute-price gamma sub-block), and one `beta` row.

The report distinguishes whether the printed final-period snapshot came
from the smoothed or the filtered path (`tvOut.smoothed`), since a caller
may have set `tvpCtl.smooth = 0` to skip the RTS smoother.

## Examples

```gauss
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
call printQuaidsTVP(tvOut);
```

## Source

`quaidstvpfit.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md), [printQuaidsTrend](printQuaidsTrend.md)
(the closest existing precedent this printer's layout follows)
