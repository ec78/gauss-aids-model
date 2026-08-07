# printQuaidsRobust

## Purpose

Prints the reduced-form coefficient table with robust/cluster-robust
standard errors captured in a `quaidsRobustOut` structure.

## Format

```gauss
call printQuaidsRobust(rOut);
```

## Parameters

- `rOut` (*`quaidsRobustOut` structure*) - the result of
  [quaidsRobustFit](quaidsRobustFit.md).

## Returns

Nothing (prints to the console): the number of observations/clusters,
then the reduced-form coefficient table (one column per independently-
estimated good) with standard errors.

## Remarks

Separated from [quaidsRobustFit](quaidsRobustFit.md) so the fit can be
run silently and printed only when wanted -- mirrors the
`quaidsFit()`/`printQuaids()` and `quaidsCurvatureFit()`/
`printQuaidsCurvature()` splits.

## Examples

```gauss
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, clusterId);
call printQuaidsRobust(rOut);
```

## Source

`quaidsrobust.src`

## See Also

[quaidsRobustFit](quaidsRobustFit.md)
