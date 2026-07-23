# printQuaidsCurvature

## Purpose

Prints the curvature-constrained coefficient table and the Slutzky-matrix
eigenvalues at the reference point captured in a `quaidsCurvOut`
structure.

## Format

```gauss
call printQuaidsCurvature(cOut);
```

## Parameters

- `cOut` (*`quaidsCurvOut` structure*) - the result of
  [quaidsCurvatureFit](quaidsCurvatureFit.md).

## Returns

Nothing (prints to the console): the coefficient table (intercept | gamma
| beta, one column per good) with standard errors, convergence
diagnostics, and the Slutzky-matrix eigenvalues at the reference point.

## Remarks

Separated from [quaidsCurvatureFit](quaidsCurvatureFit.md) so the fit can
be run silently and printed only when wanted -- mirrors the
`quaidsFit()`/`printQuaids()` and `quaidsElasFit()`/`printQuaidsElas()`
splits.

## Examples

```gauss
struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call printQuaidsCurvature(cOut);
```

## Source

`quaidscurvature.src`

## See Also

[quaidsCurvatureFit](quaidsCurvatureFit.md)
