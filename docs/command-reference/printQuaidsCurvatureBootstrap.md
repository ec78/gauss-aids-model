# printQuaidsCurvatureBootstrap

## Purpose

Prints the curvature-constrained coefficient table with the delta-method
and bootstrap standard errors side by side, plus the bootstrap run
diagnostics, captured in a `quaidsCurvBootOut` structure.

## Format

```gauss
call printQuaidsCurvatureBootstrap(bootOut);
```

## Parameters

- `bootOut` (*`quaidsCurvBootOut` structure*) - the result of
  [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md).

## Returns

Nothing (prints to the console): the coefficient table (intercept | gamma
| beta, plus a trailing lambda row for QUAIDS) with the delta-method SE
and the bootstrap SE printed as separate labeled rows underneath each
coefficient row, plus the requested/completed/failed/attempts/seed
bootstrap diagnostics.

## Remarks

Separated from [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)
so the bootstrap can run silently and be printed only when wanted --
mirrors every other `...Fit()`/`print...()` split in this codebase. Both
standard errors are shown together deliberately: the existing delta-method
SE is never replaced, only supplemented, so a caller can see directly
where the two disagree (typically wherever the estimated Cholesky factor
sits at or near the boundary of the negative-semidefinite cone).

## Examples

```gauss
struct quaidsCurvBootOut bootOut;
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, 42);
call printQuaidsCurvatureBootstrap(bootOut);
```

## Source

`quaidscurvature.src`

## See Also

[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md),
[quaidsCurvatureFit](quaidsCurvatureFit.md)
