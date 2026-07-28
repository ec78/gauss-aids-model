# printQuaidsZero

## Purpose

Prints the first-stage probit summary (per-good zero-share fraction and
convergence) and the Shonkwiler-Yen-corrected coefficient table captured
in a `quaidsZeroOut` structure.

## Format

```gauss
call printQuaidsZero(zOut);
```

## Parameters

- `zOut` (*`quaidsZeroOut` structure*) - the result of
  [quaidsZeroFit](quaidsZeroFit.md).

## Returns

Nothing (prints to the console): the probit-stage summary table, then the
corrected coefficient table (intercept | gamma | beta | `[lambda]` | `u` |
`delta`, one column per good) with standard errors and convergence
diagnostics.

## Remarks

Separated from [quaidsZeroFit](quaidsZeroFit.md) so the fit can be run
silently and printed only when wanted -- mirrors the
`quaidsFit()`/`printQuaids()` and `quaidsCurvatureFit()`/
`printQuaidsCurvature()` splits.

## Examples

```gauss
struct quaidsZeroOut zOut;
zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);
```

## Source

`quaidszerocorrect.src`

## See Also

[quaidsZeroFit](quaidsZeroFit.md)
