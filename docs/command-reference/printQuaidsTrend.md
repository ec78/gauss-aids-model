# printQuaidsTrend

## Purpose

Prints the report captured in a `quaidsTrendOut` structure: the headline
joint trend test, a one-line plain-language read of it, and the
trend-slope coefficient table (so a user can see which specific
price/income responses show a hint of drift, not just whether one exists
anywhere in the system).

## Format

```gauss
call printQuaidsTrend(tOut);
```

## Parameters

- `tOut` (*`quaidsTrendOut` structure*) - the result of
  [quaidsTrendFit](quaidsTrendFit.md).

## Returns

Nothing (prints to the console): a reminder that this is a screening
diagnostic, not time-varying-parameter estimation; the joint
`chi2(trendDf)` test statistic and p-value; a one-line read ("some
evidence of coefficient drift" vs. "no strong evidence"); and a table of
`tOut.b1`'s trend-slope coefficients with standard errors, one row per
trend regressor (intercept-trend, then each relative-price trend, then
the expenditure trend) and one column per good.

## Remarks

Separated from [quaidsTrendFit](quaidsTrendFit.md) so the fit can be run
silently and printed only when wanted -- mirrors the
`quaidsFit()`/`printQuaids()` and `quaidsZeroFit()`/`printQuaidsZero()`
splits.

## Examples

```gauss
tOut = quaidsTrendFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsTrend(tOut);
```

## Source

`quaidstrend.src`

## See Also

[quaidsTrendFit](quaidsTrendFit.md)
