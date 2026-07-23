# printQuaidsWelfare

## Purpose

Prints the compensating/equivalent variation estimates and standard
errors captured in a `quaidsWelfareOut` structure.

## Format

```gauss
call printQuaidsWelfare(wOut);
```

## Parameters

- `wOut` (*`quaidsWelfareOut` structure*) - the result of
  [quaidsWelfareFit](quaidsWelfareFit.md).

## Returns

Nothing (prints to the console): CV and EV with standard errors, and the
implied utility levels at the base point and at the new prices with
unchanged nominal expenditure.

## Remarks

Separated from [quaidsWelfareFit](quaidsWelfareFit.md) so welfare
measures can be computed silently and printed only when wanted -- mirrors
the `quaidsElasFit()`/`printQuaidsElas()` and
`quaidsCurvatureFit()`/`printQuaidsCurvature()` splits.

## Examples

```gauss
struct quaidsWelfareOut wOut;
wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
call printQuaidsWelfare(wOut);
```

## Source

`quaidswelfare.src`

## See Also

[quaidsWelfareFit](quaidsWelfareFit.md)
