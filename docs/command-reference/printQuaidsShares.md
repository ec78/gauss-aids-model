# printQuaidsShares

## Purpose

Prints the predicted budget shares and standard errors captured in a
`quaidsSharesOut` structure.

## Format

```gauss
call printQuaidsShares(sharesOut);
```

## Parameters

- `sharesOut` (*`quaidsSharesOut` structure*) - the result of
  [quaidsSharesFit](quaidsSharesFit.md).

## Returns

Nothing (prints to the console): a table of good name, predicted share,
and standard error, one row per good.

## Remarks

Separated from [quaidsSharesFit](quaidsSharesFit.md) so the computation
can run silently and be printed only when wanted -- mirrors every other
`...Fit()`/`print...()` split in this codebase.

## Examples

```gauss
struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call printQuaidsShares(sharesOut);
```

## Source

`quaidsshares.src`

## See Also

[quaidsSharesFit](quaidsSharesFit.md)
