# printQuaidsRobustBootstrap

## Purpose

Prints the coefficient table with the base closed-form robust/cluster SE
and the bootstrap SE side by side, plus the bootstrap run diagnostics,
from a `quaidsRobustBootOut` structure.

## Format

```gauss
call printQuaidsRobustBootstrap(rbOut);
```

## Parameters

- `rbOut` (*`quaidsRobustBootOut` structure*) - the result of
  [quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md).

## Returns

Nothing (prints to the console): cluster count and bootstrap run
bookkeeping, then the coefficient table with `closed` (the base
closed-form `seRobust`) and `boot` (`seBoot`) standard errors printed
side by side per coefficient.

## Remarks

Separated from [quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md) so
the fit can be run silently and printed only when wanted -- mirrors the
`quaidsCurvatureBootstrapFit()`/`printQuaidsCurvatureBootstrap()` split.

## Examples

```gauss
struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, clusterId=clusterId, seed=42);
call printQuaidsRobustBootstrap(rbOut);
```

## Source

`quaidsrobust.src`

## See Also

[quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md)
