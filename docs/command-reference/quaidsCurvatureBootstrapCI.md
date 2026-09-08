# quaidsCurvatureBootstrapCI

## Purpose

Percentile bootstrap confidence intervals for a
[quaidsCurvatureFit](quaidsCurvatureFit.md) coefficient vector, computed
directly from an already-computed
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) result's raw
draws -- no new resampling or refitting.

## Format

```gauss
library optmt, quaids;
#include quaidscurvature.src

{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, alpha);
```

Not loaded by `library quaids;` alone -- see
[quaidsCurvatureFit](quaidsCurvatureFit.md)'s Format section and
[`docs/public-api.json`](../public-api.json)'s `optional_modules` entry.

## Parameters

- `bootOut` (*`quaidsCurvBootOut` structure*) - the result of
  [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md).
- `alpha` (*scalar*) - a value in `(0, 1)`, e.g. `0.05` for a 95%
  interval. **Required -- no default**, matching
  `quaidsCurvatureBootstrapFit()`'s own convention of never silently
  guessing an inference-affecting parameter on the caller's behalf.

## Returns

- `ciLower`, `ciUpper` - matrices the same shape as `bootOut.b`: the
  `alpha/2` and `1 - alpha/2` empirical quantiles of `bootOut.bBoot`'s
  columns, reshaped back to `bootOut.b`'s row/column layout.

## Remarks

**Percentile CIs from a small `B` are necessarily crude.** This library's
own test fixtures use `B=15`/`B=5` to bound runtime (see
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
timing notes) -- a caller who chose a small `B` there should expect
correspondingly wide, noisy intervals here, not a false sense of
precision. No minimum-`B` error is enforced; the caller is trusted to
weigh this tradeoff, same as the `B` choice itself.

## Examples

```gauss
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, 42);

{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, 0.05);
print "95% CI for the first coefficient, good 1:" ciLower[1,1] ciUpper[1,1];
```

## Source

`quaidscurvature.src`

## See Also

[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md),
[printQuaidsCurvatureBootstrap](printQuaidsCurvatureBootstrap.md),
[quaidsCurvatureFit](quaidsCurvatureFit.md)
