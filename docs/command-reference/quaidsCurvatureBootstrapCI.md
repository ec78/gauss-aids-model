# quaidsCurvatureBootstrapCI

## Purpose

Percentile bootstrap confidence intervals for a
[quaidsCurvatureFit](quaidsCurvatureFit.md) coefficient vector, computed
directly from an already-computed
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) result's raw
draws -- no new resampling or refitting.

## Format

```gauss
{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, alpha);
```

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

**A real bug found while building this, now fixed everywhere it
occurred**: GAUSS's `reshape()` fills row-major, not column-major like
`vec()` -- confirmed empirically (`reshape(vec(X), rows(X), cols(X))`
does **not** recover `X` in general). Building this proc's own ground-
truth cross-check surfaced that `quaidsCurvatureFit()`'s `cOut.se` and
`quaidsCurvatureBootstrapFit()`'s `bootOut.seBoot` had exactly this bug
since Milestone 10/15 respectively -- their individual cells were
silently scrambled relative to `cOut.b`/`bootOut.b`, invisible to the
shape/sign/finiteness checks already in place (those properties are
permutation-invariant). Both were fixed alongside this proc's own
(previously unshipped) correct implementation. See
`GOLD_STANDARD_TODO.md`'s Milestone 18 section.

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
