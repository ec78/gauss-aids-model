# quaidsTVPControlCreate

## Purpose

Returns a default `quaidsTVPControl` structure for
[quaidsTVPFit](quaidsTVPFit.md) -- mirrors
[quaidsControlCreate](quaidsControlCreate.md)'s own minimalism.

## Format

```gauss
tvpCtl = quaidsTVPControlCreate();
```

## Parameters

None.

## Returns

`tvpCtl` is a `quaidsTVPControl` structure:

- `smooth` - `1` by default: also run the RTS smoother (see
  [quaidsTVPFit](quaidsTVPFit.md)). Set to `0` to return only the
  filtered (real-time-causal) path, skipping the full-sample smoother.
- `othnam` - `0` by default: auto-generated `W1..Wn` good names. Set to
  a legacy character matrix of `n` names to override. Unlike
  `quaidsControl.othnam` (`string`-typed, defaults to `""`), this field is
  `matrix`-typed and defaults to `0` -- a `string`-typed field was found
  to reject a real character-matrix assignment (`error G0071 : Type
  mismatch`), so this field uses the same `matrix` type every output
  struct's own name field (`wnam`/`xnam`/etc.) already uses successfully.
  **A second gotcha**: a value built with `$|` (vertical string
  concatenation -- e.g. `"Food" $| "Housing" $| ...`) is still rejected
  even by this `matrix`-typed field; only the legacy character-matrix form
  `$+`/`ftocv()` produce is accepted. Coerce with a leading `0$+` first:
  `tvpCtl.othnam = 0$+someDollarPipeBuiltNames;`.
- `maxIters` - `1e5` by default (CMLMT's own default): the iteration
  cap for the `Q` maximum-likelihood step, **per CMLMT pass** (see
  Remarks). Must be a finite positive integer scalar. *(Added in 0.4.0.)*
- `maxTime` - `1e8` by default (CMLMT's own default): the wall-clock cap
  for the `Q` maximum-likelihood step, in **minutes**, **per CMLMT pass**
  (see Remarks). Must be a finite positive scalar. *(Added in 0.4.0.)*

## Remarks

Deliberately minimal: `H` and `q0` (the observation covariance and
starting state-innovation-covariance diagonal) are NOT control-struct
fields -- they are required, data-dependent arguments to
[quaidsTVPFit](quaidsTVPFit.md) itself, not tunable defaults.

**Optimizer caps.** The `maxIters`/`maxTime` defaults are CMLMT's own,
so a fit made with the defaults is identical to one made before these
fields existed -- and, like before, effectively unbounded. Set finite
values to bound run time. A converged fit does not depend on a cap it
never reached.

**Both caps apply per CMLMT pass, not per fit.** `sslib`'s `ssFitTVP()`
may run two CMLMT passes with the same control: the estimation pass and,
because `Q` is positivity-transformed, a second, covariance-only pass
started from the first pass's result. The worst-case run time is
therefore roughly **twice** each cap, plus up to one iteration of
overrun per pass (`maxTime` is checked once per CMLMT iteration). So
`maxIters = 500, maxTime = 5` bounds the MLE at about 1000 iterations or
10 minutes, not 500 or 5.

An exhausted cap is not an error. `tvOut.mleRetcode` reports the
**first (estimation) pass**: [quaidsTVPFit](quaidsTVPFit.md) returns
normally with `tvOut.mleRetcode /= 0` (`2` for the iteration cap, `11`
for the time cap), and [printQuaidsTVP](printQuaidsTVP.md) flags it. A
cap hit in the second pass is **not exposed**: `ssFitTVP()` discards
that pass's return code upstream, and `quaidsTVPFit()` does not report
the second pass's covariance. Invalid values (zero, negative,
fractional `maxIters`, non-scalar, missing, NaN or infinite) stop
`quaidsTVPFit()` with a clear error before any estimation. A struct
declared without this proc (`struct quaidsTVPControl tvpCtl;`) leaves
both fields empty, and `quaidsTVPFit()` then uses the defaults.

## Examples

```gauss
tvpCtl = quaidsTVPControlCreate();
tvpCtl.smooth = 0;     // filtered path only, skip the RTS smoother
tvpCtl.maxIters = 500; // cap the Q MLE at 500 iterations ...
tvpCtl.maxTime = 5;    // ... or 5 minutes, per CMLMT pass (worst case ~2x)

tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
if tvOut.mleRetcode /= 0;
    print "Q MLE did not converge within the caps (first pass)";
endif;
```

## Source

`quaidstvpfit.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md), [quaidsControlCreate](quaidsControlCreate.md)
