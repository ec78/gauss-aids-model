# quaidsTVPDiagnostics

## Purpose

Summarizes one completed [quaidsTVPFit](quaidsTVPFit.md) result without
refitting the model. It checks optimizer success, fitted/starting `Q` validity,
and whether any fitted `Q` element is small relative to the caller's own start.
It is part of the optional TVP module and uses the same `sslib` v2.1.2 setup
and include chain documented by [quaidsTVPFit](quaidsTVPFit.md).

## Format

```gauss
dOut = quaidsTVPDiagnostics(tvOut, scaleRatio=1e-3);
```

## Parameters

- `tvOut` (*`quaidsTVPOut` structure*) - result returned by
  [quaidsTVPFit](quaidsTVPFit.md).
- `scaleRatio` (*scalar*, default `1e-3`) - strictly positive threshold below
  one. An element cautions when `Qfit/q0 < scaleRatio`; equality does not.

## Returns

`dOut` is a `quaidsTVPDiagOut` structure:

- `status`, `statusLabel` - `0`/`OK`, `1`/`FIT-FAILED`,
  `2`/`INVALID-Q`, or `5`/`SCALE-CAUTION`. Precedence is `FIT-FAILED`,
  `INVALID-Q`, `SCALE-CAUTION`, then `OK`.
- `mleRetcode`, `fitOk` - optimizer return code and its `retcode == 0` flag.
- `qValid` - one only when `q0` and `Qfit` are conformable, finite, strictly
  positive column vectors.
- `scaleRatio`, `scaleCaution`, `nScaleCollapse` - applied threshold, any-hit
  flag, and number of fitted variances below it.
- `minQOverQ0`, `maxQOverQ0` - range of the elementwise fitted/start ratios.
  When `Q` is invalid, `scaleCaution`, `nScaleCollapse`, `minQOverQ0`, and
  `maxQOverQ0` are missing; the validated `scaleRatio` remains populated.

## Remarks

This is the recommended first diagnostic after every TVP fit. It performs no
optimization or Kalman-filter call, so it is cheap enough for normal examples
and workflows. `FIT-FAILED` remains the top-level status when the optimizer did
not converge, even if another field also shows invalid `Q`.

`SCALE-CAUTION` is not a statistical boundary test. It says only that one or
more fitted process variances became small relative to the caller-chosen `q0`.
It does not prove that a variance is zero, that the model is identified, or
that the fitted point is globally optimal.

## Examples

```gauss
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
dOut = quaidsTVPDiagnostics(tvOut);

print dOut.statusLabel;
if not dOut.fitOk;
    print "TVP optimizer did not converge; do not treat this as a fitted result.";
elseif not dOut.qValid;
    print "The fitted or starting Q vector is invalid.";
elseif dOut.scaleCaution;
    print "One or more fitted Q values are small relative to q0.";
endif;
```

## Source

`quaidstvpfit.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md),
[quaidsTVPMultiStartDiagnostics](quaidsTVPMultiStartDiagnostics.md),
[printQuaidsTVP](printQuaidsTVP.md)
