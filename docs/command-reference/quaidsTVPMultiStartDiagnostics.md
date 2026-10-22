# quaidsTVPMultiStartDiagnostics

## Purpose

Runs three complete TVP `Q`-MLE fits from a fixed `0.1x`, `1x`, and `10x`
starting grid, then compares their common-scale fixed-`Q` filter likelihoods.
This opt-in diagnostic reports observed start stability and start-relative
scale caution; it does not claim global optimality.
It is part of the optional TVP module and requires the same `sslib` v2.1.2
setup and include chain documented by [quaidsTVPFit](quaidsTVPFit.md).

## Format

```gauss
msOut = quaidsTVPMultiStartDiagnostics(w, prices, totexp, H, q0, tvpCtl,
    scaleRatio=1e-3);
```

## Parameters

- `w`, `prices`, `totexp`, `H`, `q0`, `tvpCtl` - the same inputs as
  [quaidsTVPFit](quaidsTVPFit.md). `q0` must be a finite, strictly positive
  `k_states x 1` vector.
- `scaleRatio` (*scalar*, default `1e-3`) - strictly positive threshold below
  one. Each arm compares its fitted `Q` with that arm's scaled start.

Only `tvpCtl.maxIters` and `tvpCtl.maxTime` affect these fits. State paths and
printed reports are not returned.

## Returns

`msOut` is a `quaidsTVPMultiStartOut` structure. Its three-row fields are in
fixed `0.1x`, `1x`, `10x` order:

- `status`, `statusLabel` - `0`/`START-STABLE-OBSERVED`, `1`/`FIT-FAILED`,
  `2`/`INVALID-Q`, `3`/`INVALID-LOGLIKELIHOOD`, `4`/`START-SENSITIVE`, or
  `5`/`SCALE-CAUTION`. Precedence is `FIT-FAILED`, `INVALID-Q`,
  `INVALID-LOGLIKELIHOOD`, `START-SENSITIVE`, `SCALE-CAUTION`, then
  `START-STABLE-OBSERVED`.
- `multipliers`, `retcodes`, `fitOk`, `qValid`, `llValid`, `scaleCaution`,
  `loglikelihoods`, `nScaleCollapse` - per-arm results.
- `qfits` - `k_states x 3` fitted `Q` matrix, one column per arm.
- `scaleRatio`, `anyScaleCaution` - applied threshold and aggregate caution.
- `bestStartIndex`, `bestMultiplier`, `bestLL`, `maxLLGap` - best observed
  valid arm and maximum likelihood gap. Exact ties select the first arm.
- `tolerance` - fixed likelihood tolerance, `1e-6*nobs`.
- `startStable` - one only when all fits converge, all `Q` and likelihood
  results are valid, and `maxLLGap <= tolerance`.

An invalid-Q arm leaves its `qfits` column, likelihood, collapse count, and
scale-caution field missing. An invalid returned likelihood is normalized to
missing rather than exposed as a raw nonfinite value. `bestStartIndex`,
`bestMultiplier`, `bestLL`, and
`maxLLGap` are missing unless every arm has valid `Q` and likelihood output.
`anyScaleCaution` is one if any valid arm cautions, zero only when all three
`Q` vectors are valid and none cautions, and otherwise missing. `startStable`
is always a populated zero/one scalar.

## Remarks

This call performs three sequential complete MLE fits. Its worst-case optimizer
budget is therefore approximately three times the per-fit caps, plus the
documented possible one-iteration overruns. Set finite caps before calling it;
do not add it casually to routine smoke tests or interactive workflows.

Likelihood agreement does not require fitted `Q` vectors to agree. Materially
different `Q` values can lie along a flat likelihood ridge. Consequently,
`START-STABLE-OBSERVED` means only that these three starts agreed within the
documented tolerance. It is not proof of global optimality, identification, or
valid inference. `SCALE-CAUTION` may coexist with `startStable == 1`: the
likelihoods can agree while fitted process variances shrink relative to their
own starts.

## Examples

```gauss
tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 5;

msOut = quaidsTVPMultiStartDiagnostics(w, prices, totexp, H, q0, tvpCtl);
print msOut.statusLabel;
print msOut.multipliers~msOut.retcodes~msOut.loglikelihoods;

if msOut.startStable;
    print "Observed likelihood agreement across the fixed start grid.";
endif;
if not scalmiss(msOut.anyScaleCaution);
    if msOut.anyScaleCaution;
        print "At least one valid arm has start-relative Q scale caution.";
    endif;
endif;
```

## Source

`quaidstvpfit.src`

## See Also

[quaidsTVPDiagnostics](quaidsTVPDiagnostics.md),
[quaidsTVPFit](quaidsTVPFit.md),
[quaidsTVPControlCreate](quaidsTVPControlCreate.md)
