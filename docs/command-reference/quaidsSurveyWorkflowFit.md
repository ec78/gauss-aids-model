# quaidsSurveyWorkflowFit

## Purpose

Runs the applied workflow with `weight` fitting the estimator itself
(Milestone 26), and recomputes predicted shares and elasticities at a
sampling-weighted evaluation point.

## Format

```gauss
struct quaidsWorkflowOut wfOut;
wfOut = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, clusterId, weight);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl`, `clusterId` - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `weight` - `Tx1` nonnegative sampling-weight vector. The total weight must
  be positive. Zero weights are allowed. **Milestone 26: this weight now
  does double duty** -- it is forwarded into
  [quaidsWorkflowFit](quaidsWorkflowFit.md)'s own optional `weight`
  argument (so the estimator itself, not just the evaluation point, is
  fit under this weighting), *and* it still computes the weighted
  evaluation point exactly as before. This is a real, deliberate behavior
  change from the original (Milestone 25) version of this proc, which
  left the estimator unweighted -- see Remarks.

## Returns

Returns the same `quaidsWorkflowOut` structure as
[quaidsWorkflowFit](quaidsWorkflowFit.md), with these differences:

- `bestB`/`bestV` (and every other core-fit field) reflect the WEIGHTED
  estimator -- exactly what a direct
  `quaidsWorkflowFit(..., weight)` call would produce.
  `weighted`/`weightSum`/`effN` are also filled from that weighted fit.
- `evalIntcpt`, `evalPrices`, and `evalTotexp` are the sampling-weighted
  means of the intercept block, prices, and total expenditure.
- `shares`, `incomeElas`, `priceElas`, `compPriceElas`, and their classical
  delta-method standard errors are recomputed at that weighted evaluation
  point, using the weighted `bestB`/`bestV` above.
- If robust/cluster-robust post-estimation is available, the robust
  shares/elasticity standard errors are also recomputed at the weighted
  evaluation point.
- Survey metadata fields are filled: `surveyWeighted`,
  `surveyWeightValid`, `surveyWeightSum`, `surveyWeightNPositive`,
  `surveyWeightMin`, and `surveyWeightMax`.

## Remarks

**Milestone 26 behavior change**: earlier releases of this proc
(Milestone 25) left the underlying `quaidsFit()` estimator unweighted --
`weight` only ever affected the post-estimation evaluation point. Since
Milestone 26, `weight` is forwarded into
[quaidsWorkflowFit](quaidsWorkflowFit.md)'s own new optional `weight`
argument too, so the estimator itself is now genuinely fit under this
weighting. This closes a gap the Milestone 25 release explicitly flagged
as incomplete. If you have existing code relying on this proc leaving
`bestB` identical to the unweighted `quaidsWorkflowFit()` call, that
assumption no longer holds for a non-uniform `weight` -- a uniform weight
(e.g. `ones(nobs,1)`) still reproduces the unweighted fit exactly.

Use this when household- or person-level microdata should both fit a
sampling-weighted estimator and report population-representative
post-estimation summaries at the same weighted point. Full design-based
estimation (formal strata as a concept distinct from clustering,
replicate-weight/BRR/jackknife variance, finite-population correction) is
still a future roadmap item -- see `weight`'s own
[quaidsFit](quaidsFit.md) documentation for the current scope of what
"weighted" means here.

Invalid weights fail fast with a clear diagnostic.

## Examples

```gauss
struct quaidsWorkflowOut wfSurvey;
wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0, sampwt);

if wfSurvey.postValid;
    print wfSurvey.evalTotexp;
    print wfSurvey.shares;
    print wfSurvey.incomeElas;
endif;
```

## Source

`quaidssurvey.src`

## See Also

[quaidsWorkflowFit](quaidsWorkflowFit.md),
[quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md),
[quaidsSharesFit](quaidsSharesFit.md), [quaidsElasFit](quaidsElasFit.md),
[quaidsRobustCovariance](quaidsRobustCovariance.md)
