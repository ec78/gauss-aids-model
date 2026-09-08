# quaidsSurveyWorkflowFit

## Purpose

Runs the applied workflow with `weight` fitting the estimator itself, and
recomputes predicted shares and elasticities at a sampling-weighted
evaluation point.

## Format

```gauss
wfOut = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, weight);
wfOut = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, weight, clusterId=householdId);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl` - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `weight` (*Tx1 vector, required*) - nonnegative sampling-weight vector.
  The total weight must be positive. Zero weights are allowed. Required
  with no default -- this is the entire purpose of this proc, so it is
  not keyword-callable, matching this project's "never silently guess an
  inference-affecting parameter" precedent (`B` in the bootstrap procs,
  `replicateWeights`/`scaleFactor` in
  [quaidsReplicateWeightFit](quaidsReplicateWeightFit.md)). This weight
  does double duty: it is forwarded into
  [quaidsWorkflowFit](quaidsWorkflowFit.md)'s own optional `weight`
  argument (so the estimator itself, not just the evaluation point, is
  fit under this weighting), *and* it also computes the weighted
  evaluation point -- see Remarks.
- `clusterId` (*OPTIONAL keyword argument, default `0`*) - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).

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

`weight` fits the underlying `quaidsFit()` estimator itself (forwarded
into [quaidsWorkflowFit](quaidsWorkflowFit.md)'s own optional `weight`
argument), not just the post-estimation evaluation point -- `bestB` is
not identical to an unweighted `quaidsWorkflowFit()` call for a
non-uniform `weight`. A uniform weight (e.g. `ones(nobs,1)`) reproduces
the unweighted fit exactly.

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
wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, sampwt);

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
