# quaidsSurveyWorkflowFit

## Purpose

Runs the applied workflow and recomputes predicted shares and elasticities at
a sampling-weighted evaluation point.

## Format

```gauss
struct quaidsWorkflowOut wfOut;
wfOut = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, clusterId, weight);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl`, `clusterId` - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `weight` - `Tx1` nonnegative sampling-weight vector. The total weight must
  be positive. Zero weights are allowed.

## Returns

Returns the same `quaidsWorkflowOut` structure as
[quaidsWorkflowFit](quaidsWorkflowFit.md), with these differences:

- `evalIntcpt`, `evalPrices`, and `evalTotexp` are the sampling-weighted
  means of the intercept block, prices, and total expenditure.
- `shares`, `incomeElas`, `priceElas`, `compPriceElas`, and their classical
  delta-method standard errors are recomputed at that weighted evaluation
  point.
- If robust/cluster-robust post-estimation is available, the robust
  shares/elasticity standard errors are also recomputed at the weighted
  evaluation point.
- Survey metadata fields are filled: `surveyWeighted`,
  `surveyWeightValid`, `surveyWeightSum`, `surveyWeightNPositive`,
  `surveyWeightMin`, and `surveyWeightMax`.

## Remarks

This is a workflow/post-estimation helper, not a survey-weighted estimator.
The underlying coefficient estimates still come from `quaidsFit()` through
`quaidsWorkflowFit()`, using the original unweighted moment conditions. The
sampling weights affect the representative point at which shares and
elasticities are evaluated.

Use this when household- or person-level microdata should report
population-representative post-estimation summaries, while keeping the
estimator itself unchanged. Full design-based estimation, strata, replicate
weights, and survey-weighted covariance formulas are future roadmap items.

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
