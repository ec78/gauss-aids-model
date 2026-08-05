# quaidsWorkflowScenarioFit

## Purpose

Runs the default applied workflow with [quaidsWorkflowFit](quaidsWorkflowFit.md)
and adds exact compensating/equivalent variation for one explicit price-change
scenario.

## Format

```gauss
struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    intcptPt, pricesPt0, pricesPt1, totexpPt0);
wfOut = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    intcptPt, pricesPt0, pricesPt1, totexpPt0, clusterId=householdId, weight=myWeight);   // Milestone 28
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl` - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `intcptPt` - intercept/control vector for the welfare scenario, including
  the leading constant.
- `pricesPt0` - base log-price vector.
- `pricesPt1` - new log-price vector.
- `totexpPt0` - base log total expenditure. Nominal expenditure is held fixed
  at this value when computing CV/EV.
- `clusterId` (*OPTIONAL keyword argument, default `0`*) - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `weight` (*OPTIONAL keyword argument, default `0`*) - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md). **Milestone 28 note**: this
  proc's required scenario arguments (`intcptPt`, `pricesPt0`, `pricesPt1`,
  `totexpPt0`) were moved ahead of `clusterId`/`weight` in the parameter
  list, since GAUSS requires every required (non-defaulted) parameter to
  precede any keyword-defaulted one -- `clusterId` sat earlier in the
  signature before this conversion.

## Returns

`wfOut` is the same `quaidsWorkflowOut` structure returned by
[quaidsWorkflowFit](quaidsWorkflowFit.md), including its compact preflight
diagnostic summary fields. Welfare scenario fields are filled when the base
fit converges:

- `welfareValid`
- `cv`, `ev`
- `seCV`, `seEV`
- `welfareRobustValid`, `seCVRobust`, `seEVRobust`
- `u0`, `u1`
- `scenarioIntcpt`, `scenarioPrices0`, `scenarioPrices1`,
  `scenarioTotexp0`

If the base fit does not converge, the workflow returns the core fit fields but
leaves `postValid`, `robustValid`, `postRobustValid`, `welfareValid`, and
`welfareRobustValid` equal to `0`.

## Remarks

This proc preserves `quaidsWorkflowFit()`'s shorter signature and makes the
welfare scenario explicit. The welfare fields are exactly the output of:

```gauss
wOut = quaidsWelfareFit(wfOut.bestB, wfOut.bestV,
    intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
```

When `wfOut.postRobustValid == 1`, robust/cluster-robust welfare standard
errors are also computed as:

```gauss
wOutR = quaidsWelfareFit(wfOut.robustBestB, wfOut.robustBestV,
    intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
```

## Example

```gauss
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);

pricesPt1 = wfOut.evalPrices;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);

if wfScenario.welfareValid;
    print wfScenario.cv wfScenario.seCV;
    print wfScenario.ev wfScenario.seEV;
endif;

if wfScenario.welfareRobustValid;
    print wfScenario.cv wfScenario.seCVRobust;
    print wfScenario.ev wfScenario.seEVRobust;
endif;
```

## Source

`quaidsworkflow.src`

## See Also

[quaidsWorkflowFit](quaidsWorkflowFit.md),
[quaidsWelfareFit](quaidsWelfareFit.md),
[quaidsRobustCovariance](quaidsRobustCovariance.md)
