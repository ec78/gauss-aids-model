# quaidsWorkflowScenarioFit

## Purpose

Runs the default applied workflow with [quaidsWorkflowFit](quaidsWorkflowFit.md)
and adds exact compensating/equivalent variation for one explicit price-change
scenario.

## Format

```gauss
struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    clusterId, intcptPt, pricesPt0, pricesPt1, totexpPt0);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl`, `clusterId` - same as
  [quaidsWorkflowFit](quaidsWorkflowFit.md).
- `intcptPt` - intercept/control vector for the welfare scenario, including
  the leading constant.
- `pricesPt0` - base log-price vector.
- `pricesPt1` - new log-price vector.
- `totexpPt0` - base log total expenditure. Nominal expenditure is held fixed
  at this value when computing CV/EV.

## Returns

`wfOut` is the same `quaidsWorkflowOut` structure returned by
[quaidsWorkflowFit](quaidsWorkflowFit.md), with welfare scenario fields filled
when the base fit converges:

- `welfareValid`
- `cv`, `ev`
- `seCV`, `seEV`
- `u0`, `u1`
- `scenarioIntcpt`, `scenarioPrices0`, `scenarioPrices1`,
  `scenarioTotexp0`

If the base fit does not converge, the workflow returns the core fit fields but
leaves `postValid`, `robustValid`, and `welfareValid` equal to `0`.

## Remarks

This proc preserves `quaidsWorkflowFit()`'s shorter signature and makes the
welfare scenario explicit. The welfare fields are exactly the output of:

```gauss
wOut = quaidsWelfareFit(wfOut.bestB, wfOut.bestV,
    intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
```

## Example

```gauss
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);

pricesPt1 = wfOut.evalPrices;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    0, wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);

if wfScenario.welfareValid;
    print wfScenario.cv wfScenario.seCV;
    print wfScenario.ev wfScenario.seEV;
endif;
```

## Source

`quaidsworkflow.src`

## See Also

[quaidsWorkflowFit](quaidsWorkflowFit.md), [quaidsWelfareFit](quaidsWelfareFit.md)
