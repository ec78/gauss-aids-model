# quaidsWorkflowFit

## Purpose

Runs the default applied workflow in one silent call: fit the demand system
with [quaidsFit](quaidsFit.md), evaluate predicted shares and elasticities at
the sample mean, and compute robust or cluster-robust standard errors when the
base fit converges.

## Format

```gauss
struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, clusterId);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl` - same as
  [quaidsFit](quaidsFit.md).
- `clusterId` - same as [quaidsRobustFit](quaidsRobustFit.md): scalar `0`
  for heteroskedasticity-robust standard errors, or a `Tx1` cluster-label
  vector.

## Returns

`wfOut` is a `quaidsWorkflowOut` structure containing:

- Core fit metadata and coefficient blocks copied from `quaidsOut`:
  `model`, `converged`, `iterations`, `finalErr`, `b`, `v`, `bS`, `vS`,
  `bestB`, `bestV`, `symValid`, and `overidValid`.
- Restriction/model-choice summaries:
  `symStat`, `symPval`, `symDf`, `symReject05`, `overidFstat`,
  `overidPvf`, `overidDf`, `overidReject05`, `quadraticValid`,
  `quadraticStat`, `quadraticPval`, `quadraticDf`, `quadraticReject05`.
  The quadratic summary is filled only for unconstrained QUAIDS fits
  (`aCtl.linear=0`, `aCtl.homogenous=0`), matching
  [quaidsQuadraticTest](quaidsQuadraticTest.md)'s own contract.
- Evaluation point fields: `evalIntcpt`, `evalPrices`, `evalTotexp`, all set
  to the sample mean used for post-estimation.
- Mean-point predicted shares: `shares`, `sharesSE`, `sharesV`.
- Mean-point elasticities: `incomeElas`, `incomeElasSE`, `priceElas`,
  `priceElasSE`, `compPriceElas`, `compPriceElasSE`.
- Robust/cluster-robust outputs: `robustB`, `robustV`, `robustSE`,
  `robustNClusters`.
- Welfare scenario outputs are present in the struct but left missing by
  `quaidsWorkflowFit()`: `welfareValid`, `cv`, `ev`, `seCV`, `seEV`,
  `u0`, `u1`, `scenarioIntcpt`, `scenarioPrices0`, `scenarioPrices1`,
  `scenarioTotexp0`. Use
  [quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md) to fill them.
- Flags: `postValid` and `robustValid`, equal to `1` only when the base
  `quaidsFit()` converged and the post-estimation products were computed.

## Remarks

This is a workflow convenience layer, not a new estimator. It deliberately
composes existing public procs so the numerical results match the explicit
manual sequence:

1. `qOut = quaidsFit(...)`
2. `sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, mean point, aCtl)`
3. `elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, mean point, aCtl)`
4. `rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, clusterId)`

If `qOut.converged /= 1`, the workflow returns the core fit fields but leaves
post-estimation fields missing and sets `postValid = robustValid = 0`.

For a one-call workflow that also computes CV/EV for a price-change scenario,
use [quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md).

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);

if wfOut.postValid;
    print wfOut.shares;
    print wfOut.incomeElas;
endif;
```

## Source

`quaidsworkflow.src`

## See Also

[quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md),
[quaidsFit](quaidsFit.md), [quaidsSharesFit](quaidsSharesFit.md),
[quaidsElasFit](quaidsElasFit.md), [quaidsRobustFit](quaidsRobustFit.md)
