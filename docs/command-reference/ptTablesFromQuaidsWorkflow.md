# ptTablesFromQuaidsWorkflow

## Purpose

Builds a `pubtable` table bundle from a `quaidsWorkflowOut` object. Optional
-- requires the `pubtable` package and `src/pubtable_quaids.src`.

## Format

```gauss
struct ptTable tbls;
tbls = ptTablesFromQuaidsWorkflow(wfOut);
```

## Parameters

- `wfOut` - a `quaidsWorkflowOut` returned by
  [quaidsWorkflowFit](quaidsWorkflowFit.md) or
  [quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md), with
  `wfOut.postValid == 1`.

## Returns

`tbls` is a `struct ptTable` array:

1. Predicted budget shares.
2. Income elasticities.
3. Uncompensated (Marshallian) price elasticities.
4. Compensated (Hicksian) price elasticities.
5. Welfare scenario table, only when `wfOut.welfareValid == 1`.

## Remarks

This adapter reports the workflow's applied outputs. For coefficient tables,
use [ptFromQuaids](ptFromQuaids.md) on a `quaidsOut` fit.

## Example

```gauss
library pubtable, quaids;
#include quaids.sdf
#include pubtable_quaids.src

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);
wfTbls = ptTablesFromQuaidsWorkflow(wfOut);

call ptExport(wfTbls[1], "shares.md");
call ptExport(wfTbls[2], "income_elasticities.md");
```

## Source

`pubtable_quaids.src`

## See Also

[quaidsWorkflowFit](quaidsWorkflowFit.md),
[quaidsWorkflowScenarioFit](quaidsWorkflowScenarioFit.md),
[ptTablesFromQuaidsElas](ptTablesFromQuaidsElas.md)
