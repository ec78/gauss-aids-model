# printQuaidsReplicateWeight

## Purpose

Prints the coefficient table with replicate-weight standard errors and the
replicate run diagnostics (requested/completed/failed) from a
`quaidsReplicateOut` structure.

## Format

```gauss
call printQuaidsReplicateWeight(rOut);
```

## Parameters

- `rOut` (*`quaidsReplicateOut` structure*) - from
  [quaidsReplicateWeightFit](quaidsReplicateWeightFit.md).

## Returns

None -- prints to the console only.

## Remarks

Separated from [quaidsReplicateWeightFit](quaidsReplicateWeightFit.md) so
the fit itself stays silent, matching this library's own silent-Fit-proc
convention. The report header includes `rOut.method` verbatim (e.g.
`"REPLICATE-WEIGHT (JK1) STANDARD ERRORS"`), so the printed report always
discloses which design the caller declared -- purely for readability;
`method` has no effect on the underlying computation.

## Examples

```gauss
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl,
    replicateWeights, scaleFactorJK1, weight=surveyWeight, method="JK1");

call printQuaidsReplicateWeight(rOut);
```

## Source

`quaidsreplicate.src`

## See Also

[quaidsReplicateWeightFit](quaidsReplicateWeightFit.md)
