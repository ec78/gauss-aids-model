# quaidsPreflight

## Purpose

Runs an estimator-free preflight diagnostic pass on the raw QUAIDS/AIDS
inputs. It checks dimensions, finite values, share adding-up, zero/negative
shares, price/expenditure/instrument variation, first-stage IV strength,
basic design invertibility, cluster counts, and a simple convergence-risk
flag.

## Format

```gauss
pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl);
pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, clusterId=householdId, weight=myWeight);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr`, `aCtl` - same shapes and
  meaning as [quaidsFit](quaidsFit.md).
- `clusterId` (*OPTIONAL keyword argument, default `0`*) - scalar `0` for
  the heteroskedasticity-robust case, or a `Tx1` cluster-label vector.
- `weight` (*OPTIONAL keyword argument, default `0`*) - scalar `0`, or a
  `Tx1` vector, validated the same way [quaidsFit](quaidsFit.md)
  validates its own `weight` (finite, nonnegative, positive sum) and
  surfaced via `weightValid`/`weightSum`/`effN` below. An invalid weight
  is a hard preflight error, the same tier as an invalid `clusterId`.
  Scalar `0` is the only scalar sentinel; scalar nonzero weights are
  invalid.
  Milestone 28 converted both `clusterId` and `weight` from required
  positional arguments (a `0` sentinel meant "not supplied") to genuine
  GAUSS keyword-defaulted parameters -- omit either one entirely, or
  supply by name in any order, instead of always typing `, 0, 0`.

## Returns

`pOut` is a `quaidsPreflightOut` structure. Key fields include:

- `ok`, `nErrors`, `nWarnings`.
- Data-shape and finite-value flags: `dimensionsOk`, `finiteOk`,
  `nonFiniteCount`.
- Share diagnostics: `shareAddOk`, `maxShareSumDev`, `zeroShareCount`,
  `negativeShareCount`.
- Variation diagnostics: `priceStd`, `minPriceStd`, `totexpStd`,
  `instrStd`, `minInstrStd`, and low-variation flags.
- Design/IV diagnostics: `designCols`, `designDf`, `designInvOk`,
  `ivDiagnosticsValid`, `ivFstat`, `ivPvf`, `weakIV`.
- Cluster diagnostics: `clusterValid`, `nClusters`, `minClusterSize`,
  `singletonClusters`, `clusterWarning`.
- Weight diagnostics (Milestone 26): `weightValid`, `weightSum` (sum of
  the raw weight vector, or `nobs` when unweighted), `effN` (Kish's
  effective sample size, `(sum w)^2 / sum(w^2)`; also `nobs` when
  unweighted).
- `convergenceRisk`: `0` low, `1` elevated, `2` high.

## Remarks

`quaidsPreflight()` is silent and does not fit the demand system. It is meant
to be run before [quaidsFit](quaidsFit.md), especially in applied scripts
where a bad input matrix or fragile design would otherwise surface later as
a numerical failure.

Hard errors for `ok` include dimension/control problems, non-finite data,
share adding-up failures, negative budget shares, singular first-stage
design, or invalid cluster labels. Warnings include zero shares, low
variation, weak first-stage IV diagnostics, few/singleton clusters, and
elevated convergence risk.

The first-stage F statistic comes from the same internal first-stage helper
used by `quaidsFit()`. `weakIV` uses the conventional simple rule
`min(ivFstat) < 10`; it is a screening diagnostic, not a substitute for a
full weak-IV analysis.

## Examples

```gauss
aCtl = quaidsControlCreate();

pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, clusterId=householdId);

if not pOut.ok;
    call printQuaidsPreflight(pOut);
    end;
endif;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Source

`quaidsdiagnostics.src`

## See Also

[printQuaidsPreflight](printQuaidsPreflight.md),
[quaidsFit](quaidsFit.md), [quaidsWorkflowFit](quaidsWorkflowFit.md),
[quaidsRobustFit](quaidsRobustFit.md)
