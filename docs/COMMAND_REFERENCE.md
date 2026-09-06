# QUAIDS Command Reference

This command reference follows the standard GAUSS documentation pattern:
each user-facing procedure has a page with `Purpose`, `Format`,
`Parameters`, `Returns`, `Remarks`, `Examples`, `Source`, and `See Also`
sections.

## User Guides

- [Usage guide](USAGE_GUIDE.md)
- [Methodology notes](METHODOLOGY_NOTES.md)
- [Feature support matrix](FEATURE_SUPPORT_MATRIX.md)

## Control Struct

- [quaidsControlCreate](command-reference/quaidsControlCreate.md)
- [getDefaultQuaidsControl](command-reference/getDefaultQuaidsControl.md)
- [quaidsSetHomogeneity](command-reference/quaidsSetHomogeneity.md)
- [quaidsGetHomogeneity](command-reference/quaidsGetHomogeneity.md)

## Estimation

- [quaidsFit](command-reference/quaidsFit.md)
- [printQuaids](command-reference/printQuaids.md)
- [quaids](command-reference/quaids.md)
- [quaidsFull](command-reference/quaidsFull.md)
- [quaidsWorkflowFit](command-reference/quaidsWorkflowFit.md)
- [quaidsWorkflowScenarioFit](command-reference/quaidsWorkflowScenarioFit.md)
- [quaidsSurveyWorkflowFit](command-reference/quaidsSurveyWorkflowFit.md)

## Hypothesis Tests

- [quaidsHomogeneityTest](command-reference/quaidsHomogeneityTest.md)
- [quaidsJointTest](command-reference/quaidsJointTest.md)
- [quaidsQuadraticTest](command-reference/quaidsQuadraticTest.md)

## Elasticities and Diagnostics

- [quaidsPreflight](command-reference/quaidsPreflight.md)
- [printQuaidsPreflight](command-reference/printQuaidsPreflight.md)
- [quaidsElasFit](command-reference/quaidsElasFit.md)
- [printQuaidsElas](command-reference/printQuaidsElas.md)
- [quaidsElas](command-reference/quaidsElas.md)
- [quaidsSharesFit](command-reference/quaidsSharesFit.md)
- [printQuaidsShares](command-reference/printQuaidsShares.md)
- [quaidsSlutzky](command-reference/quaidsSlutzky.md)

## Welfare Analysis

- [quaidsWelfareFit](command-reference/quaidsWelfareFit.md)
- [printQuaidsWelfare](command-reference/printQuaidsWelfare.md)

## Zero Budget Share Correction

- [quaidsZeroFit](command-reference/quaidsZeroFit.md)
- [printQuaidsZero](command-reference/printQuaidsZero.md)

## Robust and Cluster-Robust Standard Errors

- [quaidsRobustFit](command-reference/quaidsRobustFit.md)
- [quaidsRobustCovariance](command-reference/quaidsRobustCovariance.md)
- [printQuaidsRobust](command-reference/printQuaidsRobust.md)
- [quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md)
- [quaidsRobustBootstrapCovariance](command-reference/quaidsRobustBootstrapCovariance.md)
- [printQuaidsRobustBootstrap](command-reference/printQuaidsRobustBootstrap.md)

## Replicate-Weight Standard Errors

- [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)
- [printQuaidsReplicateWeight](command-reference/printQuaidsReplicateWeight.md)

## Curvature Imposition (optional, requires `optmt`)

`src/quaidscurvature.src` is **not** listed in `package.json`'s `src`
array and is not loaded by `library quaids;` -- it has a hard compile-time
dependency on the `optmt` package's struct types, and core estimation
(`quaidsFit`, elasticities, welfare, zero-share correction, robust/
replicate-weight standard errors, the applied workflow) needs no external
package at all. A caller who wants curvature imposition `#include`s
`src/quaidscurvature.src` directly, after loading both `quaids` and
`optmt`:

```gauss
library optmt, quaids;
#include quaidscurvature.src
```

See [`docs/public-api.json`](public-api.json)'s `optional_modules` entry
and the
[Imposing Curvature section of the usage guide](USAGE_GUIDE.md#imposing-curvature-diewert-wales).

- [quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)
- [printQuaidsCurvature](command-reference/printQuaidsCurvature.md)
- [quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md)
- [printQuaidsCurvatureBootstrap](command-reference/printQuaidsCurvatureBootstrap.md)
- [quaidsCurvatureBootstrapCI](command-reference/quaidsCurvatureBootstrapCI.md)

## Reporting (optional, requires `pubtable`)

`src/pubtable_quaids.src` is **not** listed in `package.json`'s `src`
array and is not loaded by `library quaids;` -- it has a hard compile-time
dependency on the `pubtable` package's struct types. A caller who wants
these procs `#include`s `src/pubtable_quaids.src` directly, after loading
both `quaids` and `pubtable`. See
[pubtable_quaids.src's own header comment](../src/pubtable_quaids.src) and
the [Reporting section of the usage guide](USAGE_GUIDE.md#reporting-pubtable).

- [ptModelFromQuaids](command-reference/ptModelFromQuaids.md)
- [ptFromQuaids](command-reference/ptFromQuaids.md)
- [ptModelFromQuaidsElas](command-reference/ptModelFromQuaidsElas.md)
- [ptFromQuaidsElas](command-reference/ptFromQuaidsElas.md)
- [ptTablesFromQuaidsElas](command-reference/ptTablesFromQuaidsElas.md)
- [ptTablesFromQuaidsWorkflow](command-reference/ptTablesFromQuaidsWorkflow.md)
- [ptFromQuaidsFamily](command-reference/ptFromQuaidsFamily.md)

## Compatibility (Deprecated)

Retained callable through the `0.x` series (not removed before `1.0.0`)
but not part of the supported API for new code -- see the
[compatibility policy](../README.md#compatibility-policy) and
[`docs/public-api.json`](public-api.json).

- [quaidsElas_](command-reference/quaidsElas_.md) -- use
  [quaidsElasFit](command-reference/quaidsElasFit.md) instead.
