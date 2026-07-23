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

## Estimation

- [quaidsFit](command-reference/quaidsFit.md)
- [printQuaids](command-reference/printQuaids.md)
- [quaids](command-reference/quaids.md)
- [quaidsFull](command-reference/quaidsFull.md)

## Hypothesis Tests

- [quaidsHomogeneityTest](command-reference/quaidsHomogeneityTest.md)
- [quaidsJointTest](command-reference/quaidsJointTest.md)

## Elasticities and Diagnostics

- [quaidsElasFit](command-reference/quaidsElasFit.md)
- [printQuaidsElas](command-reference/printQuaidsElas.md)
- [quaidsElas](command-reference/quaidsElas.md)
- [quaidsSlutzky](command-reference/quaidsSlutzky.md)

## Curvature Imposition (requires `optmt`)

`src/quaidscurvature.src` is listed in `package.json`'s `src` array (real,
required public API) but has a hard compile-time dependency on GAUSS's
`optmt` package (`package.json`'s `deps` array lists `optmt` accordingly)
-- `library quaids;` requires `optmt` installed and loaded too.

- [quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)
- [printQuaidsCurvature](command-reference/printQuaidsCurvature.md)

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
- [ptFromQuaidsFamily](command-reference/ptFromQuaidsFamily.md)
