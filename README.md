# GAUSS QUAIDS Library

A [GAUSS](https://www.aptech.com) application package for estimating
**Almost Ideal Demand System** models: linearized AIDS (Stone price index),
iterated AIDS (nonlinear translog price index), and **QUAIDS** (Banks,
Blundell & Lewbel 1997 quadratic-log-expenditure extension), with
instrumental-variables treatment of (endogenous) log total expenditure.
Estimation imposes homogeneity and/or Slutzky symmetry via iterated FGLS
with cross-equation restrictions applied through a minimum-distance
reparametrization.

Use cases: consumer demand estimation, welfare analysis, elasticity
calculation, testing demand-theory restrictions (homogeneity, symmetry,
overidentification).

This library is **pre-alpha** (package version `0.5.0`). The estimation
core, hypothesis tests, elasticities, diagnostics, dataframe entry point,
`pubtable` export, and release tooling are complete and tested; full
documentation (this doc set) and a final integration gate are in progress
-- see `GOLD_STANDARD_TODO.md` for the release roadmap.

## Requirements

- GAUSS 26 or later.
- No external GAUSS packages are required for estimation. The optional
  `pubtable` package (LaTeX/Markdown/CSV/RTF/HTML/XLSX table export) is
  needed only if you use `src/pubtable_quaids.src` -- see
  [Reporting](#reporting-optional-pubtable) below.

## Installation

Install the release zip in GAUSS using **Tools > Install Application**,
then load the library:

```gauss
library quaids;
```

Or build and install from this repo directly:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact
```

This builds `quaids <version>.zip`, verifies it, and installs it to your
GAUSS package directory (`<GaussHome>\pkgs\quaids` by default). See
[scripts/run_release_verification.ps1](scripts/run_release_verification.ps1)
for options.

## Quick Start

```gauss
library quaids;

// w: TxN budget shares. prices: TxN log prices (absolute). totexp: Tx1 log
// total expenditure (treated as endogenous). instr: TxH instruments.
// intcpt: TxK extra intercept-shifter variables, or 0 for none.
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;          // 0 = QUAIDS, 1 = LA-AIDS
aCtl.maxiter = 100;       // 1 = one-step Stone-index LA-AIDS
aCtl.homogenous = 1;      // impose homogeneity (and test/report symmetry)

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call printQuaidsElas(elasOut);

call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);
```

`quaidsFit()` is silent (no console output) and returns a `quaidsOut`
struct. For the original console-report behavior in one call, use the
backward-compatible wrapper instead:

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
```

Named GAUSS dataframes can be used instead of assembling matrices by hand:

```gauss
data = loadd("mydata.csv");
qOut = quaidsFull(data, shareVars, priceVars, "totexp", "instr", extraVars, aCtl);
```

## Main Features

- LA-AIDS (`aCtl.maxiter = 1`, Stone price index), iterated AIDS
  (`aCtl.linear = 1`, `aCtl.maxiter > 1`, nonlinear translog price index),
  and QUAIDS (`aCtl.linear = 0`) from a single estimator, `quaidsFit`.
- Instrumental-variables (control-function) treatment of log total
  expenditure is always applied -- `quaidsFit` always requires an `instr`
  argument.
- Homogeneity and Slutzky symmetry imposed via iterated FGLS with
  cross-equation restrictions (minimum-distance reparametrization), or left
  unconstrained (`aCtl.homogenous = 0`) for hypothesis testing.
- Standalone Wald tests for homogeneity (`quaidsHomogeneityTest`) and joint
  homogeneity+symmetry (`quaidsJointTest`) against an unconstrained fit,
  plus a built-in symmetry-given-homogeneity test and overidentification
  test as part of `quaidsFit` itself.
- Income and price (Marshallian/Hicksian) elasticities with delta-method
  standard errors at any evaluation point (`quaidsElasFit`) -- not just the
  sample mean/quartiles.
- Slutzky negativity diagnostic (`quaidsSlutzky`), evaluated
  observation-by-observation.
- Dataframe/column-name entry point (`quaidsFull`) alongside the matrix API.
- Optional LaTeX/Markdown/CSV/RTF/HTML/XLSX export via the `pubtable`
  package (`src/pubtable_quaids.src`).
- Cross-implementation validated against an independent R (`micEconAids`)
  reference on published data, and against known-true synthetic DGPs across
  all six (LA-AIDS/AIDS/QUAIDS) x (with/without-IV) combinations.

## Documentation

- [Command reference](docs/COMMAND_REFERENCE.md): one page per public
  procedure, with purpose, format, parameters, returns, remarks, examples,
  source, and related commands.
- [Usage guide](docs/USAGE_GUIDE.md): choosing an API, model-choice
  switches, IV requirements, elasticities/diagnostics workflow, `pubtable`
  export.
- [Methodology notes](docs/METHODOLOGY_NOTES.md): the estimator itself --
  iterated linearized/nonlinear FGLS with cross-equation homogeneity/
  symmetry restrictions via minimum distance, citing Deaton & Muellbauer
  (1980) and Banks, Blundell & Lewbel (1997).
- [Feature support matrix](docs/FEATURE_SUPPORT_MATRIX.md): LA-AIDS vs.
  iterated AIDS vs. QUAIDS support for IV, hypothesis tests, elasticities,
  diagnostics, and export.
- [CLAUDE.md](CLAUDE.md): detailed context file for AI coding assistants
  (and human contributors) working on this repository -- design decisions,
  real bugs found and fixed, GAUSS-specific gotchas.

## Reporting (optional, `pubtable`)

```gauss
library pubtable, quaids;
#include pubtable_quaids.src   // not in package.json's src array -- see docs/COMMAND_REFERENCE.md

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);
call ptExport(coefTbl, "results.tex");   // .tex/.md/.csv/.rtf/.html/.xlsx by extension

struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);  // 3x1: income, uncompensated, compensated
call ptExportAll(elasTbls, "elasticities");
```

Requires the [pubtable](https://github.com/aptech/gauss_table_creator)
package installed separately. See `examples/pubtable_export_example.e` for
a full runnable example.

## Examples

The `examples/` directory contains runnable GAUSS programs:

| File | Description |
| --- | --- |
| `quaids_example.e` | End-to-end synthetic-data workflow: fit, print, eyeball-compare to true parameters |
| `pubtable_export_example.e` | Export a coefficient table and elasticity tables to LaTeX/Markdown/CSV |

## Testing

Run the source-tree test suite (7 test files, 150 checks, plus a
package-manifest consistency check):

```powershell
powershell -ExecutionPolicy Bypass -File tests\run_source_tests.ps1
```

After building and installing the package, verify the installed public API:

```gauss
tgauss -b -x tests/package_public_api.e
```

Build and verify a release artifact:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build_package.ps1 -Force
```

Run the full release gate, including artifact installation into the GAUSS
package directory:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_release_verification.ps1 -InstallArtifact
```

## Citation

See [CITATION.cff](CITATION.cff) for citation metadata, and cite the
underlying methodology:

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312-326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527-539.

## License

The package license is listed in `package.json` (MIT). See
[LICENSE](LICENSE).

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312-326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527-539.
- Blanciforti, L., Green, R., King, G. (1986). *U.S. Consumer Behavior Over
  the Postwar Period: An Almost Ideal Demand System Analysis*. Giannini
  Foundation Monograph No. 40. Used for published-data cross-validation --
  see `tests/fixtures/published/SOURCE.md`.

## Author

[Eric Clower](mailto:eric.clower78@gmail.com)
