# Changelog

All notable changes to this project are documented here. This project is
pre-alpha and does not yet follow strict semantic versioning guarantees
(see `GOLD_STANDARD_TODO.md` for the release roadmap); version numbers
below match `package.json` at the time each milestone landed.

## 0.7.0 - 2026-07-23

### Added
- `quaidsWelfareFit()`/`printQuaidsWelfare()` (`src/quaidswelfare.src`):
  exact compensating variation (CV) and equivalent variation (EV) for a
  price change, holding nominal expenditure fixed, with delta-method
  standard errors. Works for LA-AIDS, iterated AIDS, and QUAIDS alike --
  a closed-form evaluation of the already-fitted expenditure function at
  two points, no new estimation and no new package dependency (unlike
  Milestone 10's curvature imposition).
- `quaidsWelfareOut` struct (`src/quaids.sdf`).
- `tests/quaids_welfare_test.e` (20 checks): exact zero-price-change
  identity, exact round-trip inverse-function identity, first-order
  (Marshallian) limiting-case consistency, and sign-agreement between CV
  and EV -- exercised on both a QUAIDS and an AIDS/LA-AIDS fit from the
  existing `_quaidsSyntheticDGP` fixture.

### Notes
- The QUAIDS expenditure-function inversion (Banks, Blundell & Lewbel
  1997) was verified by Shephard's-lemma cross-check against the already-
  validated QUAIDS share equation before implementation, after an initial
  from-memory attempt failed that same check (a sign/inversion error,
  caught before any code was written). See `docs/METHODOLOGY_NOTES.md`.

## 0.6.0 - 2026-07-22

### Added
- `quaidsCurvatureFit()`/`printQuaidsCurvature()` (`src/quaidscurvature.src`):
  local curvature (Slutzky negative semidefiniteness) imposition for
  LA-AIDS/AIDS via the Diewert-Wales (1987) Cholesky reparametrization,
  evaluated at the sample mean. Built on GAUSS's `optmt` package (now a
  real package dependency -- `package.json`'s `deps` array is no longer
  empty) for the small (`vech(A)`-only) profiled nonlinear IV step; the
  remaining coefficients are exactly identified by OLS once the
  reparametrized `gamma` is substituted in, reusing the same
  `moment`/`solpd` primitives as `quaidsFit()` itself.
- `quaidsCurvOut` struct (`src/quaids.sdf`).
- `tests/quaidsfixtures.src`: `_quaidsCurvatureSyntheticDGP()`, a
  synthetic dataset whose true gamma is curvature-consistent at its own
  sample mean by construction, found by direct seed screening (most seeds
  make the underlying nonlinear system numerically unstable -- the same
  kind of iteration sensitivity already documented for the main estimator).
- `tests/quaids_curvature_test.e` (17 checks): recovery, exact adding-up/
  homogeneity/symmetry, near-exact negative-semidefiniteness at the
  reference point, and a non-vacuousness check confirming the
  unconstrained fit genuinely violates curvature on this fixture.
- QUAIDS (the quadratic extension) curvature imposition is explicitly
  deferred, not silently absent -- see `docs/FEATURE_SUPPORT_MATRIX.md`
  and `GOLD_STANDARD_TODO.md`'s Milestone 10 section.

### Known limitation
- Standard errors from `quaidsCurvatureFit()` are a simplified,
  homoskedastic-NLS delta-method approximation, not a full SUR/GMM
  sandwich, and are unreliable specifically when the estimated Cholesky
  factor has boundary (near-zero) entries -- a known, documented
  complication of Cholesky-based negative-semidefinite-cone estimation,
  not unique to this implementation. Point estimates and the exact
  curvature property are unaffected.

## 0.5.0 - 2026-07-20

### Added
- `src/pubtable_quaids.src`: optional adapter onto the `pubtable` package
  for LaTeX/Markdown/CSV/RTF/HTML/XLSX export. `ptModelFromQuaids`/
  `ptFromQuaids` build coefficient tables (one comparison column per good)
  from `quaidsOut`; `ptModelFromQuaidsElas`/`ptFromQuaidsElas`/
  `ptTablesFromQuaidsElas` build income/uncompensated/compensated
  elasticity tables from `quaidsElasOut`; `ptFromQuaidsFamily` dispatches
  on either struct type. Guarded by `#ifDef QUAIDS_SDF_INCLUDED`, following
  the same pattern as `pubtable`'s own bundled `pubtable_qardl.src`, but
  kept inside this repo rather than the installed `pubtable` package (see
  `GOLD_STANDARD_TODO.md`'s Milestone 6 section for why).
- `examples/pubtable_export_example.e`: exports a coefficient table and all
  three elasticity tables to `.tex`/`.md`/`.csv`.
- `tests/quaids_pubtable_test.e` (30 checks): exact numeric parity between
  `pubtable`'s `ptModel.estimates`/`stdErrors` and the `quaidsOut`/
  `quaidsElasOut` values they're built from, plus an end-to-end export
  smoke test that reads the generated files back.
- `#ifndef`/`#define QUAIDS_SDF_INCLUDED`/`#endif` include guard on
  `src/quaids.sdf`, so `pubtable_quaids.src` has a symbol to detect.
- Package build/release tooling (`scripts/build_lcg.ps1`,
  `scripts/build_package.ps1`, `scripts/verify_release_artifact.ps1`,
  `scripts/run_release_verification.ps1`, `tests/verify_package_manifest.ps1`,
  `tests/run_source_tests.ps1`) and an installed-package public API gate
  (`tests/package_public_api.e`). The package is now actually installed and
  loadable via `library quaids;`.
- Full documentation set: `README.md`, `docs/COMMAND_REFERENCE.md` plus 18
  `docs/command-reference/*.md` pages, `docs/USAGE_GUIDE.md`,
  `docs/METHODOLOGY_NOTES.md`, `docs/FEATURE_SUPPORT_MATRIX.md`.
- Published-data validation extended to iterated AIDS (`aCtl.linear=1,
  aCtl.maxiter>1`) against R's `micEconAids::aidsEst(method="IL", ...)`,
  alongside the existing LA-AIDS check
  (`tests/quaids_published_validation_test.e`, now 19 checks).

### Changed
- `examples/quaids_example.e` and `examples/pubtable_export_example.e` now
  use `library quaids;` against the installed package, matching
  README.md/USAGE_GUIDE.md's documented usage pattern, instead of
  `#include`-ing the source tree directly (only meaningful once the
  package became installable).

### Fixed
- `tests/package_public_api.e` now actually calls `printQuaids()` and
  `quaidsElas()` -- both are required, `package.json`-listed public procs
  that the installed-package gate exercised only indirectly (via their
  split `quaidsFit`/`quaidsElasFit`/`printQuaidsElas` components) until the
  final integration pass found the gap.

No public API in `src/` changed across this tooling/documentation/testing
work, so the package version did not bump for it -- see
`GOLD_STANDARD_TODO.md`'s Milestones 7-9 for the full writeups.

## 0.4.0 - 2026-07-20

### Added
- `quaidsElasFit(b, v, intcpt, prices, totexp, aCtl)`: silent,
  struct-returning elasticity computation (`quaidsElasOut`: point estimates
  plus delta-method standard errors) at any evaluation point.
- `printQuaidsElas(elasOut)`: separated console printer for a
  `quaidsElasOut`.
- `tests/quaids_elasticities_test.e` (17 checks): parity between
  `quaidsElasFit()`/`printQuaidsElas()` and the pre-split `quaidsElas_()`,
  plus three exact algebraic identities (Engel aggregation, Cournot
  aggregation, elasticity homogeneity) checked at an out-of-sample
  observation and a synthetic counterfactual price scenario.

### Changed
- `quaidsElas(b, v, intcpt, prices, totexp, aCtl)` is now a thin wrapper
  (`quaidsElasFit()` then `printQuaidsElas()`); its signature and printed
  output are unchanged.

## 0.3.0 - 2026-07-20

### Added
- `src/quaidstests.src`: `quaidsHomogeneityTest(qOut)` and
  `quaidsJointTest(qOut)`, standalone Wald chi-squared tests against an
  unconstrained (`aCtl.homogenous == 0`) `quaidsOut` fit.
- `tests/quaids_hypothesis_tests_test.e` (19 checks): size and power for
  both new tests, a power check for the existing symmetry-given-homogeneity
  test, and the first exercise in this repo's history of the
  overidentification test (`ninst > nu`).

## 0.2.0 - 2026-07-20

### Added
- `src/quaids.src` split into `quaidsiv.src` (`_quaidsIVFirstStage()`),
  `quaidselas.src` (elasticities), `quaidsslutzky.src` (Slutzky negativity
  diagnostic), and `quaids.src` (core estimation, printing, legacy
  wrapper).
- `quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars,
  aCtl)`: dataframe/column-name entry point (`src/quaidsformula.src`).
- `tests/quaids_formula_parity_test.e` (17 checks): `quaidsFull()` vs.
  `quaidsFit()` parity, including the `extraVars == 0` path.
- `tests/quaidsfixtures.src`: shared synthetic 5-good AIDS/QUAIDS data
  generator used by the validation test suite below.
- `tests/quaids_synthetic_validation_test.e` (22 checks): recovers true DGP
  parameters within a documented tolerance across LA-AIDS/iterated-AIDS/
  QUAIDS x with/without-IV.
- `tests/quaids_published_validation_test.e` (11 checks): cross-checks
  `quaidsFit()` against an independent R (`micEconAids`) reference on
  published `Blanciforti86` food-demand data.

### Fixed
- Stone price index starting value (`quaidsFit()`'s "STARTING VALUE"
  block): was computed on a partially-relative, partially-absolute price
  matrix, silently producing wrong results for every `aCtl.maxiter==1`
  (LA-AIDS) call using default starting values, and a bad (though
  self-correcting via iteration) starting point for `aCtl.maxiter>1`
  calls. Found via the published-data cross-check against R (~5x
  discrepancy), confirmed by isolating the formula in Python. After the
  fix, GAUSS matches the R reference to within ~0.021 max absolute
  difference (was ~5x off before).
- Two dead-code bugs in `quaidsFit()`'s `intcpt == 0` branch (uninitialized
  `xnam`; a `G0071 Type mismatch` assigning a native string literal into a
  `matrix`-typed struct field), surfaced by the first-ever exercise of that
  code path in `quaids_formula_parity_test.e`.

## 0.1.0 - 2026-07-20

### Added
- Initial `quaids` package scaffold: renamed from `aids` (Milestone 0;
  QUAIDS is the more general model actually implemented, and `aids` risked
  colliding as a bare identifier). Dead code removed (`rankControl`/
  `latentControl` structs and their unused constructors). `package.json`,
  `LICENSE` (MIT), `CITATION.cff`, `.gitignore`, `CLAUDE.md` added.
- `quaidsFit(w, intcpt, prices, totexp, instr, aCtl)`: silent,
  struct-returning estimation core (`quaidsOut`, ~75 fields), split out of
  the original combined estimate-and-print `quaids()`/`aids()` proc.
- `printQuaids(qOut)`: separated console estimation-report printer.
- `tests/quaids_schema_test.e` (34 checks): `quaidsOut` field values/
  shapes, that `quaidsFit()` prints nothing, and that the legacy `quaids()`
  wrapper's returned matrices are byte-identical to the struct fields
  they're drawn from.

### Changed
- `quaids(w, intcpt, prices, totexp, instr, aCtl)` (formerly `aids()`) is
  now a thin wrapper: `quaidsFit()` then `printQuaids()` then the legacy
  elasticities-at-four-points/descriptive-statistics/Slutzky report.
  Signature and printed output verified byte-for-byte identical to the
  pre-split code at the time of the split.
