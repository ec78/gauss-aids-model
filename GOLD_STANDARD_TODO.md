# GAUSS AIDS Library Gold Standard Roadmap

Status date: 2026-07-19

This is the release-readiness checklist and roadmap for turning this repository
into the reference GAUSS implementation of the Almost Ideal Demand System (AIDS)
model family: LA-AIDS, iterated (nonlinear price index) AIDS, and QUAIDS, with
instrumental-variables treatment of total expenditure. It follows the same
structure as `GOLD_STANDARD_TODO.md` in the `gauss-qardl` repository so the two
libraries stay consistent to maintain and to use.

## Current Status Snapshot

The repository is pre-alpha: two commits, four files, no package structure.

- `aids_rev.src` — one ~1,300-line proc, `aids()`, that does everything:
  IV first-stage regression for log total expenditure, iterated FGLS-style
  estimation of the AIDS/QUAIDS share system, homogeneity imposition,
  symmetry-constrained re-estimation via minimum distance, an
  overidentification test, absolute-price-effect recovery, elasticities at
  four fixed evaluation points (mean/Q1/Q2/Q3) with delta-method standard
  errors, descriptive statistics, and a Slutzky-negativity eigenvalue
  diagnostic. Estimation and printing are interleaved throughout — there is
  no way to run the model and get back a result object without triggering
  ~10 pages of console output.
- `aidsutil.src` — `aidsControlCreate()` (used), plus `rankControlCreate()`
  and `latentControlCreate()`, which build control structures that nothing in
  this repo calls. Dead code, likely copied from a template.
- `aid_model.sdf` — three structs: `aidsControl` (used), `rankControl` and
  `latentControl` (both dead, matching the unused constructors above).
- `aids_example.e` — one synthetic 5-good dataset with parameters chosen to
  satisfy homogeneity/symmetry by construction, run through `aids()` with
  eyeball comparison of printed estimates to true values. No assertions, no
  pass/fail signal, not runnable as a regression test.
- No `package.json`, no `docs/`, no `tests/`, no `examples/` directory, no
  `README`, `CHANGELOG`, `CITATION`, `LICENSE` decision, or `CLAUDE.md`.

## What GAUSS Already Provides — Do Not Duplicate

Checked against the installed GAUSS 26 runtime (`src/`) and installed
packages (`pkgs/`), and against `aptech/gauss-llm-reference`.

- **No built-in SUR / systems-of-equations estimator and no existing
  AIDS/QUAIDS/demand-system implementation anywhere in the GAUSS runtime or
  shipped packages.** The iterated, cross-equation-restricted (homogeneity +
  symmetry) FGLS core in `aids_rev.src` is genuinely new functionality — it
  is this library's reason to exist. Keep and harden it; do not look for a
  built-in replacement.
- **`gmmFitIV` / `gmm.sdf` / `gmm_est.src` / `gmm_hac.src` /
  `gmm_weight_mat.src`** — a full single-equation IV-GMM estimator with
  robust/HAC weighting and formula-string/dataset support. The hand-rolled
  first-stage 2SLS block at the top of `aids()` (moment-matrix IV regression
  of log total expenditure on instruments) is a candidate to route through
  `gmmFitIV` instead of hand-coded `moment()`/`invpd()` algebra — but this
  needs an explicit evaluation (Milestone 2), because the multi-equation
  system needs raw residuals and moment blocks threaded into the share-system
  normal equations in a specific layout that `gmmFitIV`'s single-equation
  output may not expose directly. Decide, don't assume.
- **`pubtable` (`pkgs/pubtable`)** — a complete publication-table engine:
  LaTeX/HTML/RTF/CSV/XLS/Markdown export, model-comparison tables,
  significance stars, GOF rows. It already ships an adapter pattern for other
  in-house libraries (`pubtable_qardl.src`, `pubtable_cmlmt.src`,
  `pubtable_maxlikmt.src`, `pubtable_tsmt.src`, `pubtable_optmt.src`) guarded
  by `#ifDef <LIB>_SDF_INCLUDED`. Build `pubtable_aids.src` the same way
  instead of hand-rolling any LaTeX/CSV/Markdown export. All of the current
  `printfm()`-based console tables in `aids_rev.src` should eventually be
  reachable through `pubtable` model objects as well as the current console
  form.
- **Dataframes and formula strings** — `loadd()`, `asdf()`, formula parsing
  (`"w1 + w2 + w3 ~ p1 + p2 + p3 + totexp"`-style), `getColNames()`,
  `getColTypes()`, `selif()`/`packr()` for missing-data handling. The current
  API is matrix-only and positional (`w, intcpt, prices, totexp, instr,
  aCtl`). Add a formula/dataframe entry point on top of the matrix core,
  mirroring `qardl`'s `applyQARDLFormula()` / `qardlFull(..., formula=...)`
  pattern, rather than writing a bespoke formula parser.
- **Core primitives already used correctly** in `aids_rev.src` — `moment()`,
  `invpd()`, `solpd()`, `design()`/`vech()`/`xpnd()` for the symmetric-matrix
  restriction algebra, `eigh()` for the Slutzky check, `cdfchic`/`cdftc`/
  `cdffc`/`cdfnc`, `printfm()`, `quantile()`. Keep using these; no rewrite
  needed here.
- **No built-in curvature/negativity imposition.** GAUSS has generic NLP
  solvers (`sqpSolveMT`/`cmlmt`/`optmt`) but nothing AIDS-specific. Imposing
  local curvature (Diewert-Wales Cholesky reparametrization) would be new
  work built on top of `optmt`/`cmlmt`, not a duplicate of anything existing.
  Treat as a P2 feature, not a release blocker.

## Target Model Coverage

- **LA-AIDS** (linearized, Stone price index, one-step): already reachable
  via `aCtl.maxiter = 1`, but not documented, named, or tested as a distinct,
  supported entry point.
- **Iterated AIDS** (nonlinear translog price index, iterated FGLS): the main
  loop already does this when `aCtl.linear = 1`.
- **QUAIDS** (Banks-Blundell-Lewbel 1997 quadratic log-expenditure term): the
  `lx2`/`lambda` block already implements this when `aCtl.linear = 0`, but it
  is an implicit branch of one giant proc, not a named, independently tested
  model choice.
- **Restriction levels**: unconstrained, homogeneity-constrained, and
  homogeneity+symmetry-constrained estimation are all already computed in one
  pass. Needs: a standalone homogeneity test (currently only visible
  qualitatively as "does the reference absolute-price effect look like
  zero"), formalized alongside the existing symmetry-given-homogeneity Wald
  test and the existing overidentification test.
- **IV treatment of total expenditure**: already present and a genuine
  differentiator versus textbook AIDS code most researchers pass around.
  Keep, document explicitly, and add weak-instrument guidance around the
  first-stage F-statistic that is already computed.
- **Elasticities**: income, uncompensated price, compensated price, with
  delta-method standard errors — already implemented but hardcoded to four
  evaluation points (mean, Q1, median, Q3). Generalize to arbitrary
  user-supplied points/covariate profiles.
- **Curvature imposition** (P2, optional): local negativity imposition as an
  opt-in estimation mode, with the existing Slutzky-eigenvalue diagnostic
  kept as the always-available post hoc check.

## Roadmap

Each milestone should exit with source tests, examples, and docs updated
together — no milestone is "done" with code alone.

### Milestone 0 — Repository Hygiene

- [ ] Remove or justify `rankControl`/`latentControl` and their constructors
  in `aid_model.sdf`/`aidsutil.src` — dead code from an unrelated template.
- [ ] Decide the package name/public proc prefix (`aids`, `aidsQuaids`, or a
  neutral `aidsystem`-style prefix) before any public API is locked in —
  check for collisions with common variable names like `aids` used as a
  local control-structure field name today (`aCtl.aids`).
- [ ] Add `.gitignore`, `LICENSE` decision, `CITATION.cff`.
- [ ] Move the current root-level `.e`/`.src`/`.sdf` files into `src/` and
  `examples/` to match `dccelib`/`qardl` layout.

### Milestone 1 — API and Output Schema Baseline

Goal: make `aids()` callable without side-effect printing and return a
predictable structure, before anything else is built on top of it.

- [ ] Split estimation from printing: `aids()` should return a struct
  (`aidsOut`) with parameter estimates/covariances/fit stats/residuals/model
  metadata; printing becomes a separate `printAIDS(aidsOut)` call, following
  the `qardl`/`printQARDL` split.
- [ ] Define `aidsOut` fields: model family (LA-AIDS/AIDS/QUAIDS), homogeneity
  and symmetry flags, `n` goods, sample size, instrument list, alpha/gamma/
  beta/lambda blocks and their covariances, residuals, first-stage IV
  diagnostics, log-det-sigma fit criterion, and the raw `b`/`v` matrices for
  backward compatibility.
- [ ] Add `getDefaultAidsControl()` alongside the existing
  `aidsControlCreate()` if naming should align with `qardl`'s
  `getDefault...Control` convention — otherwise document why the existing
  name stays.
- [ ] Add schema tests asserting `aidsOut` field names/shapes.
- [ ] Keep the current positional call signature working; do not break
  `aids_example.e`-style calls while restructuring.

### Milestone 2 — Modular Source Split + Formula/Dataframe Entry Point

- [ ] Split `aids_rev.src` into focused files: core estimation, IV first
  stage, homogeneity/symmetry testing, elasticities, Slutzky diagnostics,
  printing — mirroring how `qardl.src`/`icmean.src`/`wtestlrb.src`/etc. are
  separated in `gauss-qardl`.
- [ ] Evaluate routing the IV first-stage regression through `gmmFitIV`
  instead of the hand-rolled moment-matrix 2SLS block; document the decision
  either way.
- [ ] Add a formula/dataframe entry point (e.g. `aidsFull(data, formula,
  instFormula, aCtl)`) built on `loadd()`/`asdf()`, matching the
  `applyQARDLFormula` pattern instead of a bespoke parser.
- [ ] Add formula-vs-matrix parity tests.

### Milestone 3 — Validation Fixture and Benchmark Harness

- [ ] Add deterministic synthetic fixtures (expand on `aids_example.e` with
  actual pass/fail assertions instead of eyeballed printouts), covering
  LA-AIDS, iterated AIDS, and QUAIDS, each with and without IV.
- [ ] Identify a published-replication target: Deaton & Muellbauer (1980)
  original AIDS UK data, Banks-Blundell-Lewbel (1997) QUAIDS, or a commonly
  used teaching dataset (e.g. the food-expenditure data distributed with
  Poi's Stata `quaids`/`aidsills`) — confirm redistribution rights before
  committing any external data file.
- [ ] Add a cross-implementation comparison against an existing open
  implementation where licensing allows (R `micEconAids`/`easyNCC`, Stata
  `quaids`) for at least coefficient and elasticity parity.
- [ ] Document tolerance policy for deterministic expected-output tests.

### Milestone 4 — Hypothesis Testing Completeness

- [ ] Add a standalone homogeneity Wald/LR test (unrestricted vs.
  homogeneity-constrained), not just the qualitative "near-zero reference
  price effect" signal.
- [ ] Formalize the existing symmetry-given-homogeneity test and the existing
  overidentification test as independently callable, independently tested
  procs.
- [ ] Add a joint homogeneity+symmetry test.
- [ ] Document degrees of freedom and asymptotic assumptions for each test.

### Milestone 5 — Elasticities and Diagnostics Generalization

- [ ] Generalize `elas()`/`elas_()` to accept arbitrary evaluation points
  (currently hardcoded to mean/Q1/median/Q3 inside `aids()`).
- [ ] Keep the Slutzky-eigenvalue negativity diagnostic as the default
  always-on check.
- [ ] Scope and, if justified, implement optional local curvature imposition
  (Diewert-Wales Cholesky reparametrization) as an opt-in estimation mode
  built on `optmt`/`cmlmt` — P2, not a release blocker.

### Milestone 6 — Reporting via `pubtable`

- [ ] Add `pubtable_aids.src` following the `pubtable_qardl.src` pattern
  (`#ifDef AIDS_SDF_INCLUDED` guard, `ptModelFromAids`, `ptFromAids`).
- [ ] Route elasticity tables and coefficient tables through `pubtable`
  model objects in addition to the console `printfm()` output.
- [ ] Add LaTeX/Markdown/CSV export examples.

### Milestone 7 — Package Build and Release Tooling

- [ ] Add `package.json` (name, version, `src` array, deps, keywords,
  license) matching `dccelib`/`qardl` conventions.
- [ ] Add `tests/run_source_tests.ps1`-style runner and a package-manifest
  consistency check.
- [ ] Add repeatable package build/install scripts and an installed-package
  smoke test (`library aids;` public API gate), matching `qardl`'s
  `tests/package_public_api.e`.
- [ ] Add `CHANGELOG.md` and version the first tagged release.

### Milestone 8 — Documentation

- [ ] `README.md` as the front door (install, quick start, model choices).
- [ ] `docs/COMMAND_REFERENCE.md` plus one `docs/command-reference/*.md` page
  per public proc.
- [ ] `docs/USAGE_GUIDE.md` covering LA-AIDS vs. iterated AIDS vs. QUAIDS,
  with/without IV, formula vs. matrix API.
- [ ] `docs/METHODOLOGY_NOTES.md` documenting the exact estimator (iterated
  linearized/nonlinear FGLS with cross-equation homogeneity/symmetry
  restrictions imposed via minimum distance), citing Deaton & Muellbauer
  (1980) and Banks, Blundell & Lewbel (1997).
- [ ] `docs/FEATURE_SUPPORT_MATRIX.md` (model family x
  diagnostics/tests/elasticities/export/IV support).
- [ ] `CLAUDE.md` context file for future Claude Code sessions, matching the
  one in `gauss-qardl`.

### Milestone 9 — Final Gold Standard Integration Gate

- [ ] Source tests pass; installed-package tests pass; all examples run
  against the installed package.
- [ ] Package exports match public docs.
- [ ] At least one published or independently reproduced validation exists
  per model family (LA-AIDS, iterated AIDS, QUAIDS).
- [ ] Remaining unsupported features (e.g. curvature imposition, if deferred)
  are explicitly documented rather than silently absent.
- [ ] Package artifact, metadata, changelog, and docs have matching version
  numbers.

## Definition of Done for a Gold Standard Release

- [ ] `aids()` (and formula-based `aidsFull()`) return structured output with
  no forced console printing.
- [ ] LA-AIDS, iterated AIDS, and QUAIDS are each documented, tested, and
  independently callable model choices, with and without IV.
- [ ] Homogeneity, symmetry, and overidentification each have standalone,
  tested hypothesis-test procs.
- [ ] Elasticities are computable at arbitrary evaluation points with
  delta-method standard errors.
- [ ] Slutzky negativity diagnostics ship by default; curvature imposition is
  either implemented or explicitly documented as future work.
- [ ] `pubtable_aids.src` provides LaTeX/Markdown/CSV export.
- [ ] Package builds, installs, and passes an installed-package public API
  test, matching the `qardl`/`dccelib` release process.
- [ ] Full doc set (`README`, command reference, usage guide, methodology
  notes, feature support matrix, `CLAUDE.md`) exists and is synchronized with
  the code.
