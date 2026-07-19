# GAUSS AIDS Library Gold Standard Roadmap

Status date: 2026-07-19

This is the release-readiness checklist and roadmap for turning this repository
into the reference GAUSS implementation of the Almost Ideal Demand System (AIDS)
model family: LA-AIDS, iterated (nonlinear price index) AIDS, and QUAIDS, with
instrumental-variables treatment of total expenditure. It follows the same
structure as `GOLD_STANDARD_TODO.md` in the `gauss-qardl` repository so the two
libraries stay consistent to maintain and to use.

## Current Status Snapshot

The repository is pre-alpha. **Milestones 0 (repository hygiene) and 1
(API/output-schema baseline) are complete** as of 2026-07-19:

- Milestone 0: dead code removed, files moved into `src/`/`examples/`,
  package/proc naming decided (`quaids`), license decided (MIT).
  `package.json`, `LICENSE`, `CITATION.cff`, `.gitignore`, and `CLAUDE.md`
  exist at the repo root.
- Milestone 1: estimation split from printing. `quaidsFit()` is a new,
  silent, struct-returning entry point (`quaidsOut`, ~75 fields). `quaids()`
  remains the original call (unchanged signature, unchanged printed
  behavior, verified byte-for-byte against the pre-Milestone-1 baseline). A
  first automated test, `tests/quaids_schema_test.e`, checks `quaidsOut`
  field shapes/values, silence, and wrapper/struct consistency.

`docs/` still does not exist — that is Milestone 8. A fuller `tests/`
harness (installed-package tests, published/deterministic numerical
fixtures) is Milestones 3 and 7.

- `src/quaids.src` (formerly `aids_rev.src`) — one ~1,300-line proc, `quaids()`
  (formerly `aids()`), that does everything:
  IV first-stage regression for log total expenditure, iterated FGLS-style
  estimation of the AIDS/QUAIDS share system, homogeneity imposition,
  symmetry-constrained re-estimation via minimum distance, an
  overidentification test, absolute-price-effect recovery, elasticities at
  four fixed evaluation points (mean/Q1/Q2/Q3) with delta-method standard
  errors, descriptive statistics, and a Slutzky-negativity eigenvalue
  diagnostic. Estimation and printing are interleaved throughout — there is
  no way to run the model and get back a result object without triggering
  ~10 pages of console output.
- `src/quaidsutil.src` (formerly `aidsutil.src`) — `quaidsControlCreate()`
  (formerly `aidsControlCreate()`, used), plus formerly also
  `rankControlCreate()` and `latentControlCreate()`, which built control
  structures that nothing in this repo called. **Removed at Milestone 0** —
  dead code, likely copied from a template.
- `src/quaids.sdf` (formerly `aid_model.sdf`) — formerly three structs:
  `aidsControl` (used, renamed to `quaidsControl`), plus `rankControl` and
  `latentControl` (both dead, matching the unused constructors above).
  **`rankControl`/`latentControl` removed at Milestone 0.**
- `examples/quaids_example.e` (formerly `aids_example.e`) — one synthetic
  5-good dataset with parameters chosen to satisfy homogeneity/symmetry by
  construction, run through `quaids()` with
  eyeball comparison of printed estimates to true values. No assertions, no
  pass/fail signal, not runnable as a regression test.
- Still missing: `README.md`, `CHANGELOG.md`, `docs/`, `tests/` — see
  Milestones 7–8.

## What GAUSS Already Provides — Do Not Duplicate

Checked against the installed GAUSS 26 runtime (`src/`) and installed
packages (`pkgs/`), and against `aptech/gauss-llm-reference`.

- **No built-in SUR / systems-of-equations estimator and no existing
  AIDS/QUAIDS/demand-system implementation anywhere in the GAUSS runtime or
  shipped packages.** The iterated, cross-equation-restricted (homogeneity +
  symmetry) FGLS core in `src/quaids.src` is genuinely new functionality — it
  is this library's reason to exist. Keep and harden it; do not look for a
  built-in replacement.
- **`gmmFitIV` / `gmm.sdf` / `gmm_est.src` / `gmm_hac.src` /
  `gmm_weight_mat.src`** — a full single-equation IV-GMM estimator with
  robust/HAC weighting and formula-string/dataset support. The hand-rolled
  first-stage 2SLS block at the top of `quaids()` (moment-matrix IV regression
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
  by `#ifDef <LIB>_SDF_INCLUDED`. Build `pubtable_quaids.src` the same way
  instead of hand-rolling any LaTeX/CSV/Markdown export. All of the current
  `printfm()`-based console tables in `src/quaids.src` should eventually be
  reachable through `pubtable` model objects as well as the current console
  form.
- **Dataframes and formula strings** — `loadd()`, `asdf()`, formula parsing
  (`"w1 + w2 + w3 ~ p1 + p2 + p3 + totexp"`-style), `getColNames()`,
  `getColTypes()`, `selif()`/`packr()` for missing-data handling. The current
  API is matrix-only and positional (`w, intcpt, prices, totexp, instr,
  aCtl`). Add a formula/dataframe entry point on top of the matrix core,
  mirroring `qardl`'s `applyQARDLFormula()` / `qardlFull(..., formula=...)`
  pattern, rather than writing a bespoke formula parser.
- **Core primitives already used correctly** in `src/quaids.src` — `moment()`,
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

### Milestone 0 — Repository Hygiene — COMPLETE (2026-07-19)

- [x] Remove or justify `rankControl`/`latentControl` and their constructors
  in `aid_model.sdf`/`aidsutil.src` — dead code from an unrelated template.
  Removed both structs and both constructors; `quaidsControl`/
  `quaidsControlCreate()` are the only surviving struct/constructor.
- [x] Decide the package name/public proc prefix. **Decision: `quaids`.**
  QUAIDS is the more general model actually implemented (linear AIDS is a
  special case), and `aids` was judged too likely to collide/confuse as a
  bare identifier. All public procs renamed: `aids()` -> `quaids()`,
  `slutzky()` -> `quaidsSlutzky()`, `elas()` -> `quaidsElas()`, `elas_()` ->
  `quaidsElas_()`, `aidsControlCreate()` -> `quaidsControlCreate()`. The
  `aidsControl` struct is renamed to `quaidsControl`. The now-vestigial
  `aCtl.aids` field is unchanged (still unread; see Milestone 1) and no
  longer collides with the package name.
- [x] Add `.gitignore`, `LICENSE` decision, `CITATION.cff`.
  **License decision: MIT** (matches `gauss-qardl`), copyright Eric Clower.
- [x] Move the current root-level `.e`/`.src`/`.sdf` files into `src/` and
  `examples/` to match `dccelib`/`qardl` layout, renaming files to match the
  new `quaids` prefix: `aid_model.sdf` -> `src/quaids.sdf`, `aidsutil.src` ->
  `src/quaidsutil.src`, `aids_rev.src` -> `src/quaids.src`, `aids_example.e`
  -> `examples/quaids_example.e`.
- [x] Verify `examples/quaids_example.e` still runs correctly against the
  renamed/relocated source via `tgauss -b -x` after the rename, and that the
  renamed/relocated codebase compiles cleanly and produces byte-identical
  output to the pre-Milestone-0 baseline (see Milestone 0 Verification
  below).

#### Milestone 0 Verification

Two independent checks were run with GAUSS 26.1.1's console runner
(`tgauss.exe -b -x`), since this repository has no test harness yet
(Milestone 3/7):

1. **Compile check** — `#include`d `src/quaids.sdf`, `src/quaidsutil.src`,
   and `src/quaids.src` together from a fresh `tgauss -b -x -e "..."`
   invocation with no other statements. **Result: PASS** — compiled and ran
   with no parse/link errors, confirming the renamed procs/structs
   (`quaids()`, `quaidsSlutzky()`, `quaidsElas()`, `quaidsElas_()`,
   `quaidsControl`, `quaidsControlCreate()`) all resolve correctly with no
   duplicate or dangling identifiers left over from the rename.
2. **Behavioral parity check** — ran the pre-Milestone-0 code
   (`aids_rev.src`/`aidsutil.src`/`aid_model.sdf`/`aids_example.e`, checked
   out from commit `ac5d924` into a scratch directory) and the
   post-Milestone-0 code (`src/quaids.src`/`src/quaidsutil.src`/
   `src/quaids.sdf`/`examples/quaids_example.e`) side by side, both via
   `tgauss -b -x`, both against the same fixed-seed (`seed = 11`) synthetic
   5-good dataset generated inside the example script. Diffed the two
   captured `output file=out` result files byte-for-byte.
   **Result: PASS** — `diff` reported zero differences (656/656 lines
   identical) across the full instrumental-regression tables, iteration log,
   homogeneity/symmetry-constrained coefficient tables, overidentification
   test, elasticities at all four evaluation points, descriptive statistics,
   and Slutzky-eigenvalue diagnostic. This confirms Milestone 0 was a pure
   rename/reorganization with no change to numerical behavior.

**Standard applied**: Milestone 0 is repository hygiene only (renaming,
moving files, deleting dead code, adding metadata) — no estimation logic
should change. The pass bar was therefore "compiles and runs" (compile
check) plus "identical numerical output on identical input" (parity check),
rather than a correctness check against an external/published benchmark,
which is out of scope until Milestone 3.

### Milestone 1 — API and Output Schema Baseline — COMPLETE (2026-07-19)

Goal: make `quaids()` callable without side-effect printing and return a
predictable structure, before anything else is built on top of it.

**Design decision, stated explicitly since it reads as a deviation from the
literal wording below**: `quaids()` itself was *not* changed to return a
struct instead of its four legacy matrices, because that would have broken
`{ b1, v1, b2, v2 } = quaids(...)`-style calls (the "keep the current
positional call signature working" requirement, immediately below).
Instead, a new proc `quaidsFit()` was added as the silent, struct-returning
primary entry point, and `quaids()` became a thin backward-compatible
wrapper around it (calls `quaidsFit()`, calls `printQuaids()`, reproduces
the legacy elasticities/descriptive-stats/Slutzky report, returns the same
four matrices). This mirrors the `predictQARDL`/`predictARDL` legacy-wrapper
pattern already used in `gauss-qardl`.

- [x] Split estimation from printing: `quaidsFit()` returns a struct
  (`quaidsOut`) with parameter estimates/covariances/fit stats/residuals/
  model metadata; printing is `printQuaids(quaidsOut)`, a separate proc,
  following the `qardl`/`printQARDL` split. `printQuaids()` covers the
  estimation-stage report (IV first stage, iteration summary,
  homogeneity-constrained table, overidentification test, symmetry test +
  symmetry-constrained table). Elasticities/descriptive-stats/Slutzky remain
  separate, explicitly-callable reports (`quaidsElas()`, `quaidsSlutzky()`),
  not part of `quaidsOut`, since Milestone 5 plans to generalize
  elasticities to arbitrary evaluation points rather than bake a fixed set
  into the struct now.
- [x] Define `quaidsOut` fields: model family (`qOut.model`,
  `"LA-AIDS"|"AIDS"|"QUAIDS"`), homogeneity/symmetry validity flags
  (`homogenous`, `symValid`), `n` goods, sample size (`nobs`), instrument
  list (`ninst`, `znam`), coefficient/covariance blocks for every stage
  (IV first-stage `iv*`, homogeneity-constrained `homog*`, overidentification
  `overid*`, symmetry test `sym*`, symmetry-constrained `symc*`), residuals
  (`u`), first-stage IV diagnostics (`ivRsq`/`ivFstat`/`ivPvf`/...),
  log-det-sigma fit criteria (`homogCrit`, `symcCrit`), and the raw
  `b`/`v`/`bS`/`vS` matrices for backward compatibility (exactly what
  `quaids()` returns as `b1`/`v1`/`b2`/`v2`). ~75 fields total; see
  `src/quaids.sdf` and the field-by-field notes in `CLAUDE.md`.
- [x] Add `getDefaultQuaidsControl()` alongside the existing
  `quaidsControlCreate()`, aligned with `qardl`'s `getDefault...Control`
  convention (a thin alias). Also upgraded both to GAUSS structure-inference
  return typing (`proc (struct quaidsControl) = ...`).
- [x] Add schema tests asserting `quaidsOut` field names/shapes:
  `tests/quaids_schema_test.e`, 34 checks (see Milestone 1 Verification).
- [x] Keep the current positional call signature working: `quaids()`'s
  signature, return values, and full printed console report are unchanged —
  verified byte-for-byte identical to the pre-Milestone-1 code on the same
  fixed-seed synthetic dataset (see Milestone 1 Verification).

Also folded in while touching `quaidsControl`: removed the `stone`, `aids`,
and `varname` fields flagged as dead-but-kept at Milestone 0 (confirmed
still unread; no external consumers of this pre-alpha struct exist yet, so
removing was safe and was re-verified against the same parity/schema tests).

#### Milestone 1 Verification

Same standard as Milestone 0 (compiles/runs, plus deterministic parity
where behavior must not change), extended with a real schema/unit test
since there is now a struct-shaped output to assert against.

1. **Compile check** — as at Milestone 0, `#include`d all three source files
   and ran with no other statements. **PASS.**
2. **Behavioral parity check (`quaids()` legacy wrapper)** — reran the exact
   Milestone 0 procedure: same fixed-seed (`seed = 11`) synthetic dataset,
   `tgauss -b -x examples/quaids_example.e`, diffed the captured
   `output file=out` result against the untouched Milestone-0-era baseline
   (itself already diffed against the original pre-Milestone-0 code).
   **PASS — zero differences across all 656 lines**, including the
   per-iteration convergence log, which `quaidsFit()` cannot print directly
   (it doesn't print at all) but reconstructs faithfully from a stored
   `qOut.iterHistory` matrix so `printQuaids()` reproduces it unchanged. This
   was re-run again after removing the three dead `quaidsControl` fields,
   with the same zero-diff result.
3. **New: `tests/quaids_schema_test.e`** — 34 checks against a `quaidsFit()`
   call on the same synthetic dataset:
   - A "silence window" check: nothing prints between calling `quaidsFit()`
     and it returning.
   - Metadata/shape checks: `model`, `homogenous`/`linear` echoing,
     dimension fields (`nobs`/`n`/`nint`/`ninst`/`nu`), name-vector lengths,
     `intcptFull`/`u` shapes.
   - Iteration bookkeeping: this particular fixture (`err=.001`,
     `maxiter=100`) does not actually converge within the iteration cap
     (confirmed against the printed log: iteration 99, err≈20.9) — the test
     asserts `quaidsFit()` honestly reports `converged == 0` and
     `iterations == maxiter`, rather than assuming convergence.
   - Per-stage shape checks: `ivB`, `homogB`, `symcB` dimensions;
     `overidValid == 0` (correctly, since `ninst == nu` here — exactly
     identified); `symValid == 1`; `symDf`/`symStat`/`symPval` sanity.
   - Final-output shape checks: `b`/`bS`/`bestB` all reshape to `(ng+1) x n`
     after absolute-price recovery, as expected when `homogenous == 1`.
   - **Cross-check**: calls the legacy `quaids()` wrapper on the identical
     inputs and asserts its four returned matrices are *exactly* equal
     (`maxc(maxc(abs(...))) == 0`, not just close) to `qOut.b`/`qOut.v`/
     `qOut.bS`/`qOut.vS` — i.e., the wrapper and the new struct-returning
     core are provably drawing from the same numbers, not two independently
     drifting code paths.
   **PASS — all 34 checks pass** (`SCHEMA TEST: ALL 34 CHECKS PASSED`).

**Standard applied**: Milestone 1 is a structure-preserving refactor of a
correctness-sensitive estimator, so the bar was deliberately strict —
`quaids()`'s output had to remain byte-identical (not just "close enough"),
and the new struct's fields had to be checked both for internal consistency
(shapes, flags) and for consistency with the legacy code path they replace
(the exact-equality cross-check), not merely "did it run."

**One real bug caught and fixed during this work**: the first `quaidsOut`
draft declared the name-vector fields (`xnam`, `wnam`, `znam`, `unam`,
`enam`) as `string array`, which produced a hard `G0071 Type mismatch`
runtime error — GAUSS's classic `0$+"X"$+ftocv(...)` idiom produces a
character matrix, not the newer native `string array` type. Fixed by
declaring those fields `matrix`. Documented in `CLAUDE.md` so it isn't
rediscovered the hard way again.

**A pre-existing anomaly identified, not fixed**: while relocating the
symmetry-constrained table's print logic, noticed that the original code
pairs the *homogeneity-stage* `b`'s IV-residual-block point estimates with
the *symmetry-stage*'s standard errors/t/p-values in that one table row.
This was carried over unchanged (correct per the "structure-preserving, not
correctness-fixing" scope of this milestone) and flagged in `CLAUDE.md` for
review during Milestone 3/4 validation.

### Milestone 2 — Modular Source Split + Formula/Dataframe Entry Point

- [ ] Split `src/quaids.src` into focused files: core estimation, IV first
  stage, homogeneity/symmetry testing, elasticities, Slutzky diagnostics,
  printing — mirroring how `qardl.src`/`icmean.src`/`wtestlrb.src`/etc. are
  separated in `gauss-qardl`.
- [ ] Evaluate routing the IV first-stage regression through `gmmFitIV`
  instead of the hand-rolled moment-matrix 2SLS block; document the decision
  either way.
- [ ] Add a formula/dataframe entry point (e.g. `quaidsFull(data, formula,
  instFormula, aCtl)`) built on `loadd()`/`asdf()`, matching the
  `applyQARDLFormula` pattern instead of a bespoke parser.
- [ ] Add formula-vs-matrix parity tests.

### Milestone 3 — Validation Fixture and Benchmark Harness

- [ ] Add deterministic synthetic fixtures (expand on `quaids_example.e` with
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
  (currently hardcoded to mean/Q1/median/Q3 inside `quaids()`).
- [ ] Keep the Slutzky-eigenvalue negativity diagnostic as the default
  always-on check.
- [ ] Scope and, if justified, implement optional local curvature imposition
  (Diewert-Wales Cholesky reparametrization) as an opt-in estimation mode
  built on `optmt`/`cmlmt` — P2, not a release blocker.

### Milestone 6 — Reporting via `pubtable`

- [ ] Add `pubtable_quaids.src` following the `pubtable_qardl.src` pattern
  (`#ifDef QUAIDS_SDF_INCLUDED` guard, `ptModelFromQuaids`, `ptFromQuaids`).
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

- [ ] `quaids()` (and formula-based `quaidsFull()`) return structured output with
  no forced console printing.
- [ ] LA-AIDS, iterated AIDS, and QUAIDS are each documented, tested, and
  independently callable model choices, with and without IV.
- [ ] Homogeneity, symmetry, and overidentification each have standalone,
  tested hypothesis-test procs.
- [ ] Elasticities are computable at arbitrary evaluation points with
  delta-method standard errors.
- [ ] Slutzky negativity diagnostics ship by default; curvature imposition is
  either implemented or explicitly documented as future work.
- [ ] `pubtable_quaids.src` provides LaTeX/Markdown/CSV export.
- [ ] Package builds, installs, and passes an installed-package public API
  test, matching the `qardl`/`dccelib` release process.
- [ ] Full doc set (`README`, command reference, usage guide, methodology
  notes, feature support matrix, `CLAUDE.md`) exists and is synchronized with
  the code.
