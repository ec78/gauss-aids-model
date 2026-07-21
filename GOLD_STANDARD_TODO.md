# GAUSS AIDS Library Gold Standard Roadmap

Status date: 2026-07-20

This is the release-readiness checklist and roadmap for turning this repository
into the reference GAUSS implementation of the Almost Ideal Demand System (AIDS)
model family: LA-AIDS, iterated (nonlinear price index) AIDS, and QUAIDS, with
instrumental-variables treatment of total expenditure. It follows the same
structure as `GOLD_STANDARD_TODO.md` in the `gauss-qardl` repository so the two
libraries stay consistent to maintain and to use.

## Current Status Snapshot

The repository is pre-alpha, package version `0.5.0`. **Milestones 0
(repository hygiene), 1 (API/output-schema baseline), 2 (modular source
split + dataframe entry point), 3 (validation fixtures), 4 (hypothesis
testing completeness), 5 (elasticities/diagnostics generalization), 6
(reporting via `pubtable`), and 7 (package build and release tooling) are
all complete** as of 2026-07-20:

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
- Milestone 2: `src/quaids.src` split into `quaidsiv.src` (IV first stage),
  `quaidselas.src` (elasticities), `quaidsslutzky.src` (Slutzky diagnostic),
  and `quaids.src` (core estimation + printing + legacy wrapper) — see the
  Milestone 2 entry below for why the rest of the estimation core stayed
  together. `gmmFitIV` evaluated and explicitly not adopted (documented
  reasons in `CLAUDE.md`). New dataframe entry point `quaidsFull()`
  (`quaidsformula.src`), verified against the matrix API by a new
  `tests/quaids_formula_parity_test.e` (17 checks) — which in the process
  caught and fixed two real pre-existing bugs in a previously-untested code
  path (`intcpt == 0`).
- Milestone 3: deterministic synthetic recovery fixtures across all 6
  model/endogeneity combinations (`tests/quaids_synthetic_validation_test.e`,
  22 checks), which surfaced a real numerical-reliability finding (the
  iterative estimator fails to converge cleanly for roughly half of random
  seeds in this DGP family). Tolerance policy documented. Published
  replication against `Blanciforti86` (Blanciforti, Green & King 1986;
  bundled in R's `micEconAids`), committed with repo-owner approval, plus a
  cross-implementation comparison against R (close agreement, ~0.021 max
  abs difference) and Python (broadly consistent, documented caveat on one
  equation) — `tests/quaids_published_validation_test.e`, 11 checks. This
  comparison **found and fixed a real correctness bug**: `quaidsFit()`'s
  Stone-index starting value used a mutated (partly relative, partly
  absolute) price matrix, silently producing wrong results for every
  `aCtl.maxiter==1` (LA-AIDS) call with default starting values. See the
  Milestone 3 entry below for the full writeup — this fix changes numerical
  output relative to every prior milestone's frozen baseline, intentionally
  and correctly.
- Milestone 4: two new standalone Wald tests, `quaidsHomogeneityTest()` and
  `quaidsJointTest()` (`src/quaidstests.src`), both validated for size and
  power (not just "it runs"). Strengthened the existing symmetry-given-
  homogeneity test with a power check, and exercised the overidentification
  test for the first time ever in this repo (every prior fixture was
  exactly identified, `ninst==nu`) — `tests/quaids_hypothesis_tests_test.e`,
  19 checks. The first implementation attempt at the homogeneity test was
  wrong (misread the `quaidsFit()` docstring) and was caught immediately by
  the size check rejecting a true null — see the Milestone 4 entry below.
- Milestone 5: split `quaidsElas()` into `quaidsElasFit()` (silent,
  struct-returning) and `printQuaidsElas()`, mirroring the Milestone 1
  `quaidsFit()`/`printQuaids()` split — verified byte-for-byte identical
  printed output. Validated correctness away from the four standard
  evaluation points (mean/Q1/median/Q3) using three *exact* algebraic
  identities (Engel aggregation, Cournot aggregation, elasticity
  homogeneity), not tolerance-based approximations —
  `tests/quaids_elasticities_test.e`, 17 checks. `quaidsSlutzky()` needed
  no change (already accepted arbitrary points). Curvature imposition
  scoped and explicitly deferred — even the reference implementation used
  for Milestone 3 validation only diagnoses curvature, doesn't impose it.
- Milestone 6: `src/pubtable_quaids.src`, an optional `pubtable` adapter
  (`ptModelFromQuaids`/`ptFromQuaids` for coefficient tables,
  `ptModelFromQuaidsElas`/`ptFromQuaidsElas`/`ptTablesFromQuaidsElas` for
  elasticity tables, `ptFromQuaidsFamily` dispatcher) enabling
  LaTeX/Markdown/CSV export, following the `pubtable_qardl.src` pattern —
  see the Milestone 6 entry below for why it lives in this repo's own
  `src/` rather than physically bundled inside the installed `pubtable`
  package the way `pubtable_qardl.src` is. Added a `#ifndef`/`#define`
  include guard to `src/quaids.sdf` (`QUAIDS_SDF_INCLUDED`) so the adapter
  has a symbol to detect. Found and worked around a real, empirically-
  verified type-system gap (no builtin converts a legacy character-matrix
  name vector to `pubtable`'s native `string`/`string array` fields) and a
  real row-count bug (the IV control-function residual coefficient rows
  were missing from the first draft's row names) — both caught by running
  the adapter against a real `quaidsFit()` result, not by inspection.
  `tests/quaids_pubtable_test.e`, 30 checks, including an end-to-end export
  smoke test that reads the generated files back.
- Milestone 7: package build/release tooling (`scripts/build_lcg.ps1`,
  `scripts/build_package.ps1`, `scripts/verify_release_artifact.ps1`,
  `scripts/run_release_verification.ps1`, `tests/verify_package_manifest.ps1`,
  `tests/run_source_tests.ps1`) and an installed-package public API gate
  (`tests/package_public_api.e`, run via `library quaids;` against a real
  install at `c:\gauss26\pkgs\quaids`), adapted from `gauss-qardl`.
  `CHANGELOG.md` added, reconstructing the 0.1.0-0.5.0 history. Found and
  fixed 3 real pre-existing/newly-introduced bugs by actually building,
  installing, and running the package (not by re-reading the scripts) —
  see the Milestone 7 entry below, including a genuinely pre-existing
  dead-but-accidentally-live `quantile()` duplicate in `src/quaids.src`
  dating from before Milestone 0, invisible to every `#include`-based test
  but not to `library`-based loading.

`docs/` still does not exist — that is Milestone 8.

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

### Milestone 2 — Modular Source Split + Formula/Dataframe Entry Point — COMPLETE (2026-07-19)

**Scoping decision, stated explicitly**: the checklist below names six
things to split into separate files ("core estimation, IV first stage,
homogeneity/symmetry testing, elasticities, Slutzky diagnostics, printing").
Elasticities, Slutzky, printing, and the IV first stage were split out.
Homogeneity/symmetry testing was **not** split out of `quaidsFit()`'s core,
because unlike the IV first stage it is not cleanly separable: the starting
values, iteration loop, variance computation, overidentification test, and
symmetry test/symmetry-constrained stage all share heavily mutated
intermediate state (`m`, `gg`, `gw`, `ng`, and iteration-final `_beta`/
`lambda`/`lx`/`b_p`/`lx2`, which the variance and overidentification-test
formulas both need). Forcing a split now would mean either a long, brittle
inter-proc parameter list or inventing an "in-progress fit state" struct —
real maintainability value, but not urgent on pre-alpha code with no other
consumers, and safer to attempt after Milestone 3 gives this a validation
harness to catch regressions in a deeper refactor. This matches
`gauss-qardl`'s own stated "Roadmap Rules": build fixture infrastructure
before a broad rewrite, prefer small releasable increments. Full reasoning
in `CLAUDE.md`.

- [x] Split `src/quaids.src` into focused files: **done for IV first stage**
  (`quaidsiv.src`, private `_quaidsIVFirstStage()`), **elasticities**
  (`quaidselas.src`), and **Slutzky diagnostics** (`quaidsslutzky.src`).
  **Deferred** for "core estimation, homogeneity/symmetry testing" — see
  scoping decision above; these plus printing (`printQuaids()`) and the
  legacy wrapper (`quaids()`) remain in `quaids.src` as one cohesive unit,
  matching how `gauss-qardl`'s own `qardl.src` bundles `qardl()`/
  `qardlECM()`/`plotQARDL()`/etc. together while giving genuinely distinct
  features (`icmean.src`, `wtestlrb.src`, `ardlbounds.src`, `qirf.src`)
  their own files.
- [x] Evaluate routing the IV first-stage regression through `gmmFitIV`.
  **Decision: not adopted.** `gmmFitIV` is single-equation and, under
  `"onestep"`/`"unadj"` settings, mathematically identical to the classical
  2SLS already computed — no accuracy difference — but it does not expose
  the raw `zzi`/moment-matrix building blocks that `quaidsFit()`'s
  downstream system covariance and overidentification-test formulas need in
  a specific layout. Adopting it would add a package dependency for zero
  net simplification. Full reasoning in `CLAUDE.md`.
- [x] Add a dataframe entry point: `quaidsFull(data, shareVars, priceVars,
  totexpVar, instrVars, extraVars, aCtl)` (`src/quaidsformula.src`).
  **Not** a `"y ~ x1 + x2"` formula string, by design: AIDS/QUAIDS is a
  multi-equation system (N shares against N parallel prices) with no
  natural single-equation-formula representation. Column-name string
  arrays, matched by position, are the fit instead — documented in
  `quaidsformula.src`'s header and in `CLAUDE.md`.
- [x] Add formula-vs-matrix parity tests: `tests/quaids_formula_parity_test.e`,
  17 checks, including the `extraVars == 0` path — see Milestone 2
  Verification below.

#### Milestone 2 Verification

1. **Compile check** — all seven source files together. **PASS.**
2. **Behavioral parity (`quaids()` legacy wrapper)** — same procedure as
   Milestones 0/1: fixed-seed synthetic dataset, diffed against the
   untouched Milestone-0-era baseline. **PASS — zero differences**, run
   after the file split and again after the two bug fixes below.
3. **Schema test** (`tests/quaids_schema_test.e`, 34 checks) — re-run after
   the file split. **PASS.**
4. **New: formula parity test** (`tests/quaids_formula_parity_test.e`, 17
   checks) — builds the identical synthetic dataset as both plain matrices
   and a named-column dataframe (`asDF`/`dfaddcol`), estimates both ways,
   and asserts exact equality on `u`, `ivB`, `homogB`, `homogV`, `symcB`,
   `symStat`, and the final `b`/`v`/`bS`/`vS`, plus the `extraVars == 0`
   path against `quaidsFit(..., intcpt=0, ...)`. **PASS — 17/17**, after
   two real bug fixes (below).

**Two pre-existing bugs found and fixed**, both in the `intcpt == 0` branch
of `quaidsFit()`'s name-setup block — a branch that existed unchanged since
the original `aids_rev.src` but had never been exercised by any test or
example in this repo, because the synthetic fixture always passes a
non-zero `intcpt`:

1. `xnam` was read uninitialized (`G0152`) when `intcpt == 0`, since it was
   only assigned in the `else` branch.
2. After fixing #1 with a native string-literal assignment, assigning that
   into the `matrix`-typed `qOut.xnam` field threw `G0071 Type mismatch` —
   same class of bug as the Milestone 1 `string array` vs. `matrix` issue,
   different specific cause (native `string` literal vs. the legacy
   `0$+"X"`-built character-matrix type). Fixed with the `0$+` idiom.

Both fixes are scoped to the previously-dead branch; the already-tested
`else` branch was untouched, and re-running the Milestone 0/1 parity and
schema tests after the fix confirmed zero change to already-verified
behavior. See `CLAUDE.md` for the full GAUSS type-system notes this
surfaced (`type()` codes for matrix/string/string array, dataframe column
selections being plain type-6 matrices, and where GAUSS does vs. doesn't
coerce between the legacy character-matrix and native string types).

**Standard applied**: same as Milestones 0/1 — behavior-preserving changes
(file moves, the IV extraction) had to produce zero-diff output; genuinely
new code (`quaidsFull()`) had to be verified against the existing matrix API
on identical data, not just "ran without error." The two bugs found here
are a direct product of that standard: a parity test that actually exercises
a previously-dead code path is what caught them, not code review.

### Milestone 3 — Validation Fixture and Benchmark Harness — COMPLETE (2026-07-20)

- [x] Add deterministic synthetic fixtures (expand on `quaids_example.e`
  with actual pass/fail assertions instead of eyeballed printouts),
  covering LA-AIDS, iterated AIDS, and QUAIDS, each with and without IV.
  `tests/quaidsfixtures.src` (`_quaidsSyntheticDGP()`, a shared 5-good
  homogeneity/adding-up-true-by-construction generator, parameterized by
  quadratic-term and endogeneity switches) plus
  `tests/quaids_synthetic_validation_test.e` (22 checks across all 6
  model/endogeneity combinations). See "Seed sensitivity finding" below —
  this surfaced a real numerical-reliability issue in the iterative
  estimator, not just a green checkmark.
- [x] Identify a published-replication target, with explicit repo-owner
  approval before committing external data (given 2026-07-20). Used
  `Blanciforti86` (annual U.S. food-consumption data, 1947–1978, 4 food
  groups) bundled in the R package `micEconAids` (Arne Henningsen, GPL ≥ 2
  on CRAN), sourced from Blanciforti, Green & King (1986), *U.S. Consumer
  Behavior Over the Postwar Period: An Almost Ideal Demand System
  Analysis*, Giannini Foundation Monograph No. 40. Committed as
  `tests/fixtures/published/blanciforti86_food32.csv`, attribution and
  license note in `tests/fixtures/published/SOURCE.md`.
- [x] Add a cross-implementation comparison against R and Python (given
  explicit go-ahead to install both, 2026-07-20). R: installed
  `micEconAids` (CRAN, binary packages, no compilation needed) and ran
  `aidsEst(..., instNames=...)` (3SLS) with the same identification
  strategy GAUSS uses (instrument `log(xFood)` with `log(xAgg)`).
  **Result: close agreement, max abs difference ≈0.021** across
  alpha/beta/gamma, after fixing a real bug this comparison surfaced (see
  below). Python: hand-coded (no comparably-established Python AIDS
  package exists) independent replica of the same specification;
  broadly consistent but with larger residual differences on one equation,
  attributed to the from-scratch replica rather than to GAUSS since R
  agrees closely with GAUSS on exactly the coefficients where the Python
  replica diverges most (see `tests/fixtures/published/
  python_reference_check.py`'s header for the full accounting). Committed
  as `tests/quaids_published_validation_test.e` (11 checks, R numbers used
  as the hard assertion target; Python kept as documented supplementary
  evidence, not an assertion source, given the above caveat). Reference
  scripts (`generate_r_reference.R`, `python_reference_check.py`) are
  committed alongside the data for reproducibility.
- [x] Document tolerance policy for deterministic expected-output tests —
  see "Tolerance Policy" below.

#### Real bug found and fixed: Stone-index starting value used the wrong price matrix

The published-data cross-check (GAUSS vs. R, both on `Blanciforti86`, same
identification strategy) initially disagreed by an order of magnitude —
`beta` off by roughly 5x, not the few-percent gap expected between two
different-but-valid IV algorithms. Root cause, in `quaidsFit()`'s
"STARTING VALUE" block (`src/quaids.src`): `stone = prices*meanc(w)` was
applied to `prices` *after* it had already been converted to relative form
in columns `1:n-1` (each minus the reference good's price) while column
`n` stayed absolute — a mutation done earlier, for the homogeneity
reparametrization. Weighting that mixed relative/absolute matrix by mean
shares does not compute the Stone price index; algebraically (verified by
direct derivation, not just observation) it computes
`StandardStoneIndex − ln(p_n)·(1 − meanShare_n)`, a distortion that tracks
the reference good's own price trend rather than a valid deflator.

**Impact**: since `aCtl.maxiter == 1` (LA-AIDS, the Stone-index one-step
model) never iterates past this starting value, the distorted deflator
*was* the final answer for every LA-AIDS call using default (`aCtl.b0==0`)
starting values — not an occasional glitch. For `aCtl.maxiter > 1`
(iterated AIDS/QUAIDS), this only supplied a bad *starting point*; the
correct nonlinear `a(p)` formula used inside the iteration loop itself was
never affected. This likely also explains part of why the Milestone-3
seed-sensitivity probe (below) found so many non-converging seeds: a
materially-wrong starting point makes convergence failure more likely, on
top of whatever intrinsic conditioning issues remain.

**Fix**: reconstruct absolute prices before computing `stone`
(`(prices[.,1:n-1] + prices[.,n])~prices[.,n]`) rather than changing
anything about the relative-price convention used elsewhere in the proc —
confirmed correct by an isolated Python check (patch the formula, rerun,
watch the gap with R close from ~5x-off to matching within normal
cross-implementation noise) before touching `src/quaids.src`.

**This changes numerical output** for `aCtl.maxiter==1` calls (and, to a
lesser extent, shifts the iteration path — not necessarily the converged
answer — for `aCtl.maxiter>1` calls) relative to every prior milestone's
frozen baseline. That is expected and correct: those baselines were
captured from the original, buggy `aids_rev.src`. **The Milestone 0/1/2
"verified byte-for-byte identical" claims elsewhere in this document remain
true as historical statements about those specific milestones (structure-
preserving refactors of the code as it existed then) — they are not claims
that current output matches that old baseline anymore, and it
intentionally no longer does for `aCtl.maxiter==1`.** `examples/
quaids_example.e` (seed 11, `aCtl.maxiter=100`) was already one of the
non-converging seeds identified below even before this fix, so its
iteration path and final numbers changed too; this was not re-tuned to a
better-behaved seed as part of this fix, since doing so wasn't necessary to
validate the fix itself.

#### Seed sensitivity finding (numerical reliability)

While calibrating the synthetic fixtures' tolerances, a multi-seed probe
(8 seeds, `tobs=3000`, `aCtl.err=.0001`, `aCtl.maxiter=100`) found that the
**iterative estimator (QUAIDS and iterated linear AIDS) fails to converge,
or converges to numerically nonsensical estimates (errors of magnitude
200–2500 against true parameters of magnitude ~0.1–2), for roughly half of
random seeds** in this DGP family — independent of whether the model has a
quadratic term or genuine endogeneity (the same seeds failed or succeeded
across both `quadratic`/`endogenous` settings, holding `prices` fixed for a
given seed, which points at price-draw-dependent conditioning of the
iteration rather than something specific to QUAIDS or IV). The synthetic
fixtures deliberately use `seed=204`, one of the seeds confirmed to converge
cleanly across all six model/endogeneity combinations, and that choice is
called out in `quaids_synthetic_validation_test.e`'s own comments — this is
not silently cherry-picked. This finding itself is valuable Milestone-3
output: it identifies a real robustness gap (the iterated FGLS estimator
has no globally-convergent guarantee and no fallback/damping for bad
starting points) worth a dedicated numerical-reliability pass once more of
the roadmap is built out — likely as part of hardening work analogous to
`gauss-qardl`'s Milestone 13. Not fixed here beyond the one confirmed
starting-value bug above; a full numerical-reliability pass (e.g. damped
iteration, multiple starting points, convergence diagnostics surfaced to
the caller) is out of scope for a validation milestone.

#### Seed sensitivity finding (numerical reliability)

While calibrating the synthetic fixtures' tolerances, a multi-seed probe
(8 seeds, `tobs=3000`, `aCtl.err=.0001`, `aCtl.maxiter=100`) found that the
**iterative estimator (QUAIDS and iterated linear AIDS) fails to converge,
or converges to numerically nonsensical estimates (errors of magnitude
200–2500 against true parameters of magnitude ~0.1–2), for roughly half of
random seeds** in this DGP family — independent of whether the model has a
quadratic term or genuine endogeneity (the same seeds failed or succeeded
across both `quadratic`/`endogenous` settings, holding `prices` fixed for a
given seed, which points at price-draw-dependent conditioning of the
iteration rather than something specific to QUAIDS or IV). The synthetic
fixtures deliberately use `seed=204`, one of the seeds confirmed to converge
cleanly across all six model/endogeneity combinations, and that choice is
called out in `quaids_synthetic_validation_test.e`'s own comments — this is
not silently cherry-picked. This finding itself is valuable Milestone-3
output: it identifies a real robustness gap (the iterated FGLS estimator
has no globally-convergent guarantee and no fallback/damping for bad
starting points) worth a dedicated numerical-reliability pass once more of
the roadmap is built out — likely as part of hardening work analogous to
`gauss-qardl`'s Milestone 13. Not fixed here; fixing it would mean changing
estimation numerics, out of scope for a validation milestone.

#### Tolerance Policy

Deterministic synthetic-fixture tests in this repo (`tests/*_test.e`)
follow two different standards depending on what's being checked:

- **Structure-preserving refactors** (Milestones 0–2: renames, file moves,
  the `quaidsFit()`/`printQuaids()`/`quaids()` split, the IV-stage
  extraction) are held to **exact/byte-identical equality** against a
  frozen prior-behavior baseline. Nothing here should introduce numerical
  drift, so any difference at all is a fail. This is what
  `examples/quaids_example.e`'s parity check and
  `quaids_formula_parity_test.e`'s cross-API checks do
  (`maxc(maxc(abs(a-b))) == 0`). This standard applies to *structural*
  changes, not correctness fixes: the confirmed Stone-index bug fix
  documented above intentionally breaks byte-parity against the
  pre-fix/pre-Milestone-0 baseline for `aCtl.maxiter==1` calls, because the
  old baseline was itself wrong. After a correctness fix, future structural
  refactors should be checked for byte-parity against a freshly captured
  *post-fix* baseline, not the old one.
- **Statistical recovery fixtures** (Milestone 3 synthetic validation) use
  an **absolute-error tolerance calibrated against an observed multi-seed
  probe, not a guessed number** — see `quaids_synthetic_validation_test.e`.
  Structural coefficient rows: `0.10`. The IV-residual coefficient row
  (consistently ~10x noisier than every other row across every
  model/endogeneity combination tested — see the seed-sensitivity finding
  above): `0.50`. LA-AIDS (Stone index, one-step): `1.20` throughout,
  looser because Stone-index approximation bias is a real, expected
  property of that method, not test slack. Sample size (`tobs=3000`) and
  seed (`204`) are fixed and documented, not swept — this is a regression
  guard against the *implementation* drifting from correct behavior on a
  known-good case, not a power analysis of the estimator.
- Tolerances are **not** meant to certify general-purpose statistical
  accuracy across arbitrary data — that is what the deferred
  published-replication and cross-implementation comparisons above are for.

### Milestone 4 — Hypothesis Testing Completeness — COMPLETE (2026-07-20)

- [x] Add a standalone homogeneity Wald/LR test (unrestricted vs.
  homogeneity-constrained), not just the qualitative "near-zero reference
  price effect" signal. New `quaidsHomogeneityTest(qOut)`
  (`src/quaidstests.src`), a Wald test on an unconstrained
  (`aCtl.homogenous=0`) fit's recovered gamma matrix: `df = n-1`. Validated
  for both **size** (fails to reject, p=0.63, on a homogeneity-true-by-
  construction fixture) and **power** (rejects, p≈0, on a fixture with a
  deliberate, clean homogeneity violation injected into the observed
  shares) — see "New hypothesis tests" below for how the exact restriction
  vector/covariance were derived and why the first attempt was wrong.
- [x] Add a joint homogeneity+symmetry test. New `quaidsJointTest(qOut)`
  (`src/quaidstests.src`), same unconstrained-fit input, `df = (n-1) +
  (n-1)(n-2)/2`. Also validated for size and power.
- [x] Formalize the existing symmetry-given-homogeneity test and the
  existing overidentification test as independently tested (not
  necessarily independently *callable* as separate procs — see the scoping
  note below, which extends Milestone 2's coupling finding).
  `tests/quaids_hypothesis_tests_test.e` adds: a **power** check for the
  symmetry-given-homogeneity test (previously only implicitly
  size-checked; now confirmed to reject on a deliberately asymmetric-but-
  homogeneous fixture, p≈0), and the **first-ever exercise of the
  overidentification test** in this repo's history — every fixture through
  Milestone 3 used exactly-identified instruments (`ninst==nu`), so
  `qOut.overidValid` was always `0` and that branch had literally never
  run. A new 2-instrument fixture confirms it runs, has the right shape/df,
  and doesn't spuriously reject when both instruments are valid by
  construction.
- [x] Document degrees of freedom and asymptotic assumptions for each test
  — see "New hypothesis tests" below and the header comment of
  `src/quaidstests.src`.

#### Scoping note: why the existing tests weren't extracted into standalone procs

Milestone 2 declined to split `quaidsFit()`'s overidentification-test and
symmetry-test computations out of the main estimation proc, because they
depend on heavily mutated intermediate state from the iteration/variance
stages (`Ji`, `Sgma`, `S`, `O`, `D`, `zzi`, `m1`, and iteration-final
`_beta`/`lambda`/`lx`/`b_p`/`lx2`) that would otherwise need threading
through a long, brittle parameter list. That coupling is unchanged by this
milestone (nothing about the estimation core moved), so the same reasoning
applies here: those two tests remain computed inline as part of
`quaidsFit()`, exposed as `qOut` fields, not as separate callable procs.
The two genuinely *new* tests (`quaidsHomogeneityTest`, `quaidsJointTest`)
**are** standalone, independently callable procs — they only need the
*already-finished* `qOut.b`/`qOut.v` from an unconstrained fit, not any
of that intermediate iteration state, so there was no coupling problem to
work around for them.

#### New hypothesis tests: derivation and validation

Both new tests are Wald tests built the same way: extract a linear
combination `L` such that `L'*vec(qOut.b)` is the vector of restrictions
being tested, use `L'*qOut.v*L` as its covariance (`qOut.b`/`qOut.v` here
must come from an **unconstrained**, `aCtl.homogenous=0` fit — the final,
absolute-price-form gamma matrix and its covariance), and compute
`stat = r'*inv(V)*r ~ chi2(df)`.

- **Homogeneity** (`df = n-1`): the restriction is `sum_j gamma_ij = 0` for
  each independently-estimated equation `i = 1..n-1` (equation `n` is
  recovered from the others via adding-up and contributes no new
  information to a Wald test, matching the `n1 = n-1` convention used
  throughout this library).
- **Joint** (`df = (n-1) + (n-1)(n-2)/2`): homogeneity's `n-1` restrictions
  plus symmetry restrictions `gamma_ij = gamma_ji` for `i<j`,
  `i,j = 1..n-1`. Note a **symmetric** gamma matrix with adding-up already
  imposed automatically satisfies homogeneity too (row `i` sum = `sum_j
  gamma_ij` = `sum_j gamma_ji` [symmetry] = column `i` sum = 0 by
  adding-up) — so there is no separate "symmetry only, given adding-up, on
  a fully unconstrained fit" test offered; test joint, or use the existing
  symmetry-given-homogeneity test if homogeneity itself is not in
  question.

**The first implementation attempt was wrong**, caught by validating size
before trusting the formula: the docstring for `quaidsFit()` describes an
internal *reparametrized* representation (used mid-computation, before the
"recovers absolute price effects" step) where a row sum stands in for one
of the raw price coefficients. Reading that as describing `qOut.b`
directly (i.e., assuming a single row already *was* the row-sum statistic)
gave a homogeneity test that rejected with `p ≈ 0` on a fixture where
homogeneity is true by construction — a dead giveaway of a wrong formula,
not a real finding. Re-deriving from the actual recovery code (`qOut.b`
holds the raw, non-reparametrized gamma post-recovery) and summing
*across* the price rows for each equation, rather than reading a single
row, fixed it: `p = 0.63` on the true-null fixture, `p ≈ 0` on a
deliberately-violated one. Constructing that violated fixture also took
two attempts — a hand-built asymmetric gamma matrix turned out to be
internally inconsistent (the price aggregator `a_p` uses a quadratic form
`p'*Γ*p`, which is mathematically invariant to using `Γ` or `Γ'`, so an
asymmetric "true" gamma silently degenerates to its own symmetrized
average inside `a_p` while the share equations still see the full
asymmetric version — not wrong, just not what was intended, and reflected
in bad recovery). Injecting a clean, explicit violation directly into the
observed shares instead (`w[.,1] += c*prices[.,1]`, `w[.,2] -=
c*prices[.,1]`, preserving adding-up exactly) sidestepped that and gave a
fixture that worked as intended. Both dead ends are a demonstration of why
size *and* power checks matter — a formula that merely "runs" proves
nothing.

### Milestone 5 — Elasticities and Diagnostics Generalization — COMPLETE (2026-07-20)

- [x] Generalize `elas()`/`elas_()` to accept arbitrary evaluation points
  (currently hardcoded to mean/Q1/median/Q3 inside `quaids()`).
  **Re-scoped based on what the code actually did**: `quaidsElas_()`
  already accepted any point as an argument (`intcpt`/`prices`/`totexp` are
  point values, not sample statistics) — the real gap was that
  `quaidsElas()` mixed computation with printing, and `quaids()` only ever
  called it at four fixed points, so there was no clean, silent,
  struct-returning way to ask for elasticities anywhere else. Split
  `quaidsElas()` into `quaidsElasFit()` (silent, returns `quaidsElasOut`:
  point estimates + delta-method standard errors, no printing) and
  `printQuaidsElas()` (the separated printer), mirroring the
  `quaidsFit()`/`printQuaids()` split from Milestone 1. `quaidsElas()`
  itself is now a thin wrapper (`quaidsElasFit()` then
  `printQuaidsElas()`) — unchanged signature, **verified byte-for-byte**
  identical printed output on `examples/quaids_example.e`.
  `tests/quaids_elasticities_test.e` (17 checks) validates the split
  (parity against `quaidsElas_()` directly) and, more importantly,
  validates *correctness away from the four standard points* using three
  **exact** algebraic identities (Engel aggregation, Cournot aggregation,
  elasticity homogeneity — consequences of adding-up/homogeneity holding,
  not approximations) at a real out-of-sample observation and at a fully
  synthetic counterfactual price scenario. All held to floating-point
  precision (~1e-16).
- [x] Keep the Slutzky-eigenvalue negativity diagnostic as the default
  always-on check. Unchanged — `quaidsSlutzky()` already accepted an
  arbitrary `intcpt`/`prices`/`totexp` sample (any number of rows) as
  input, so it was already general in the same sense the elasticities
  functions needed to become; no code change was needed here, and
  `quaids()` still calls it unconditionally. Confirmed unaffected by the
  full regression suite.
- [x] Scope and, if justified, implement optional local curvature imposition
  (Diewert-Wales Cholesky reparametrization) as an opt-in estimation mode
  built on `optmt`/`cmlmt` — P2, not a release blocker. **Scoped, not
  implemented.** Even the reference implementation used for Milestone 3's
  validation, R's `micEconAids`, only offers curvature *diagnosis*
  (`aidsMono()`, `aidsConcav()` — check monotonicity/concavity post-hoc,
  matching what `quaidsSlutzky()` already does here) and does not offer
  curvature *imposition* as a constrained-estimation mode. That weakens the
  case for this GAUSS library needing to leap ahead of the reference
  implementation on a P2, "if justified" item. Revisit if a concrete use
  case emerges; not pursued now.

### Milestone 6 — Reporting via `pubtable` — COMPLETE (2026-07-20)

- [x] Add `pubtable_quaids.src` following the `pubtable_qardl.src` pattern
  (`#ifDef QUAIDS_SDF_INCLUDED` guard, `ptModelFromQuaids`, `ptFromQuaids`).
- [x] Route elasticity tables and coefficient tables through `pubtable`
  model objects in addition to the console `printfm()` output.
- [x] Add LaTeX/Markdown/CSV export examples.

**Where the adapter lives, and why that diverges from the `pubtable_qardl.src`
precedent**: `pubtable_qardl.src` is not part of the `gauss-qardl` git repo
at all — it is bundled *inside the installed pubtable package itself*
(`c:\gauss26\pkgs\pubtable\src\pubtable_qardl.src`, listed in pubtable's own
`package.json` `src` array alongside its Aptech-authored `cmlmt`/`maxlikmt`/
`optmt`/`tsmt` adapters). Physically matching that precedent for
`pubtable_quaids.src` would have meant writing into a shared, installed
package outside this repo's git history — a change that would affect every
other project on this machine that loads `pubtable`, not something
reversible by a normal `git` operation in `gauss-aids-model`. Given the
choice (explicitly asked of the repo owner rather than assumed), the adapter
instead lives at `src/pubtable_quaids.src`, self-contained and git-tracked
in this repo, following the same "no file self-includes another" convention
already used by every other file in `src/` (see the "GAUSS language
conventions" section of `CLAUDE.md`) — a caller `#include`s `quaids.sdf`,
`quaids.src` (or at least `quaidselas.src`), `pubtable.sdf`/`pubtable.src`
(or `library pubtable;`), and finally `pubtable_quaids.src`, in that order.
It matches the *pattern* (naming, `#ifDef QUAIDS_SDF_INCLUDED` guard,
`ptModelFromX`/`ptFromX` function shapes, `dynargsGet`) exactly, just not
the physical file location. Not added to `package.json`'s `src` array,
since (unlike every other file there) it has a hard compile-time dependency
on `pubtable.sdf`'s `ptModel`/`ptTable` struct types (the proc return-type
annotations `proc (struct ptModel) = ...` are outside the `#ifDef` guard,
so they need those types declared regardless of whether `QUAIDS_SDF_INCLUDED`
is defined) — adding it to `src` would make `pubtable` a hard dependency for
the whole package to even compile, contradicting `package.json`'s empty
`deps` array. This mirrors how `tests/fixtures/published/generate_r_reference.R`/
`python_reference_check.py` are real, working, documented parts of this
repo without being part of the installable package.

**A required include-guard addition**: `pubtable_quaids.src`'s `#ifDef
QUAIDS_SDF_INCLUDED` branches need that symbol defined by something —
`quaids.sdf` did not previously guard itself the way `qardl.sdf` does
(`#ifndef QARDL_SDF_INCLUDED` / `#define QARDL_SDF_INCLUDED` / ... /
`#endif`). Added the same guard to `src/quaids.sdf` (harmless for every
existing caller, since none of them included it twice) specifically so
`pubtable_quaids.src` has a symbol to test.

**A real, empirically-verified type-system finding — not just an analogy
to the Milestone 1 lesson**: `quaidsOut`/`quaidsElasOut` name vectors
(`xnam`, `wnam`, `unam`, ...) are legacy character matrices (GAUSS type 6),
but `pubtable`'s `ptModel.termNames`/`ptTable.rowNames`/`colNames` are
natively typed `string array` (type 15), and `ptModel.name`/`ptTable.title`
are scalar native `string` (type 13). CLAUDE.md already documented that
assigning a bare native `string` into a `matrix`-typed struct field throws
`G0071 Type mismatch`; this milestone needed the *reverse* direction (a
legacy char-matrix into a `string array`/`string`-typed field), which is a
different code path and was not previously verified — so it was tested
directly with `tgauss`, not assumed. Confirmed empirically: direct
assignment throws `G0071` there too. **There is no dedicated conversion
builtin** (`strtrim()` was tried first, since it looked like the obvious
candidate — it errors with "Invalid argument type" on a legacy char-matrix
input). The working idiom, also confirmed empirically: concatenating a
legacy char-matrix with a native string via `$|` forces element-wise
conversion to a native string array (`cm $| ""`, then slice off the trailing
blank row: `(cm $| "")[1:rows(cm)]`); indexing a *single* row of the result
with a scalar index (`sa[1]`, not `sa[1:1]`) further demotes it to a true
scalar `string` (type 13) for `ptModel.name`/`ptTable.title`. Both idioms
are `_ptQuaidsToStrArray()`/`_ptQuaidsToStr()` in `pubtable_quaids.src`.

**A real row-count bug found by running against a real fit, not by
inspection**: the first draft of `ptModelFromQuaids()` built row names from
only the "structurally interpretable" blocks (`qOut.xnam` | `GAMMA_`+
`qOut.wnam` | `BETA_LX` | `LAMBDA_LX2`), mirroring how `quaidsElas_()` reads
`b` — but `quaidsElas_()` only *reads* the first `1+nint+n+nendog` rows of
`b` and silently ignores the rest; it does not mean `qOut.bestB` has only
that many rows. `qOut.bestB`/`qOut.bS` always have `nu` additional trailing
rows for the IV control-function residual coefficient(s) (`qOut.unam`,
always `nu >= 1` since `quaidsFit()` always treats log total expenditure as
endogenous — see `qOut.ng = 1+nint+n1+nendog+nu` in `src/quaids.src`).
Running the adapter against a real `quaidsFit()` QUAIDS/IV result
immediately threw `ptModelSetNames: termNames must contain 10 labels` (9
built vs. 10 actual) — caught by testing against real output, not by
re-reading the formula. Fixed by appending `qOut.unam`-derived names for
the trailing block; a coefficient *report* arguably should show the fitted
residual coefficients anyway, since `printQuaids()`'s own console table
already reports them (the "Residuals of instrumental regressions" row).
Re-verified against LA-AIDS (`linear=1`), QUAIDS (`linear=0`), and both
`homogenous=1`/`homogenous=0` fits to confirm the row-count formula holds
across all four combinations, not just the one that first exposed the bug.

**What was built**:
- `src/pubtable_quaids.src`: `ptModelFromQuaids(name, qOut, eqIdx)` /
  `ptFromQuaids(qOut)` (one coefficient column per good via
  `ptModelCompare`, mirroring `pubtable_qardl.src`'s per-quantile
  `ptFromQardl`) for `quaidsOut`; `ptModelFromQuaidsElas(name, elasOut)` /
  `ptFromQuaidsElas(elasOut)` (income elasticities) and
  `ptTablesFromQuaidsElas(elasOut)` (3-table bundle: income elasticities,
  uncompensated/Marshallian price elasticities, compensated/Hicksian price
  elasticities — the latter two built directly as `ptTable`s via the
  internal `_ptQuaidsElasMatrixTable()` helper, since an `n x n` matrix
  with a value row and an `(se)` row per good does not fit `ptModel`'s
  single-coefficient-vector shape) for `quaidsElasOut`; and
  `ptFromQuaidsFamily(x)` (an `isStructType`-based dispatcher, mirroring
  `ptFromArdlFamily`).
- `examples/pubtable_export_example.e`: manual, eyeball-comparison example
  (no assertions, matching `quaids_example.e`'s own style) exporting a
  coefficient table and all three elasticity tables to `.tex`/`.md`/`.csv`.
- `tests/quaids_pubtable_test.e` (30 checks): exact numeric parity between
  `ptModel.estimates`/`stdErrors` and the `qOut.bestB`/`qOut.bestV`/
  `elasOut.er`/`elasOut.ser` values they're built from (not just "it
  runs"); row/column shape and title checks; the dispatcher tested against
  both struct types; and an end-to-end export smoke test that writes real
  `.tex`/`.md`/`.csv` files, reads them back, and checks their content
  (booktabs `\begin{tabular}`, Markdown pipes, CSV commas, and specific
  row-label substrings) — not just that `ptExport()` returns without
  error.
- `.gitignore`: added the specific generated export filenames from both
  the example and the test (not a blanket `*.tex`/`*.md`/`*.csv` pattern,
  which would have also hidden `CLAUDE.md`/`GOLD_STANDARD_TODO.md` and the
  committed `tests/fixtures/published/blanciforti86_food32.csv`).

**Requires `pubtable` installed** (this machine has it at
`c:\gauss26\pkgs\pubtable`) to run the new example/test — same "installed
locally, not a repo dependency" status as R/Python for the Milestone 3
published-data cross-checks. `package.json`'s `deps` stays empty; this is
optional, activated-by-inclusion functionality, not a hard package
dependency.

### Milestone 7 — Package Build and Release Tooling — COMPLETE (2026-07-20)

- [x] `package.json` (name, version, `src` array, deps, keywords, license)
  matching `dccelib`/`qardl` conventions — already existed since Milestone
  0; re-verified consistent.
- [x] `tests/run_source_tests.ps1`-style runner and a package-manifest
  consistency check (`tests/verify_package_manifest.ps1`).
- [x] Repeatable package build/install scripts (`scripts/build_lcg.ps1`,
  `scripts/build_package.ps1`, `scripts/verify_release_artifact.ps1`,
  `scripts/run_release_verification.ps1`) and an installed-package smoke
  test (`library quaids;` public API gate, `tests/package_public_api.e`),
  matching `qardl`'s `tests/package_public_api.e`.
- [x] `CHANGELOG.md`, reconstructing the 0.1.0-0.5.0 version history from
  `GOLD_STANDARD_TODO.md`'s own milestone records. **Not tagging a git
  release as part of this milestone** — nothing in this repo has been
  committed yet (every milestone's changes, 0 through 7, are staged but
  uncommitted, per the "never commit unless asked" policy this whole
  engagement has followed), and a tag requires a commit to point at.
  Version-numbering infrastructure (`CHANGELOG.md`, `package.json`) is in
  place; the actual `git commit`/`git tag` is a repo-owner decision for
  whenever they choose to commit this work.

Scripts adapted from `gauss-qardl`'s `scripts/`/`tests/` tooling
(`build_lcg.ps1`, `build_package.ps1`, `verify_release_artifact.ps1`,
`run_release_verification.ps1`, `verify_package_manifest.ps1`,
`run_source_tests.ps1`, `package_public_api.e`), scaled down for this
repo's smaller current scope:

- `verify_release_artifact.ps1`'s `requiredEntries` list omits
  `README.md`/`docs/COMMAND_REFERENCE.md` (Milestone 8, not built yet) —
  add them once Milestone 8 lands, noted in that script's own header
  comment.
- `run_release_verification.ps1` omits `gauss-qardl`'s separate
  new-model-benchmark/validation-benchmark/examples-smoke scripts, since
  the equivalent validation already lives inside `run_source_tests.ps1`'s
  7 tgauss test files (in particular `quaids_synthetic_validation_test.e`
  and `quaids_published_validation_test.e`).
- `run_source_tests.ps1` checks this repo's own `PASS`/`FAIL`-line and
  `ALL N CHECKS PASSED`/`N CHECKS FAILED` convention (documented in
  `CLAUDE.md`'s "Testing status" section as more reliable than `tgauss`'s
  process exit code for this harness), not just GAUSS-level compile/
  execute error text the way `gauss-qardl`'s version does.
- `tests/package_public_api.e` builds its own small inline synthetic
  dataset (the same DGP shape as `examples/quaids_example.e`) rather than
  reusing `tests/quaidsfixtures.src`'s private `_quaidsSyntheticDGP()` --
  that helper is tests/-only source, not part of the installed package,
  and the point of this test is to exercise exactly what an
  installed-package consumer actually has available.
- `pubtable_quaids.src` is deliberately not exercised by
  `package_public_api.e`: it is not in `package.json`'s `src` array (see
  Milestone 6), so `library quaids;` does not load it.

**No version bump for this milestone**: build/release tooling and
`CHANGELOG.md` don't change GAUSS public API surface (no new/changed
procs in `src/`), so per this repo's established policy (version bumps
are keyed to public API surface changes, not every milestone — see
Milestone 3's R/Python reference scripts, which also didn't bump the
version), the package stays at `0.5.0`.

**Real bugs found by actually building, installing, and running the
package — not by re-reading the scripts**, following this repo's
established validation standard (never trust a derived formula or script
without running it against a real case first):

1. **`build_package.ps1`'s cleanup step deleted the entire staged
   `examples/`/`tests/` directories**, not just the generated run
   artifacts it was meant to strip. `Get-ChildItem -LiteralPath <dir>
   -Include <patterns> -Recurse` silently ignores `-Include` when the base
   path isn't itself a wildcard — a known PowerShell footgun — so
   `-Recurse` returned every file with no effective filter, and
   `Remove-Item -Force` deleted all of them. Caught immediately by
   `verify_release_artifact.ps1` failing with "missing required entry:
   examples/quaids_example.e" the first time the script actually ran.
   Fixed by removing named files by explicit literal path (no
   `-Include`/`-Recurse` combination) plus a `-Filter` (not `-Include`)
   pass for `*.log`, which does not have this bug.
2. **`build_lcg.ps1`'s proc-detection regex only matched one of the three
   GAUSS proc-declaration forms actually used in this codebase**:
   `proc (struct X) = name(...)` (matched), but not `proc N = name(...)`
   (bare digit return count, e.g. `proc 3 = quaidsElas_(...)`) or
   `proc name(...)` (no return spec at all, e.g. `proc quantile(x, s);`).
   This silently dropped `quaidsSlutzky`, `_quaidsIVFirstStage`, the
   legacy `quaids()` wrapper, `printQuaids`, `quaidsElas_`,
   `printQuaidsElas`, `quaidsElas`, and a private `quantile` helper from
   the generated `.lcg` catalog — invisible via every source-tree
   `#include`-based test (which doesn't go through the catalog at all),
   surfaced only as `Undefined symbol` errors from `library quaids;` when
   `tests/package_public_api.e` actually called them against a real
   install. Fixed by extending the regex to match all three forms;
   verified the regenerated catalog lists every proc in every `src/` file.
3. **A genuinely pre-existing, previously-invisible bug in
   `src/quaids.src` itself**: a private `quantile(x, s)` helper
   (duplicating GAUSS's builtin `quantile()`) that the original author
   clearly intended to delete was left accidentally live, because GAUSS
   comments do not nest and the wrapping "delete this" comment's own
   doc-header used an inner `/**...**/`-style block whose closing marker
   closed the *outer* comment early. This has been silently live and
   silently duplicating the builtin's behavior at its 3 call sites (in the
   legacy `quaids()` wrapper's elasticities-at-four-points block) since
   before Milestone 0 — invisible via `#include`-based compilation (which
   just locally shadows the builtin name within that compile unit), but a
   real GAUSS builtin name cannot be redefined via `library`-based lazy
   loading, so it surfaced as `Undefined symbol: 'quantile'` resolving
   `quaids.src`'s *own* proc definition once the catalog fix above let
   `library quaids;` actually try to load it. Fixed by deleting the dead
   code outright (matching the original author's evident intent) rather
   than just repairing the comment nesting — the 3 call sites now
   correctly resolve to the GAUSS builtin, the same as
   `quaidsslutzky.src`'s identical `quantile()` calls already did. Full
   regression suite (all 7 source-tree tests) re-run afterward to confirm
   no numeric-output regression outside the (expected, and previously
   untested) legacy quartile-elasticities block.
4. A minor self-inflicted bug while writing bug #3's own explanatory
   comment: GAUSS's lexer does not tolerate literal `"` characters inside
   a `/* ... */` block comment (an odd count breaks it with `error G0097
   String not closed`, even though the characters are inside a comment,
   not a string) — caught immediately by re-running the test suite, fixed
   by rephrasing the comment without quote characters.

**A real install location decision, asked of the repo owner rather than
assumed**: fully validating `tests/package_public_api.e` requires
`library quaids;` to resolve, which means installing the built package
into GAUSS's real, shared package directory (`c:\gauss26\pkgs\quaids`) --
not touching any existing package, but still writing outside this git
repo into shared machine state, the same category of decision as the
`pubtable_quaids.src` location question at Milestone 6. Asked; repo owner
chose a real install. `c:\gauss26\pkgs\quaids` now exists as a working
installed copy of the package (verified via `library quaids;` +
`tests/package_public_api.e`), alongside `qardl` and `pubtable`.

**Verification**: `scripts/run_release_verification.ps1 -BuildArtifact
-ForceArtifact -InstallArtifact` run end to end -- all 7 source-tree tests
(150 checks) pass, the release `.zip` is built and verified (contains
every file `package.json`'s `src` array promises plus the root files/dirs
this repo currently ships), the artifact is installed to
`c:\gauss26\pkgs\quaids`, and `tests/package_public_api.e` passes against
that real installed copy via `library quaids;` (`quaidsControlCreate`/
`getDefaultQuaidsControl`, `quaidsFit`/`quaids`, `quaidsFull`,
`quaidsElasFit`/`quaidsElas`/`printQuaidsElas`, `quaidsSlutzky`,
`quaidsHomogeneityTest`/`quaidsJointTest`).

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
