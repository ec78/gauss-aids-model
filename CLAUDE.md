# CLAUDE.md — GAUSS QUAIDS Library

Context file for Claude Code sessions working on this repository.

## What this library does

Estimates **Almost Ideal Demand System** models: linearized AIDS (Stone price
index), iterated AIDS (nonlinear translog price index), and **QUAIDS**
(Banks, Blundell & Lewbel 1997 quadratic-log-expenditure extension), with
optional instrumental-variables treatment of (endogenous) log total
expenditure. Estimation imposes homogeneity and/or Slutzky symmetry via
iterated FGLS with cross-equation restrictions applied through a
minimum-distance reparametrization. Use cases: consumer demand estimation,
welfare analysis, elasticity calculation, testing demand-theory restrictions.

The library is **pre-alpha** (package version `0.5.0`) and is not yet
packaged as an installable GAUSS application package (`library quaids;` does
not work yet). See `GOLD_STANDARD_TODO.md` for the full roadmap — this file
is the quick-orientation companion to it, and should be kept synchronized
with it.

**Naming**: the package and its public procs use a `quaids` prefix (decided
at Milestone 0), even though the estimator also covers plain linear AIDS —
QUAIDS is the more general model actually implemented, and `aids` was judged
too likely to collide/confuse as a bare identifier. "AIDS"/"Almost Ideal
Demand System" remains the correct term for the model family in docs, papers,
and comments; only the GAUSS identifier prefix changed.

## Repository layout (post-Milestone-9)

```
src/
  quaids.sdf        # Struct definitions: quaidsControl, quaidsOut,
                    #   quaidsElasOut. Guarded by #ifndef/#define
                    #   QUAIDS_SDF_INCLUDED / #endif (added Milestone 6) so
                    #   pubtable_quaids.src has a symbol to detect.
  quaidsutil.src    # quaidsControlCreate() / getDefaultQuaidsControl().
  quaidsiv.src      # _quaidsIVFirstStage() -- private helper, the
                    #   instrumental-variables first-stage regression of log
                    #   total expenditure. The one internal phase of
                    #   quaidsFit() that was cleanly separable; see the
                    #   Milestone 2 scoping note below for why the rest of
                    #   quaidsFit() (starting values, iteration, variance,
                    #   overidentification test, symmetry test) was not
                    #   further split.
  quaidselas.src    # quaidsElas_() (silent, low-level), quaidsElasFit()
                    #   (silent, struct-returning: point estimates +
                    #   standard errors), printQuaidsElas() (the separated
                    #   printer), quaidsElas() (backward-compatible
                    #   wrapper: fit then print) -- elasticities at a
                    #   point. See "Milestone 5: elasticities
                    #   generalization" below.
  quaidsslutzky.src # quaidsSlutzky() -- Slutzky negativity diagnostic.
  quaids.src        # quaidsFit() (silent, struct-returning estimation core;
                    #   calls _quaidsIVFirstStage()), printQuaids() (the
                    #   separated estimation-report printer), quaids()
                    #   (backward-compatible wrapper: fits, prints,
                    #   reproduces the legacy elasticities/descriptive
                    #   stats/Slutzky report via quaidsElas()/
                    #   quaidsSlutzky(), returns the 4 legacy matrices).
  quaidsformula.src # quaidsFull() -- dataframe/column-name entry point;
                    #   selects w/intcpt/prices/totexp/instr from a
                    #   dataframe by column name and calls quaidsFit().
  quaidstests.src   # quaidsHomogeneityTest(), quaidsJointTest() -- standalone
                    #   Wald tests, operating on an already-computed
                    #   unconstrained (aCtl.homogenous=0) quaidsOut. See
                    #   "Milestone 4: new hypothesis tests" below.
  pubtable_quaids.src # Optional pubtable adapter -- ptModelFromQuaids()/
                    #   ptFromQuaids() (coefficient tables),
                    #   ptModelFromQuaidsElas()/ptFromQuaidsElas()/
                    #   ptTablesFromQuaidsElas() (elasticity tables),
                    #   ptFromQuaidsFamily() dispatcher. NOT in
                    #   package.json's src array or self-included by any
                    #   other src/ file -- has a hard compile-time
                    #   dependency on pubtable.sdf's ptModel/ptTable struct
                    #   types, so a caller must #include pubtable.sdf/
                    #   pubtable.src (or `library pubtable;`) before this
                    #   file, same as quaids.sdf/quaids.src. See
                    #   "Milestone 6: reporting via pubtable" below.
examples/
  quaids_example.e  # One synthetic 5-good dataset (homogeneity/symmetry true
                    #   by construction), run through quaids() with
                    #   eyeballed comparison of printed estimates to true
                    #   parameters. Not an automated test — no assertions.
                    #   Uses `library quaids;` against the installed
                    #   package (Milestone 9 -- previously #included the
                    #   source tree directly via ../src/..., from before
                    #   Milestone 7 made the package installable; switched
                    #   once README.md/USAGE_GUIDE.md started documenting
                    #   `library quaids;` as the primary usage pattern, so
                    #   the examples actually demonstrate what the docs
                    #   promise, matching gauss-qardl's own examples/,
                    #   which are all `library qardl;`-based). Source-tree
                    #   #include-based testing still lives on in tests/
                    #   (except tests/package_public_api.e, which is also
                    #   library-based, by design -- see "Testing status"
                    #   below).
  pubtable_export_example.e # Milestone 6: same style, but exports a
                    #   quaidsFit() coefficient table and a
                    #   quaidsElasFit() elasticity report to
                    #   LaTeX/Markdown/CSV via pubtable_quaids.src.
                    #   Requires the pubtable package installed. Uses
                    #   `library quaids, pubtable;` (Milestone 9) plus a
                    #   bare `#include quaids.sdf` -- required so that
                    #   pubtable_quaids.src's #ifDef QUAIDS_SDF_INCLUDED
                    #   guard is active; `library quaids;` alone lazily
                    #   loads procs on demand and does not run quaids.sdf's
                    #   #define (confirmed empirically), matching how
                    #   pubtable's own bundled pubtable_qardl.src documents
                    #   the identical requirement for qardl.sdf. Still
                    #   #includes ../src/pubtable_quaids.src by relative
                    #   path, since that adapter is not part of the
                    #   installed quaids package (see "Milestone 6" below).
tests/
  quaids_schema_test.e         # Milestone 1: asserts quaidsOut field
                    #   values/shapes, that quaidsFit() prints nothing, and
                    #   that the legacy quaids() wrapper's returned matrices
                    #   are byte-identical to the struct fields they're
                    #   drawn from.
  quaids_formula_parity_test.e # Milestone 2: builds the same synthetic
                    #   dataset as both plain matrices and a named-column
                    #   dataframe, and asserts quaidsFull(dataframe, ...)
                    #   matches quaidsFit(matrices...) exactly, including
                    #   the extraVars == 0 ("no extra intercept shifters")
                    #   path.
  quaidsfixtures.src            # Milestone 3: _quaidsSyntheticDGP(), a
                    #   shared 5-good synthetic-data generator (not part of
                    #   the public src/ API) parameterized by whether the
                    #   true model has a quadratic term and whether total
                    #   expenditure is genuinely endogenous. Returns the
                    #   true parameters pre-stacked in qOut.bS's row order.
  quaids_synthetic_validation_test.e # Milestone 3: 22 checks recovering
                    #   true DGP parameters within a documented tolerance
                    #   across all 6 (LA-AIDS/AIDS/QUAIDS) x (with/without
                    #   IV) combinations. See the "Seed sensitivity" note
                    #   below for why this uses one specific, documented
                    #   seed rather than an arbitrary one.
  quaids_published_validation_test.e # Milestone 3: fits quaidsFit() on
                    #   real published data (Blanciforti86 food-demand
                    #   data) and checks against an independent R
                    #   reference. This is the test that caught the
                    #   Stone-index starting-value bug -- see "Milestone 3:
                    #   real bug found and fixed" below. Milestone 9:
                    #   extended to also validate iterated AIDS against R
                    #   aidsEst(method="IL", ...) -- 19 checks total (was
                    #   11) -- see "Milestone 9" below.
  fixtures/published/
    blanciforti86_food32.csv  # The data itself; see SOURCE.md for
                    #   attribution/license.
    SOURCE.md                 # Attribution, license note, repo-owner
                    #   approval record.
    generate_r_reference.R    # Reproduces the R numbers hardcoded into
                    #   quaids_published_validation_test.e. Requires R +
                    #   the CRAN package micEconAids.
    python_reference_check.py # Independent from-scratch Python replica;
                    #   supplementary evidence, not the assertion source
                    #   (see its header for why). Requires numpy/pandas.
  quaids_hypothesis_tests_test.e # Milestone 4: 19 checks -- size AND power
                    #   for quaidsHomogeneityTest()/quaidsJointTest(), a
                    #   power check for the existing symmetry-given-
                    #   homogeneity test, and the first-ever exercise of
                    #   the overidentification test (every prior fixture
                    #   was exactly identified, ninst==nu).
  quaids_elasticities_test.e   # Milestone 5: 17 checks -- parity between
                    #   quaidsElasFit()/printQuaidsElas() and the
                    #   pre-split quaidsElas_(), plus three EXACT
                    #   algebraic identities (Engel aggregation, Cournot
                    #   aggregation, elasticity homogeneity) checked at a
                    #   real out-of-sample observation and a synthetic
                    #   counterfactual price scenario -- not just the four
                    #   points quaids() has always used.
  quaids_pubtable_test.e       # Milestone 6: 30 checks -- exact numeric
                    #   parity between pubtable ptModel.estimates/
                    #   stdErrors and the qOut.bestB/qOut.bestV/
                    #   elasOut.er/ser values they're built from, shape/
                    #   title checks, the ptFromQuaidsFamily dispatcher,
                    #   and an end-to-end export smoke test that writes
                    #   real .tex/.md/.csv files and reads them back.
                    #   Requires the pubtable package installed.
                    #   All test/fixture files run from tests/ (or
                    #   tests/fixtures/published/ for the R/Python
                    #   scripts) as the working directory.
  package_public_api.e   # Milestone 7: installed-package release gate --
                    #   `library quaids;` (not #include) against a real
                    #   install, exercising quaidsControlCreate/
                    #   getDefaultQuaidsControl, quaidsFit/printQuaids/
                    #   quaids, quaidsFull, quaidsElasFit/quaidsElas/
                    #   printQuaidsElas, quaidsSlutzky,
                    #   quaidsHomogeneityTest/quaidsJointTest. Builds its
                    #   own small inline synthetic dataset rather than
                    #   reusing quaidsfixtures.src (tests/-only, not part
                    #   of the installed package). Run after
                    #   scripts/run_release_verification.ps1
                    #   -InstallArtifact. See "Milestone 7: package build
                    #   and release tooling" below. Milestone 9: found and
                    #   fixed a real gap here -- printQuaids() and
                    #   quaidsElas() were never actually called (only
                    #   their split components were), so a load-order bug
                    #   or stale .lcg entry specific to either one could
                    #   have passed this gate undetected. See "Milestone 9"
                    #   below.
  run_source_tests.ps1    # Milestone 7: runs verify_package_manifest.ps1
                    #   then all 7 tgauss test files above, checking this
                    #   repo's own PASS/FAIL-line convention (not just
                    #   tgauss's exit code -- see "Testing status" below).
  verify_package_manifest.ps1  # Milestone 7: package.json src array vs.
                    #   actual src/ directory consistency (no dupes,
                    #   nothing missing, nothing unlisted except the
                    #   documented pubtable_quaids.src allowlist entry).
scripts/
  build_lcg.ps1     # Milestone 7: writes lib/<name>.lcg, the plain-text
                    #   symbol-location catalog GAUSS's `library`
                    #   mechanism reads (verified against the real,
                    #   installed qardl.lcg/pubtable.lcg format -- not a
                    #   stub). Run against an INSTALLED copy of the
                    #   package (e.g. the staging dir build_package.ps1 +
                    #   Expand-Archive produce), not this repo directly.
  build_package.ps1 # Milestone 7: stages package.json plus whichever of
                    #   README.md/CHANGELOG.md/CITATION.cff/CITATION.md/
                    #   LICENSE/llms.txt exist, plus src/docs/examples/
                    #   scripts/tests, strips generated run artifacts, and
                    #   zips "<name> <version>.zip" in the repo root (already
                    #   gitignored via the pre-existing `quaids *.zip`
                    #   pattern from Milestone 0).
  verify_release_artifact.ps1  # Milestone 7: checks a built .zip's name/
                    #   CHANGELOG.md entry match package.json's version, no
                    #   stale artifacts sit in the repo root, and the
                    #   archive contains every file package.json's src
                    #   array promises plus this repo's current root
                    #   files/dirs -- README.md and all four docs/*.md
                    #   files added to the required-entries list at
                    #   Milestone 8, once they existed.
  run_release_verification.ps1 # Milestone 7: orchestrator -- source tests,
                    #   build/verify the release artifact, optionally
                    #   install it into a GAUSS package directory
                    #   (defaults to <GaussHome>/pkgs, i.e. for real,
                    #   alongside every other package on the machine, only
                    #   with -InstallArtifact), then the installed-package
                    #   public API gate.
docs/
  COMMAND_REFERENCE.md  # Milestone 8: index of every public proc, grouped
                    #   (control struct, estimation, hypothesis tests,
                    #   elasticities/diagnostics, optional pubtable
                    #   reporting), linking to one command-reference/*.md
                    #   page each. Cross-checked against the actual source
                    #   by tests/verify_package_manifest.ps1 (every
                    #   documented proc must exist in src/; every linked
                    #   page must exist).
  USAGE_GUIDE.md    # Milestone 8: choosing an API (quaidsFit vs. quaids
                    #   vs. quaidsFull), the LA-AIDS/iterated-AIDS/QUAIDS
                    #   switch table, why IV is always required, symmetry/
                    #   homogeneity/overidentification workflow,
                    #   elasticities at arbitrary points, pubtable
                    #   reporting, current limitations.
  METHODOLOGY_NOTES.md  # Milestone 8: the estimator itself -- Stone vs.
                    #   translog price index, the QUAIDS quadratic term,
                    #   the IV control-function approach (and how it
                    #   differs from R micEconAids's 3SLS), the full
                    #   estimation algorithm phase-by-phase, elasticity/
                    #   Slutzky formulas, citing Deaton & Muellbauer (1980)
                    #   and Banks, Blundell & Lewbel (1997).
  FEATURE_SUPPORT_MATRIX.md  # Milestone 8: LA-AIDS x iterated-AIDS x
                    #   QUAIDS support for IV, hypothesis tests,
                    #   elasticities, Slutzky, curvature (not supported,
                    #   documented not silently absent), dataframe API,
                    #   pubtable export, synthetic/published validation.
  command-reference/  # Milestone 8: one *.md page per public proc (18
                    #   pages) -- Purpose/Format/Parameters/Returns/
                    #   Remarks/Examples/Source/See Also, matching
                    #   gauss-qardl's page template. Covers every proc in
                    #   quaids.sdf's load-bearing src/ files plus the
                    #   optional pubtable_quaids.src adapter's 6 procs
                    #   (documented despite being outside package.json's
                    #   src array, since they're real public API surface).
package.json      # GAUSS package manifest (name: quaids, version: 0.5.0,
                  #   license: MIT). pubtable_quaids.src deliberately not
                  #   listed in its src array -- see "Milestone 6" below.
LICENSE           # MIT, copyright Eric Clower.
CITATION.cff      # Citation metadata; cites Deaton & Muellbauer (1980) and
                  #   Banks, Blundell & Lewbel (1997).
CHANGELOG.md      # Milestone 7: version history 0.1.0-0.5.0, reconstructed
                  #   from this file's own milestone records (Keep a
                  #   Changelog style). No git tag cut yet -- nothing in
                  #   this repo has been committed as of Milestone 7.
README.md         # Milestone 8: front door -- install (Tools > Install
                  #   Application, or the scripts/ build+install
                  #   one-liner), quick start, model-choice summary,
                  #   feature list, links into docs/, testing/release
                  #   commands.
.gitignore        # Compiled .gcg artifacts, tmp/, .claude/, packaged zips,
                  #   generated `output file=...` run artifacts.
GOLD_STANDARD_TODO.md  # Living roadmap: release blockers, milestones,
                  #   definition of done. Read this before any nontrivial
                  #   change and update it as milestones close.
```

All nine milestones are complete: 0 (repo hygiene), 1 (API/output-schema
baseline), 2 (modular source split + dataframe entry point), 3 (validation
fixtures, including published-data cross-implementation validation), 4
(hypothesis testing completeness), 5 (elasticities/diagnostics
generalization), 6 (reporting via `pubtable`), 7 (package build and
release tooling), 8 (documentation), and 9 (final gold standard
integration gate).

**The package is now actually installed** at `c:\gauss26\pkgs\quaids`
(Milestone 7), alongside `qardl` and `pubtable` on this machine --
`library quaids;` works. Rebuild/reinstall after any `src/` change with
`powershell scripts\run_release_verification.ps1 -BuildArtifact
-ForceArtifact -InstallArtifact` (run from the repo root).

**pubtable is installed in this environment** at `c:\gauss26\pkgs\pubtable`
(package version `1.0.0`) and is what `src/pubtable_quaids.src` (Milestone
6) targets. Like R/Python below, it is not a repo dependency — nothing in
`src/` or the core test suite requires it; only
`examples/pubtable_export_example.e` and `tests/quaids_pubtable_test.e` do.

**R and Python are installed in this environment** (as of 2026-07-20, with
explicit repo-owner approval) for cross-implementation validation: R 4.5.0
at `C:\Program Files\R\R-4.5.0\bin\` (not on `PATH` — invoke `Rscript.exe`
by full path) with the CRAN package `micEconAids` (and its dependencies
`micEcon`, `systemfit`) installed to the user library; Python 3.12 with
`numpy`/`pandas`/`scipy` installed via `pip`. Neither is a repo dependency
(nothing in `src/` or the test suite requires them to run) — they exist
only to regenerate/extend the published-validation reference numbers in
`tests/fixtures/published/`.

**Milestone 2 scoping note**: `quaidsFit()`'s starting-value construction,
iteration loop, variance computation, overidentification test, and symmetry
test/symmetry-constrained stage were deliberately *not* further split into
separate procs/files. They share heavily mutated intermediate state across
phases (`m`, `gg`, `gw`, `ng`, and iteration-loop-local `_beta`/`lambda`/
`lx`/`b_p`/`lx2` all flow into the variance and overidentification-test
formulas) — splitting them now would mean either passing a long, brittle
parameter list between new procs or bundling "in-progress fit state" into a
new struct, for a maintainability win that's real but not urgent on
pre-alpha code with no other consumers yet. Revisit after Milestone 3's
validation/fixture harness exists to safety-net a deeper refactor of this
kind — see `GOLD_STANDARD_TODO.md`'s "Roadmap Rules" for why that ordering
matters. The IV first stage (`quaidsiv.src`) was extracted because it
*is* cleanly separable: it only reads `intcpt`/`prices`/`instr`/`totexp`/`n`
and its outputs (`prices` converted to relative, `u`, `zzi`, `m1`, and the
IV diagnostics) are either consumed once downstream or not at all after this
phase.

## The `quaidsFit()` / `printQuaids()` / `quaids()` split

**`quaidsFit()` is the primary, silent, struct-returning entry point.** It
does 100% of the estimation with zero printing and returns a `quaidsOut`
struct (defined in `src/quaids.sdf`, ~75 fields, grouped by phase: metadata,
IV first-stage diagnostics, homogeneity-constrained stage, overidentification
test, symmetry test + symmetry-constrained stage, and the final
absolute-price-form `b`/`v`/`bS`/`vS`).

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

**`printQuaids(qOut)`** reproduces the estimation-stage console report (IV
first-stage table, iteration summary, homogeneity-constrained coefficient
table, overidentification test, symmetry test + symmetry-constrained table)
from a `quaidsOut` struct alone. It does *not* print elasticities,
descriptive statistics, or the Slutzky diagnostic — those are separate,
explicitly-callable reports (`quaidsElas()`, `quaidsSlutzky()`), since
Milestone 5 plans to generalize elasticities to arbitrary evaluation points
rather than bake a fixed set into the struct.

**`quaids()` is the original, backward-compatible call** — unchanged
signature and unchanged printed behavior:

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, struct quaidsControl aCtl);
```

It calls `quaidsFit()`, calls `printQuaids(qOut)`, then reproduces the
legacy elasticities-at-four-points / descriptive-statistics / Slutzky report
exactly as the pre-Milestone-1 `aids()` proc did, then returns
`(qOut.b, qOut.v, qOut.bS, qOut.vS)`. **Verified byte-for-byte at the time**:
running `quaids()` post-refactor against the same fixed-seed synthetic
dataset used at Milestone 0 produced output identical to the pre-Milestone-1
code, including the per-iteration convergence log (reproduced from a stored
`qOut.iterHistory` matrix rather than printed live during iteration, since
`quaidsFit()` cannot print, but reproduced with matching iteration
numbers/error values so the printed table was unchanged). **This byte-parity
claim is historical, not current**: Milestone 3 found and fixed a real bug
in the Stone-index starting-value computation (see "Milestone 3: real bug
found and fixed" below), which intentionally changes numerical output for
`aCtl.maxiter==1` calls, and shifts the iteration path for `aCtl.maxiter>1`
calls, relative to this old baseline. The *refactor* (Milestone 1's
estimation/printing split) was and remains behavior-preserving with respect
to the code as it existed at that time; a later, separate, intentional
correctness fix is not a violation of that.

- `w` — `TxN` budget shares.
- `intcpt` — `TxK` extra intercept-shifter variables (demographics etc.), or
  `0` for none.
- `prices` — `TxN` **log** prices (absolute, not relative — the proc converts
  internally).
- `totexp` — `Tx1` log total expenditure (treated as endogenous).
- `instr` — `TxH` instruments for log total expenditure.
- `aCtl` — `quaidsControl` struct (see below).

`quaids()`'s (`b1`/`v1`/`b2`/`v2`) return semantics are unchanged: if
`aCtl.homogenous == 1`, `b1`/`v1` are homogeneity-constrained
estimates/covariance and `b2`/`v2` are homogeneity+symmetry-constrained. If
`aCtl.homogenous == 0`, `b1`/`v1` are the unconstrained (reparametrized)
estimates and `b2`/`v2` are `0`. Same values are available struct-side as
`qOut.b`/`qOut.v`/`qOut.bS`/`qOut.vS`; `qOut.bestB`/`qOut.bestV` always hold
"whichever is the most-constrained estimate actually fit" (symmetric if
homogeneity was imposed, else the recovered unconstrained fit) — this is
what elasticities/Slutzky are evaluated against, matching what the original
code's reused `b_s`/`v_s` locals held at that point regardless of branch.

**A pre-existing anomaly, preserved not fixed**: in the symmetry-constrained
table's "Residuals of instrumental regressions" row, the original code pairs
point estimates from the *homogeneity-stage* `b` with standard errors/t/p
from the *symmetry-stage* fit (both stages leave the IV-residual coefficient
block numerically identical, but the printed SE/t/p come from the
symmetry-constrained covariance). This looks like it could be an original
authoring inconsistency. It was carried over unchanged (`qOut.homogB` paired
with `qOut.symcSE`/`qOut.symcT`/`qOut.symcPvt` in `printQuaids()`) since
Milestone 1 is a structure-preserving refactor, not a correctness pass —
flag for review during Milestone 3/4 validation work.

## `quaidsFull()` — dataframe entry point

```gauss
struct quaidsOut qOut;
qOut = quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars, aCtl);
```

`data` is an already-loaded dataframe (`data = loadd(fname);` or built with
`asDF()`/`dfaddcol()` — `quaidsFull()` does not load files itself).
`shareVars`/`priceVars` are string arrays of column names, **matched by
position** (`shareVars[i]` and `priceVars[i]` must be the same good — there
is no name-matching magic). `totexpVar` is a single column-name string.
`instrVars` is a string array. `extraVars` is a string array, or the scalar
`0` for "no extra intercept shifters" (matching `quaidsFit()`'s `intcpt ==
0` convention — checked via `type(extraVars)`, not value equality, since
comparing a string array to `0` with `==` is itself a type-mismatch trap).
Internally it just selects the named columns and calls `quaidsFit()` — no
duplicated estimation logic. Verified numerically identical to the
equivalent `quaidsFit(matrices...)` call, including the `extraVars == 0`
path, by `tests/quaids_formula_parity_test.e` (17 checks).

### GAUSS type-system notes learned while building this

- `type()` codes actually seen on GAUSS 26.1.1: `6` = plain matrix (and, it
  turns out, **also** a dataframe and a dataframe column selection —
  `data[., "col"]` and `data[., stringArray]` both return plain type-`6`
  matrices, not some distinct "dataframe" type; dataframe-ness doesn't leak
  into downstream numeric code). `13` = native `string`. `15` = native
  `string array`. Legacy `0$+"X"$+ftocv(...)`-built character matrices are
  *also* type `6` (matrix) — this is why `quaidsOut`'s name-vector fields
  are declared `matrix`, not `string array` (see below).
- Struct-field assignment is strict about matching one of these types
  exactly; ordinary expressions are more forgiving. `"CONSTANT"|xnam`
  (vertical-concatenating a native string literal onto a legacy char-matrix)
  silently coerces and works fine as an *expression*. Directly assigning a
  bare native string literal to a `matrix`-typed struct field
  (`qOut.xnam = "CONSTANT";`) does **not** coerce — it throws `G0071 Type
  mismatch`. Use the `0$+` prefix idiom (`0$+"CONSTANT"`) to get the
  legacy-char-matrix form explicitly before assigning into a `matrix`
  struct field.

### Two pre-existing bugs found (and fixed) by the formula parity test

`tests/quaids_formula_parity_test.e` exercises `intcpt == 0` (no extra
intercept-shifter columns) — a code path that existed unchanged from the
original `aids_rev.src` but had **never actually been run** by anything in
this repo before, because the synthetic example always passes a non-zero
`intcpt`. Writing a real caller for that path immediately surfaced two bugs,
both in `quaidsFit()`'s name-setup block (`src/quaids.src`, just after the
`intcpt == 0` branch):

1. `xnam` was never assigned inside the `if intcpt == 0;` branch (only in
   `else;`), so the following unconditional `xnam = "CONSTANT"|xnam;` read
   an uninitialized variable — `G0152 Variable not initialized`.
2. After fixing that by assigning `xnam = "CONSTANT";` in the `if` branch,
   assigning it into `qOut.xnam` (`matrix` type) threw `G0071 Type
   mismatch`, per the type-system note above — a native string literal
   isn't the same type as the char-matrix `xnam` holds in the `else`
   branch. Fixed with `xnam = 0$+"CONSTANT";`.

Both fixes are scoped to the previously-dead branch only; the `else` branch
(exercised by every prior test) was not touched, and re-running the
Milestone 0/1 parity and schema tests after the fix confirmed zero change to
already-verified behavior. Lesson: a code path with no caller is a code path
with no evidence it works, regardless of how long it's sat there unchanged.

## Milestone 3: synthetic validation findings

`tests/quaids_synthetic_validation_test.e` fits `quaidsFit()` against
`tests/quaidsfixtures.src`'s `_quaidsSyntheticDGP()` (a 5-good dataset with
homogeneity/adding-up true by construction, parameterized by whether the
true model has a quadratic term and whether total expenditure is genuinely
endogenous) and checks recovered parameters against the true ones, not just
that the code ran. Full tolerance rationale is in `GOLD_STANDARD_TODO.md`'s
Milestone 3 "Tolerance Policy" section — summary: structural coefficients
within `0.10`, the IV-residual coefficient row within `0.50` (it's
consistently the noisiest row by roughly 10x — a control-function term on
an estimated regressor, not a "deep" structural parameter), LA-AIDS within
`1.20` throughout (Stone-index approximation bias is a real property of
that method).

### Seed sensitivity — a real numerical-reliability finding

Calibrating those tolerances required an 8-seed probe (`tobs=3000`,
`aCtl.err=.0001`, `aCtl.maxiter=100`) that turned up something worth
flagging prominently: **the iterative estimator (QUAIDS and iterated linear
AIDS) fails to converge, or "converges" to wildly wrong estimates (errors
of magnitude 200–2500 against true parameters of magnitude ~0.1–2), for
roughly half of random seeds** in this DGP family. The failure pattern
tracked with the *price* draw for a given seed, not with whether the model
had a quadratic term or genuine endogeneity — the same seeds failed (or
succeeded) whether `quadratic`/`endogenous` were on or off, holding
`prices` fixed. This points at the iteration's conditioning being sensitive
to the specific price data realization, not at a bug specific to QUAIDS or
IV handling.

The synthetic fixtures use `seed=204`, one of the seeds confirmed (by that
probe) to converge cleanly across all six model/endogeneity combinations —
documented as such in the test file's comments, not silently cherry-picked.
**Partially explained, not fully fixed**: the Stone-index starting-value
bug (below) likely explains *some* of this non-convergence (a materially
wrong starting point makes convergence failure more likely), and was fixed.
But re-running the same 8-seed probe after the fix would be needed to know
how much of the non-convergence rate that actually accounts for — not done,
since it wasn't necessary to validate the fix itself (the published-data
comparison already did that directly). `quaidsFit()`'s iteration still has
no globally-convergent guarantee, no damping, and no fallback for a bad
starting point. Worth a dedicated numerical-reliability pass (analogous to
`gauss-qardl`'s Milestone 13) once more of the roadmap exists to build on.

### Milestone 3: real bug found and fixed — Stone-index starting value

The published-data cross-check below initially disagreed with R by roughly
5x on `beta` — far beyond the few-percent gap expected between two
different-but-valid IV algorithms. Root cause, in `quaidsFit()`'s "STARTING
VALUE" block (`src/quaids.src`): `stone = prices*meanc(w)` was applied to
`prices` *after* it had already been converted to relative form in columns
`1:n-1` (each minus the reference good's price) while column `n` stayed
absolute. Weighting that mixed matrix by mean shares does not compute the
Stone price index — algebraically (verified by direct derivation, then
confirmed empirically by patching the formula in an isolated Python check
and watching the gap with R close) it computes
`StandardStoneIndex − ln(p_n)·(1 − meanShare_n)`, a distortion tracking the
reference good's own price trend.

**Impact**: `aCtl.maxiter==1` (LA-AIDS) never iterates past this starting
value, so the distorted deflator *was* the final answer for every LA-AIDS
call using default starting values (`aCtl.b0==0`) — not an occasional
glitch, every time. For `aCtl.maxiter>1`, this was only a bad *starting
point*; the iteration loop's own `a_p` formula is correct and unaffected.

**Fix** (`src/quaids.src`, "STARTING VALUE" block): reconstruct absolute
prices before computing `stone`:
```gauss
stone = (prices[., 1:n-1] + prices[., n])~prices[., n];
stone = stone*meanc(w);
```
rather than changing the relative-price convention used elsewhere in the
proc. Confirmed this fix (not something else) closed the gap: after
patching, GAUSS's estimates matched R's independent `micEconAids`
3SLS-with-instrument reference to within ~0.021 max absolute difference
(see below) instead of ~5x off.

**This changes numerical output** for `aCtl.maxiter==1` calls, and shifts
(not necessarily worsens) the iteration path for `aCtl.maxiter>1` calls,
relative to every prior milestone's frozen baseline — expected and correct,
since those baselines were captured from the original, buggy code.
`examples/quaids_example.e` (seed 11) was already one of the non-converging
seeds identified above even pre-fix; its output changed too and was not
re-tuned to a better-behaved seed as part of this fix.

### Published-data validation against R and Python

`tests/quaids_published_validation_test.e` fits `quaidsFit()`
(`aCtl.linear=1, aCtl.maxiter=1` — LA-AIDS) on `Blanciforti86` (annual U.S.
food-consumption data, 1947–1978, 4 food groups — see
`tests/fixtures/published/SOURCE.md`), instrumenting `log(xFood)` (total
food expenditure) with `log(xAgg)` (total aggregate expenditure, `corr ≈
0.97` in logs — a strong, genuinely informative instrument, not a
near-degenerate one: first-stage `R² ≈ 0.998`).

- **R** (`tests/fixtures/published/generate_r_reference.R`,
  `micEconAids::aidsEst(..., instNames=...)`, dispatches to 3SLS via
  `systemfit`): **max abs difference from GAUSS ≈ 0.021** across
  alpha/beta/gamma. This is the hard assertion target in the test
  (tolerance `0.05`, real headroom above the observed gap). The residual
  ~0.02 is attributable to R's 3SLS and GAUSS's control-function
  (residual-inclusion) approach being different, both valid, IV algorithms
  for the same model — not expected to match bit-for-bit.
- **Python** (`tests/fixtures/published/python_reference_check.py`,
  hand-coded from the Deaton-Muellbauer equations — no comparably
  established Python AIDS package exists the way R has `micEconAids`):
  broadly consistent (one equation matches GAUSS almost exactly) but with
  larger residual differences on another equation than either the R
  comparison or an earlier, simpler Python no-IV replica showed. Since R —
  an independently authored, widely-used implementation — agrees closely
  with GAUSS on exactly the coefficients where this Python script diverges
  most, that divergence is attributed to the from-scratch replica (most
  likely in how it forms the IV-residual-augmented GLS weighting), not to
  GAUSS. Kept in the repo for transparency and as a documented starting
  point; **not used as a pass/fail assertion source**, only R is.
- Both scripts are runnable standalone from `tests/fixtures/published/`
  (`Rscript generate_r_reference.R`, `python python_reference_check.py`)
  to regenerate or extend these numbers.

## Milestone 4: new hypothesis tests

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);  // aCtl.homogenous = 0
{ stat, pval, df } = quaidsHomogeneityTest(qOut);
{ statJ, pvalJ, dfJ } = quaidsJointTest(qOut);
```

Both (`src/quaidstests.src`) are Wald chi2 tests and **require an
unconstrained fit** (`qOut.homogenous == 0` — they error clearly if not).
They read `qOut.b`/`qOut.v` (the final, absolute-price-form unconstrained
gamma matrix and its covariance) and build a restriction vector/covariance
via a selection matrix `L`: `stat = (L'*vec(b))' * inv(L'*V*L) * (L'*vec(b))`.

- **`quaidsHomogeneityTest`**, `df = n-1`: tests `sum_j gamma_ij = 0` for
  each independently-estimated equation (equation `n` is recovered via
  adding-up and adds no information to a Wald test).
- **`quaidsJointTest`**, `df = (n-1) + (n-1)(n-2)/2`: homogeneity's
  restrictions plus symmetry (`gamma_ij = gamma_ji`, `i<j`, `i,j=1..n-1`).
  A symmetric gamma with adding-up already implies homogeneity (row `i`
  sum = column `i` sum by symmetry = 0 by adding-up) — there's no separate
  "symmetry given adding-up alone" test on an unconstrained fit; use this
  joint test, or the existing symmetry-given-homogeneity test if
  homogeneity itself isn't in question.

Both are validated for size *and* power in `tests/quaids_hypothesis_tests_test.e`
(19 checks) — full derivation, both dead ends hit while getting there, and
why size+power (not just "it runs") is the standard applied, are in
`GOLD_STANDARD_TODO.md`'s Milestone 4 section. Short version: the first
implementation read `qOut.b`'s row layout wrong (mistook an internal
pre-recovery reparametrization, described in `quaidsFit()`'s own
docstring, for the final post-recovery form) and rejected a true null with
`p≈0` — caught immediately by running the size check before trusting the
formula, not by code review.

Same milestone also **exercised the overidentification test for the first
time in this repo's history**: every fixture through Milestone 3 used
exactly-identified instruments (`ninst==nu`), so `qOut.overidValid` was
always `0` and that branch of `quaidsFit()` had never actually run. A new
2-instrument fixture in the same test file confirms it runs, has the
right shape (`overidGamma` is `ninst x n`) and df (`ninst-nu`), and
doesn't spuriously reject when both instruments are valid by construction.

**Why the existing symmetry/overID tests weren't extracted into their own
procs** (unlike the two new ones above): they depend on heavily mutated
intermediate iteration/variance state (`Ji`, `Sgma`, `S`, `O`, `D`, `zzi`,
`m1`, iteration-final `_beta`/`lambda`/`lx`/`b_p`/`lx2`) — the same coupling
Milestone 2 already declined to untangle for the estimation core itself,
unchanged here since nothing about the estimation core moved. The two new
tests avoid this entirely: they only need the *finished* `qOut.b`/`qOut.v`,
not any of that intermediate state.

## Milestone 5: elasticities generalization

```gauss
struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPoint, pricesPoint, totexpPoint, aCtl);
call printQuaidsElas(elasOut);   // optional -- omit for a silent call
```

**The roadmap's framing turned out to be slightly off, worth recording**:
it described `quaidsElas()`/`quaidsElas_()` as "hardcoded to mean/Q1/
median/Q3," but `quaidsElas_()` already took `intcpt`/`prices`/`totexp` as
*point arguments*, not sample statistics — it always could evaluate
anywhere. The actual gap was that `quaidsElas()` mixed computation with
printing (no way to get a result object silently) and that `quaids()`
(the legacy wrapper) only ever *called* it at four fixed points. Fixed by
splitting `quaidsElas()` the same way Milestone 1 split `quaidsFit()`:

- **`quaidsElasFit(b, v, intcpt, prices, totexp, aCtl)`** — silent, returns
  `struct quaidsElasOut` (`er`/`ep`/`epc` point estimates, `ser`/`sep`/`sepc`
  raw numeric delta-method standard errors — not the pre-formatted
  `"(0.123)"` display strings the old code produced inline, which stay in
  the printer where they belong).
- **`printQuaidsElas(elasOut)`** — the separated printer.
- **`quaidsElas(...)`** — unchanged signature, now a thin wrapper
  (`quaidsElasFit()` then `printQuaidsElas()`). **Verified byte-for-byte**
  identical printed output against the pre-split version on
  `examples/quaids_example.e`.
- `quaidsElas_()` itself is untouched — it was already the right shape.

**Validating "arbitrary point" actually works, not just "compiles"**:
`tests/quaids_elasticities_test.e` checks three **exact** algebraic
identities that any valid AIDS/QUAIDS elasticity set must satisfy given
adding-up/homogeneity (which this estimator always imposes by
construction) — these are consequences of the functional form, not
separately-estimated statistical quantities, so they should hold to
floating-point precision, not just "within tolerance":

- Engel aggregation: `sum_i(w_i * er_i) = 1`
- Cournot aggregation: `sum_i(w_i * ep_ij) + w_j = 0`, for each price `j`
- Elasticity homogeneity: `sum_j(ep_ij) + er_i = 0`, for each good `i`

Checked at a real out-of-sample observation and at a fully synthetic
counterfactual (a hypothetical 20% price increase on one good, evaluated
at the sample-mean point otherwise) — both held to ~1e-16. A test proc
bug surfaced immediately while writing this (`struct quaidsElasOut
elasOut` must be declared in the parameter list for field access to work,
same as any GAUSS struct-typed parameter — omitting it gives a `G0008
Syntax error` on the first `.field` access, not a silent failure).

**`quaidsSlutzky()` needed no change**: it already accepts an arbitrary
`intcpt`/`prices`/`totexp` sample (any number of rows) as input, so "keep
the Slutzky diagnostic general" was already true going into this
milestone.

**Curvature imposition (Diewert-Wales Cholesky reparametrization) — scoped,
not implemented**: the roadmap listed this as P2/"if justified." Even
`micEconAids` (the reference implementation used for Milestone 3's
cross-check) only offers curvature *diagnosis* (`aidsMono()`,
`aidsConcav()` — post-hoc checks, matching what `quaidsSlutzky()` already
does here), not *imposition* as a constrained-estimation mode. That's real
evidence against this library needing to leap ahead of the reference
implementation on a deliberately-optional item — revisit only if a
concrete use case shows up.

### `quaidsControl` fields (`src/quaids.sdf` / `quaidsControlCreate()`)

| Field | Default | Meaning |
|---|---|---|
| `linear` | `0` | `1` = LA-AIDS (linear); `0` = QUAIDS (quadratic log-expenditure term) |
| `maxiter` | `50` | `1` = one-step linearized AIDS with Stone price index; `>1` = iterate |
| `homogenous` | `1` | `1` = impose homogeneity (and test/report symmetry); `0` = unconstrained |
| `alpha0` | `0` | Fixed value of the translog price-index intercept `alpha_0` |
| `err` | `.0001` | Relative parameter-change convergence tolerance |
| `othnam` | `""` | Optional alternate variable names for printed output |
| `b0` | `0` | Optional user-supplied starting values; `0` = use linearized-AIDS starting values |

The `stone`, `aids`, and `varname` fields flagged as dead at Milestone 0 were
removed from `quaidsControl` at Milestone 1 (also dropped from
`quaidsControlCreate()`), since grep confirmed no read anywhere in
`quaids.src` and there are no external consumers of this pre-alpha struct yet.
Also added: `getDefaultQuaidsControl()` (alias for `quaidsControlCreate()`,
naming-convention parity with `gauss-qardl`'s `getDefault...Control` procs),
and structure-inference return typing (`proc (struct quaidsControl) =
quaidsControlCreate();`) so callers no longer need to pre-declare
`struct quaidsControl aCtl;` before assignment.

## Milestone 6: reporting via pubtable

```gauss
library pubtable;
#include quaids.sdf
#include quaidsutil.src
#include quaidsiv.src
#include quaidselas.src
#include quaidsslutzky.src
#include quaids.src
#include pubtable_quaids.src

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);            // one comparison column per good
call ptExport(coefTbl, "results.tex");   // .tex/.md/.csv/.rtf/.html/.xlsx by extension

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);  // 3x1: income, uncompensated, compensated
call ptExport(elasTbls[1], "income_elasticities.md");
```

`src/pubtable_quaids.src` is an optional adapter onto the `pubtable`
package (`c:\gauss26\pkgs\pubtable`), following the same pattern as
`pubtable`'s own bundled `pubtable_qardl.src`
(`c:\gauss26\pkgs\pubtable\src\pubtable_qardl.src`): every proc is guarded
by `#ifDef QUAIDS_SDF_INCLUDED` with a clear `_library_missing_error()` stub
otherwise, using the `(...)`/`dynargsGet()` variadic-argument idiom for the
struct-typed inputs.

**Location deliberately diverges from the `pubtable_qardl.src` precedent**:
that file lives *inside the installed pubtable package itself*, not inside
`gauss-qardl`'s own repo — pubtable ships adapters for both Aptech's own
estimators (`cmlmt`/`maxlikmt`/`optmt`/`tsmt`) and third-party ones
(`qardl`) bundled together in its own `src/`. Matching that physically for
`pubtable_quaids.src` would mean writing into a shared, installed package
outside this repo's git history, affecting every other project on this
machine that loads `pubtable` — asked of the repo owner explicitly rather
than assumed (see `GOLD_STANDARD_TODO.md`'s Milestone 6 section for the
full reasoning). Decision: `src/pubtable_quaids.src` stays inside this
repo, self-contained and git-tracked, matching every other file in `src/`
(none of which `#include` their own dependencies — the caller does). It is
**not** listed in `package.json`'s `src` array, since (unlike every other
file there) its `proc (struct ptModel) = ...`/`proc (struct ptTable) = ...`
return-type annotations sit outside the `#ifDef` guard and need
`pubtable.sdf`'s struct types declared unconditionally — adding it to
`src` would make `pubtable` a hard compile-time dependency for the whole
package, contradicting `package.json`'s empty `deps` array.

**Public procs**:
- `ptModelFromQuaids(name, qOut, eqIdx)` / `ptFromQuaids(qOut)` — one
  equation's (good's) coefficient column from `qOut.bestB`/`qOut.bestV`
  (the same "most-constrained estimate actually fit" elasticities/Slutzky
  are evaluated against), or an `n`-good comparison table via
  `ptModelCompare` (mirrors `pubtable_qardl.src`'s per-quantile
  `ptFromQardl`). Row order: intercept block (`qOut.xnam`) | price/gamma
  block (`GAMMA_`+good name, one row per good) | `BETA_LX` |
  `LAMBDA_LX2` (QUAIDS only) | one row per IV control-function residual
  term (`qOut.unam`, `qOut.nu` rows — always present, since `quaidsFit()`
  always treats log total expenditure as endogenous).
- `ptModelFromQuaidsElas(name, elasOut)` / `ptFromQuaidsElas(elasOut)` —
  income elasticities (`elasOut.er`/`elasOut.ser`) as a single-column
  `ptModel`/`ptTable`.
- `ptTablesFromQuaidsElas(elasOut)` — 3x1 `struct ptTable` array (income
  elasticities, uncompensated/Marshallian price elasticities, compensated/
  Hicksian price elasticities), suitable for `ptExportAll()`. The two
  price-elasticity tables are `n x n` matrices with a value row and an
  `(se)` row per good — built directly via the internal
  `_ptQuaidsElasMatrixTable()` helper rather than through `ptModel`, since
  that shape doesn't fit `ptModel`'s single-coefficient-vector design.
- `ptFromQuaidsFamily(x)` — `isStructType`-based dispatcher (`quaidsOut` ->
  `ptFromQuaids`, `quaidsElasOut` -> `ptFromQuaidsElas`), mirroring
  `pubtable_qardl.src`'s `ptFromArdlFamily`.

**A real type-system finding, empirically verified (not assumed from the
Milestone 1 lesson)**: `quaidsOut`/`quaidsElasOut` name vectors (`xnam`,
`wnam`, `unam`, ...) are legacy character matrices (GAUSS type 6); `pubtable`
struct fields for names/titles are natively typed `string array` (type 15)
or scalar `string` (type 13). Direct assignment throws `G0071 Type
mismatch` in *this* direction too (matrix into string-typed field), not
just the matrix-field direction CLAUDE.md already documented. Tested
several candidate conversions directly with `tgauss` rather than guessing:
`strtrim()` looked like the obvious builtin and instead errors ("Invalid
argument type") on a legacy char-matrix input — **there is no dedicated
conversion builtin**. The idiom that does work: concatenating the char-matrix
with a native string via `$|` forces element-wise conversion to a native
string array (`cm $| ""`, slicing off the trailing blank row for a full
array; a *scalar*-indexed single row of the result, `sa[1]` and not
`sa[1:1]`, further demotes to a true scalar `string`). Encapsulated as
`_ptQuaidsToStrArray()`/`_ptQuaidsToStr()` in `pubtable_quaids.src`.

**A real row-count bug, found by running against a real fit**: the first
draft of `ptModelFromQuaids()` built row names only from the blocks
`quaidsElas_()` reads (`qOut.xnam` | `GAMMA_`+`qOut.wnam` | `BETA_LX` |
`LAMBDA_LX2`) — but `quaidsElas_()` only *reads* the first
`1+nint+n+nendog` rows of `b` and silently ignores the rest; `qOut.bestB`
itself always has `qOut.nu` additional trailing rows for the IV
control-function residual coefficient(s). Running the adapter against a
real `quaidsFit()` QUAIDS/IV result immediately threw
`ptModelSetNames: termNames must contain 10 labels` (9 built vs. 10
actual) — caught by execution, not by re-deriving the formula on paper.
Fixed by appending `qOut.unam`-derived names for the trailing block, which
also makes the table more complete: `printQuaids()`'s own console report
already shows these as the "Residuals of instrumental regressions" row, so
a pubtable coefficient report should too. Re-verified against LA-AIDS
(`linear=1`) and QUAIDS (`linear=0`), both `homogenous=1` and
`homogenous=0`, to confirm the row-count formula holds across all four
combinations.

## Milestone 7: package build and release tooling

```powershell
# From the repo root, after any src/ change:
powershell scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact

# Or step by step:
powershell tests\run_source_tests.ps1                       # manifest check + all 7 tgauss tests
powershell scripts\build_package.ps1 -Force                 # writes "quaids <version>.zip", self-verifies
powershell scripts\run_release_verification.ps1 -InstallArtifact -SkipInstalledPackageTest  # installs to <GaussHome>/pkgs/quaids
tgauss -b -x package_public_api.e                            # from tests/, exercises `library quaids;`
```

Adapted from `gauss-qardl`'s `scripts/`/`tests/` release tooling, scaled
down for this repo's current scope — see `GOLD_STANDARD_TODO.md`'s
Milestone 7 section for the full list of what was intentionally omitted
(separate benchmark/examples-smoke scripts; this repo's equivalent
validation already lives in `run_source_tests.ps1`'s 7 tgauss files) and
what's deferred to Milestone 8 (`verify_release_artifact.ps1`'s
`requiredEntries` list doesn't require `README.md`/
`docs/COMMAND_REFERENCE.md` yet, since neither exists).

**The package is now really installed** at `c:\gauss26\pkgs\quaids`,
verified via `library quaids;` — not just staged/zipped. This was an
explicit repo-owner decision (asked via the same "writing outside this git
repo into shared machine state" reasoning as the Milestone 6 adapter-
location question), since fully validating the release pipeline requires
`library quaids;` to actually resolve, which only happens against GAUSS's
real, shared package directory.

**Three real bugs found by actually building, installing, and running the
package — not by re-reading the scripts**, the same "never trust a
derived formula/script without running it" standard this whole project has
followed since Milestone 3:

1. **`build_package.ps1`'s cleanup step deleted the entire staged
   `examples/`/`tests/` directories.** `Get-ChildItem -LiteralPath <dir>
   -Include <patterns> -Recurse` silently ignores `-Include` when the base
   path isn't itself a wildcard (a known PowerShell footgun) — `-Recurse`
   returned every file with no effective filter, and `Remove-Item -Force`
   deleted all of them. `verify_release_artifact.ps1` caught it
   immediately: "missing required entry: examples/quaids_example.e" on the
   first real run. Fixed by removing named generated files by explicit
   literal path, plus a `-Filter` (not `-Include`) pass for `*.log`, which
   doesn't have this bug.
2. **`build_lcg.ps1`'s proc-detection regex only matched one of the three
   GAUSS proc-declaration forms this codebase actually uses.** It matched
   `proc (struct X) = name(...)`, but not the bare-digit form `proc N =
   name(...)` (e.g. `proc 3 = quaidsElas_(...)`) or the no-return-spec
   form `proc name(...)` (e.g. `proc quantile(x, s);`). This silently
   dropped `quaidsSlutzky`, `_quaidsIVFirstStage`, `quaids()`,
   `printQuaids`, `quaidsElas_`, `printQuaidsElas`, `quaidsElas`, and a
   private `quantile` helper from the generated `.lcg` catalog — invisible
   to every `#include`-based test (which never goes through the catalog),
   surfaced only as `Undefined symbol` errors from `library quaids;` when
   `tests/package_public_api.e` called them for real. Fixed the regex to
   match all three forms; verified the regenerated catalog lists every
   proc in every `src/` file.
3. **A genuinely pre-existing bug in `src/quaids.src`, dating from before
   Milestone 0**: a private `quantile(x, s)` helper duplicating GAUSS's
   builtin `quantile()`. The original author clearly intended to delete
   it — it sat inside a comment block — but GAUSS comments do not nest,
   and the wrapping comment's own doc-header used an inner
   `/**...**/`-style block whose closing marker closed the *outer*
   "delete this" comment early, leaving the entire proc live and
   uncommented. Silently duplicated the builtin's behavior at its 3 call
   sites (`quaids()`'s legacy elasticities-at-four-points block) for as
   long as this repo has existed, invisible via `#include`-based
   compilation (which just locally shadows the builtin name) — but a real
   GAUSS builtin name cannot be redefined via `library`-based lazy
   loading, so it surfaced as `Undefined symbol: 'quantile'` resolving
   `quaids.src`'s *own* proc definition, once bug #2's fix let `library
   quaids;` actually try to load it. Fixed by deleting the dead code
   outright (the original author's evident intent), not just repairing
   the comment nesting — the 3 call sites now resolve to the GAUSS
   builtin, same as `quaidsslutzky.src`'s identical calls always have.
   Full 7-test regression suite re-run afterward; no assertion anywhere
   depended on the old local implementation's exact output (only the
   legacy, unasserted quartile-elasticities console block uses it), so no
   test needed updating.

Also hit and fixed a self-inflicted issue while writing bug #3's own
explanatory comment: GAUSS's lexer does not tolerate literal `"`
characters inside a `/* ... */` block comment — an odd count throws
`error G0097 String not closed`, even though the characters are inside a
comment, not a string. Caught immediately by re-running the test suite;
fixed by rephrasing without quote characters. Worth remembering when
writing GAUSS comments generally, not just here.

**No version bump for this milestone**: build/release tooling doesn't
change GAUSS public API surface, so per this repo's established policy
(version bumps track public API changes, not every milestone — see
Milestone 3's R/Python scripts) the package stays at `0.5.0`. Not tagging
a git release either: nothing in this repo has been committed yet as of
Milestone 7 (every milestone's changes are staged but uncommitted, per the
"never commit unless asked" policy), and a tag needs a commit to point at
— `CHANGELOG.md`/`package.json` versioning infrastructure is in place for
whenever the repo owner chooses to commit and tag.

## Milestone 8: documentation

Full doc set added: `README.md`, `docs/COMMAND_REFERENCE.md` plus 18
`docs/command-reference/*.md` pages (one per public proc, Purpose/Format/
Parameters/Returns/Remarks/Examples/Source/See Also template, matching
`gauss-qardl`'s page style), `docs/USAGE_GUIDE.md`, `docs/METHODOLOGY_NOTES.md`,
`docs/FEATURE_SUPPORT_MATRIX.md`. `CLAUDE.md` itself was already complete
and kept synchronized with every milestone since Milestone 0 — nothing
further needed there.

**Coverage includes the optional `pubtable` adapter**: all 12 procs across
`quaids.sdf`'s load-bearing `src/` files, plus all 6 procs in
`src/pubtable_quaids.src` (not in `package.json`'s `src` array, but real,
working, documented public API surface — the same reasoning `README.md`
and `USAGE_GUIDE.md` already apply to the optional `pubtable` export
path).

**Followed through on two Milestone-7-deferred TODOs rather than leaving
them stale**, both explicitly flagged in their own header comments at the
time:

- `tests/verify_package_manifest.ps1` now cross-checks
  `docs/COMMAND_REFERENCE.md` against the actual source the same way
  `gauss-qardl`'s version does: every documented proc must be defined
  somewhere in `src/` (using the same intentionally-unlisted allowlist the
  `package.json`-`src`-array check already has, so `pubtable_quaids.src`'s
  documented procs don't spuriously fail), and every linked
  command-reference page must exist. **Verified this actually catches a
  real problem, not just that it runs** — the same "never trust a check
  without testing it" standard as everywhere else in this project:
  temporarily renamed one documented link to a nonexistent proc name,
  confirmed the script threw a clear, correctly-identifying error, then
  reverted.
- `scripts/verify_release_artifact.ps1`'s `requiredEntries` list now
  includes `README.md` and all four top-level `docs/*.md` files (previously
  deliberately omitted, per that script's own Milestone 7 comment, since
  they didn't exist yet).

**Verification**: rebuilt the release artifact (`scripts/build_package.ps1
-Force`), confirmed `docs/` and its 18 command-reference pages are
actually present inside the `.zip` (not just referenced by the build
script's directory list), then re-ran the full pipeline
(`scripts/run_release_verification.ps1 -InstallArtifact`) end to end — all
7 source-tree test files (150 checks), the extended manifest check, a real
reinstall to `c:\gauss26\pkgs\quaids`, and `tests/package_public_api.e`
against that install all still pass.

**No version bump**: documentation doesn't change GAUSS public API
surface, so per this repo's established policy the package stays at
`0.5.0`.

## Milestone 9: final gold standard integration gate

A genuine gate, not a rubber stamp: running it top-to-bottom surfaced
three real, previously-undetected gaps, each found by actually exercising
the system (not by re-reading it) — the same standard applied at every
prior milestone.

**Examples didn't match what the docs promised.** README.md,
`docs/USAGE_GUIDE.md`, and every command-reference page (all written at
Milestone 8) document `library quaids;` as the primary usage pattern,
since the package became genuinely installable at Milestone 7. But
`examples/quaids_example.e` and `examples/pubtable_export_example.e`
still `#include`d the source tree directly — a leftover from before
Milestone 7. Fixed: both now use `library quaids;`.
`pubtable_export_example.e` additionally needs a bare `#include
quaids.sdf` — confirmed empirically (not assumed) that `library quaids;`
alone does not activate `pubtable_quaids.src`'s `#ifDef
QUAIDS_SDF_INCLUDED` guard, since `library` lazily compiles individual
procs on demand rather than eagerly running `quaids.sdf`'s `#define`;
explicitly including the `.sdf` (resolved via the installed package's
search path) does activate it — matching exactly what `pubtable`'s own
bundled `pubtable_qardl.src` documents as required for `qardl.sdf` ("just
include qardl.sdf first"). Both examples re-verified against the freshly
rebuilt/reinstalled package.

**The installed-package gate didn't exercise two of its own required
procs.** Auditing "package exports match public docs" line by line found
`tests/package_public_api.e` never actually called `printQuaids()` or
`quaidsElas()` — both real, required, `package.json`-listed public procs,
previously only exercised indirectly via their split components
(`quaidsFit`, `quaidsElasFit`+`printQuaidsElas`). A load-order bug or
stale `.lcg` entry specific to either one could have passed this gate
undetected. Fixed by adding direct calls to both; re-ran the full pipeline
to confirm.

**Published-data validation only covered one of three model families.**
Extended `tests/quaids_published_validation_test.e` (11 -> 19 checks) to
also validate iterated AIDS (`aCtl.linear=1, aCtl.maxiter>1`) against R
`micEconAids::aidsEst(method="IL", ...)` — the Iterated Linear Least
Squares Estimator (Blundell & Robin 1999), which uses the same
starting-value/iteration structure `quaidsFit()` does (LA-AIDS starting
point, then iterate with the translog price index). Observed max absolute
difference ~0.11 (vs. ~0.021 for the existing LA-AIDS/3SLS check),
tolerance set to `0.15`. The wider gap has a real, understood cause, not
approximation slop: `micEconAids`'s `method="IL"` does not support
instrumental variables — confirmed by direct testing, not assumed:
combining `method="IL"` with `instNames` **segfaults** R's `aidsEst`
rather than erroring cleanly. So the R reference here is SUR-only, while
GAUSS's iterated fit always instruments log total expenditure — the
comparison spans both a different algorithm and an IV-vs-no-IV
difference. **QUAIDS has no independent reference implementation
available**: `micEconAids` does not implement a quadratic log-expenditure
term at all, and no other comparably-established QUAIDS implementation
exists (consistent with the Milestone 3 search that produced only a
from-scratch Python replica, itself kept as supplementary evidence, not a
QUAIDS reference). QUAIDS's validation remains the known-true
synthetic-DGP recovery in `tests/quaids_synthetic_validation_test.e` — a
real, non-circular check, but a different, weaker tier of evidence than
cross-implementation agreement on real published data. Documented
explicitly in `docs/FEATURE_SUPPORT_MATRIX.md` rather than silently
equated with the other two model families.

**Two stale Definition-of-Done aspirations reconciled against the actual,
deliberate, already-shipped design**, in `GOLD_STANDARD_TODO.md`: (1) "LA-
AIDS, iterated AIDS, and QUAIDS... with and without IV" assumed an
exogenous-total-expenditure mode that was deliberately never built —
`instr` has been a required argument since Milestone 2/3 by design, not
oversight. (2) "Homogeneity, symmetry, and overidentification each have
standalone... procs" is literally false for two of the three: only
homogeneity and joint homogeneity+symmetry are separately-callable procs
(`quaidsHomogeneityTest`/`quaidsJointTest`); symmetry-given-homogeneity
and overidentification are automatically computed and reported as part of
`quaidsFit()` itself (`qOut.symStat`/`qOut.symPval`,
`qOut.overidValid`/...), a deliberate Milestone 4 scoping decision. Both
aspirations were corrected to match what was actually built and tested,
rather than left as checkboxes nothing could ever satisfy.

**Version consistency verified directly**: `package.json` (`0.5.0`), the
rebuilt artifact filename (`quaids 0.5.0.zip`), the installed copy's
`package.json` (`0.5.0`), and `CHANGELOG.md`'s `## 0.5.0` entry all agree
— read each one directly rather than trusting
`verify_release_artifact.ps1`'s own pass alone, since this is the final
gate. `CHANGELOG.md`'s `0.5.0` entry gained `### Changed`/`### Fixed`
sections documenting Milestones 7-9's tooling/docs/test work (no version
bump -- no `src/` public API changed across any of it).

**Full pipeline re-run end to end after every fix above**: all 7
source-tree test files (158 checks total — schema 34, formula parity 17,
synthetic validation 22, published validation 19, hypothesis tests 19,
elasticities 17, pubtable adapter 30), the extended manifest check, a
real rebuild/reinstall to `c:\gauss26\pkgs\quaids`,
`tests/package_public_api.e` against that install, and both examples run
directly against the installed package.

## What GAUSS already provides — do not duplicate

Full detail and evaluation status is in `GOLD_STANDARD_TODO.md` under "What
GAUSS Already Provides." Summary:

- **No built-in SUR/systems-of-equations estimator and no existing
  AIDS/QUAIDS/demand-system implementation** anywhere in the GAUSS runtime or
  shipped packages (`pkgs/`). The cross-equation-restricted iterated FGLS
  core here is genuinely new — this is the library's reason to exist.
- **`gmmFitIV`** (`gmm.sdf`/`gmm_est.src`/`gmm_hac.src`) — evaluated at
  Milestone 2, **decision: not adopted**. It's a single-equation IV-GMM
  estimator (`y` is `N x 1`); under its `"onestep"`/`"unadj"` settings it's
  mathematically identical to the classical 2SLS the hand-rolled first stage
  already computes, so adopting it wouldn't change any estimates. It also
  doesn't expose `zzi` (`inv(Z'Z)`) or the full `[Z, endog]` moment matrix in
  the layout `quaidsFit()`'s later system-wide covariance and
  overidentification-test formulas need (`D*zzi*D'`, `z1z1i`, `zz1`, ...) —
  adopting it would mean unpacking a `gmmOut` struct to reconstruct those
  building blocks anyway, adding a package dependency for zero net
  simplification. The extraction that *did* happen (`quaidsiv.src`,
  Milestone 2) keeps the hand-rolled computation but gives it its own file.
- **`pubtable`** (`pkgs/pubtable`) is a complete LaTeX/HTML/RTF/CSV/XLS/
  Markdown table-export engine, already used by `gauss-qardl` via a small
  adapter bundled inside pubtable itself (`pubtable_qardl.src`, guarded by
  `#ifDef QARDL_SDF_INCLUDED`). **Adopted at Milestone 6**: this repo's own
  `src/pubtable_quaids.src` follows the same pattern (guard, `ptModelFromX`/
  `ptFromX` shape) but lives in this repo rather than physically inside the
  installed `pubtable` package — see "Milestone 6: reporting via pubtable"
  above for why.
- **`loadd()`/`asdf()`/dataframe column selection** — used at Milestone 2 for
  `quaidsFull()` (see below). **Not** a `"y ~ x1 + x2"` formula string:
  AIDS/QUAIDS is a multi-equation system (N shares against N parallel
  prices), which doesn't fit GAUSS's single-equation formula grammar. Column
  name lists (`data[., stringArrayOfNames]`), matched positionally, are the
  natural fit instead.
- **Primitives already used correctly and worth keeping**: `moment()`,
  `invpd()`, `solpd()`, `design()`/`vech()`/`xpnd()` (symmetric-matrix
  restriction algebra), `eigh()` (Slutzky check), `cdfchic`/`cdftc`/`cdffc`/
  `cdfnc`, `printfm()`, `quantile()`.

## GAUSS language conventions observed in this codebase

- **Variable naming**: short lowercase names (`w`, `n`, `nz`, `ge`, `gg`,
  `b`, `u`), matching the original author's terse econometrics-code style —
  preserve it rather than renaming to verbose identifiers.
- **Locals**: all local variables declared in one `local` statement at the
  top of each proc.
- **Struct declarations inside procs**: `struct quaidsControl aCtl;` appears
  as a formal parameter, not in the `local` list — standard GAUSS syntax.
  `struct quaidsOut qOut;` (a local, not a parameter) is declared the same
  way, separately from the `local` statement, inside `quaidsFit()`/`quaids()`.
- **Character-matrix name vectors vs. the native `string array` type**: name
  vectors built with the classic `0$+"X"$+ftocv(...)` idiom (`xnam`, `wnam`,
  `znam`, `unam`, `enam`) are legacy character matrices, not the newer
  `string array` type — struct fields holding them must be declared `matrix`,
  not `string array`, or you get a `G0071 Type mismatch` at assignment. Hit
  and fixed this exact error while building `quaidsOut` at Milestone 1. The
  *reverse* direction (legacy char-matrix into a native `string`/
  `string array`-typed field, e.g. a `pubtable` `ptModel`/`ptTable`) throws
  the same error and has no dedicated conversion builtin — see "Milestone
  6: reporting via pubtable" above for the `$|`-concatenation idiom that
  does work.
- **Loop style**: `do while ok; ... endo;` / `do while i<=n; ... i=i+1;
  endo;` — not GAUSS's newer `for` loop syntax.
- **Symmetric-restriction idiom**: `design(vec(xpnd(seqa(1,1,k*(k+1)/2))))`
  builds the selection matrix `R` such that `vec(G) = R*vech(G)` for
  symmetric `G` — this is how homogeneity/symmetry constraints get imposed
  via minimum distance. Reuse this idiom rather than re-deriving it.
- **Matrix concatenation**: `~` horizontal, `|` vertical, as usual in GAUSS.
- **Relative vs. absolute prices**: `prices` is converted to relative
  (`prices[.,1:n-1] - prices[.,n]`) near the top of `quaidsFit()`, then
  converted back to absolute before the recovery step at the end — mutating
  the input matrix in place *within `quaidsFit()`'s own local scope*. Because
  GAUSS passes matrix arguments by value, this mutation does **not** leak
  back into `quaids()`'s own `prices` local — the legacy wrapper's `prices`
  stays in its original absolute-price form throughout, which is exactly
  what lets `quaids()` reuse its own `prices`/`totexp`/`w`/`instr` arguments
  unchanged for the elasticities/descriptive-stats/Slutzky calls after
  `quaidsFit()` returns. `intcpt` is different: `quaidsFit()`'s internal
  constant-column-prepend is a real structural change (not a round-trip), so
  `quaidsFit()` returns the mutated version as `qOut.intcptFull` for
  `quaids()` to reuse.

## Testing status

Seven automated tests exist, all run from `tests/` as the working directory:

```
tgauss -b -x quaids_schema_test.e
tgauss -b -x quaids_formula_parity_test.e
tgauss -b -x quaids_synthetic_validation_test.e
tgauss -b -x quaids_published_validation_test.e
tgauss -b -x quaids_hypothesis_tests_test.e
tgauss -b -x quaids_elasticities_test.e
tgauss -b -x quaids_pubtable_test.e
```

- `quaids_schema_test.e` (Milestone 1, 34 checks): asserts `quaidsOut` field
  values/shapes, that `quaidsFit()` prints nothing between call and return,
  and that the legacy `quaids()` wrapper's four returned matrices are
  exactly (not approximately) equal to the `quaidsOut` fields they're drawn
  from.
- `quaids_formula_parity_test.e` (Milestone 2, 17 checks): asserts
  `quaidsFull(dataframe, ...)` produces numerically identical output to
  `quaidsFit(matrices...)` on the same underlying data, including the
  `extraVars == 0` path.
- `quaids_synthetic_validation_test.e` (Milestone 3, 22 checks): asserts
  `quaidsFit()` recovers true DGP parameters within a documented tolerance
  across LA-AIDS/iterated-AIDS/QUAIDS x with/without-IV. See "Milestone 3:
  synthetic validation findings" above for the tolerances and the seed
  sensitivity finding this surfaced.
- `quaids_published_validation_test.e` (Milestone 3, 19 checks): asserts
  `quaidsFit()` on real published data (`Blanciforti86`) matches an
  independent R reference within tolerance, plus adding-up/homogeneity/
  symmetry sanity checks. This is the test that caught the Stone-index
  starting-value bug — see "Milestone 3: real bug found and fixed" above.
  Milestone 9 extended it (11 -> 19 checks) to also validate iterated
  AIDS against R `aidsEst(method="IL", ...)` — see "Milestone 9" above.
- `quaids_hypothesis_tests_test.e` (Milestone 4, 19 checks): size and power
  for `quaidsHomogeneityTest()`/`quaidsJointTest()`, a power check for the
  existing symmetry-given-homogeneity test, and the first-ever exercise of
  the overidentification test. See "Milestone 4: new hypothesis tests"
  above.
- `quaids_elasticities_test.e` (Milestone 5, 17 checks): parity between
  `quaidsElasFit()`/`printQuaidsElas()` and `quaidsElas_()`, plus three
  exact algebraic identities checked at points other than the four
  standard ones. See "Milestone 5: elasticities generalization" above.
- `quaids_pubtable_test.e` (Milestone 6, 30 checks; requires the `pubtable`
  package installed): exact numeric parity between `pubtable`
  `ptModel.estimates`/`stdErrors` and the `qOut`/`elasOut` values they're
  built from, shape/title checks, the `ptFromQuaidsFamily` dispatcher, and
  an end-to-end export smoke test that writes real `.tex`/`.md`/`.csv`
  files and reads them back. See "Milestone 6: reporting via pubtable"
  above.

All seven print one `PASS`/`FAIL` line per check and a final `...: ALL N
CHECKS PASSED` (or a failure count) summary line — check that line, since
`tgauss`'s exit code is not currently a reliable pass/fail signal for this
harness. `tests/run_source_tests.ps1` (Milestone 7) runs
`verify_package_manifest.ps1` plus all 7 of these in one shot and checks
this same summary-line convention (not just GAUSS-level compile/execute
errors).

An eighth test, `tests/package_public_api.e` (Milestone 7), is different
in kind from the seven above: it loads `library quaids;` against a real
*installed* copy of the package (currently `c:\gauss26\pkgs\quaids`) rather
than `#include`-ing the source tree, so it only runs correctly after
`scripts/run_release_verification.ps1 -InstallArtifact` (or equivalent)
has installed the package — it is a release gate, not a routine
`run_source_tests.ps1` member. See "Milestone 7: package build and release
tooling" above.

`examples/quaids_example.e` remains a manual, eyeball-comparison smoke
script (no assertions) — superseded for correctness-checking purposes by
`quaids_synthetic_validation_test.e`, but kept as a simple, readable
end-to-end usage example. `examples/pubtable_export_example.e` (Milestone
6) is the same style for the reporting layer — exports a coefficient table
and elasticity tables to `.tex`/`.md`/`.csv`, requires `pubtable` installed.
Both use `library quaids;` against the installed package (Milestone 9),
matching README.md/USAGE_GUIDE.md's documented usage pattern — see
"Milestone 9" below.

To run the examples from a GAUSS 26 console/batch shell:

```
tgauss -b -x examples/quaids_example.e
tgauss -b -x examples/pubtable_export_example.e
```

`quaids_example.e` no longer depends on its working directory, since it
loads the installed package via `library quaids;` rather than `#include`.
`pubtable_export_example.e` still needs `examples/` as the working
directory (or paths adjusted), since it `#include`s
`../src/pubtable_quaids.src` by relative path — that adapter is not part
of the installed package (see "Milestone 6" above).

## Package manifest

`package.json` lists (relative to `src/`, in load order): `quaids.sdf`,
`quaidsutil.src`, `quaidsiv.src`, `quaidselas.src`, `quaidsslutzky.src`,
`quaids.src`, `quaidsformula.src`, `quaidstests.src`. `quaids.src` must load
after `quaidsiv.src`/`quaidselas.src`/`quaidsslutzky.src` since it calls
procs they define; `quaidsformula.src` must load after `quaids.src` since
`quaidsFull()` calls `quaidsFit()`. `quaidstests.src` has no load-order
dependency on the others beyond `quaids.sdf` (it only reads an
already-computed `quaidsOut`). `src/pubtable_quaids.src` (Milestone 6) is
deliberately **not** in this array — it has a hard dependency on
`pubtable.sdf`'s struct types, and adding it would make `pubtable` a hard
dependency for the whole package to compile; see "Milestone 6: reporting
via pubtable" above. **Buildable/installable as of Milestone 7**: `scripts/
build_lcg.ps1`/`build_package.ps1`/`run_release_verification.ps1` build a
release `.zip` and can install it to a real GAUSS package directory
(`c:\gauss26\pkgs\quaids` on this machine) so `library quaids;` works — see
"Milestone 7: package build and release tooling" above. If you add a new
*required* `.src` file, add it to the `"src"` array (respecting load
order), bump the version, and rebuild/reinstall (`scripts/
run_release_verification.ps1 -BuildArtifact -ForceArtifact
-InstallArtifact`) so the installed copy and its `.lcg` catalog stay in
sync with the source tree; optional/adapter files like
`pubtable_quaids.src` still warrant a version bump (real new public API
surface) but stay out of `src`.

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312–326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527–539.
- Blanciforti, L., Green, R., King, G. (1986). *U.S. Consumer Behavior Over
  the Postwar Period: An Almost Ideal Demand System Analysis*. Giannini
  Foundation Monograph No. 40. Published-replication dataset (as
  `Blanciforti86` in the R package `micEconAids`), used by
  `tests/quaids_published_validation_test.e` — see
  `tests/fixtures/published/SOURCE.md`.
- Aptech GAUSS coding conventions: https://github.com/aptech/gauss-llm-reference
- Sibling library for structure/process conventions: `gauss-qardl`
  (`CLAUDE.md`, `GOLD_STANDARD_TODO.md`, `docs/command-reference/`,
  `pubtable_qardl.src` adapter pattern).
