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

The library is **pre-alpha** (package version `0.2.0`) and is not yet
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

## Repository layout (post-Milestone-2)

```
src/
  quaids.sdf        # Struct definitions: quaidsControl and quaidsOut.
  quaidsutil.src    # quaidsControlCreate() / getDefaultQuaidsControl().
  quaidsiv.src      # _quaidsIVFirstStage() -- private helper, the
                    #   instrumental-variables first-stage regression of log
                    #   total expenditure. The one internal phase of
                    #   quaidsFit() that was cleanly separable; see the
                    #   Milestone 2 scoping note below for why the rest of
                    #   quaidsFit() (starting values, iteration, variance,
                    #   overidentification test, symmetry test) was not
                    #   further split.
  quaidselas.src    # quaidsElas_(), quaidsElas() -- elasticities at a point.
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
examples/
  quaids_example.e  # One synthetic 5-good dataset (homogeneity/symmetry true
                    #   by construction), run through quaids() with
                    #   eyeballed comparison of printed estimates to true
                    #   parameters. Not an automated test — no assertions.
                    #   Includes src/*.sdf and src/*.src directly via
                    #   relative paths (../src/...), matching gauss-qardl's
                    #   source-tree testing convention, since there is no
                    #   installable package build yet.
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
                    #   real bug found and fixed" below.
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
                    #   All test/fixture files run from tests/ (or
                    #   tests/fixtures/published/ for the R/Python
                    #   scripts) as the working directory.
package.json      # GAUSS package manifest (name: quaids, version: 0.2.0,
                  #   license: MIT).
LICENSE           # MIT, copyright Eric Clower.
CITATION.cff      # Citation metadata; cites Deaton & Muellbauer (1980) and
                  #   Banks, Blundell & Lewbel (1997).
.gitignore        # Compiled .gcg artifacts, tmp/, .claude/, packaged zips,
                  #   generated `output file=...` run artifacts.
GOLD_STANDARD_TODO.md  # Living roadmap: release blockers, milestones,
                  #   definition of done. Read this before any nontrivial
                  #   change and update it as milestones close.
```

Milestones 0 (repo hygiene), 1 (API/output-schema baseline), 2 (modular
source split + dataframe entry point), and 3 (validation fixtures,
including published-data cross-implementation validation) are all
complete. `docs/` does not exist yet — that is Milestone 8. A fuller
`tests/` harness (installed-package tests) is Milestone 7.

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
  adapter (`pubtable_qardl.src`, guarded by `#ifDef QARDL_SDF_INCLUDED`). Any
  export/reporting work should add `pubtable_quaids.src` the same way, not a
  bespoke exporter. Still not done — candidate for Milestone 6.
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
  and fixed this exact error while building `quaidsOut` at Milestone 1.
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

Four automated tests exist, all run from `tests/` as the working directory:

```
tgauss -b -x quaids_schema_test.e
tgauss -b -x quaids_formula_parity_test.e
tgauss -b -x quaids_synthetic_validation_test.e
tgauss -b -x quaids_published_validation_test.e
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
- `quaids_published_validation_test.e` (Milestone 3, 11 checks): asserts
  `quaidsFit()` on real published data (`Blanciforti86`) matches an
  independent R reference within tolerance, plus adding-up/homogeneity/
  symmetry sanity checks. This is the test that caught the Stone-index
  starting-value bug — see "Milestone 3: real bug found and fixed" above.

All four print one `PASS`/`FAIL` line per check and a final `...: ALL N
CHECKS PASSED` (or a failure count) summary line — check that line, since
`tgauss`'s exit code is not currently a reliable pass/fail signal for this
harness.

`examples/quaids_example.e` remains a manual, eyeball-comparison smoke
script (no assertions) — superseded for correctness-checking purposes by
`quaids_synthetic_validation_test.e`, but kept as a simple, readable
end-to-end usage example.

To run the example from a GAUSS 26 console/batch shell:

```
tgauss -b -x examples/quaids_example.e
```

(must be run with `examples/` as the working directory, or with paths
adjusted, since the example includes source files via `../src/...`).

## Package manifest

`package.json` lists (relative to `src/`, in load order): `quaids.sdf`,
`quaidsutil.src`, `quaidsiv.src`, `quaidselas.src`, `quaidsslutzky.src`,
`quaids.src`, `quaidsformula.src`. `quaids.src` must load after
`quaidsiv.src`/`quaidselas.src`/`quaidsslutzky.src` since it calls procs
they define; `quaidsformula.src` must load after `quaids.src` since
`quaidsFull()` calls `quaidsFit()`. The package is not yet
buildable/installable as a GAUSS package (`.lcg`) — that is Milestone 7. If
you add a new `.src` file, add it to the `"src"` array (respecting load
order) and bump the version.

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
