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

The library is **pre-alpha** (package version `0.1.0`) and is not yet
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

## Repository layout (post-Milestone-1)

```
src/
  quaids.sdf      # Struct definitions: quaidsControl and quaidsOut.
  quaidsutil.src  # quaidsControlCreate() / getDefaultQuaidsControl().
  quaids.src      # quaidsFit() (silent, struct-returning estimation core),
                  #   printQuaids() (the separated estimation-report
                  #   printer), quaids() (backward-compatible wrapper: fits,
                  #   prints, reproduces the legacy elasticities/descriptive
                  #   stats/Slutzky report, returns the 4 legacy matrices),
                  #   plus quaidsSlutzky(), quaidsElas_(), quaidsElas().
examples/
  quaids_example.e  # One synthetic 5-good dataset (homogeneity/symmetry true
                  #   by construction), run through quaids() with eyeballed
                  #   comparison of printed estimates to true parameters.
                  #   Not an automated test — no assertions. Includes
                  #   src/*.sdf and src/*.src directly via relative paths
                  #   (../src/...), matching gauss-qardl's source-tree
                  #   testing convention, since there is no installable
                  #   package build yet.
tests/
  quaids_schema_test.e  # Milestone 1 schema test: asserts quaidsOut field
                  #   values/shapes, that quaidsFit() prints nothing, and
                  #   that the legacy quaids() wrapper's returned matrices
                  #   are byte-identical to the struct fields they're drawn
                  #   from. Run from tests/ as the working directory.
package.json      # GAUSS package manifest (name: quaids, version: 0.1.0,
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

Milestones 0 (repo hygiene) and 1 (API/output-schema baseline) are complete.
`docs/` does not exist yet — that is Milestone 8. A fuller `tests/` harness
(installed-package tests, more fixtures) is Milestone 7.

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
`(qOut.b, qOut.v, qOut.bS, qOut.vS)`. **Verified byte-for-byte**: running
`quaids()` post-refactor against the same fixed-seed synthetic dataset used
at Milestone 0 produces output identical to the pre-Milestone-1 code,
including the per-iteration convergence log (reproduced from a stored
`qOut.iterHistory` matrix rather than printed live during iteration, since
`quaidsFit()` cannot print, but reproduced with matching iteration
numbers/error values so the printed table is unchanged).

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
- **`gmmFitIV`** (`gmm.sdf`/`gmm_est.src`/`gmm_hac.src`) is a full
  single-equation IV-GMM estimator with robust/HAC weighting. It is a
  candidate to replace the hand-rolled first-stage 2SLS block in `quaids()`
  (the `moment()`/`invpd()`/`solpd()` code near the top of the proc), but
  this needs an explicit evaluation, not an assumption — the system estimator
  needs raw residuals/moment blocks in a specific layout.
- **`pubtable`** (`pkgs/pubtable`) is a complete LaTeX/HTML/RTF/CSV/XLS/
  Markdown table-export engine, already used by `gauss-qardl` via a small
  adapter (`pubtable_qardl.src`, guarded by `#ifDef QARDL_SDF_INCLUDED`). Any
  export/reporting work should add `pubtable_quaids.src` the same way, not a
  bespoke exporter.
- **`loadd()`/`asdf()`/formula strings** should back a future dataframe/
  formula entry point (à la `qardl`'s `applyQARDLFormula`), instead of a
  hand-written formula parser.
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

`tests/quaids_schema_test.e` is the first automated test: it asserts
`quaidsOut` field values/shapes, that `quaidsFit()` prints nothing between
call and return, and that the legacy `quaids()` wrapper's four returned
matrices are exactly (not approximately) equal to the `quaidsOut` fields
they're drawn from. Run it from `tests/` as the working directory:

```
tgauss -b -x quaids_schema_test.e
```

It prints one `PASS`/`FAIL` line per check and a final
`SCHEMA TEST: ALL N CHECKS PASSED` (or a failure count) summary line — check
that line, since `tgauss`'s exit code is not currently a reliable pass/fail
signal for this harness.

`examples/quaids_example.e` remains a manual, eyeball-comparison smoke
script (no assertions) — see `GOLD_STANDARD_TODO.md` Milestone 3 for the plan
to add deterministic, published-benchmark fixtures with real assertions
beyond the current schema/shape-level checks.

To run the example from a GAUSS 26 console/batch shell:

```
tgauss -b -x examples/quaids_example.e
```

(must be run with `examples/` as the working directory, or with paths
adjusted, since the example includes source files via `../src/...`).

## Package manifest

`package.json` lists `quaids.sdf`, `quaidsutil.src`, `quaids.src` (relative
to `src/`) in load order (structs first). The package is not yet
buildable/installable as a GAUSS package (`.lcg`) — that is Milestone 7. If
you add a new `.src` file, add it to the `"src"` array and bump the version.

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312–326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527–539.
- Aptech GAUSS coding conventions: https://github.com/aptech/gauss-llm-reference
- Sibling library for structure/process conventions: `gauss-qardl`
  (`CLAUDE.md`, `GOLD_STANDARD_TODO.md`, `docs/command-reference/`,
  `pubtable_qardl.src` adapter pattern).
