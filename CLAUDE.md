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

## Repository layout (post-Milestone-0)

```
src/
  quaids.sdf      # Struct definition: quaidsControl (renamed from
                  #   aidsControl). The dead rankControl/latentControl
                  #   structs from a copy-pasted template were removed here.
  quaidsutil.src  # quaidsControlCreate() (renamed from aidsControlCreate()).
                  #   The dead rankControlCreate()/latentControlCreate()
                  #   constructors were removed here.
  quaids.src      # The entire estimation library: quaids() (renamed from
                  #   aids()), plus three helper procs it calls:
                  #   quaidsSlutzky(), quaidsElas_(), quaidsElas() (renamed
                  #   from slutzky(), elas_(), elas()).
examples/
  quaids_example.e  # One synthetic 5-good dataset (homogeneity/symmetry true
                  #   by construction), run through quaids() with eyeballed
                  #   comparison of printed estimates to true parameters.
                  #   Not an automated test — no assertions. Includes
                  #   src/*.sdf and src/*.src directly via relative paths
                  #   (../src/...), matching gauss-qardl's source-tree
                  #   testing convention, since there is no installable
                  #   package build yet.
package.json      # GAUSS package manifest (name: quaids, version: 0.1.0,
                  #   license: MIT).
LICENSE           # MIT, copyright Eric Clower.
CITATION.cff      # Citation metadata; cites Deaton & Muellbauer (1980) and
                  #   Banks, Blundell & Lewbel (1997).
.gitignore        # Compiled .gcg artifacts, tmp/, .claude/, packaged zips.
GOLD_STANDARD_TODO.md  # Living roadmap: release blockers, milestones,
                  #   definition of done. Read this before any nontrivial
                  #   change and update it as milestones close.
```

Milestone 0 (repo hygiene) is complete: dead code removed, `src/`/`examples/`
split done, package/proc naming decided (`quaids`), license decided (MIT).
`docs/` and `tests/` directories do not exist yet — those are Milestones 3
and 8.

## The `quaids()` proc

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, struct quaidsControl aCtl);
```

- `w` — `TxN` budget shares.
- `intcpt` — `TxK` extra intercept-shifter variables (demographics etc.), or
  `0` for none.
- `prices` — `TxN` **log** prices (absolute, not relative — the proc converts
  internally).
- `totexp` — `Tx1` log total expenditure (treated as endogenous).
- `instr` — `TxH` instruments for log total expenditure.
- `aCtl` — `quaidsControl` struct (see below).

Returns, if `aCtl.homogenous == 1`: `b1`/`v1` = homogeneity-constrained
estimates/covariance, `b2`/`v2` = homogeneity+symmetry-constrained
estimates/covariance. If `aCtl.homogenous == 0`: `b1`/`v1` are the
unconstrained (reparametrized) estimates and `b2`/`v2` are `0`.

**This single call also prints ~10 pages of console output** (first-stage IV
regression tables, iteration log, coefficient tables, an overidentification
test if `ninst > nu`, a symmetry-given-homogeneity test, elasticities at four
fixed points, descriptive statistics, and a Slutzky-negativity eigenvalue
summary). There is currently no way to get `quaids()`'s results without that
output — splitting estimation from printing is Milestone 1.

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
| `stone` | `0` | Present in the struct; not read by `quaids.src` today — audit before relying on it |
| `aids` | `1` | Present in the struct; not read by `quaids.src` today — audit before relying on it |
| `varname` | `"serie"` | Present in the struct; not read by `quaids.src` today — audit before relying on it |

`stone`, `aids`, and `varname` are set by `quaidsControlCreate()` but not
referenced anywhere in `quaids.src`'s estimation logic — confirmed by grep at
Milestone 0. Left in place deliberately (Milestone 0 scope was dead
*procs*/*structs*, not unused *fields* on an otherwise-live struct); resolve
as part of the Milestone 1 output-schema/control-structure audit.

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
- **Loop style**: `do while ok; ... endo;` / `do while i<=n; ... i=i+1;
  endo;` — not GAUSS's newer `for` loop syntax.
- **Symmetric-restriction idiom**: `design(vec(xpnd(seqa(1,1,k*(k+1)/2))))`
  builds the selection matrix `R` such that `vec(G) = R*vech(G)` for
  symmetric `G` — this is how homogeneity/symmetry constraints get imposed
  via minimum distance. Reuse this idiom rather than re-deriving it.
- **Matrix concatenation**: `~` horizontal, `|` vertical, as usual in GAUSS.
- **Relative vs. absolute prices**: `prices` is converted to relative
  (`prices[.,1:n-1] - prices[.,n]`) near the top of `quaids()` for
  estimation, then converted back to absolute before the elasticities
  section — mutating the input matrix in place. Be careful with this if
  refactoring; it is a real footgun for anyone reading only part of the proc.

## Testing status

There is no automated test suite yet. `examples/quaids_example.e` is a
manual, eyeball-comparison smoke script (verified to still run correctly
after the Milestone 0 rename — see `GOLD_STANDARD_TODO.md` for the run log),
not a CI-style test — see `GOLD_STANDARD_TODO.md` Milestone 3 for the plan to
add deterministic fixtures with real assertions.

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
