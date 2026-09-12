# CLAUDE.md — GAUSS QUAIDS Library

Durable orientation for Claude Code sessions working on this repository.
This file should stay short and stable. For **current work state**
(what's in progress, recent decisions, open issues), read
`PROJECT_STATUS.md` instead — do not look for it here.

**Context efficiency:** Read or inspect only the files needed for the current task unless broader repository review is necessary. Avoid dumping large file contents, logs, or test output into the conversation when a concise summary or targeted excerpt is sufficient.

## What this project is

A GAUSS application package (`library quaids;`) estimating Almost Ideal
Demand System models (LA-AIDS, iterated AIDS, QUAIDS) with IV-treated
total expenditure, homogeneity/symmetry restrictions, elasticities,
welfare measures, robust/cluster/bootstrap/replicate-weight inference,
and an applied workflow layer. Full feature list, quick start, and model
support tiers: **`README.md`** (the canonical user-facing description —
don't duplicate it here). Methodology: `docs/METHODOLOGY_NOTES.md`.

Package name `quaids` (not `aids`) — decided early to avoid a
too-generic bare identifier; "AIDS"/"Almost Ideal Demand System" is still
the correct term for the model family in prose.

## Repository layout

```
src/                    # Installed package source (see package.json's
                         #   src array for load order). One .src file per
                         #   estimator/feature area:
  quaids.sdf             #   struct definitions (quaidsControl, quaidsOut, ...)
  quaidsutil.src          quaidsControlCreate()
  quaidsiv.src            IV first-stage helper (private)
  quaidszerocorrect.src   quaidsZeroFit() -- Shonkwiler-Yen zero-share correction
  quaidselas.src          quaidsElasFit()/quaidsElas()/printQuaidsElas()
  quaidsshares.src        quaidsSharesFit()/printQuaidsShares()
  quaidsslutzky.src       quaidsSlutzky()
  quaids.src              quaidsFit()/printQuaids()/quaids() -- estimation core
  quaidsformula.src       quaidsFull() -- dataframe entry point
  quaidstests.src         quaidsHomogeneityTest/quaidsJointTest/quaidsQuadraticTest
  quaidswelfare.src       quaidsWelfareFit() -- CV/EV
  quaidsrobust.src        quaidsRobustFit()/quaidsRobustBootstrapFit()
  quaidsdiagnostics.src   quaidsPreflight()
  quaidsworkflow.src      quaidsWorkflowFit()/quaidsWorkflowScenarioFit()
  quaidssurvey.src        quaidsSurveyWorkflowFit()
  quaidsreplicate.src     quaidsReplicateWeightFit()
  quaidstrend.src         quaidsTrendFit() -- TVP screening diagnostic
  quaidstvp.src           TVP-AIDS state-vector plumbing (private, WIP -- see PROJECT_STATUS.md)
  quaidscurvature.src     quaidsCurvatureFit() -- OPTIONAL, requires optmt,
                         #   deliberately NOT in package.json's src array
                         #   (see "Optional modules" below)
  pubtable_quaids.src     OPTIONAL pubtable export adapter, also not in src array
examples/                # 14-file numbered, read-in-order example suite
                         #   (00_real_data_quickstart.e .. 13_pubtable_reporting.e)
                         #   + example_data.src (shared synthetic-data generator).
                         #   See examples/README.md.
tests/
  quaids*_test.e          # One assertion-based test file per feature area.
  quaidsfixtures.src      #   shared synthetic-data generators (private, not
                         #   part of the installed package)
  guard_error_cases/      #   standalone scripts confirming specific bad-input
                         #   guards fail loudly and clearly
  fixtures/published/     #   real published data (Blanciforti86) + R/Python
                         #   cross-validation reference scripts
  run_source_tests.ps1    #   runs the manifest/API/docs-consistency checks
                         #   plus every tgauss test file
  run_examples_smoke.ps1  #   runs every examples/*.e file for real
  verify_package_manifest.ps1
docs/
  COMMAND_REFERENCE.md, USAGE_GUIDE.md, METHODOLOGY_NOTES.md,
  DATA_PREPARATION_GUIDE.md, TROUBLESHOOTING_GUIDE.md,
  FEATURE_SUPPORT_MATRIX.md, public-api.json (machine-checked API contract)
  command-reference/*.md  # one page per public proc (49 pages)
scripts/
  build_lcg.ps1, build_package.ps1, verify_release_artifact.ps1,
  run_release_verification.ps1, run_release_gate.ps1 (single go/no-go
  release command), verify_public_api.ps1, verify_docs_consistency.ps1,
  verify_docs_quality.ps1
.github/workflows/tests.yml   # self-hosted CI, push-to-master only (see below)
dev/
  GOLD_STANDARD_TODO.md   # Living roadmap + full milestone-by-milestone
                         #   engineering history/decision log. Consult this
                         #   for "why was X built this way" on anything
                         #   already shipped -- do NOT duplicate its content
                         #   here or in PROJECT_STATUS.md.
  PUBLIC_RELEASE_ROADMAP.md  # Phased plan for the public release (mostly complete)
package.json             # GAUSS package manifest (src array = installed
                         #   package's load order; deps array = required
                         #   external GAUSS packages)
CHANGELOG.md, CITATION.cff, README.md, CONTRIBUTING.md, SUPPORT.md
PROJECT_STATUS.md        # Current work state -- read this, not chat history
```

### Optional modules (not in `package.json`'s `src` array)

`quaidscurvature.src` (needs `optmt`) and `pubtable_quaids.src` (needs
`pubtable`) are real, documented, tested public API, but are excluded
from the installed package's lazy-load catalog so that neither `optmt`
nor `pubtable` becomes a hard dependency for core estimation. A caller
loads them explicitly:
```gauss
library optmt, quaids;
#include quaidscurvature.src
```
Both files are listed in `tests/verify_package_manifest.ps1`'s
`intentionallyUnlisted` allowlist. **Adding any other new required
`.src` file**: add it to `package.json`'s `src` array (respecting load
order — a file that calls procs in another file must load after it),
bump the version, and rebuild/reinstall.

## Development environment

- **GAUSS 26** at `C:\gauss26` (`tgauss.exe` at `C:\gauss26\tgauss.exe`).
  Installed packages relevant here: `optmt`, `pubtable`, and this
  package itself (`quaids`, at `C:\gauss26\pkgs\quaids`, rebuilt via
  the release scripts below — not by hand-editing that directory).
- **`sslib` (gauss-state-space)** — a separate Aptech-licensed
  state-space/Kalman-filter package this project's TVP-AIDS work builds
  on, installed at `C:\gauss26\pkgs\sslib`. Sourced from a specific
  commit of the separate `gauss-state-space` repo (not built by this
  repo's own scripts) — check `PROJECT_STATUS.md` for which commit and
  whether it's still present before assuming it's available or current.
- **`tsmt` package shadowing on this machine**: `gauss.cfg`'s
  `extra_lib_path` resolves `$(PACKAGEDIR)\*\lib` alphabetically, and
  `pkgs\timeseries\lib\tsmt.lcg` (an older/incomplete catalog also named
  `tsmt.lcg`) shadows the real `pkgs\tsmt\lib\tsmt.lcg` — `library tsmt;`
  (directly, or transitively via `library sslib;`) then fails to resolve
  real TSMT procs (`cusum`, `constrain_stationary`, etc.) with plain
  "Undefined symbol" errors that give no hint shadowing is the cause.
  Fix per-invocation without touching the shared `gauss.cfg`: set env
  var `GAUSS26_CFG` to a directory holding a copy of `gauss.cfg` whose
  `extra_lib_path` lists `$(PACKAGEDIR)\tsmt\lib` explicitly before the
  `*` wildcard.
- **R 4.5.0** (`C:\Program Files\R\R-4.5.0\bin\Rscript.exe`, package
  `micEconAids`) and **Python 3.12** (numpy/pandas/scipy) are installed
  only to regenerate the published-data cross-validation reference
  numbers in `tests/fixtures/published/` — neither is a runtime
  dependency of anything in `src/`.
- This repo is **public** on GitHub. CI runs on a **self-hosted**
  runner (GAUSS is licensed, not available on GitHub-hosted runners),
  triggered on `push` to `master` only — never `pull_request`, since a
  self-hosted runner on a public repo is a real fork/PR code-execution
  risk under that trigger.

## Key commands

```powershell
# Fast source-tree test suite (manifest/API/docs checks + all tgauss tests)
powershell -ExecutionPolicy Bypass -File tests\run_source_tests.ps1
# Add -SkipBootstrap / -SkipCurvature / -SkipPubtable to skip slower or
# optional-package-gated groups.

# Run every examples/*.e file for real
powershell -ExecutionPolicy Bypass -File tests\run_examples_smoke.ps1

# Rebuild + reinstall the package after any src/ change, then verify
# the installed public API (library quaids; against a real install)
powershell -ExecutionPolicy Bypass -File scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact

# Single release go/no-go gate (full suite incl. bootstrap, build,
# install, installed-API test, all example smoke tests, convergence
# sweep, checksum/record)
powershell -ExecutionPolicy Bypass -File scripts\run_release_gate.ps1
```

A single test file directly: `tgauss -b -x <file>.e` from `tests/` (or
`examples/`) as the working directory.

## Coding conventions

- **Variable naming**: short lowercase names (`w`, `n`, `gg`, `b`, `u`),
  matching the original author's terse econometrics-code style.
- All locals for a proc declared in one `local` statement at the top.
- Struct-returning procs are declared `proc (struct TypeName) = name(...);`
  with no variable name in the return slot, so callers can skip
  pre-declaring the target (`qOut = quaidsFit(...);` works with no prior
  `struct quaidsOut qOut;`) — both via `#include` and via `library
  quaids;`. Do this for every new struct-returning proc. Caveats: a
  variable's inferred type can't be retyped by a later call returning a
  *different* struct type in the same scope (`error G0504`); inference
  does not propagate through a plain struct-to-struct copy (`error
  G0008`); inside a proc body `struct T var;` is never optional
  regardless of inference (it doubles as the local-variable declaration).
- **Symmetric-restriction idiom**: `design(vec(xpnd(seqa(1,1,k*(k+1)/2))))`
  builds the selection matrix `R` such that `vec(G) = R*vech(G)` for
  symmetric `G` — reuse this for any new homogeneity/symmetry-style
  minimum-distance restriction rather than re-deriving it.
- **Relative vs. absolute prices**: `quaidsFit()` converts `prices` to
  relative form internally, then back to absolute before its final
  recovery step. GAUSS passes matrices by value, so this never leaks
  into a caller's own `prices` local.
- **Don't touch already-shipped, tested estimation core without a
  strong reason.** The default pattern for new functionality is a new
  sibling `.src` file/proc (e.g. `quaidsSharesFit()`, `quaidsZeroFit()`),
  not editing `quaids.src` itself — `quaidsFit()`'s iteration/variance/
  restriction machinery shares heavily mutated intermediate state across
  phases, and splitting or editing it is a real risk best taken only
  when unavoidable (e.g. adding sampling weights had to touch it, since
  weighting threads through the whole estimation core).
- **Version bump policy**: bump `package.json`'s version on any real new
  or changed *public* API surface (new/changed proc signature, new
  struct field, a deliberate behavior-changing bugfix to shipped output)
  — not for pure tooling, docs, or example changes. Every change gets a
  `CHANGELOG.md` entry regardless.
- Loop style is `do while ...; ... endo;`, not GAUSS's `for` loop.

### Known GAUSS-26 language gotchas (verified in this codebase)

- `{a, b, c}` matrix literal: commas separate **rows**, not columns.
- `gamma`, `quantile` and other builtin/function names are reserved —
  can't be reused as local variable names.
- `sign()` is not a GAUSS builtin — use `.>` comparisons.
- Legacy `$+` character-matrix concatenation truncates each cell to 8
  characters — use `printQuaidsElas.src`'s pattern (value and `(SE)` on
  separate printed rows) instead.
- GAUSS identifiers are **case-insensitive** — `K` and `k` collide.
- `reshape()` fills **row-major**; `reshape(v, rows(X), cols(X))` does
  **not** invert `vec(X)` in general — use
  `reshape(v, cols(X), rows(X))'` instead.
- `trap 1` does **not** catch every failure mode — a call-arity mismatch
  (e.g. `eighv()` returning fewer values than expected on some degenerate
  inputs) or certain internal engine errors (`glm()`, deep indexing
  errors) abort the whole job even inside a `trap 1` block. Pre-check the
  specific known trigger condition before the call instead of relying on
  `trap`.
- A struct-returning proc call **can** be safely wrapped in
  `trap 1,1; ...; trap oldtrp,1; if scalmiss(x); ...`, but there is no
  way to `scalmiss()` a whole struct directly — check specific fields
  the caller reads afterward.
- `library`-based lazy loading does not activate a `#define` (e.g.
  `quaids.sdf`'s `#ifndef QUAIDS_SDF_INCLUDED` guard) for a file outside
  the package unless that file `#include`s the `.sdf` explicitly —
  struct availability cannot be assumed from load order alone under
  `library`.
- A literal `"` (odd count) or a bare `*` right after `/` inside `/* */`
  comment *text* can break the comment lexer (`error G0097`/`G0562`) —
  even though it's inside a comment, not code.
- Keyword-defaulted parameters and `...` (dynargs) cannot coexist in one
  proc; required parameters must precede any defaulted ones. For a
  `library`-loaded (not `#include`d) proc, `build_lcg.ps1`'s catalog
  needs the `proc (...) = name(a, b=default);` declaration on **one
  line** — a wrapped multi-line signature can silently lose its
  `: keywords` catalog tag and then reject a keyword-omitted call that
  works fine under `#include`.
- GAUSS's `run "file.e";` does not return control to the calling script
  — do not chain multiple `run` statements expecting sequential execution.
- `library`-based lazy loading resolves a plain global-variable
  declaration (e.g. `struct ssControl _ssActiveCtl;` or a bare top-level
  assignment) in one file only if *some other symbol from that same
  file* has already been referenced — a hand-built `.lcg` catalog (e.g.
  `build_lcg.ps1`) only catalogs `proc`/`struct`-type-definition entries,
  not plain globals, so a proc in file B that reads a global declared in
  file A fails with "Undefined symbol" if file A was never independently
  triggered to load first. Reference a real proc from the defining file
  before the dependent one, or `#include` both directly, instead of
  relying on `library` lazy-loading order. Confirmed in `sslib`
  (`sstvp.src`'s TVP filters depend on globals declared in `ssmain.src`).

## Testing expectations

- Every `tests/quaids*_test.e` file prints one `PASS`/`FAIL` line per
  check and a final `ALL N CHECKS PASSED` (or `N CHECKS FAILED`) summary
  line. **Check that line** — `tgauss`'s process exit code is not a
  reliable pass/fail signal for this harness.
- `tests/run_source_tests.ps1` is the standard local gate; it also runs
  `verify_package_manifest.ps1`, `verify_public_api.ps1`,
  `verify_docs_consistency.ps1`, and `verify_docs_quality.ps1` first.
- New estimation logic needs two independent checks, not one: (1) a
  synthetic fixture with a known true answer (exact recovery where the
  math allows it, e.g. noiseless data; documented tolerance otherwise),
  and (2) either a cross-implementation check (R/Python on real
  published data) or an internal-consistency check against an
  already-correct sibling proc. Never trust a single check, and when a
  derived formula/reference disagrees with the code, verify by hand
  before assuming the code is wrong.
- A new test-only fixture belongs in `tests/quaidsfixtures.src`, not
  duplicated per test file, unless it's genuinely one-off.
- `tests/guard_error_cases/*.e` confirm one specific bad-input guard
  fails with a specific, clear diagnostic — add one here for any new
  required-input validation.
- `tests/package_public_api.e` / `package_public_api_core_only.e` are
  release gates, not part of the routine suite — they run `library
  quaids;` against a real *installed* copy, so they only make sense
  after a rebuild/install.

## Project-wide constraints

- **Never commit or push without being explicitly asked.**
- Don't mutate the shared, installed GAUSS package directory
  (`C:\gauss26\pkgs\quaids`) casually — only via the release scripts,
  and only when a rebuild is actually needed for what you're testing.
- If another Claude session might be working in a different GAUSS
  package's repo on this same machine at the same time, be aware shared
  machine state (installed packages, `gauss.cfg`) is a real collision
  risk — see `PROJECT_STATUS.md` if this is currently relevant.
- Keep `CLAUDE.md` and `PROJECT_STATUS.md` in sync with reality as you
  work: update `PROJECT_STATUS.md` at meaningful checkpoints (a
  milestone/stage completed, before ending a session); only add to this
  file when something is genuinely durable across many future sessions
  (a permanent convention, a permanent environment fact, a recurring
  language gotcha) — not project history, which belongs in
  `dev/GOLD_STANDARD_TODO.md`.
