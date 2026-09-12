# PROJECT_STATUS.md

Current work state for this repository. Read this (plus `git log`/`git
status`) at the start of a new session instead of relying on prior chat
history. See `CLAUDE.md` for durable project knowledge and
`dev/GOLD_STANDARD_TODO.md` for the full historical decision log.

_Last updated: 2026-09-12_

## Current Objective

Building **TVP-AIDS** (time-varying-parameter AIDS via a Kalman filter),
a repo-owner-requested extension beyond the now-largely-complete public
release roadmap. Staged as Stage 0–6 (see `dev/GOLD_STANDARD_TODO.md`'s
"TVP-AIDS initiative" section for the full plan). Stages 0 and 1 are
complete and committed. **Stage 2** (wire `sslib`'s `kalmanFilterTVP()`/
`kalmanFilterDiffuseTVP()` into Stage 1's state-vector construction, with
a caller-supplied fixed `Q`/`H`) has **not been started** — no Stage 2
code exists yet, but its blocking dependency question (was `sslib`
installed and usable?) is now resolved, so Stage 2 can begin.

## Completed Work

- **Public release roadmap, Phases 0–5**: package taken from internal
  pre-alpha to a documented public alpha (`0.1.0` → `0.2.0`). Split
  `optmt`-dependent curvature imposition out of the core installed
  package (opt-in adapter, matching the existing `pubtable` adapter
  pattern); defined support tiers per estimator/feature (README's
  "Model & Feature Support Tiers"); reconciled doc contradictions and
  added automated doc-quality gates; added a real-data quickstart
  example, data prep and troubleshooting guides; added
  `CONTRIBUTING.md`/`SUPPORT.md` and `scripts/run_release_gate.ps1`.
- **TVP-AIDS Stage 0** (`src/quaidstrend.src`, `quaidsTrendFit()`): a
  cheap, one-shot LA-AIDS diagnostic for linear coefficient drift, ahead
  of the full Kalman-filter effort. Shipped in `0.2.0`.
- **TVP-AIDS Stage 1** (`src/quaidstvp.src`, private helpers only — no
  public API): homogeneity+symmetry-respecting state-vector construction
  (`_quaidsTVPGammaIndex`, `_quaidsTVPBuildZ`, `_quaidsTVPStoneIndex`,
  `_quaidsTVPStateToB`), validated by a noiseless synthetic-recovery
  test. Committed as `3e3264b`. See `dev/GOLD_STANDARD_TODO.md`'s
  TVP-AIDS section for the full design writeup.
- **`sslib` installed and verified** at `C:\gauss26\pkgs\sslib`: a real
  copy (not a junction) built via `git archive` of `gauss-state-space`'s
  committed `main` HEAD (`9132c35`) — deliberately not that repo's live
  working tree, which has unrelated uncommitted work in progress
  (`src/ssstructural.src`). Catalog built with this repo's own
  `scripts/build_lcg.ps1` (confirmed genuinely package-agnostic).
  Confirmed working with a real `tgauss` run (see Known Issues for the
  two environment gotchas this required working around, now also in
  `CLAUDE.md`).
- **Context-management restructuring**: `CLAUDE.md` cut from ~5,100 lines
  down to a durable-only orientation file; this `PROJECT_STATUS.md`
  created for current-state tracking. Committed as `e090225`.

## Decisions

- **Stone price index, not the full translog index**, for TVP-AIDS's
  first pass — the translog index makes the AIDS share equation bilinear
  in the state once alpha/gamma/beta are all time-varying, breaking a
  standard linear-Gaussian Kalman filter. Approved by the repo owner. A
  relinearized full-AIDS version is legitimate future work, not started.
- **Symmetry imposed as a hard shared-state-element constraint** (gamma_ij
  and gamma_ji are literally the same state) — deliberately different
  from `quaidsFit()`'s own GLS-projection approach to the same
  restriction. Both confirmed legitimate via a noiseless recovery test;
  not expected to match exactly on real data.
- **`sslib` (gauss-state-space) is the Kalman-filter dependency**, not
  `tsmt`'s own `kalmanFilter()` (documented `k_endog>1` limitation).
- **`sslib` version pinning decided in principle** ("pin to a specific
  commit") but **no mechanism built yet** — `package.json`'s `deps` array
  is a bare string list with no room for a commit hash. Currently
  installed from `9132c35`; this should be the pin once a mechanism
  exists. Blocks a clean Stage 6.

## Tests / Validation

- `tests/run_source_tests.ps1` (33 files, no flags skipped) passed clean
  as of the Stage 1 commit — includes `tests/quaidstvp_test.e` (56
  checks: exact noiseless-recovery + loose real-data plausibility vs.
  `quaidsFit()`'s `bestB`).
- `sslib` install verified directly with `tgauss` (not just file
  presence): `library cmlmt, tsmt, sslib;` then referencing
  `ssControlCreate()` before `kalmanFilterTVP`/`kalmanFilterDiffuseTVP`
  resolves and compiles cleanly under the `GAUSS26_CFG` override
  described in Known Issues.
- Full release gate (`scripts/run_release_gate.ps1`) has **not** been
  re-run since the Phase 5 release work — not required for Stage 1/the
  sslib install (no public API surface changed), but worth running
  before any future version bump/release.

## Known Issues

- **`tsmt` package shadowing + the `library`/cross-file-global lazy-load
  gotcha** — both now documented as durable environment/language facts
  in `CLAUDE.md` (Development environment / Known GAUSS-26 language
  gotchas) rather than here, since they'll recur for any future session
  touching `sslib`. In short: any `tgauss` invocation using `tsmt` or
  `sslib` needs the `GAUSS26_CFG` override described there, and Stage 2
  code should call `ssControlCreate()` before any `sstvp.src` proc.
- `gauss-state-space`'s documented collision risk is real and live, not
  hypothetical: its `main` working tree currently has an uncommitted,
  in-progress change to `src/ssstructural.src` (someone else's
  analytic-gradient work, unrelated to anything Stage 2 needs). Re-check
  its `git status` before reading from it again, and never install
  `sslib` from its live working tree — use `git archive` of a specific
  commit, as this session did.
- **`README.md`'s prose still says "public alpha (package version
  `0.1.0`)"** (line ~16) — actual current version is `0.2.0`. Stale
  reference, not yet fixed.
- Iterated AIDS and QUAIDS (the estimator's own nonlinear iteration, not
  the TVP work) remain **experimental** per README's support-tier table
  — a known, documented estimator property, not a bug.

## Next Steps

1. **Stage 2**: wire `sslib`'s `kalmanFilterTVP()`/`kalmanFilterDiffuseTVP()`
   into Stage 1's `_quaidsTVPBuildZ()` output, with a caller-supplied
   (not yet estimated) fixed `Q`/`H`. Establish `sslib` as a real
   `package.json` dependency and decide/build the commit-pinning
   mechanism (see Decisions) — pin to `9132c35` unless a newer commit is
   deliberately chosen.
2. **Stage 3**: hyperparameter MLE via `sslib`'s `ssFitTVP()`.
3. **Stage 4**: the TVP smoother (`gauss-state-space`'s
   `ssKalmanSmoothTVP()`, already built in that repo).
4. **Stage 5**: `quaidsTVPElasFit()`.
5. **Stage 6**: printer/docs/example/packaging, version bump.

## Handoff Notes

- Working tree: `master` up to date with `origin/master` at `e090225`;
  this file has uncommitted edits from this session (not committed —
  not asked to).
- No Stage 2 code written yet — this session was scoped to resolving the
  `sslib` dependency question only, which is now done.
- Update this file (not `CLAUDE.md`, not chat history) at the end of a
  meaningful unit of work or before starting a fresh session. Only
  promote something to `CLAUDE.md` if it will still be true and relevant
  many sessions from now.
