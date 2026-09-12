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
complete and committed. **Stage 2** (wiring `sslib`'s
`kalmanFilterTVP()`/`kalmanFilterDiffuseTVP()` into Stage 1's state-vector
construction, with a caller-supplied fixed `Q`/`H`) was just starting when
this status file was created — no Stage 2 code has been written yet.

## Completed Work

- **Public release roadmap, Phases 0–5**: package taken from internal
  pre-alpha to a documented public alpha (`0.1.0` → `0.2.0`). Split
  `optmt`-dependent curvature imposition out of the core installed
  package (opt-in adapter, matching the existing `pubtable` adapter
  pattern); defined support tiers per estimator/feature (README's
  "Model & Feature Support Tiers"); reconciled doc contradictions and
  added automated doc-quality gates (`scripts/verify_docs_consistency.ps1`,
  `verify_docs_quality.ps1`); added a real-data quickstart example, data
  prep and troubleshooting guides; added `CONTRIBUTING.md`/`SUPPORT.md`
  and a single release go/no-go script (`scripts/run_release_gate.ps1`).
- **TVP-AIDS Stage 0** (`src/quaidstrend.src`, `quaidsTrendFit()`): a
  cheap, one-shot LA-AIDS diagnostic for linear coefficient drift, ahead
  of the full Kalman-filter effort. Shipped in `0.2.0`.
- **TVP-AIDS Stage 1** (`src/quaidstvp.src`, private helpers only — no
  public API): the homogeneity+symmetry-respecting state-vector
  construction (`_quaidsTVPGammaIndex`, `_quaidsTVPBuildZ`,
  `_quaidsTVPStoneIndex`, `_quaidsTVPStateToB`), validated by a noiseless
  synthetic-recovery test. Committed as `3e3264b`, pushed to
  `origin/master`. See `CLAUDE.md`'s repository-layout note and
  `dev/GOLD_STANDARD_TODO.md`'s TVP-AIDS section for the full design
  writeup (the Stone-index-vs-translog course correction, the two-valid-
  ways-to-impose-symmetry finding, real bugs found along the way).

## Decisions

- **Stone price index, not the full translog index**, for TVP-AIDS's
  first pass. The translog index makes the AIDS share equation bilinear
  in the state once alpha/gamma/beta are all time-varying, which breaks
  a standard linear-Gaussian Kalman filter. Approved by the repo owner
  after this was found and disclosed mid-Stage-1. A relinearized
  full-AIDS version is legitimate future work, not started.
- **Symmetry imposed as a hard shared-state-element constraint** (gamma_ij
  and gamma_ji are literally the same state, not two coefficients
  reconciled after the fact) — deliberately different from
  `quaidsFit()`'s own GLS-projection approach to the same restriction.
  Confirmed via a noiseless recovery test that both are legitimate and
  not expected to match exactly on real data.
- **`sslib` (gauss-state-space) is the Kalman-filter dependency**, not
  `tsmt`'s own `kalmanFilter()` (documented `k_endog>1` limitation) —
  decided early in the TVP-AIDS research phase.
- Dependency **version pinning for `sslib` was decided in principle**
  ("pin to a specific commit/tag") but **no mechanism has been designed
  or built yet** — GAUSS's `package.json` `deps` array is a bare string
  list with no room for a commit hash. This blocks a clean Stage 6.

## Tests / Validation

- `tests/run_source_tests.ps1` (33 files, no flags skipped) passed clean
  as of the Stage 1 commit — includes the new `tests/quaidstvp_test.e`
  (56 checks: exact noiseless-recovery + loose real-data plausibility
  vs. `quaidsFit()`'s `bestB`).
- A real coupling bug was found and fixed while wiring Stage 1 in:
  `tests/quaidsfixtures.src`'s new TVP fixture originally called a
  `src/quaidstvp.src` proc directly, which broke every
  `guard_error_cases/*.e` script (none of them load `quaidstvp.src`).
  Fixed by inlining the needed logic into the fixture instead.
- Full release gate (`scripts/run_release_gate.ps1`, build+install+
  examples smoke) has **not** been re-run since the Phase 5 release work
  — not required for Stage 1 (no public API surface changed), but worth
  running before any future version bump/release.

## Known Issues

- **`sslib` is not currently installed** at `C:\gauss26\pkgs\sslib` (or
  anywhere else checked) in this environment. Stage 2 cannot proceed
  until it's reinstalled — check whether a private staging copy still
  exists, or reinstall from the `gauss-state-space` repo, before writing
  Stage 2 code.
- **`README.md`'s prose still says "public alpha (package version
  `0.1.0`)"** (line ~16) — the actual current version (`package.json`,
  `CITATION.cff`, `CHANGELOG.md`) is `0.2.0`, bumped for Stage 0's
  `quaidsTrendFit()`. Stale reference, not yet fixed.
- Iterated AIDS and QUAIDS (the estimator's own nonlinear iteration, not
  the TVP work) remain **experimental** per README's support-tier table
  — high documented convergence-failure rates. This is a known,
  documented property of the estimator, not a bug to fix.
- A separate Claude session may be working on `gauss-state-space` in its
  own repo on this same machine — that repo's own history (not tracked
  here) documented a real shared-working-directory collision earlier in
  this project's life. If picking up `sslib` work, check
  `gauss-state-space`'s own status before assuming its working directory
  is exclusively yours.

## Next Steps

1. Resolve the `sslib` install-location issue above before writing any
   Stage 2 code.
2. **Stage 2**: wire `sslib`'s `kalmanFilterTVP()`/`kalmanFilterDiffuseTVP()`
   into Stage 1's `_quaidsTVPBuildZ()` output, with a caller-supplied
   (not yet estimated) fixed `Q`/`H`. Establish `sslib` as a real
   `package.json` dependency and decide/build the commit-pinning
   mechanism (see Decisions above).
3. **Stage 3**: hyperparameter MLE via `sslib`'s `ssFitTVP()`.
4. **Stage 4**: the TVP smoother (`gauss-state-space`'s
   `ssKalmanSmoothTVP()`, already built in that repo per earlier work).
5. **Stage 5**: `quaidsTVPElasFit()`.
6. **Stage 6**: printer/docs/example/packaging, version bump.

## Handoff Notes

- Working tree is clean; `master` is up to date with `origin/master` at
  commit `3e3264b`.
- This session just restructured `CLAUDE.md` (previously ~5100 lines of
  full milestone-by-milestone history) down to a short durable-only
  orientation file, and created this file, per the repo owner's explicit
  request to reduce reliance on chat history for project context.
  `dev/GOLD_STANDARD_TODO.md` (the pre-existing "living roadmap") was
  left untouched — it already holds the detailed historical record
  `CLAUDE.md` used to duplicate.
- Going forward: update this file (not `CLAUDE.md`, not chat history) at
  the end of a meaningful unit of work or before starting a fresh
  session. Only promote something to `CLAUDE.md` if it will still be
  true and relevant many sessions from now.
