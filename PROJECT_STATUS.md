# PROJECT_STATUS.md

Current work state for this repository. Read this (plus `git log`/`git
status`) at the start of a new session instead of relying on prior chat
history. See `CLAUDE.md` for durable project knowledge and
`dev/GOLD_STANDARD_TODO.md` for the full historical decision log.

_Last updated: 2026-09-13_

## Current Objective

Building **TVP-AIDS** (time-varying-parameter AIDS via a Kalman filter),
a repo-owner-requested extension beyond the now-largely-complete public
release roadmap. Staged as Stage 0–6 (see `dev/GOLD_STANDARD_TODO.md`'s
"TVP-AIDS initiative" section for the full plan). Stages 0, 1, and 2 are
now **complete, committed (`086c97a`), and pushed to `origin/master`**.
**Stage 3** (hyperparameter MLE via `sslib`'s `ssFitTVP()`) has not been
started. Nothing functional remains open on Stage 2 itself; only the
open items already listed under Next Steps (sslib install
durability/commit-pinning mechanism) carry forward into general
TVP-AIDS upkeep, not Stage 2 specifically.

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
- **TVP-AIDS Stage 2** (`src/quaidstvpkalman.src`, new private file --
  `_quaidsTVPBuildModel()`, `_quaidsTVPKalmanFit()`, plus an internal
  `_quaidsTVPReplicateConstant()` helper; committed as `086c97a`, pushed):
  wires
  sslib's `kalmanFilterDiffuseTVP()`/`kalmanFilterTVP()` to Stage 1's
  `_quaidsTVPBuildZ()` output, with a caller-supplied fixed `Q`/`H` and a
  standard random-walk state transition (`T=I`, `c=0`, `R=I`). Deliberately
  a SEPARATE file from `src/quaidstvp.src` (Stage 1) -- putting Stage 2
  code directly into quaidstvp.src was tried first and broke Stage 1's own
  sslib-free compilation/test (`quaidstvp_test.e` started failing with
  "Undefined structure 'tvpModel'"), so it was reverted in favor of a
  second file, mirroring the quaidscurvature.src/pubtable_quaids.src
  optional-adapter pattern exactly. Validated against a real `tgauss` run
  (not just compiled): with `Q`/`H` forced to ~0 and the diffuse filter,
  the final-period filtered state recovers Stage 1's own noiseless
  synthetic true state to floating-point precision (~8e-17 max abs diff),
  confirming the filter-wiring is correct independently of Stage 1's own
  OLS-based exact-recovery check. New test `tests/quaidstvp_kalman_test.e`
  (8 checks) plus two new guard-error cases
  (`tvp_bad_Q_shape.e`/`tvp_bad_H_shape.e`); wired into
  `run_source_tests.ps1` behind a new `-SkipTVPKalman` flag (CI passes it
  -- see Decisions). No version bump (no public API surface -- every new
  proc is private, `_`-prefixed).
  - **Environment finding, corrected this session**: despite the prior
    session's own record that `sslib` was "installed and verified" at
    `C:\gauss26\pkgs\sslib`, that directory was actually ABSENT at the
    start of this session -- confirmed directly (`Get-ChildItem` on the
    real path, sandbox-disabled), not assumed from the stale doc.
    Reinstalled via the same documented method (git archive of
    gauss-state-space's pinned commit `9132c35`, NOT its live working
    tree, which was re-confirmed dirty again this session -- same
    contributor's in-progress `ssstructural.src`/`ssmain.src` changes as
    before) + this repo's own `build_lcg.ps1`. Cause of the disappearance
    is unknown (not investigated -- out of scope); treat `sslib`'s
    presence at that path as NOT durable across sessions until a better
    mechanism exists (see Next Steps).
  - A durable, repo-tracked `GAUSS26_CFG` override config now lives at
    `tests/gauss26_cfg_override/gauss.cfg` (a copy of `C:\gauss26\gauss.cfg`
    with `tsmt\lib` added explicitly to `extra_lib_path`, ahead of the `*`
    wildcard) -- `run_source_tests.ps1` points `GAUSS26_CFG` at it only for
    the sslib-dependent child processes (`quaidstvp_kalman_test.e` and its
    two guard cases), leaving every other test's environment untouched.
    Supersedes the ad hoc user-profile-`%TEMP%`-based copy used earlier
    in this session and in the prior session's own verification (which
    would not be visible to the self-hosted CI runner's separate service
    account). **Before committing, this copy was found to carry a live
    `fred_api_key` value from the real `C:\gauss26\gauss.cfg`** (unrelated
    to the override's actual purpose) — blanked before the commit that
    added it; confirmed via `git diff --cached | grep -iE
    "api[_-]?key|secret|password|token|license|hostid"` that nothing else
    in the diff matched. Re-run that same check if this file is ever
    regenerated by copying the real `gauss.cfg` again.
  - Committed as `086c97a` and pushed to `origin/master` — this will have
    triggered the self-hosted push-triggered CI run
    (`.github/workflows/tests.yml`, with `-SkipTVPKalman` per its own
    updated comment); not confirmed green from inside this session (no
    CI-status tool available here) — check it directly
    (`gh run list`/the Actions tab) at the start of the next session if
    not already known to have passed.

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
- **`sslib` stays OUT of `package.json`'s `deps` array, and
  `quaidstvpkalman.src` stays unlisted in its `src` array** — same
  reasoning that already keeps `optmt`/`pubtable` out of `deps` and
  `quaidscurvature.src`/`pubtable_quaids.src` out of `src`: `deps` is read
  as "hard requirement to even install/compile the core package," not "a
  dependency of one of its optional adapters," and listing
  `quaidstvpkalman.src` would make `sslib` exactly that. The `9132c35` pin
  is recorded only in this file and `quaidstvpkalman.src`'s own header
  comment, not in package.json, pending a real pinning mechanism.
- **Stage 2 code split across two files, not one** —
  `src/quaidstvp.src` (Stage 1, no sslib dependency) and the new
  `src/quaidstvpkalman.src` (Stage 2, hard sslib dependency). Tried as one
  file first; broke Stage 1's own sslib-free compilation immediately
  (confirmed via a real failing test run, not predicted), so reverted to
  two files before anything was committed. `quaidstvp.src` itself was
  restored to be byte-identical to its Stage 1 commit (`git diff` confirms
  no changes survived in that file).
- **`sslib`-dependent tests gated behind a new `-SkipTVPKalman` flag**,
  passed by `.github/workflows/tests.yml`'s push-triggered CI run (the
  same treatment as `-SkipCurvature`/`-SkipPubtable`, but for a stronger
  reason: `sslib` isn't a package.json dependency at all, and was found
  genuinely MISSING from this machine's own `C:\gauss26\pkgs` once already
  this initiative -- its presence is not yet a safe assumption for an
  unattended CI run the way optmt/pubtable's is).

## Tests / Validation

- `tests/run_source_tests.ps1` (33 files, no flags skipped) passed clean
  as of the Stage 1 commit — includes `tests/quaidstvp_test.e` (56
  checks: exact noiseless-recovery + loose real-data plausibility vs.
  `quaidsFit()`'s `bestB`).
- `tests/run_source_tests.ps1 -SkipBootstrap` (this machine's routine
  local gate) passed clean this session with the new
  `tests/quaidstvp_kalman_test.e` (8 checks) and its two new guard cases
  included — confirmed BOTH with and without `-SkipTVPKalman` (the latter
  correctly excludes all three new sslib-dependent scripts and leaves
  everything else, including the untouched `quaidstvp_test.e`, passing).
  Not yet re-run with `-SkipBootstrap` absent (full local gate) this
  session — nothing in that group touches TVP-AIDS, low risk, but hasn't
  been re-confirmed since Stage 1.
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

1. **Stage 3**: hyperparameter MLE via `sslib`'s `ssFitTVP()`. Decide/build
   a real commit-pinning mechanism for `sslib` at that point if it becomes
   more pressing (see Decisions — currently just documented, not
   mechanized, and deliberately not a `package.json` `deps`/`src` entry).
   Consider whether `sslib`'s disappearance-and-reinstall this session
   warrants a more durable install step (a script, not a one-off manual
   `git archive`) before relying on it further.
2. **Stage 4**: the TVP smoother (`gauss-state-space`'s
   `ssKalmanSmoothTVP()`, already built in that repo).
3. **Stage 5**: `quaidsTVPElasFit()`.
4. **Stage 6**: printer/docs/example/packaging, version bump.
5. (Housekeeping, not blocking) Confirm the push-triggered CI run for
   `086c97a` passed — not verified from inside the session that pushed it
   (see Handoff Notes).

## Handoff Notes

- Working tree: `master` clean, up to date with `origin/master` at
  `086c97a` (Stage 2, pushed this session — two commits ahead of this
  file's previous note at `e090225`: `cb960f9` sslib-install, `086c97a`
  Stage 2). Nothing uncommitted.
- `086c97a`'s push triggers the self-hosted, push-only CI workflow
  (`.github/workflows/tests.yml`, `-SkipBootstrap -SkipTVPKalman`) — this
  session pushed but had no way to check the run's result afterward (no
  CI-status tool available). **Check it before trusting master is green**
  if that matters for whatever comes next (e.g. before building on top of
  Stage 2, or before a release).
- `sslib` was found MISSING from `C:\gauss26\pkgs\sslib` at the start of
  this session despite the prior session's own record that it was
  installed and verified there — reinstalled the same way (git archive of
  `gauss-state-space`'s pinned `9132c35`, live working tree re-confirmed
  dirty again). Treat its presence there as not durable across sessions;
  re-verify with `Get-ChildItem C:\gauss26\pkgs\sslib` (not just trusting
  this file) before relying on it, exactly as this session had to.
- Full Stage 2 functional validation (diffuse filter exact-recovers
  Stage 1's noiseless true state to ~8e-17) was done via ad hoc scratch
  scripts before being formalized into the committed-to-repo test file —
  the scratch scripts themselves were NOT kept (session scratchpad, not
  part of this repo).
- A live secret (`fred_api_key`) was caught in `tests/gauss26_cfg_override/
  gauss.cfg` before commit and blanked — see the Completed Work entry
  above. Worth remembering if that file is ever regenerated.
- Update this file (not `CLAUDE.md`, not chat history) at the end of a
  meaningful unit of work or before starting a fresh session. Only
  promote something to `CLAUDE.md` if it will still be true and relevant
  many sessions from now.
