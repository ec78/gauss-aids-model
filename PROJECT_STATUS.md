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
"TVP-AIDS initiative" section for the full plan). Stages 0–3 are now
**complete, committed, and pushed to `origin/master`** (Stages 0–2 as
`086c97a`/doc-sync `8ed8ea3`; **Stage 3 as `a36c7a6`**, CI confirmed
`success` via `gh run list`). **Stage 4** (the TVP smoother,
`ssKalmanSmoothTVP()`) has not been started. Stage 3 is still Q-only MLE
(H stays caller-fixed) — see Decisions for why, a repo-owner-approved
scope call made explicitly before Stage 3 was written, not assumed.

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
- **TVP-AIDS Stage 3** (`src/quaidstvpmle.src`, new private file --
  `_quaidsTVPQUpdate()`, `_quaidsTVPMLEFit()`): hyperparameter MLE for
  the diagonal state innovation covariance Q, via `sslib`'s `ssFitTVP()`.
  Scope decided explicitly with repo-owner sign-off before writing code
  (see Decisions): Q ONLY (H stays caller-fixed, same as Stage 2) and Q
  DIAGONAL (not a full covariance) -- both to sidestep a documented
  state-space variance-identification risk, not arbitrary
  simplifications. Positivity enforced via `sslib`'s own
  `ssControl.positive_vars` squaring transform (`sstransformParms`), not
  a hand-rolled log-variance transform. Depends directly on Stage 2's own
  `_quaidsTVPBuildModel()`/`_quaidsTVPReplicateConstant()` (a real
  proc-level dependency, unlike Stage 2's own deliberately Stage-1
  decoupled design) -- a caller must `#include quaidstvpkalman.src`
  before this file.
  - **Real, non-obvious API gotcha found and documented** (now also in
    CLAUDE.md's language-gotchas list): `ssFitTVP()`'s own `y` parameter
    is **nobs x k_endog** (this library's usual Txn convention), the
    OPPOSITE of `_quaidsTVPKalmanFit()`'s `k_endog x nobs` -- `ssFitTVP()`
    transposes internally. Passing the transposed form compiles fine and
    fails deep inside `kalmanFilterDiffuseTVP` with a generic "Matrix
    dimensions are incompatible", not a clear argument-shape error at the
    call boundary. Confirmed directly against a real failure before
    fixing, not assumed from either proc's doc comment (which are
    themselves easy to misread by analogy from Stage 2's own contract).
  - New fixture `_quaidsTVPDynamicSyntheticDGP()` added to
    `tests/quaidsfixtures.src` -- unlike Stage 1's
    `_quaidsTVPStaticSyntheticDGP()` (noiseless, time-invariant), this
    simulates a genuine random-walk state with KNOWN diagonal Q and KNOWN
    H, i.e. the actual correctly-specified DGP Stage 3's MLE assumes, so
    recovering the true Q is a meaningful correctness check.
  - **Finding, empirically confirmed (not guessed)**: individual
    per-state-element Q diagonal entries are only loosely identified even
    at `tobs=500` (worst single element off by ~54% of its own true value
    at the test's seed), while the AGGREGATE (mean across all `k_states`
    elements) is much better identified (~7% off, same seed) -- a milder,
    within-Q version of the same general state-space
    variance-identification phenomenon that motivated fixing H in the
    first place (see Decisions). `tests/quaidstvp_mle_test.e`'s own
    synthetic-recovery check is therefore on the MEAN of the fitted Q
    diagonal, not each element individually -- a deliberate, documented
    choice, not a weakened test.
  - New test `tests/quaidstvp_mle_test.e` (7 checks: convergence,
    positivity/no-blowup sanity, mean-Q recovery, an EXACT
    internal-consistency check against Stage 2's independently-validated
    `_quaidsTVPKalmanFit()` at the fitted Q -- confirmed to match to
    floating-point precision -- and a loose final-state plausibility
    check) plus four new guard-error cases
    (`tvp_mle_bad_q0_shape.e`/`tvp_mle_nonpositive_q0.e`/
    `tvp_mle_bad_H_shape.e`/`tvp_mle_bad_y_shape.e`); wired into
    `run_source_tests.ps1` under the SAME `-SkipTVPKalman` flag as Stage
    2 (not a new flag -- identical underlying reason to skip). No version
    bump (no public API surface -- every new proc is private,
    `_`-prefixed).
  - `src/quaidstvpmle.src` added to `verify_package_manifest.ps1`'s
    `intentionallyUnlisted` allowlist, same reasoning as
    `quaidstvpkalman.src`.
  - **Mid-session incident, found and fixed before Stage 3 work could
    proceed**: `quaidstvp_kalman_test.e` (Stage 2, previously passing)
    started failing with `error G0159 : Wrong number of parameters
    'init_diffTVP' expected 2 arguments, received 1` at the very start of
    this session's work, despite NO local change to this repo. Root
    cause: `C:\gauss26\pkgs\sslib\src\sstvp.src`/`ssstructural.src` (the
    SHARED installed package directory) had been modified that same
    morning (confirmed via `Get-ChildItem` timestamps, not assumed) by
    something other than this session. `ListAgents` found a concurrent
    session (`gauss-state-space-ea`) actively working in the
    `gauss-state-space` repo; messaged them directly rather than touching
    shared state unilaterally. They confirmed: their own repo working
    tree was clean/unrelated (not mid-edit on the installed copy
    themselves), and the signature change (`init_diffTVP`/
    `init_stationaryTVP` gaining a required `stationary_states` second
    parameter, unused for the pure-diffuse case, added so every init
    function shares one dispatch signature for a new "mixed" init mode)
    is from upstream commit `ae921ce`, intentional and finished, not WIP
    -- so someone/something else had refreshed the shared install to a
    newer commit than this repo's own documented `9132c35` pin, between
    last session's end and this session's start. Fixed by updating
    `_quaidsTVPKalmanFit()`'s call site (`src/quaidstvpkalman.src`) to
    `init_diffTVP(tvpm, 0)` (the second arg is genuinely unused for the
    diffuse case, confirmed by reading `sstvp.src`'s own comment) --
    re-ran `quaidstvp_kalman_test.e` standalone afterward to confirm all
    8 checks pass again before continuing to Stage 3. See Known Issues
    for the now-stale `9132c35` pin this incident exposed, and Decisions
    for CLAUDE.md's new durable gotcha about shared-package-directory
    collision risk.
  - Committed as `a36c7a6` and pushed to `origin/master` (user explicitly
    asked for the commit+push). The self-hosted push-triggered CI run
    completed with `success` (confirmed via `gh run list` from inside
    this session, unlike the prior two commits where that had to wait
    for a follow-up session) -- run id `34784706305`.

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
  is a bare string list with no room for a commit hash. **The
  previously-documented `9132c35` pin is now KNOWN STALE**, not just
  theoretically at risk: this session found the actual installed copy at
  `C:\gauss26\pkgs\sslib` includes at least upstream commit `ae921ce`
  (confirmed by the `init_diffTVP` arity break/fix — see Completed Work's
  Stage 3 entry), and `gauss-state-space-ea`'s own concurrent session
  reported their repo (clean, matching origin) is at `7d5ed72`, LIKELY
  (not independently confirmed) close to what's actually installed. Blocks
  a clean Stage 6 even more concretely now than before.
- **`sslib` stays OUT of `package.json`'s `deps` array, and
  `quaidstvpkalman.src`/`quaidstvpmle.src` stay unlisted in its `src`
  array** — same reasoning that already keeps `optmt`/`pubtable` out of
  `deps` and `quaidscurvature.src`/`pubtable_quaids.src` out of `src`:
  `deps` is read as "hard requirement to even install/compile the core
  package," not "a dependency of one of its optional adapters," and
  listing either file would make `sslib` exactly that. The stale
  `9132c35` pin is recorded only in this file and
  `quaidstvpkalman.src`'s own header comment, not in package.json,
  pending a real pinning mechanism (see above — now a real, not
  hypothetical, gap).
- **TVP-AIDS Stage 3 estimates Q ONLY via MLE; H stays caller-supplied
  and FIXED** (repo-owner sign-off, given explicitly this session before
  any Stage 3 code was written) — `sslib`'s own `test/sstvpfit.inc`
  header documents that jointly estimating both Q and H via unconstrained
  MLE hit the classic state-space variance-identification problem and
  never converged even given hundreds of iterations on a univariate
  local-level model. A future stage could revisit joint estimation (e.g.
  profiling H from a static `quaidsFit()` residual covariance as a
  smarter starting point) if it becomes a real need — not attempted.
- **Stage 3's Q is DIAGONAL, not a full covariance** (same sign-off) —
  one free variance per state element via `sslib`'s own
  `positive_vars`-squaring transform, not a Cholesky-parameterized full
  covariance. Standard TVP-VAR/TVP-AIDS simplifying assumption
  (independent per-state-element random-walk innovations), and the far
  cheaper/safer choice given the identification risk above (a full
  covariance would be `k_states*(k_states+1)/2` free parameters).
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
  unattended CI run the way optmt/pubtable's is). Stage 3's own
  sslib-dependent test/guard cases reuse this SAME flag rather than a new
  one -- identical underlying reason to skip.

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
- `tests/run_source_tests.ps1 -SkipBootstrap` (this machine's routine
  local gate) re-run clean this session with Stage 3's new
  `tests/quaidstvp_mle_test.e` (7 checks) and its four new guard cases
  included, AND with the `init_diffTVP` arity fix applied to
  `quaidstvpkalman.src` — full suite (all `guard_error_cases`, all
  `gaussTests` including `quaidstvp_kalman_test.e` and
  `quaidstvp_mle_test.e`) reported `run_source_tests.ps1: PASS`. Not yet
  re-run WITHOUT `-SkipBootstrap` (full local gate) this session.

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
- **The collision risk is not limited to that repo's own working tree —
  the SHARED INSTALLED COPY at `C:\gauss26\pkgs\sslib` itself changed
  mid-session** (this session, not a hypothetical): `sstvp.src`/
  `ssstructural.src` there were modified the same morning by something
  other than this session, silently breaking already-committed Stage 2
  code (`init_diffTVP` arity). Now fixed (see Completed Work), and a new
  durable CLAUDE.md gotcha records the general pattern, but the
  `9132c35` pin is now confirmed stale (see Decisions) and there is still
  no mechanism preventing this from recurring. `ListAgents` found the
  likely source (`gauss-state-space-ea`, a concurrent session); after
  investigating on their end, they clarified: the drift past `9132c35`
  predates their own session (a PREVIOUS session pushed `ae921ce`/
  `bba9467` to `gauss-state-space`'s `origin/main`, now at `7d5ed72`, and
  something refreshed the shared install from that before this session
  started); they themselves temporarily copied `src`/`test` into the
  shared install today for their own testing but reverted it back to
  match `origin/main` exactly (verified via diff) before this incident
  was even raised. So the actual actor that FIRST refreshed the shared
  install past `9132c35` remains unidentified, but the install's current
  content is understood: consistent with `gauss-state-space`
  `origin/main` around `7d5ed72` as of 2026-09-13 — this should become
  the new reference point once a real pinning mechanism exists (see
  Decisions), not `9132c35`. They also flagged pushing a further new
  commit (`1b82f62`, additive-only -- analytic-gradient support for a
  correlated/non-diagonal H in the exact-diffuse filter, previously an
  error case) after this incident, NOT yet confirmed present in the
  installed copy -- low risk (additive, no signature changes, their own
  43-file suite passed) but worth knowing if something TVP-related
  behaves unexpectedly in a future session.
- **`README.md`'s prose still says "public alpha (package version
  `0.1.0`)"** (line ~16) — actual current version is `0.2.0`. Stale
  reference, not yet fixed.
- Iterated AIDS and QUAIDS (the estimator's own nonlinear iteration, not
  the TVP work) remain **experimental** per README's support-tier table
  — a known, documented estimator property, not a bug.

## Next Steps

1. Decide/build a real commit-pinning mechanism for `sslib` — a
   confirmed-real gap (the `9132c35` pin is stale; see Decisions/Known
   Issues), not just a theoretical one. Consider whether the pin should
   live somewhere more durable than this file + a header comment, given
   it has now silently drifted at least once without anyone noticing
   until a test broke. Current best reference point if this is tackled:
   `gauss-state-space` `origin/main` was at `7d5ed72` as of 2026-09-13,
   plus a further additive commit `1b82f62` not yet confirmed installed
   (see Known Issues).
2. **Stage 4**: the TVP smoother (`gauss-state-space`'s
   `ssKalmanSmoothTVP()`, already built in that repo).
3. **Stage 5**: `quaidsTVPElasFit()`.
4. **Stage 6**: printer/docs/example/packaging, version bump.

## Handoff Notes

- Working tree: `master` clean, up to date with `origin/master` at
  `a36c7a6` (Stage 3, committed and pushed this session at the user's
  explicit request — three commits ahead of this file's previous note at
  `086c97a`: `8ed8ea3` doc-sync, then `a36c7a6` Stage 3 itself, bundling
  `src/quaidstvpmle.src`, the `init_diffTVP` arity fix to
  `src/quaidstvpkalman.src`, `tests/quaidstvp_mle_test.e`, four new guard
  cases, the `_quaidsTVPDynamicSyntheticDGP()` fixture,
  `run_source_tests.ps1`/`verify_package_manifest.ps1` wiring, and two
  new CLAUDE.md gotchas). Nothing uncommitted.
- `gh run list` confirmed this session: CI runs for `086c97a`, `8ed8ea3`,
  AND `a36c7a6` all completed with `success`. `086c97a`/`8ed8ea3` (like
  every push so far) ran with `-SkipBootstrap -SkipTVPKalman`, so they
  never actually exercised any sslib-dependent test — a reminder that CI
  green here means "the non-sslib suite passed," not "the sslib path was
  validated," given `-SkipTVPKalman` is always passed by CI. The Stage
  3/sslib path (`quaidstvp_kalman_test.e`, `quaidstvp_mle_test.e`, all
  six `tvp_*`/`tvp_mle_*` guard cases) is validated only by this
  session's own local `run_source_tests.ps1 -SkipBootstrap` run (passed
  clean, no flags related to TVP skipped), not by CI.
- `sslib`'s presence at `C:\gauss26\pkgs\sslib` was reconfirmed at this
  session's start (real files, `lib/sslib.lcg` catalog present) — NOT
  missing this time, unlike last session's finding. But its CONTENT had
  silently drifted past the documented `9132c35` pin (see Known Issues'
  new entry) — presence alone is no longer sufficient reassurance;
  arity/signature drift is now a demonstrated real risk on top of the
  already-documented disappearance risk. Re-verify both before relying on
  it in a future session.
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
