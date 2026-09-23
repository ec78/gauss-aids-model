# PROJECT_STATUS.md

Current work state for this repository. Read this (plus `git log`/`git
status`) at the start of a new session instead of relying on prior chat
history. See `CLAUDE.md` for durable project knowledge and
`dev/GOLD_STANDARD_TODO.md` for the full historical decision log.

_Last updated: 2026-09-23_

## Current Objective

Building **TVP-AIDS** (time-varying-parameter AIDS via a Kalman filter),
a repo-owner-requested extension beyond the now-largely-complete public
release roadmap. Staged as Stage 0–6 (see `dev/GOLD_STANDARD_TODO.md`'s
"TVP-AIDS initiative" section for the full plan). **Stages 0–5 complete,
committed, and pushed to `origin/master`** (Stages 0–2 as `086c97a`/
doc-sync `8ed8ea3`; Stage 3 as `a36c7a6`; Stage 4 as `e3ef801`/doc-sync
`3d9f3d7`; Stage 5 as `c58c511`). CI confirmed `success` for every one of
these via `gh run list` (Stage 5's own run id `34841421508`).

**Stage 6 (the final, last stage) is functionally COMPLETE and fully
verified, but NOT yet committed** — per this session's own standing
instruction, confirm with the user before committing even though the
work itself is done. The full release-verification pipeline
(`scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact
-InstallArtifact`) is green end to end: every source test, every guard
case, the build, the install, the installed-package public API test, and
all 15 example smoke tests (including the new TVP-AIDS one) all pass with
zero failures. See Completed Work's Stage 6 entry for the full account,
including several real bugs found and fixed only once real `sslib`
access became available (an array-typed struct-field sentinel, a
`string`-vs-`matrix` struct field type mismatch, a `$|`-vs-`$+`
character-matrix type mismatch found in TWO separate places, a `diag()`
vs `diagrv()` mistake, an MLE-hangs-at-`n1=4` scale limit, and a
`$+`-broadcast printer bug) — none of these were catchable by
`#include`-based testing alone, only by actually building, installing,
and running against a real `sslib` install and the installed package.
`sslib` itself is now repinned to the tagged, stable `v1.0.0` (`ad15626`)
release of `gauss-state-space`, not a commit hash on a moving branch —
see Decisions for why. Stage 3 is still Q-only MLE (H stays caller-fixed)
— see Decisions for why, a repo-owner-approved scope call made explicitly
before Stage 3 was written, not assumed.

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
- **TVP-AIDS Stage 4** (`src/quaidstvpsmooth.src`, new private file --
  `_quaidsTVPSmoothFit()`): the fixed-interval (Rauch-Tung-Striebel)
  smoother, via `sslib`'s `ssKalmanSmoothTVP()`, turning Stage 2's/Stage
  3's filtered (real-time-causal) state path into a full-sample smoothed
  one. Unlike Stage 3, has NO real proc-level dependency on either prior
  stage's own procs -- only on `struct tvpModel`/`struct kalmanResult`
  (sslib) -- so a caller feeds it whichever (tvpm, rslt) pair it already
  has (Stage 2's own tvpm + `_quaidsTVPKalmanFit()` result, or Stage 3's
  `sOut.tvpFinal`/`sOut.kfResults` -- NOT the pre-fit tvpm, which still
  holds the STARTING Q). Deliberately does not accept a raw data matrix
  in place of `rslt` (unlike `ssKalmanSmoothTVP()` itself) -- that branch
  re-filters via the ordinary non-diffuse filter using tvpm.a_0/tvpm.p_0,
  meaningless zeros defaults under this codebase's diffuse-only design.
  - Read `sslib`'s own `ssKalmanSmoothTVP()` source directly from the
    installed copy (`C:\gauss26\pkgs\sslib\src\sstvp.src`) before writing
    any code, confirming its exact `{ a_TS, p_TS } =
    ssKalmanSmoothTVP(tvpm, rsltOrY)` contract and both accepted
    second-argument forms directly from source, not by analogy from
    Stage 2/3's own conventions.
  - **Pre-check before touching the shared install, then a real finding
    confirmed AFTER this stage's own commit/push**: noticed
    `C:\gauss26\pkgs\sslib\src\sskalman.src`/`sstvp.src` had mtimes ~1.5
    hours newer than the rest of that directory at session start.
    Messaged `gauss-state-space-ea` (a concurrent session, found via
    `ListAgents`) to ask before proceeding, per last session's own
    "coordinate rather than guess" precedent, and proceeded on
    independent verification (a real `tgauss` run confirming the
    installed `ssKalmanSmoothTVP()` matches its own doc comment exactly,
    plus `gauss-state-space`'s own working tree clean at `1b82f62`) while
    awaiting a reply -- this stage's own commit (`e3ef801`) went in before
    one arrived. `gauss-state-space-ea` replied shortly after: the newer
    mtimes actually carried OLDER content -- those two files matched
    `7d5ed72` (one commit behind `1b82f62`), missing the correlated-H
    analytic-gradient path (`_cholDerivLinv` and
    `kalmanFilterDiffuseTVPGrad`'s correlated-H fix) entirely, while every
    other installed `src/*.src` file already matched `1b82f62`. `
    ssKalmanSmoothTVP()` itself was unaffected (confirmed independently
    this session, and re-confirmed by their own diff), so Stage 4's own
    work is NOT impacted -- but this would have silently bitten a future
    stage differentiating through a correlated (non-diagonal) H. Likely
    cause per their account: last session's own documented workaround
    (temporarily copying `src`/`test` into the shared install for testing,
    then reverting to `origin/main`) reverted those two files but the
    "refresh the install to the new push" step afterward was simply
    missed, not another session's active work. They copied the current
    `1b82f62` `src/sskalman.src`/`src/sstvp.src` into the installed copy
    and verified a byte-for-byte match (`diff --strip-trailing-cr`) --
    the shared install should now be genuinely consistent with `1b82f62`
    end to end. Re-verify file hashes/mtimes again if anything TVP-related
    looks stale in a future session; don't assume this one incident is
    the last word on drift (see Known Issues/Next Steps on the still-
    unbuilt pinning mechanism).
  - **Real, empirically-confirmed finding**: the RTS "smoothed variance
    <= filtered variance" tightening property can genuinely fail by a
    small amount (~0.01 absolute, one state element, one period) during
    the exact-diffuse initialization burn-in -- `ssKalmanSmoothTVP()` runs
    the ordinary RTS backward recursion throughout (per its own doc
    comment), not a specialized diffuse-smoother algorithm a still-diffuse
    filtered covariance would technically call for. A property of
    `sslib`'s own implementation, not a bug in this wrapper -- confirmed
    by locating the exact violating period (t=3 of a k_states=7/n1=2
    model needing `ceil(7/2)=4` periods to de-diffuse) and by an
    independent exact-match check against sslib's own established
    time-invariant `ssKalmanSmooth()` fed the identical filtered input
    (matches to floating-point precision, including at the violating
    period -- confirming the "violation" is inherent to the shared
    backward-recursion logic, not something specific to the TVP wrapper).
    `tests/quaidstvp_smooth_test.e`'s own tightening-property check is
    scoped to `period >= ceil(k_states/n1)` accordingly, documented
    inline, not silently weakened.
  - New test `tests/quaidstvp_smooth_test.e` (8 checks -- dimensions; the
    exact final-period state/covariance invariant; RTS tightening outside
    the diffuse burn-in; an EXACT internal-consistency check against
    sslib's own `ssKalmanSmooth()`, matching at every period; a loose
    mean-absolute-error plausibility check against the DGP's true state)
    plus two new guard cases (`tvp_smooth_bad_state_rows.e`/
    `tvp_smooth_bad_state_cols.e`), reusing Stage 2's existing
    `-SkipTVPKalman` flag. No version bump (no public API -- the new proc
    is private). `src/quaidstvpsmooth.src` added to
    `verify_package_manifest.ps1`'s `intentionallyUnlisted` allowlist.
  - **Secondary finding, fixed locally only**: GAUSS's `print`, given a
    single bare `string`-typed (type 6) expression with no leading string
    *literal* in the same statement, misformats it as numeric garbage
    (confirmed directly: `print stringVar;` and
    `print ftocv(x,w,d) $+ "suffix";` both reproduce it; a leading literal,
    even `""`, fixes it). This affects the "N CHECKS FAILED" branch of
    the shared PASS/FAIL summary idiom used across most
    `tests/quaids*_test.e` files -- latent because no committed test has
    actually failed in practice. Confirmed NOT a
    `run_source_tests.ps1` false-negative risk: its pass/fail detection
    keys off the ABSENCE of `"ALL \d+ CHECKS PASSED"`, so a garbled count
    still correctly fails the test -- cosmetic only. Fixed in this
    stage's own new test file; NOT swept across the ~30 other existing
    test files (out of scope this session) -- now a durable CLAUDE.md
    gotcha; worth a repo-wide sweep if it ever actually bites during real
    debugging.
  - `run_source_tests.ps1 -SkipBootstrap` (full suite, no TVP flags
    skipped) re-run clean this session with Stage 4's new test and guard
    cases included: `run_source_tests.ps1: PASS`.
  - Committed as `e3ef801`; doc-sync as `3d9f3d7`.
- **TVP-AIDS Stage 5** (`_quaidsTVPStateToFullB()` added to
  `src/quaidstvp.src`; new private file `src/quaidstvpelas.src` --
  `_quaidsTVPElasFit()`): elasticities at one period's (filtered OR
  smoothed) TVP-AIDS state, reusing `src/quaidselas.src`'s existing
  `_quaidsElas()` per the plan doc rather than hand-rolling new
  elasticity math. NO `sslib` dependency at all (unlike Stages 2-4) --
  `_quaidsTVPElasFit()` takes a state as a plain `k_states x 1` vector,
  not any `sslib` struct, so "filtered vs. smoothed" is entirely the
  caller's choice of which column of Stage 2/3's `rslt.filtered_state` or
  Stage 4's `a_TS` to pass in. Both new/changed procs stay
  `_`-prefixed/private, matching Stages 1-4's own disposition (Stage 6 is
  what actually publishes public API/docs/printer/example, per the plan
  doc).
  - `_quaidsTVPStateToFullB()` closes the two gaps
    `_quaidsTVPStateToB()`'s own header already flagged: (1) equation n's
    own adding-up-implied coefficients (alpha_n = 1 - sum(others), beta_n
    = -sum(others)) and (2) relative-to-absolute-price gamma conversion,
    via homogeneity's own row-sum-zero identity applied per row (INCLUDING
    row n itself, recovered via symmetry from column n) -- NOT via a
    separately-imposed adding-up identity on gamma, which falls out
    automatically once symmetry+homogeneity both hold (confirmed by a
    direct regression-guard check on the recovered output, not just
    trusted from the derivation). Reused this same derivation conceptually
    from `quaids.src`'s own "RECOVERS ABSOLUTE PRICE EFFECTS FROM
    RELATIVE" block (read directly, per this session's own design
    question) but written as fresh, much simpler direct formulas rather
    than replicating that block's general-case reshape/kron machinery
    (which also handles quaidsFit()'s own separate `ng`/nonlinear-`u`-block
    dimensions this reduced state doesn't have) -- `quaids.src` was
    consulted, not edited, per this session's own explicit instruction.
  - `_quaidsTVPElasFit()` hardcodes `intcpt = 1` (Stage 1's own nint=0
    scope decision -- no extra intercept shifters) and requires
    `aCtl.linear == 1` (guarded explicitly; the recovered `b` has no
    lambda row for `_quaidsElas()` to read if `aCtl.linear` is left at
    `quaidsControlCreate()`'s own default 0) -- a caller with no other
    reason to build a `quaidsControl` just calls `quaidsControlCreate()`
    then sets `aCtl.linear = 1`. Point elasticities only, deliberately --
    no delta-method SEs (unlike `quaidsElasFit()`), matching the plan
    doc's own choice of `_quaidsElas()` (not `quaidsElasFit()`) as the
    sibling to reuse; propagating the reduced state's own covariance
    through the (linear) recovery step is real, tractable future work, not
    attempted this stage.
  - New test `tests/quaidstvp_elas_test.e` (13 checks, NO `library`
    statement / no `-SkipTVPKalman` gating needed -- genuinely no `sslib`
    dependency): builds a full n x n TRUE absolute-price gamma directly
    (symmetric, exact-zero-row-sum by double-centering a random symmetric
    matrix -- an INDEPENDENT construction, not `_quaidsTVPStateToFullB()`'s
    own recovery formula) plus true full alpha/beta (equation n's values
    set by hand via the adding-up identities, not by calling library code),
    derives the n1-equation reduced system as a trivial submatrix/
    subvector extraction, generates noiseless data, recovers `stateHat` via
    plain pooled OLS (same pattern as Stage 1's own test), and checks the
    recovered `bFull` matches the independently-built true full system to
    floating-point precision -- a real correctness claim about the
    recovery logic, not a tautology. Plus a regression guard confirming
    homogeneity/symmetry/adding-up all hold exactly on the recovered
    output, and an exact internal-consistency check
    (`_quaidsTVPElasFit()`'s output vs. a direct `_quaidsElas()` call on
    the same recovered `bFull`) satisfying CLAUDE.md's "two independent
    checks" testing requirement. Three new guard cases
    (`tvp_elas_bad_linear.e`/`tvp_elas_bad_state_length.e`/
    `tvp_elas_bad_prices_length.e`), also NOT gated behind
    `-SkipTVPKalman`. `src/quaidstvpelas.src` added to
    `verify_package_manifest.ps1`'s `intentionallyUnlisted` allowlist.
  - **Real language-gotcha finding, now in CLAUDE.md**: `mSym` (any
    case) as an assignment TARGET is a reserved identifier in GAUSS 26 --
    `mSym = (a + a')/2;` fails with a generic `error G0008 : Syntax error
    '= (a + a')/2'` (pointing at the RHS, not the reserved name) plus a
    cascading spurious error on the next line. Found while writing this
    stage's own test fixture (`mSym` was the first natural name for a
    symmetric-matrix intermediate); bisected by varying only the
    assignment-target name in an isolated `tgauss` repro before assuming
    the expression syntax itself was at fault. Renamed to `symMat` in the
    committed test.
  - `run_source_tests.ps1 -SkipBootstrap` (no TVP flags skipped, the full
    local gate) re-run clean this session with Stage 5's new test and
    three new guard cases included.
  - Committed as `c58c511` and pushed to `origin/master` (user explicitly
    asked for the commit+push). CI confirmed `success` via `gh run list`
    (run id `34841421508`, ~2m27s).
- **TVP-AIDS Stage 6 (functionally complete, NOT yet committed)**: the
  final stage — publishes real public API, docs, example, packaging, the
  `sslib` pinning mechanism, and the version bump. A repo-owner design
  plan was written and approved (via a plan-mode review) before any code
  was written, including two explicit decisions: the public surface stays
  MINIMAL (one consolidated `quaidsTVPFit()`, not every stage exposed
  separately) and the `sslib` pin gets a real automated mechanism, not
  just better documentation.
  - Two Stage 1/5 procs promoted from private to public IN PLACE (same
    behavior, dropped leading underscore): `_quaidsTVPStateToFullB` →
    `quaidsTVPStateToFullB` (`src/quaidstvp.src`), `_quaidsTVPElasFit` →
    `quaidsTVPElasFit` (`src/quaidstvpelas.src`). Verified directly: both
    renamed sslib-FREE tests (`tests/quaidstvp_test.e`,
    `tests/quaidstvp_elas_test.e`, plus its three renamed guard cases)
    re-run clean under the new names (56 + 13 checks). A PowerShell
    `Set-Content -Encoding utf8` rewrite of `quaidstvp_elas_test.e` during
    the rename was found to have silently added a UTF-8 BOM, breaking
    `tgauss`'s own lexer (`error G0008 : Syntax error '﻿ new'`) —
    fixed via `[System.IO.File]::WriteAllText` with an explicit
    BOM-less `UTF8Encoding` instead; worth remembering for any future
    PowerShell-driven rewrite of a `.e`/`.src` file (`Write`/`Edit`
    themselves do not add a BOM, only `Set-Content -Encoding utf8` does,
    confirmed directly).
  - New file `src/quaidstvpfit.src`: `quaidsTVPFit()` (pure orchestration
    over Stages 1/3/4/5's already-tested pieces — Stone index + Z-build,
    MLE-fitted Q against a caller-fixed H via Stage 2's filter internally,
    the RTS smoother if `tvpCtl.smooth`, `quaidsTVPStateToFullB()` at the
    final period — no new estimation math), `printQuaidsTVP()`,
    `quaidsTVPControlCreate()`. New structs `quaidsTVPControl`/
    `quaidsTVPOut` in `src/quaids.sdf` (matrix/array/string fields only,
    same optmt/sslib-independent pattern `quaidsCurvOut` already uses —
    core `quaids.sdf` must keep compiling without `sslib`).
    - `quaidsTVPOut.smoothedStateCov` is `array`-typed (matches
      `filtered_state_cov`'s own type from `sslib`'s `kalmanResult`);
      the "not smoothed" sentinel uses a placeholder `arrayinit(1|1|1, 0)`,
      not a bare `0` (which would hit the same "Illegal assignment - type
      mismatch" `quaidstvpkalman.src`'s own `_quaidsTVPReplicateConstant`
      header already documents) — confirmed correct by a real passing
      check in `tests/quaidstvpfit_test.e` once `sslib` access was
      available, not left as a guess.
    - `quaidsTVPControl.othnam` is `matrix`-typed, NOT `string`-typed
      like `quaidsControl.othnam` — a real bug found via
      `examples/14_tvp_aids_estimation.e`, the FIRST place in this
      codebase's history either `othnam` field was ever assigned a real
      (non-default) value (`error G0071 : Type mismatch` on a `string`
      field). A SECOND, separate gotcha found the same way: a value
      built with `$|` (vertical string concat, e.g.
      `quaidsExampleGoodNames()`'s own construction) is rejected even by
      a `matrix`-typed field, or by `print$`/`~` horzcat, with a
      DIFFERENT error (`G0165`, same root cause) — only the legacy
      `$+`/`ftocv()` character-matrix form works; fix is a leading
      `0$+` coercion. Both are new durable `CLAUDE.md` gotchas now,
      confirmed via isolated scratch repro scripts before editing any
      real file, not guessed from the failure site alone.
    - `printQuaidsTVP()`'s coefficient/Q-diagonal-block printing
      originally used `"label: " $+ "" $+ ftocv(rowVector,w,d)` — found,
      via real output inspection (not just "did it crash"), to BROADCAST
      the label across every element instead of joining one string
      (`"alpha block: v1  alpha block: v2"`). Switching `$+` to `~` avoids
      the broadcast but then truncates any label over 8 characters (the
      already-documented legacy-character-matrix cell limit). Fixed by
      printing the label and its row of values as two SEPARATE
      `print`/`print$` calls — matching every existing multi-row printer
      in this codebase, confirmed via an isolated repro before editing
      the real printer. Also a new durable `CLAUDE.md` gotcha.
  - `sslib.pin.json` (new, repo root): now pinned to `gauss-state-space`'s
    tagged, STABLE `v1.0.0` release (`ad15626`), not a commit hash on a
    moving branch — a real, live course correction mid-session (see
    Decisions for the full story: the shared install drifted TWICE more
    even after this mechanism was built and first used, because
    `gauss-state-space-1a`'s own test runner mirrors ITS working tree into
    the shared install on every one of its own runs, making any commit
    pin inherently unstable while that session iterates on a breaking
    `v2` line). `verifiedFiles`/`lastVerified` are populated from a real
    `scripts/sync_sslib.ps1` run against this tag, not placeholders.
  - `scripts/sync_sslib.ps1` (new): `git archive` of the pinned/target
    commit from a LOCAL `gauss-state-space` checkout (never the live
    install, never that repo's own working tree) into a clean temp
    staging dir, `robocopy /MIR` into the shared install (confirmed this
    correctly PURGES any stale `lib/` directory too, since a git archive
    has no `lib/` at all — closes exactly the "source and .lcg from
    different generations" failure class `gauss-state-space-1a` flagged
    as a risk, confirmed by re-reading the script together with them, not
    just asserted), rebuild `lib/sslib.lcg` via this repo's own
    `build_lcg.ps1` (a full `Set-Content` overwrite, not an incremental
    update — also closes that same failure class, independent of the
    `/MIR` purge), recompute/rewrite `sslib.pin.json`'s `verifiedFiles`/
    `lastVerified`. Supports both resync-to-documented-pin (default) and
    deliberate repin (`-Commit`). **This mechanism already proved its own
    value live**: `scripts/verify_sslib_pin.ps1` caught real drift twice
    during this same session (see Decisions), each time BEFORE it could
    silently break a test run the way the pre-Stage-6 ad hoc process
    already had twice earlier in this initiative.
  - `scripts/verify_sslib_pin.ps1` (new): sha256-compares the installed
    copy against `sslib.pin.json`'s `verifiedFiles`; warns (default) or
    hard-fails (`-Strict`) on drift. Wired into
    `run_source_tests.ps1`'s `-SkipTVPKalman`-false branch. **Real
    PowerShell 5.1 gotcha found and fixed**: `@($obj.PSObject.
    Properties.Name)` on a TRULY EMPTY `PSCustomObject` (e.g. `{}` from
    `ConvertFrom-Json`) yields a ONE-element array containing `$null`/`""`,
    not an empty array — confirmed directly (a `Get-FileHash` call on a
    null path threw "cannot call a method on a null-valued expression"
    instead of the intended "no tracked files yet" branch being taken).
    Fixed in both this script and `sync_sslib.ps1` by checking
    `.PSObject.Properties.Count -eq 0` first, only enumerating `.Name` via
    `ForEach-Object` once confirmed non-empty.
  - `docs/public-api.json`: new `optional_modules` entry `"tvp"` (source
    is an ARRAY of six files, unlike `curvature`/`pubtable_adapter`'s one
    file each — genuinely needed, not a style choice, since the TVP-AIDS
    public surface spans that many files). Required a real code change to
    `scripts/verify_public_api.ps1` (previously assumed `module.source`
    was always a single string) — extended to accept a string OR an
    array, keyed the internal lookup by module NAME instead of `source`
    (which may no longer be a valid scalar hashtable key).
  - Version bump `0.2.0` → `0.3.0` (`package.json`, `CITATION.cff`,
    `docs/public-api.json`) + `CHANGELOG.md` entry. **Found and fixed a
    pre-existing stale reference while here**: `README.md`'s own "public
    alpha (package version `0.1.0`)" line (flagged stale — should have
    said `0.2.0` — in this file's own Known Issues since the Stage 5
    session) now correctly says `0.3.0`.
  - New `examples/14_tvp_aids_estimation.e` + `examples/README.md`/
    `README.md` updates (new "Time-Varying-Parameter Estimation (optional,
    `sslib`)" sections, examples-count bump 13→14, a new support-tier
    table row). `tests/run_examples_smoke.ps1` gained a matching
    `-SkipTVPKalman` switch and its own `GAUSS26_CFG` wiring (mirroring
    `run_source_tests.ps1`'s own pattern) for this one example.
    - **Two real bugs found only by actually running this example
      end to end**, neither catchable by `#include`-based testing alone
      since both only manifest against the installed package / real
      `sslib`: (1) `H = diag(vector)` does NOT build a diagonal matrix
      from a vector in GAUSS (that's `diagrv(eye(n), vector)`) —
      `diag()` on a vector instead degenerates to something far smaller,
      which `_quaidsTVPMLEFit`'s own shape guard correctly caught
      (`error: H must be n1 x n1`) rather than silently misbehaving;
      fixed in the example AND in the two doc pages (`quaidsTVPFit.md`,
      `README.md`) that had copied the same mistake. (2) The full 5-good
      `quaidsExampleData()` dataset (`n1=4`, `k_states=18`) makes
      `ssFitTVP()`'s CMLMT optimization pathologically slow/
      non-terminating — confirmed directly by letting a run burn ~24
      hours of real CPU time before killing it, then isolating the cause
      with short-timeout scratch scripts (still hung at `n1=4`/`tobs=200`
      within 100s; converged in ~34s at `n1=2`/`tobs=500`). Fixed by
      restricting the example to a 3-good subset (`n1=2`), matching every
      other test/example in this initiative, and documented as a real,
      confirmed scale limit (`quaidsTVPFit.md`'s Remarks, the README
      support-tier table) — not merely an untested configuration.
  - New `tests/quaidstvpfit_test.e` (16 checks, all INTERNAL CONSISTENCY —
    `quaidsTVPFit()`'s own output vs. the same already-validated Stage
    1/3/4/5 procs called directly on identical inputs — no synthetic-
    recovery-of-truth claim needed since Stage 6 adds no new math,
    including a dedicated check that a caller-supplied `tvpCtl.othnam` is
    honored, added after the type-mismatch bug above was found) + two new
    guard cases (`tvpfit_bad_H_shape.e`/`tvpfit_bad_q0_shape.e`, reusing
    `_quaidsTVPMLEFit`'s own existing error messages, no duplicate
    validation added in the public wrapper). Uses a small self-contained
    synthetic generator, NOT `tests/quaidsfixtures.src`'s existing
    `_quaidsTVPStaticSyntheticDGP`/`_quaidsTVPDynamicSyntheticDGP` — both
    of those already return pre-built relative-price/Stone-deflated
    pieces, not the raw absolute-price/totexp/full-`w` inputs
    `quaidsTVPFit()`'s own public contract actually takes. All 16 checks
    confirmed passing via a real `tgauss` run.
  - `verify_package_manifest.ps1`'s `intentionallyUnlisted` allowlist
    gained `quaidstvpfit.src` (sixth TVP-AIDS entry).
  - `CLAUDE.md` updated: the `sslib`/optional-modules notes now describe
    the new pin file/scripts as the durable mechanism (superseding the
    old "check PROJECT_STATUS.md for which commit" pointer), the
    `command-reference/*.md` page count (49 → 54), plus three new durable
    language gotchas (the `$|`-vs-`$+` character-matrix incompatibility,
    the `$+`-broadcast-across-a-multi-element-`ftocv()`-result printer
    trap, and this section's own cross-reference).
  - **Full release-verification pipeline
    (`scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact
    -InstallArtifact`) confirmed GREEN end to end**: source tests, every
    guard case, build, install, the installed-package public API test,
    and all 15 example smoke tests (including the new TVP-AIDS one) —
    zero `FAIL`/`error G\d+`/"Undefined structure" lines anywhere in the
    full output, confirmed by grepping the complete output, not just
    skimming it. This is the FIRST time in this initiative any
    `sslib`-dependent code has been validated against a real install
    AND a real package rebuild/reinstall, not just `#include`-based
    source-tree tests — and it directly surfaced 6 real, independently
    confirmed bugs (listed above) that `#include`-based testing alone had
    completely missed. **NOT committed** — working tree has all of the
    above uncommitted, per this session's own standing instruction to
    confirm before any commit even mid-stage.

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
- **`sslib` is pinned to `gauss-state-space`'s tagged `v1.0.0` release
  (`ad15626`), not a bare commit hash** — a real course correction made
  DURING Stage 6, after the newly-built pinning mechanism (`sslib.pin.json`
  + `scripts/sync_sslib.ps1`/`verify_sslib_pin.ps1`) caught the shared
  install drifting TWICE more in the same session, even after already
  being resynced once. Root cause, confirmed directly by the concurrent
  `gauss-state-space-1a` session, not guessed: their own test runner
  (`test/run-tests.ps1`) mirrors THEIR working tree into
  `C:\gauss26\pkgs\sslib` via `robocopy /MIR` on every single run of their
  own gate, so the shared install tracks whatever they're actively
  iterating on, not any commit either side pins — a genuinely different
  problem from the earlier "someone forgot to refresh the install" drift
  incidents (Stages 3-4). They are also mid-flight on a deliberately
  BREAKING `v2.0.0` line (`feature/v2-modern-api-components`, minimum
  GAUSS 26.1.4, a reshaped `ssOut` with a new nested `structural` field,
  typed-struct-return/keyword-argument public procs) — a real
  `error G0507 : Undefined structure 'ssStructuralInfo'` hit mid-session
  was diagnosed together with them as a stale-`.lcg`-vs-new-source
  generation mismatch from their OWN install-refresh timing, not a bug on
  this repo's side. `v1.0.0` is their explicit recommendation for a
  GAUSS-26.1.4-compatible, old-`ssOut`-shape target that they will not
  touch further; `scripts/sync_sslib.ps1`'s hand-built `.lcg` generator
  (a full `Set-Content` scan, not GAUSS's own `lib` command) was
  confirmed BY THEM to be strictly more robust than `lib` for `v1.0.0`
  specifically (v1's own `sstvp.src` references a global declared in
  `ssmain.src`, which trips `lib`'s own per-file-independent parser into
  silently truncating that file's index entry) — but their `.lcg` would
  need the `typed_returns`/`keywords` annotations only GAUSS's real `lib`
  command emits if this repo ever moves up to `v2` (declaration-free
  struct assignment and keyword calls are load-bearing on those
  annotations there); install `v2` via the GAUSS Package Manager or `lib`
  itself when that day comes, not `sync_sslib.ps1`'s current generator.
  They will ping this repo's session if/when `v2` reaches their `main`.
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
- `tests/run_source_tests.ps1 -SkipBootstrap`, WITHOUT `-SkipTVPKalman`
  (i.e. the full local gate, no TVP flags skipped at all) re-run clean
  this session (Stage 5) with the new `tests/quaidstvp_elas_test.e` (13
  checks) and its three new guard cases included. This run also
  re-verified Stages 2-4's own `sslib`-dependent tests still pass against
  the currently-installed `sslib` copy -- worth noting since
  `sstvp.src`/`sskalman.src` had very recent mtimes at this session's
  start (consistent with `gauss-state-space-ea`'s concurrent session being
  active), so this full-gate pass is direct re-verification, not just
  trust in a timestamp. Also confirmed standalone
  (`tgauss -b -x quaidstvp_elas_test.e`, 13/13 checks) and all three new
  guard cases fail with their expected diagnostics run individually.

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

1. **`sslib` commit-pinning mechanism: BUILT and PROVEN this session**
   (`sslib.pin.json` + `scripts/sync_sslib.ps1`/`verify_sslib_pin.ps1` —
   see Stage 6's own Completed Work entry and Decisions). Pinned to
   `gauss-state-space`'s tagged `v1.0.0` (`ad15626`), confirmed installed
   and verified clean via a real `scripts\verify_sslib_pin.ps1` run. The
   mechanism already caught real drift twice live during this session,
   working exactly as intended.
2. **Stage 5 is done** -- committed and pushed as `c58c511`, CI green.
3. **Stage 6 is functionally DONE, not yet committed.** Every piece
   (public `quaidsTVPFit()`/`printQuaidsTVP()`/`quaidsTVPControlCreate()`/
   `quaidsTVPStateToFullB()`/`quaidsTVPElasFit()`, docs, example, the
   pinning mechanism, version bump to `0.3.0`) is written AND the full
   `scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact
   -InstallArtifact` pipeline (source tests, guard cases, build, install,
   installed-package API test, all 15 example smoke tests) is confirmed
   green end to end — the real `sslib`-backed run this file's earlier
   session left as a resume step is now done, and it surfaced/fixed 6 real
   bugs in the process (see Completed Work). **Only remaining step: ask
   the user whether to commit** — per this session's own standing
   instruction to confirm before any commit even mid-stage. Suggested
   commit scope: everything currently uncommitted (see Handoff Notes)
   as one Stage 6 commit, matching how Stages 0-2/3/4/5 were each
   committed as a single unit.

## Handoff Notes

- **Working tree: `master` has substantial UNCOMMITTED, but fully
  VERIFIED, Stage 6 work** on top of `origin/master`'s `c58c511`
  (Stage 5) — do NOT assume clean/up to date. Touches:
  `src/quaidstvp.src` (rename), `src/quaidstvpelas.src` (rename), new
  `src/quaidstvpfit.src`, `src/quaids.sdf` (two new structs),
  `sslib.pin.json` (new), `scripts/sync_sslib.ps1`/`verify_sslib_pin.ps1`
  (new), `scripts/verify_public_api.ps1` (multi-file module support),
  `tests/run_source_tests.ps1`/`run_examples_smoke.ps1` (wiring),
  `tests/quaidstvpfit_test.e` (new) + two new guard cases,
  `tests/quaidstvp_elas_test.e` + its three guard cases (rename),
  `docs/public-api.json`/`docs/COMMAND_REFERENCE.md` + five new
  `docs/command-reference/*.md` pages, `README.md`/`examples/README.md`,
  new `examples/14_tvp_aids_estimation.e`,
  `package.json`/`CITATION.cff`/`CHANGELOG.md` (version bump to `0.3.0`),
  `CLAUDE.md` (sslib/optional-modules notes + three new language
  gotchas). **Not committed** — per this session's own standing
  instruction, confirm with the user before any commit even mid-stage.
  The full release-verification pipeline has been run clean against this
  exact working-tree state (see Next Steps item 3) — nothing here is
  provisional or awaiting a resume step.
- `gh run list` confirmed (as of 2026-09-16): CI runs for `086c97a`,
  `8ed8ea3`, `a36c7a6`, `e3ef801`/`3d9f3d7` (Stage 4 + doc-sync), AND
  `c58c511` (Stage 5) all completed with `success`. Every push so far
  (including Stage 5's) runs with `-SkipBootstrap -SkipTVPKalman`, so CI
  itself never actually exercises any `sslib`-dependent test — a
  reminder that CI green here means "the non-sslib suite passed," not
  "the sslib path was validated." Stages 2-4's own `sslib`-dependent
  tests (`quaidstvp_kalman_test.e`, `quaidstvp_mle_test.e`,
  `quaidstvp_smooth_test.e`, all their guard cases) are validated only by
  local `run_source_tests.ps1 -SkipBootstrap` runs (no TVP flags
  skipped), most recently Stage 5's own session (2026-09-13/14), not by
  CI. Stage 5's own new test/guard cases (`quaidstvp_elas_test.e` and its
  three guards) have NO `sslib` dependency and so ARE exercised by CI
  normally, unlike Stages 2-4's.
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
- **This session's `sslib` saga, summarized** (full account in Decisions/
  Completed Work): resynced to `2a26076` after a live coordination check;
  ran the full test suite and examples smoke suite; found a real
  `diag()`-vs-`diagrv()` bug and an `n1=4` MLE-hang scale limit in the new
  example; while fixing those, `gauss-state-space-1a`'s own test runner
  refreshed the shared install AGAIN (mirroring their own working tree,
  not any commit) and left it briefly inconsistent (`error G0507`,
  diagnosed together with them as a stale-`.lcg` timing issue, not a bug
  here); repinned to their recommended stable `v1.0.0` tag (`ad15626`)
  instead of chasing their actively-moving `v2` branch; final pipeline
  run after that repin is clean. `C:\gauss26\pkgs\sslib` was explicitly
  freed back to them at the end of this session (messaged directly) —
  check `ListAgents`/message them before assuming it's still at `v1.0.0`
  in a future session, same collision-risk caution as always.
- Update this file (not `CLAUDE.md`, not chat history) at the end of a
  meaningful unit of work or before starting a fresh session. Only
  promote something to `CLAUDE.md` if it will still be true and relevant
  many sessions from now.
