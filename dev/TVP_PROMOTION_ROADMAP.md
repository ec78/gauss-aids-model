# TVP-AIDS Promotion Roadmap

_Created: 2026-10-15. Status: **Priority 1 (active)**. P1 (D0 = option A,
v2.1.0), P2, P3, P4 and P7 implemented on 2026-10-15 in package 0.4.0, and
P5 on 2026-10-16 (see section 10). P5 found an upstream blocker: a
non-diagonal H gave a wrong diffuse likelihood in sslib v2.1.0. On
2026-10-17 (Stage 1) sslib was repinned to v2.1.1, which fixes it, and P5
case 3 became a counted full-sample check, so **P5 is closed** (section
10). P6 (bounded `n1` scale re-measurement) was done on 2026-10-17
(Stage 2, section 10). P8 remains open, so TVP-AIDS stays Experimental._

This is the durable plan for promoting TVP-AIDS (`quaidsTVPFit()` and the
Stage 1-6 machinery in `src/quaidstvp*.src`) from **Experimental** to a
supported, release-gated capability. It preserves the 2026-10-15 read-only
diagnosis (task `task_03d43536665b`) and sets the order of work, including
the decision gate that must pass before `sslib` is repinned.

Other work is **deferred and preserved, not abandoned** (see section 7).

## 1. Problem

`tests/quaidstvpfit_test.e` is the only automated check of the public
`quaidsTVPFit()`/`printQuaidsTVP()` entry point. It has been
**waived/unverified** since 2026-10-10:

- It did not finish in two earlier runs. One sat 59 minutes in its first
  sslib/CMLmt fit.
- CI and the routine local gate always pass `-SkipTVPKalman`, so no
  automated gate has ever run it.

Promotion is blocked until that test is well-posed, bounded and passing
against a deliberately chosen `sslib` baseline.

## 2. Evidence (2026-10-15 diagnosis)

All runs were read-only and hard-capped. Scripts lived in a session
scratch directory. No repo file, shared install or GAUSS configuration was
changed. Environment: GAUSS 26.1.4.5565, with `GAUSS26_CFG` pointed at
`tests/gauss26_cfg_override`.

### 2.1 Bounded test matrix (shared install as found)

| Test | Wall time | Result |
|---|---|---|
| `quaidstvp_test.e` (Stage 1, no sslib) | 1 s | ALL 56 CHECKS PASSED |
| `quaidstvp_elas_test.e` (Stage 5, no sslib) | 1 s | ALL 13 CHECKS PASSED |
| `quaidstvp_kalman_test.e` (Stage 2) | 1 s | ALL 8 CHECKS PASSED |
| `quaidstvp_mle_test.e` (Stage 3, DGP data, tobs=500) | 20 s | ALL 7 CHECKS PASSED |
| `quaidstvp_smooth_test.e` (Stage 4) | 1 s | ALL 8 CHECKS PASSED |
| `guard_error_cases/tvpfit_bad_H_shape.e` | 2 s | intended message: `_quaidsTVPMLEFit: H must be n1 x n1.` |
| `guard_error_cases/tvpfit_bad_q0_shape.e` | 1 s | intended message: `_quaidsTVPMLEFit: q0 must be k_states x 1 ...` |
| `quaidstvpfit_test.e` (Stage 6) | **600 s** | **TIMEOUT, killed.** Only check 1 (control default) had printed; still inside the first MLE. |

### 2.2 Root cause: an unbounded optimizer on an ill-posed fixture

**Unbounded optimizer (repository defect).**
- `_quaidsTVPMLEFit()` (`src/quaidstvpmle.src`) calls `ssFitTVP()` with
  CMLmt's defaults: `maxIters = 1e5` and `maxTime = 1e8` minutes
  (`C:\gauss26\pkgs\cmlmt\src\cmlmtutil.src`). In practice these never stop it.
- `quaidsTVPControl` exposes neither setting, so a caller cannot bound
  the fit.

**Measured cost.**
- About 35 ms per CMLmt iteration on the test fixture.
- A 3000-iteration capped run took about 103 s and returned retcode 2
  (iteration limit), with no convergence.
- 1e5 iterations comes to about 58 minutes for one CMLmt run. This matches
  the recorded 59-minute stall.

**Multipliers.**
- `ssFitTVP()` then runs a second, untransformed CMLmt pass from `p_final`
  (present in both sslib v1.0.0 and v2.0.0).
- The test calls the MLE four times.
- The worst case is therefore several hours.

**Ill-posed fixture (repository defect, test side).**
- The shares are random draws, not an AIDS DGP.
- `H = 0.01*eye(2)` is fixed at about one quarter of the sample share
  variance (0.0385 and 0.0341).
- `q0 = 0.001` is the flat start that `docs/command-reference/quaidsTVPFit.md`
  already warns against.
- The MLE pushes 3-4 of the 7 Q variances to the zero boundary
  (1e-5 to 1e-10), and the ML covariance fails to invert.

**Second-pass instability.** The untransformed pass warns
`CMLmt: function evaluated as complex` because variances near 0 step
negative once the positivity transform is off.

### 2.3 SSLib drift: real, but not the cause

**State of the shared install.**
- `C:\gauss26\pkgs\sslib\src\*` is byte-identical, after line-ending
  normalisation, to gauss-state-space tag **`v2.0.0` (`0326d40`)** in all
  8 files.
- Its non-source files and `package.json` are still from the v1.0.0
  install, so the manifest reports `"version": "1.0.0"` while the source
  is v2.0.0.
- `scripts/verify_sslib_pin.ps1` correctly reports drift against
  `sslib.pin.json`.

**`sslib.pin.json` is inconsistent.** `pinnedCommit` is `ad15626` (the
`v1.0.0` tag), but its `notes` text still describes an earlier repin to
`2a26076`.

**Causal test (A/B).** The same capped fit was run against `v1.0.0`
sources, `v2.0.0` sources and the shared install. Each source tree was
extracted from git into scratch and loaded with `#include`.
- Over 40 iterations, all three gave **bit-identical** objective
  trajectories and identical final Q.
- The drift therefore does not cause the hang.
- The difference is a stale manifest and pin metadata, not a numerical
  change for TVP-AIDS.

**Unexplained history.** Stage 6 recorded all 16 checks passing in one
real run. The test has a single commit (`4129201`), and v1/v2 give
identical arithmetic, so drift cannot explain the change. The most likely
explanation is that the run simply took a very long time. Its duration
was never recorded.

### 2.4 Risk: unreleased diffuse-phase fix in sslib

gauss-state-space `HEAD` commit `5c300b1` ("Fix diffuse-phase
termination") is **not in any release tag**.

**The bug, in v1.0.0 and v2.0.0.** The exact diffuse filter never ends
the diffuse phase. Rounding residue in `P_inf` can re-trigger spurious
diffuse updates, each adding about +13.8 to the log-likelihood.

**On the TVP-AIDS fixture under v1/v2.**
- The objective jumps upward between accepted iterations (for example
  -1.00 to +0.19, and -1.02 to +0.55).
- That non-smoothness sends CMLmt back and forth between its Newton and
  BFGS steps.

**Under HEAD.**
- The starting objective moves from 3.35 to 15.35.
- It then decreases smoothly and plateaus by about iteration 18.
- It **still** never converges (retcode 2 at 3000 iterations), so the fix
  alone does not unblock the test.
- Fitted Q moves materially (for example Q1 goes from 2.4e-5 to 2.0e-3).

**Consequence.** TVP-AIDS likelihoods and Q estimates produced on any
currently released sslib are contaminated for this model class. Random-walk
states with `k_states > n1` keep the diffuse phase alive for several
periods. Any numeric reference values captured before the fix would move.

## 3. Scope

**In scope:**
- bounding the MLE;
- a well-posed Stage 6 test;
- a deliberate sslib baseline and pin;
- revalidating Stages 2-6 numerics on that baseline;
- re-measuring the documented `n1` scale limit;
- a local TVP release gate;
- docs and support-tier updates.

**Out of scope for this promotion:**
- estimating H;
- full (non-diagonal) Q;
- delta-method standard errors;
- IV treatment of total expenditure;
- adding `sslib` to `package.json` `deps`, or any TVP file to its `src`
  array. `sslib` stays an optional adapter dependency (CLAUDE.md
  "Optional modules").
- changes to shipped non-TVP estimation code.

## 4. Decision gate D0: choose and release the sslib baseline (before any repin)

No repin, resync or reinstall of `C:\gauss26\pkgs\sslib` happens until
this gate is decided by the repo owner and recorded in `PROJECT_STATUS.md`
Decisions.

**Options:**
- **A (recommended).** Wait for gauss-state-space to **tag a release that
  contains `5c300b1`** (for example `v2.0.1`/`v2.1.0`), then pin that tag.
  This is the only option whose likelihood is free of the
  spurious-diffuse bias.
- **B.** Pin `v2.0.0` (`0326d40`, what the shared install already holds)
  and document the known likelihood bias as a limitation. Promotion would
  then carry that limitation in the support tier.
- **C.** Pin `v1.0.0` (`ad15626`). It has the same bias as B, plus the
  older API, so there is no advantage. Not recommended.
- Pinning an **untagged** commit (for example HEAD) is not allowed for a
  release baseline. Stage 6 moved to a tagged pin for exactly this
  reason.

**Gate criteria (all required):**
1. The chosen baseline is a **tag** in gauss-state-space, and its
   CHANGELOG states its GAUSS minimum. GAUSS 26.1.4 meets v2's floor.
2. The gauss-state-space session or owner has been consulted via
   `ListAgents`/`SendMessage` before the shared install is touched
   (CLAUDE.md, shared mutable state).
3. The install method is confirmed for v2-era catalogs:
   - `sync_sslib.ps1`'s hand-built `.lcg` does not emit the
     `typed_returns`/`keywords` annotations v2 relies on (Known Issues).
   - Either extend that script or install through the supported package
     path, then re-verify with a content diff, not mtimes.
4. The owner records the decision, the commit hash and the install method
   in `PROJECT_STATUS.md` Decisions.

## 5. Sequenced promotion plan

Run these in order. Each step's acceptance criteria must hold before the
next step starts.

**P1: Decision gate D0** (section 4).
- Acceptance: the decision is recorded.

**P2: Repin and resync** (authorized, after D0).
- Install the chosen tag.
- Correct `sslib.pin.json`: `pinnedCommit`, `verifiedFiles` hashes,
  `lastVerified`, and `notes`, removing the stale `2a26076` text.
- Acceptance:
  - `scripts\verify_sslib_pin.ps1 -Strict` exits 0;
  - the installed `package.json` version matches the tag;
  - the existing TVP tests in the 2.1 matrix still pass, apart from the
    known `quaidstvpfit_test.e` issue.

**P3: Bound the optimizer** (code; public API change, so bump the version).
- Add `maxIters` and `maxTime` fields to `quaidsTVPControl` (`src/quaids.sdf`).
  Optionally add `printIters`.
- Give them finite defaults in `quaidsTVPControlCreate()`, for example
  500 iterations and 5 minutes. Settle the exact values from measurement.
- Thread them into `ssCtl.ctl` in `_quaidsTVPMLEFit()`.
- Keep `mleRetcode` as the explicit non-convergence signal, and have
  `printQuaidsTVP()` keep flagging it.
- The second, untransformed CMLmt pass is sslib behaviour. Cap it through
  the same `ssCtl.ctl`, and file its negative-variance stepping upstream
  rather than patching sslib here.
- Acceptance:
  - a deliberately non-converging fit returns within its cap with
    `mleRetcode /= 0` instead of hanging;
  - guard case(s) for invalid cap values;
  - `CHANGELOG.md` entry, `package.json`/`CITATION.cff`/`docs/public-api.json`
    version bump, and an updated command-reference page.

**P4: Well-posed Stage 6 test.**
- Rebuild the `quaidstvpfit_test.e` fixture from a real AIDS DGP: extend
  `tests/quaidsfixtures.src` (for example a raw-input variant of
  `_quaidsTVPDynamicSyntheticDGP`) rather than using random shares.
- Take H from the DGP (or from `quaidsFit()`'s `homogSse`) and scale q0
  relative to H, as `examples/14_tvp_aids_estimation.e` does.
- Reuse one fit for checks 1, 1b and 2. Keep the `smooth = 0` fit.
- Keep all existing internal-consistency checks, and add:
  - `mleRetcode == 0`;
  - a wall-clock budget assertion (for example under 60 s per fit);
  - a capped-fit check that returns `mleRetcode /= 0`.
- Acceptance: `ALL N CHECKS PASSED` within budget on the D0 baseline,
  twice, with identical output (determinism).

**P5: Revalidate numerics on the baseline.**
- Rerun the Stage 2-6 tests and `examples/14_tvp_aids_estimation.e`.
- Re-derive any tolerances or reference values that move with the
  diffuse fix.
- Add the second independent check that CLAUDE.md requires: the
  log-likelihood at a fixed Q against an independent implementation (for
  example Python/statsmodels exact diffuse filter) on a committed fixture.
- Acceptance:
  - every TVP test passes;
  - the cross-implementation log-likelihood agrees to a documented
    tolerance;
  - no undocumented tolerance loosening.
- **Status: closed (2026-10-17).** Done for diagonal H on 2026-10-16;
  the non-diagonal-H upstream blocker was fixed in sslib v2.1.1, and all
  three fixed-Q cases now agree with statsmodels over the full sample
  (12/12). See section 10.

**P6: Re-measure the `n1` scale limit.**
- With the bounded optimizer, rerun `n1 = 3` and `n1 = 4`, each with a
  hard timeout.
- Rewrite the "MLE convergence speed degrades sharply with `k_states`"
  remark from the measured results, in `quaidsTVPFit.md` and the README
  support-tier row. The earlier 24-hour `n1 = 4` hang is probably this
  same unbounded, non-converging optimization (a hypothesis, not yet
  reproduced).
- Acceptance: the docs state measured behaviour, including retcodes, with
  no unmeasured claims.
- **Status: done (2026-10-17, Stage 2).** `n1=4` converges in 93
  iterations, deterministically. `n1=3` stops at the cap (`mleRetcode`
  11, or 2 under an iteration cap). The docs now state only these
  measurements. See section 10 and `dev/experiments/tvp_p6_scale/RESULTS.md`.

**P7: Runner hardening.**
- Add a per-test timeout to `tests/run_source_tests.ps1`'s
  `Invoke-GaussBatch`. It has none today, so a hang blocks the suite
  indefinitely.
- Report a timeout as a failure.
- Acceptance: a deliberately hanging scratch script is killed and reported
  FAIL at the timeout.

**P8: Release gate and support-tier change** (section 6).

## 6. Release acceptance criteria (TVP gate)

Every item must pass in a single session on the D0 baseline:

1. `scripts\verify_sslib_pin.ps1 -Strict` passes.
2. `tests\run_source_tests.ps1` passes **without** `-SkipTVPKalman`, with
   every TVP test under the P7 per-test timeout. This includes
   `quaidstvpfit_test.e` and all `tvp*`/`tvpfit*` guard cases.
3. `tests\run_examples_smoke.ps1` passes, including
   `14_tvp_aids_estimation.e`.
4. `scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact
   -InstallArtifact` passes, plus the installed-package API tests.
5. `scripts\run_release_gate.ps1` reports `GO`.
6. The docs are updated:
   - `CHANGELOG.md`;
   - `docs/command-reference/quaidsTVPFit.md` and
     `quaidsTVPControlCreate.md` (new control fields, the baseline's
     limitations);
   - the README support-tier row;
   - `docs/FEATURE_SUPPORT_MATRIX.md` if affected.

   `verify_docs_consistency.ps1` and `verify_docs_quality.ps1` both pass.
7. The gate run (date, sslib tag and commit, timings) is recorded in
   `PROJECT_STATUS.md`.

CI keeps `-SkipTVPKalman`, because the self-hosted runner has no `sslib`.
This TVP gate is therefore local, and its record in `PROJECT_STATUS.md`
is the evidence.

## 7. Promotion claims: safe now vs blocked

**Safe to claim now:**
- The Stage 1 state plumbing and `quaidsTVPStateToFullB()` (56/56) and
  `quaidsTVPElasFit()` (13/13) are tested and have no `sslib` dependency.
- The fixed-Q/H Kalman filter (8/8), RTS smoother (8/8) and Q-only MLE on
  DGP data (7/7, 20 s) pass against the shared install as found. That
  install's source is sslib v2.0.0 code; its manifest reports 1.0.0.
- The input guards fail clearly.
- TVP-AIDS is **Experimental**, as the README already says.

**Blocked until the gate passes:**
- that `quaidsTVPFit()` is tested or supported (its test is
  waived/unverified);
- any claim that the MLE terminates or converges in bounded time;
- numeric accuracy of TVP likelihoods or Q on sslib releases before
  v2.1.0, given the diffuse-phase bias in section 2.4. On v2.1.1 the
  fixed-Q likelihood matches statsmodels to 6.3e-13 for two diagonal-H
  cases and one non-diagonal-H case (P5, section 10);
- any result with a **non-diagonal H** on sslib v2.1.0 or earlier, which
  takes a spurious diffuse update there. Fixed in v2.1.1 (gauss-state-space
  `004c1be`, with its own regression) and covered by P5 case 3. Whether
  non-diagonal H is a supported tier of `quaidsTVPFit()` is a P8
  judgment: the evidence is one fixed-Q likelihood case, not estimated Q;
- any `n1`/`k_states` scale claim beyond what P6 measured (one dataset,
  `n1 = 3` and `4`, `dev/experiments/tvp_p6_scale/RESULTS.md`);
- reproducibility against `sslib.pin.json`, because the install currently
  fails the pin;
- any move of TVP out of Experimental.

## 8. Deferred and preserved work (not abandoned)

The items below are paused behind TVP Priority 1. Their evidence and code
stay exactly as they are.

- **AIDS direct-LM private prototype** (`dev/prototypes/aids_direct_lm/`,
  `dev/AIDS_DIRECT_LM_DESIGN.md`), including the AIDS symmetry study:
  - It is frozen, and its manifest is `symmetry/SYM_FREEZE_MANIFEST.json`.
  - Its verdict is PROCEED for a private, AIDS-only, point-estimate MVP.
  - Its owner decisions remain open.
  - Do not edit the frozen symmetry files.
- **QUAIDS direct LM** and the 10 QUAIDS no-result seeds: these await a
  separately prespecified exploratory protocol.
- **Speed**: these await a separate, controlled benchmark study. No speed
  claim is licensed by any current evidence.

## 9. Diagnosis artifacts

The scratch outputs from the 2026-10-15 diagnosis were session-local and
are not committed. Everything needed to reproduce them is described in
section 2:
- capped runs of the existing tests;
- a replica of `_quaidsTVPMLEFit` with `printIters = 1` and
  `maxIters = 40`/`3000`;
- sslib sources from `git archive ad15626`, `v2.0.0` and `HEAD`, each
  `#include`d from its own extracted directory so that nested
  `#include ss.sdf` resolves to the matching version.

## 10. Progress log

### 2026-10-15: P1-P4 and P7 (task `task_575ce63c08a6`, package 0.4.0)

**P1 / D0.** Option A was taken. gauss-state-space tagged **v2.1.0**
(`7d0e0421a008db81fb4fda0b789bdbef21ffe5f7`, GAUSS 26.1.4+), which
contains `5c300b1`. The shared install was already at v2.1.0 when this
work started; it was installed outside this repo's `sync_sslib.ps1`.
Before pinning, all 9 installed `src/` files and `package.json` were
confirmed identical to the tag after line-ending normalisation.

**P2.** `sslib.pin.json` now records the full commit, `pinnedTag` and
`pinnedVersion`, and five raw-byte hashes (`ss.sdf` was added). The
stale `2a26076` note is gone. `verify_sslib_pin.ps1 -Strict` passes. It
now also compares the installed `package.json` version with
`pinnedVersion`; a negative control with a wrong version exits 1.
`sync_sslib.ps1` keeps `pinnedTag`/`pinnedVersion` current on a repin.
Its hand-built `.lcg` limitation (section 4, criterion 3) was not
changed, because no resync was performed.

**P3.** `quaidsTVPControl.maxIters`/`maxTime` were added and threaded
into `ssCtl.ctl` in `_quaidsTVPMLEFit()` (both `ssFitTVP()` CMLMT passes).
Deviation from the plan: the **defaults are CMLMT's own (1e5, 1e8 min)**,
not finite values. The coordinator brief required backward-compatible
defaults, so existing calls are numerically unchanged. Finite caps are
opt-in and used by the test, the README and the command-reference
examples. Invalid caps fail before estimation; there are five guard
cases.

**P4.** `_quaidsTVPRawSyntheticDGP` (`tests/quaidsfixtures.src`) is a
correctly specified raw-input DGP. `quaidstvpfit_test.e` uses n = 3,
tobs = 120, seed 7, Q = (1e-5 alpha, 1e-4 gamma/beta), H = 2.5e-5 I, and
q0 = 0.1 mean(diag H). All shares lie in [0.156, 0.596].
- Result: 30/30 checks, about 6 s wall time.
- Main fit: 1.1 s, retcode 0.
- maxIters = 1: retcode 2. maxTime = 1e-6 min: retcode 11.
- Two runs gave identical output apart from the printed timings.
- Fits: the main fit, the smooth = 0 fit, one direct reference MLE, and
  two capped fits. The main fit is reused for checks 1-3.

Fixture finding: an earlier scratch draft added a level offset to
`totexp` but not to the simulation regressor. That misspecification
inflated the alpha Q and pushed the beta Q to the zero boundary. With
regressors matched exactly, 15/15 configuration and seed combinations
converged in 1-2.6 s.

**P7.** `run_source_tests.ps1` gained `-TestTimeoutSeconds` (default
1800, 0 disables). On expiry the `tgauss` tree is killed (`taskkill /T`)
and the test or guard is reported FAIL. A deliberately hanging scratch
script was killed at 10.3 s with `TimedOut = True` and left no stray
process.

**Verification (2026-10-15):**
- `run_source_tests.ps1`, with no skip flags and the TVP group included, passed: 31 test files (762 checks) and 38/38 guard cases in 158.8 min, with no timeouts.
- `verify_sslib_pin.ps1 -Strict` passed.
- The manifest, public-API, docs-consistency and docs-quality checks passed, as did `git diff --check`.
- The 0.4.0 artifact was built and verified from a scratch copy of the tree.
- `examples/14_tvp_aids_estimation.e` converged in 19.6 s.

**Still open:**
- **P5:** independent log-likelihood cross-check, plus revalidation of
  example 14 tolerances.
- **P6:** re-measure the n1 = 3/4 scale limit with caps.
- **P8:** the full release gate.

### 2026-10-15: review follow-up (task `task_6c990bb79fcc`)

An independent review of `ccd4066` returned PASS WITH FINDINGS. The
actionable findings are addressed in one follow-up commit, still package
0.4.0:
- **Per-pass cap semantics** are documented everywhere: the source
  comments, README, both command-reference pages, CHANGELOG and status.
  `ssFitTVP()` may run two CMLMT passes with one control, so the worst
  case is about 2x each cap plus up to one iteration of overrun per
  pass. `mleRetcode` is the first pass's code, and a second-pass cap hit
  is not exposed upstream.
- **Example 14** sets `maxIters = 500` and `maxTime = 5`. It converged
  against the 0.4.0 source in about 17 s, and still converged with
  `maxIters = 50`.
- **Timeouts** now also cover `run_examples_smoke.ps1` and the
  installed-package tests in `run_release_verification.ps1`, threaded
  from `run_release_gate.ps1`. `tests/verify_runner_timeouts.ps1` kills a
  GAUSS job hanging on a grandchild process through each runner's
  committed `Invoke-GaussBatch`; it passed 37/37 in about 70 s.
- **Guards:** missing and NaN `maxIters`/`maxTime` were added, for nine
  cap guards in all.
- **Pin:** all 9 v2.1.0 src files plus `package.json` are hashed, and
  unpinned installed src files are flagged. `lib/sslib.lcg` stays out of
  scope, because its determinism is unproven.
- **Release-gate blocker:** the stale `quaids 0.3.0.zip` (SHA-256
  `c5ce5e5d...d52b791`, matching the 0.3.0 release record) was moved to
  `C:\Users\eclow\Documents\quaids-artifact-backups\`.

Not addressed here, by design: the upstream complex/G0518 risk below.

**New residual risk.** sslib's untransformed second CMLMT pass warns
`function evaluated as complex` when a fitted Q is near zero; example 14
on v2.1.0 shows it even though the fit converges. In one scratch
misspecified configuration it aborted with `G0518 Not implemented for
complex` in `cmlmthess`, which `trap` cannot catch. This should be filed
upstream (section 5, P3 note); it is not patched here.

### 2026-10-16: P5 (task `task_8c1e86fde0fd`, still package 0.4.0)

Tests, fixtures and docs only. No `src/` change, no version change, and
no install. sslib and gauss-state-space are untouched.

**Environment.**
- GAUSS 26.1.4.5565 with `GAUSS26_CFG=tests\gauss26_cfg_override`.
- sslib v2.1.0 (`7d0e042`); `verify_sslib_pin.ps1 -Strict` passes.
- Validation only: Python 3.12.0, numpy 2.5.1 and statsmodels 0.15.0.

**Independent cross-check.** `tests/quaidstvp_loglik_xcheck_test.e`
reads the committed fixture `tests/fixtures/tvp_loglik/`. Its README
records provenance, independence, the convention and the commands.
- The data are the P4 DGP: n1 = 2, k_states = 7, tobs = 120, seed 7,
  written with `%.17g`. Regenerating it is byte-identical.
- There are three fixed (diag Q, H) cases.
- The reference is statsmodels' exact diffuse `KalmanFilter`, with Z_t
  rebuilt independently in Python and `loglikelihood_burn = 0`.
  `--check` reproduces it to 0.0.

| Case | H | GAUSS llf | statsmodels llf | max per-period diff | total diff |
|---|---|---|---|---|---|
| 1 (DGP Q) | 2.5e-5 I | 690.6901320251270 | 690.6901320251268 | 4.97e-14 | 2.27e-13 |
| 2 (q0 start) | 2.5e-5 I | 281.2743884578091 | 281.2743884578097 | 3.97e-13 | 6.25e-13 |
| 3 (heterogeneous Q) | [3e-5 1e-5; 1e-5 2e-5] | 688.0306020430203 | 674.2112678504466 | 14.19 (period 4) | 13.82 |

- **Tolerance: 1e-9 absolute**, per period and on the total. That is
  more than three orders above the measured rounding, and about nine
  below a real logic difference (`ln 2pi` = 0.92 per period; a spurious
  diffuse update is about 14). No existing tolerance or reference moved.
- Check (c) holds for every case: sslib's `_tvpLFNImpl()` with
  `transformed = 0`, run through `_quaidsTVPQUpdate()`, equals the filter
  path exactly. This is a **wiring check** of the untransformed
  second-pass path of `ssFitTVP()`, compared with the same sslib filter.
  It is not an independent numerical oracle, and it does not exercise
  the transformed, maximizing first pass. The statsmodels cases are the
  numerical check.
- Result: 11/11 checks pass. Two runs gave identical output. (13/13
  after the audit follow-up below.)

**Case 3 is an upstream defect, not a convention difference.** The
coordinator decided to keep it as a characterization only.
- In sslib v2.1.0's `kalmanFilterDiffuseTVP()`, a non-diagonal H is
  decorrelated by `inv(chol(H)')` row scaling, and then the absolute
  `F_inf > 1e-12` element test is applied.
- At period 4, element 1 removes the last diffuse direction. Element 2's
  rounding residue, inflated by `1/sqrt(H)`, is 2.3e-11 in the numpy
  replica. It passes the test and takes a spurious diffuse update.
- `sslib_nondiag_h_repro.py` reproduces this. The same recursion with
  LDL decorrelation, or a 1e-10 or 1e-8 threshold, matches statsmodels
  to 1.7e-13.
- It needs the diffuse phase to end mid-period. An even n1
  (`k_states mod n1 != 0`) is a common, sufficient way for that to
  happen, not a necessary one: a rank-deficient `Z_t` period can do it
  too. Even then it fires only if the residue's rounding passes 1e-12.
- The test asserts periods 1-3 and prints the full gap as
  `KNOWN-UPSTREAM-DEFECT`. (Superseded by the audit follow-up below,
  which pins the defect instead of leaving it uncounted.)
- **Remaining blocker:** an upstream gauss-state-space fix, a repin, and
  promotion of case 3 to a counted check. The reference column is already
  committed.

**Revalidation (2026-10-16).**
- `tests\run_source_tests.ps1 -TestTimeoutSeconds 1800` was run with no
  skip flags. Result: `run_source_tests.ps1: PASS`, with 32 test files
  (773 checks), 42/42 guard cases (22 of them TVP), 0 FAIL and no
  timeouts. It finished at 06:11 local time.
- TVP results: Stage 1 56/56, Stage 5 13/13, Stage 2 8/8, Stage 3 7/7,
  Stage 4 8/8, Stage 6 30/30, and the P5 cross-check 11/11 (13/13
  after the audit follow-up below).
- No reference value or tolerance needed to move.
- Example 14 was run against the 0.4.0 source through
  `run_examples_smoke.ps1 -ExamplesDir <scratch> -ExampleList
  14_tvp_aids_estimation.e -TestTimeoutSeconds 900`. The scratch copy
  uses absolute `#include`s, because installed quaids is still 0.3.0.
  - Result: PASS in 19 s. Both fits CONVERGED (23 iterations).
  - It still shows the known `function evaluated as complex` warning and
    "covariance failed to invert".

**Still open:**
- P6: the n1 = 3/4 scale re-measurement;
- P8: the release gate;
- the non-diagonal-H upstream defect above;
- the complex/G0518 upstream risk.

### 2026-10-16: P5 audit follow-up (still package 0.4.0)

An independent read-only audit of the P5 commit passed it with findings.
This follow-up changes tests and docs only: no `src/`, version, install,
reference value or upstream change.
- **Case 3 is now pinned, not left uncounted.** The test counts two
  characterization checks:
  - the first period off by more than 1e-9 is 4;
  - the GAUSS-minus-statsmodels total is `13.8193341925737059` within
    1e-9 (GAUSS 26.1.4 with sslib v2.1.0, printed with `%.16f`).
  - Any change to the upstream baseline fails the test. If the defect
    disappears, the test prints an explicit `FAIL`: re-review the fix,
    then update this characterization and the promotion gate.
  - These checks are still not validation. The full sample with a
    non-diagonal H stays unvalidated.
- **Check (c) is described accurately.** It is a wiring check of the
  untransformed (`transformed = 0`) second-pass path, compared with the
  same sslib filter. It is not the maximizing objective and not a
  numerical oracle.
- **The even-n1 wording is narrowed.** An even n1 is a common,
  sufficient way to trigger the defect, not a necessary one.
- A literal CR byte that corrupted `tests\run_source_tests.ps1` in the
  P5 revalidation notes is removed.
- **Negative controls, run on scratch copies:**
  - shifting the pinned gap by 1e-9 gives `1 CHECKS FAILED`;
  - forcing "no divergence" prints the baseline-changed `FAIL` and fails
    the first-period pin.
- **Result:** 13/13 checks pass, identical across two runs.

### 2026-10-17: Stage 1 -- sslib v2.1.1 adopted, P5 closed (still package 0.4.0)

Task `task_dd4a6d307716`. Tests, docs, comments and tooling only: no
estimation code, public API or version change.
- **Baseline.** gauss-state-space tag `v2.1.1` peels to `bf827ae`
  (`bf827aed3a72a5a3b541ce0a6a927c1665a9f992`). It is a patch release on
  v2.1.0 carrying:
  - `004c1be`: the exact diffuse filters test `F_inf * L[j,j]^2 > 1e-12`
    (LDL-decorrelation units) and count diffuse directions from
    `rank(P_inf_0)`. This fixes the P5 case 3 defect.
  - `717abb5`: `ssFit`/`ssFitTVP` convert the transformed-parameter
    covariance by the delta method instead of a second, untransformed
    CMLMT pass. `quaidsTVPFit()` (which sets `positive_vars`) therefore
    runs one CMLMT pass, and its caps bound the whole MLE.
- **Install.** The previous v2.1.0 install was copied to
  `C:\Users\eclow\Documents\sslib-backups\sslib-pre-v2.1.1-20261017-194100`
  (260 files, sha256 manifest beside it, 0 mismatches against the
  source). Then `scripts\sync_sslib.ps1 -SourceRepoPath <gauss-state-space>
  -Commit bf827ae...` installed and repinned. `verify_sslib_pin.ps1
  -Strict` passes (10 files, version 2.1.1), and all 10 tracked files
  match the tag's content after line-ending normalisation.
- **Tooling fix.** `sync_sslib.ps1` failed before doing anything under
  Windows PowerShell 5.1: its `[Parameter(Mandatory)]` attribute makes it
  an advanced script, and 5.1 leaves `$PSScriptRoot` empty while binding
  an advanced script's param defaults, so the `$RepoRoot` default threw.
  This was confirmed in an isolated two-script repro. `$RepoRoot` is now
  resolved in the script body. No other script combines `[Parameter()]`
  with that default.
- **P5 case 3.** On v2.1.1 the two pinned characterization checks failed
  exactly as designed: the gap went from 13.8193341925737059 to 0. Case 3
  is now counted like cases 1 and 2, with every period and the total
  within 1e-9. GAUSS gives 674.2112678504466 and statsmodels
  674.2112678504466, with a max per-period difference of 1.7e-13. Cases
  1 and 2 are bit-identical to v2.1.0. The committed reference and the
  tolerance are unchanged (`statsmodels_reference.py --check`: 0.000e+00).
  The cross-check now has 12 checks, down from 13, and they all pass.
- **Check (c) wording.** It is still a pure wiring check
  (`_tvpLFNImpl`, `transformed = 0`), but it is no longer described as
  `ssFitTVP()`'s second pass, which v2.1.1 removed. The caps and
  `mleRetcode` docs no longer claim two passes or a hidden second-pass
  code. Those claims are corrected in the README, the
  `quaidsTVPFit`/`quaidsTVPControlCreate` pages, `quaids.sdf`,
  `quaidstvpfit.src`, `quaidstvpmle.src` and example 14 (comments only).
- **Revalidation on v2.1.1.** Stage 1/5/4/2/3/6 tests give 56/13/8/7/8/30.
  Example 14 (scratch copy including the 0.4.0 sources, as in P5) passes
  in 11 s; both fits CONVERGED in 23 iterations. The `function evaluated
  as complex` and "covariance failed to invert" warnings seen on v2.1.0
  no longer appear.

**Still open:**
- P6: the n1 = 3/4 scale re-measurement;
- P8: the release gate. This includes the support-tier judgment for
  non-diagonal H and for the delta-method covariance, and a re-check of
  the complex/G0518 risk, whose second-pass trigger v2.1.1 removed.

### 2026-10-17: Stage 2 -- P6 bounded scale re-measurement (still package 0.4.0)

Task `task_ee20ce3a0a46`. Docs, comments and experiment evidence only: no
code, public API or version change.
- **Dependency state before and after:** `verify_sslib_pin.ps1 -Strict`
  passes (v2.1.1, `bf827ae`). `quaidstvp_loglik_xcheck_test.e` gives
  12/12 and `quaidstvpfit_test.e` gives 30/30. P5 has not regressed.
- **Harness:** `dev/experiments/tvp_p6_scale/`:
  - the template and generated `.e` inputs;
  - `run_p6.ps1`, which runs jobs sequentially with a `GAUSS26_CFG`
    override and an outer 900 s hard timeout using `taskkill /T`;
  - full logs, `summary.csv` and `RESULTS.md`.

  The setup is example 14's (`quaidsExampleData(500, 204)`, `H` from a
  static fit, `q0 = 0.1*mean(diag(H))`) with only the number of goods
  changed.
- **Results** (10 jobs, no outer timeout hit, no complex warnings):

  | `n1` | caps | runs | `mleRetcode` | iterations | CMLMT min |
  |---|---|---|---|---|---|
  | 3 | 500 x 5 | 3 | 11, 11, 11 | 457/438/437 | 5.02-5.04 |
  | 3 | 100 x 60 | 2 | 2, 2 | 101 | 1.19-1.23 |
  | 4 | 500 x 5 | 3 | 0, 0, 0 | 93 | 3.52-3.62 |
  | 4 | 100 x 60 | 2 | 0, 0 | 93 | 3.70-3.89 |

  - `n1=4`'s output is byte-identical across all 5 runs.
  - `n1=3` is identical under the iteration cap and differs under the
    time cap, where its endpoint depends on timing. CMLMT reports
    "covariance failed to invert" there.
  - Nine of `n1=4`'s 18 `Qfit` elements are under 1e-9, and one is 2.955.
    These were measured, not validated.
- **Doc changes:**
  - the README support-tier row;
  - `quaidsTVPFit.md`'s scale remark, now a measured table, plus its
    overrun wording (at most 0.036 CMLMT min past a 5-min cap);
  - example 14's header comments.

  The old "n1=4 pathologically slow" claim is withdrawn.
- **Not measured:** larger `n1 = 3` budgets, other seeds, `tobs` or `q0`
  rules, `n1 >= 5`, and `Q` accuracy.

**Still open:**
- P8: the release gate. This includes the support-tier judgment for
  non-diagonal H and for the delta-method covariance, and a re-check of
  the complex/G0518 risk.
