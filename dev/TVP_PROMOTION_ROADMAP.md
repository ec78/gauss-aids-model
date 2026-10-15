# TVP-AIDS Promotion Roadmap

_Created: 2026-10-15. Status: **Priority 1 (active)**. Nothing in this plan
has been implemented yet._

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
- numeric accuracy of TVP likelihoods or Q on released sslib, given the
  diffuse-phase bias in section 2.4;
- any `n1`/`k_states` scale claim beyond what P6 measures;
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
