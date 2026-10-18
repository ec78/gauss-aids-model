# PROJECT_STATUS.md

Current work state for this repository. Read this (plus `git log`/`git
status`) at the start of a new session instead of relying on prior chat
history. See `CLAUDE.md` for durable project knowledge and
`dev/GOLD_STANDARD_TODO.md` for the full historical decision log.

_Last updated: 2026-10-17 (Stage 3 / P8)_

## Current Objective

**PRIORITY 1 (set 2026-10-15): promote TVP-AIDS out of
Experimental. Done 2026-10-17 (P8): Supported (scoped).** Plan, evidence, decision gate and acceptance criteria:
**`dev/TVP_PROMOTION_ROADMAP.md`**.

**Progress 2026-10-15 (task `task_575ce63c08a6`, package 0.4.0, local
commit only):** P1 (D0 = option A), P2, P3, P4 and P7 are done; see
roadmap section 10.
- sslib is pinned to tag v2.1.0 (`7d0e042`), and
  `verify_sslib_pin.ps1 -Strict` passes.
- `quaidsTVPControl.maxIters`/`maxTime` were added. Their defaults are
  CMLMT's own, so existing calls are unchanged.
- `quaidstvpfit_test.e` passes 30/30 in about 6 s on a well-posed raw
  DGP, with deterministic output.
- There are five invalid-cap guards.
- `run_source_tests.ps1` has a per-test timeout
  (`-TestTimeoutSeconds`, default 1800).
- Full `run_source_tests.ps1` (TVP included, no skips) passed: 31 files, 762 checks and 38/38 guards in 158.8 min, with no timeouts.
- **Review follow-up (task `task_6c990bb79fcc`, same day, local commit
  only, still 0.4.0):** an independent review of `ccd4066` passed with
  findings, and these are now addressed:
  - The docs, source comments and changelog now say the caps apply
    **per CMLMT pass**. sslib may run two passes, so the worst case is
    about 2x each cap. `mleRetcode` is the first pass's code; a cap hit
    in the second pass is not exposed upstream.
  - `examples/14_tvp_aids_estimation.e` sets `maxIters = 500` and
    `maxTime = 5`.
  - The per-job timeout now also covers `run_examples_smoke.ps1` and the
    installed-package tests in `run_release_verification.ps1`, and is
    threaded through `run_release_gate.ps1`. The negative control
    `tests/verify_runner_timeouts.ps1` passes.
  - Four more guards cover missing and NaN caps, for nine in all.
  - The sslib pin now hashes all 9 v2.1.0 src files plus `package.json`.
    The generated `.lcg` is not hashed (documented).
  - The stale `quaids 0.3.0.zip` was moved out of the repo root (see
    Known Issues).
- **Final audit follow-up (task `task_b255e6d2f250`):** a read-only
  audit of `94c77f6` passed with four low findings, and these are now
  fixed:
  - While a timeout is active, the runners drain the output streams
    for at most 30 s instead of waiting without limit. A process that
    escaped the killed tree but still holds the stdout/stderr pipes can
    no longer stall them. `verify_runner_timeouts.ps1` covers this case.
  - The docs got a NaN in the invalid-cap list, a missing paragraph
    break, and two rewrapped long lines.
- **P5 (2026-10-16, task `task_8c1e86fde0fd`, local commit only, still
  0.4.0, tests and docs only):** an independent fixed-Q log-likelihood
  cross-check, `tests/quaidstvp_loglik_xcheck_test.e` (13/13 after the audit
  follow-up), compares
  against statsmodels 0.15.0's exact diffuse filter (numpy 2.5.1,
  Python 3.12.0) on the committed fixture `tests/fixtures/tvp_loglik/`.
  - **Diagonal H passes.** Case 1 (DGP Q) gives 690.6901320251270 in
    GAUSS vs 690.6901320251268 in statsmodels. Case 2 (q0 start) gives
    281.2743884578091 vs 281.2743884578097. The largest difference is
    4.0e-13 per period and 6.3e-13 on a total; the tolerance is 1e-9
    absolute.
  - A wiring check: `_tvpLFNImpl` with `transformed = 0`, the untransformed
    second-pass path, equals the filter path exactly. It uses the same
    sslib filter, so it is not a numerical oracle and is not the
    maximizing (transformed) pass.
  - **Non-diagonal H fails upstream.** Case 3 gives 688.0306 vs 674.2113:
    sslib's absolute `F_inf > 1e-12` test, applied after Cholesky
    scaling, takes a spurious diffuse update at period 4. This was
    hand-verified in numpy. Periods 1-3 are asserted. The defect is
    pinned by counted characterization checks (first divergence at
    period 4, total gap 13.8193341925737059 within 1e-9), so any upstream
    change fails the test. The full sample is still not validated.
  - Full `run_source_tests.ps1` (no skips, 1800 s timeout) passed:
    32 files, 773 checks, 42/42 guards and 0 FAIL. Stages 1-6 gave
    56/13/8/7/8/30.
  - Example 14 passed through the bounded `run_examples_smoke.ps1`
    (scratch copy using the 0.4.0 source includes) in 19 s, CONVERGED.
  - No reference or tolerance changes. See roadmap section 10.
- **Stage 1 (2026-10-17, task `task_dd4a6d307716`, local commit only,
  still 0.4.0; tests, docs, comments and tooling only): sslib v2.1.1
  adopted, P5 closed.** See roadmap section 10.
  - `v2.1.1` peels to `bf827ae`. Before the install, the v2.1.0 copy was
    backed up, hash-verified, to
    `C:\Users\eclow\Documents\sslib-backups\sslib-pre-v2.1.1-20261017-194100`.
    Then `sync_sslib.ps1 -Commit bf827ae...` installed and repinned, and
    `verify_sslib_pin.ps1 -Strict` passes.
  - `sync_sslib.ps1` itself was fixed: under Windows PowerShell 5.1 its
    `$RepoRoot` param default threw before doing anything, because
    `$PSScriptRoot` is empty during an advanced script's param binding.
  - **Non-diagonal H now passes.** Case 3 gives 674.2112678504466 in both
    GAUSS and statsmodels, with a max per-period difference of 1.7e-13.
    It is now a counted full-sample check, like cases 1 and 2, which are
    bit-identical to v2.1.0. The cross-check passes 12/12. The reference
    and tolerance are unchanged.
  - v2.1.1 also replaced `ssFitTVP()`'s second, untransformed CMLMT pass
    with a delta-method covariance. Docs that described two passes per
    fit, a hidden second-pass retcode, or check (c) as "the second pass"
    are corrected.
  - TVP tests give 56/13/8/7/8/30. Example 14 passes in 11 s
    (CONVERGED, 23 iterations), with no complex or "failed to invert"
    warnings.
- **Stage 2 / P6 (2026-10-17, task `task_ee20ce3a0a46`, local commit
  only, still 0.4.0; docs, comments and experiment evidence only): bounded
  `n1` scale re-measurement done.** Evidence:
  `dev/experiments/tvp_p6_scale/RESULTS.md`; also see roadmap section 10.
  - Before and after, the pin passes `-Strict`, P5 gives 12/12 and Stage
    6 gives 30/30.
  - Setup: example 14's data and setup with only the number of goods
    changed. Ten sequential jobs, each under an outer 900 s timeout,
    which was never hit.
  - **`n1=4` converges.** `mleRetcode = 0` in 93 iterations, about
    3.5-3.9 min, with byte-identical output in 5/5 runs under caps
    500 x 5 and 100 x 60.
  - **`n1=3` does not converge.** At 500 x 5 it gives `mleRetcode = 11`
    in 3/3 runs, after 437-457 iterations, with an endpoint that depends
    on timing. At 100 x 60 it gives `mleRetcode = 2` in 2/2 runs,
    deterministically. CMLMT reports "covariance failed to invert".
  - The README row, `quaidsTVPFit.md` and example 14's comments now
    state only these measurements. The old claim that `n1=4` was
    pathologically slow is withdrawn.
- **Stage 3 / P8 (2026-10-17, task `task_d49ed011d38a`, local commits
  only, still 0.4.0; docs only, no code, API or version change): the
  release gate passed, and TVP-AIDS is promoted to Supported (scoped).**
  See roadmap section 10 for the full record.
  - **Gate run 1** was on the D0 baseline, clean HEAD `d02a812`, on
    2026-10-17 from 20:39:32 to 20:45:49 local (377.1 s). Environment:
    GAUSS 26.1.4.5565 and PowerShell 5.1.26100. sslib was tag v2.1.1,
    commit `bf827aed3a72a5a3b541ce0a6a927c1665a9f992`.
    - `verify_sslib_pin.ps1 -Strict` passed (10 files, version 2.1.1).
    - `run_release_gate.ps1 -TestTimeoutSeconds 1800` reported **GO**.
      It runs `run_release_verification.ps1 -BuildArtifact
      -ForceArtifact -InstallArtifact` with no skip flags.
    - Source suite: 32 files, 774 checks, 42/42 guards (22 TVP), 0 FAIL,
      no timeouts. The TVP tests gave Stage 1 56, Stage 5 13, Kalman 8,
      MLE 7, smoother 8, `quaidstvpfit_test` 30 and cross-check 12, all
      passing. The TVP group took about 15 s.
    - Manifest, public-API (54 procs, 16 structs), docs-consistency and
      docs-quality checks passed.
    - The artifact was built and `verify_release_artifact.ps1` passed
      (63 doc pages). The installed-package `package_public_api.e` and
      `_core_only.e` tests passed. `run_examples_smoke.ps1` passed 15/15.
      In example 14, both fits CONVERGED in 23 iterations (about 0.08
      CMLMT min each).
    - No `function evaluated as complex`, G0518 or "failed to invert"
      line appeared anywhere in the gate log. That closes the complex
      re-check on v2.1.1.
    - Artifact SHA-256: `290810136E491781A7317BCC7135C0285744C3ECBCF5F4F0AB1B3FDE7F05969F`.
      Record: `release_records/0.4.0-20261017-204549.json`, gitignored.
    - The convergence sweep was captured: iterated AIDS 39/19/42% and
      QUAIDS 54.5/21.5/24%, unchanged.
    - **Warning:** the install step's `Remove-Item C:\gauss26\pkgs\quaids`
      failed ("used by another process"), as it did in the 0.3.0 gate,
      and the artifact was extracted over the old directory. A content
      diff found all 210 artifact files byte-identical in the install,
      with the generated `lib/quaids.lcg` as the only extra file.
    - This run is **not** the final certification. An independent
      review found that its artifact's 0.4.0 release notes were stale:
      they still said sslib v2.1.0, two CMLMT passes, and non-diagonal H
      unvalidated.
  - **Tier decision (coordinator-approved): Supported (scoped).** The
    scope is:
    - a fixed, diagonal H;
    - finite `maxIters`/`maxTime`;
    - a mandatory `mleRetcode == 0` check.

    There is no convergence guarantee, no standard-error claim and no
    IV/endogeneity treatment. A non-diagonal H is validated for the
    fixed-Q likelihood only and is outside the supported estimation
    scope. The Q delta-method covariance is not exposed in
    `quaidsTVPOut`, so nothing is claimed for it.
  - **Docs:**
    - CHANGELOG: Unreleased is folded into a corrected 0.4.0 section,
      and a P8 tier entry is added.
    - README: the tier row and the TVP section, including a fix to the
      stale two-pass code comment.
    - FEATURE_SUPPORT_MATRIX: a new TVP row.
    - The `quaidsTVPFit`/`quaidsTVPControlCreate` pages: scope, and
      removal of the stale second-pass text.
  - **Final gate (run 2, the P8 certification): GO.** It ran from clean
    commit `a15f190`, which holds all the P8 docs and tier changes, on
    2026-10-17 from 20:51:28 to 20:57:42 local (374.2 s; release
    verification 364.7 s). Same environment and command as run 1, with
    the same results:
    - pin `-Strict` PASS;
    - 32 files, 774 checks, 42/42 guards (22 TVP), 0 FAIL, no timeouts;
    - Stage 6 30/30 and cross-check 12/12;
    - verify_release_artifact PASS, both installed-package tests PASS;
    - examples 15/15, with example 14's two fits CONVERGED;
    - zero complex/G0518/"failed to invert" lines.

    Artifact `quaids 0.4.0.zip` SHA-256:
    `51D8E00F57FECF273C02DBBA67AD4E12F5BCE34D4CCBE324E160465B3298A65E`.
    Record: `release_records/0.4.0-20261017-205742.json`. Its
    `releaseNotes` state sslib v2.1.1, one CMLMT pass, the P5 3/3
    full-sample cross-check, the P6 measurements and the P8 Supported
    (scoped) tier. The stale two-pass and "not validated" text is gone.
    The install hit the same `Remove-Item` "in use" warning. The content
    diff matched again: 210/210 files identical, with the `.lcg` as the
    only extra.
- **Open:** the 0.4.0 `date-released` (2026-10-15, in CHANGELOG and
  CITATION.cff) is left for the owner to set at tagging. Nothing is
  pushed, tagged or published.

The rest of this section is the original 2026-10-15 read-only diagnosis
(task `task_03d43536665b`). It is kept as evidence; its "TIMEOUT",
"waived" and "drift" items are resolved by the progress above.
- **Bounded test matrix (shared install as found, hard caps):**
  - `quaidstvp_test` 56/56;
  - `quaidstvp_elas_test` 13/13;
  - `quaidstvp_kalman_test` 8/8;
  - `quaidstvp_mle_test` 7/7 (20 s);
  - `quaidstvp_smooth_test` 8/8;
  - both `tvpfit_bad_*` guards fail with their intended messages;
  - **`tests/quaidstvpfit_test.e` TIMEOUT at 600 s** (killed while still
    inside its first MLE). It remains waived/unverified, not passing.
- **Root cause (repository defect):**
  - `_quaidsTVPMLEFit()` runs CMLmt with its effectively unbounded
    defaults (`maxIters = 1e5`, `maxTime = 1e8` min).
  - `quaidsTVPControl` cannot set either.
  - At about 35 ms per iteration that is about 58 min per CMLmt run,
    matching the earlier 59-minute stall. sslib's second, untransformed
    CMLmt pass and the test's four fits multiply it.
- **Ill-posed fixture:**
  - the shares are random, not an AIDS DGP;
  - H is fixed at 0.01, about 1/4 of the share variance;
  - q0 is a flat 0.001;
  - so Q variances go to the zero boundary, and the fit hits retcode 2
    even at 3000 iterations.
- **sslib drift is real but not causal.**
  - The shared install's `src/` is byte-identical to tag `v2.0.0`
    (`0326d40`), while its `package.json` still says 1.0.0.
  - `sslib.pin.json`'s notes are stale.
  - A v1.0.0-vs-v2.0.0-vs-shared-install A/B gave bit-identical
    optimizer trajectories.
- **Risk: unreleased sslib diffuse-phase fix.**
  - gauss-state-space `5c300b1` is untagged.
  - It removes spurious diffuse updates that make the TVP likelihood
    jumpy and biased on v1 and v2, and it moves fitted Q materially.
  - It does not by itself make the fixture converge.
- **Before any repin:** decision gate D0 (roadmap section 4) must choose
  a *tagged* sslib baseline (recommended: a tag containing `5c300b1`).
- **Deferred and preserved (not abandoned):**
  - AIDS direct-LM, including the frozen symmetry study, whose evidence
    is unchanged;
  - QUAIDS direct LM / no-result seeds;
  - any speed study.

**`v0.3.0` has been RELEASED.** `scripts\run_release_gate.ps1` reported
`GO` (artifact SHA256
`c5ce5e5ded68c64c6e5bdeaad9fdaf01bd82fca384d1aa7cef202f926d52b791`;
release record at `release_records/0.3.0-20260924-093017.json`,
gitignored). Tagged `v0.3.0` on `cc22950` and pushed; GitHub release
created at https://github.com/ec78/gauss-aids-model/releases/tag/v0.3.0
with `quaids.0.3.0.zip` attached (digest matches the gate's own
checksum, confirmed). One transient, non-blocking hiccup during the
gate's install step: `Remove-Item` on `C:\gauss26\pkgs\quaids` hit "file
in use by another process" (likely a lingering handle from the prior
verification run) and left the old directory in place instead of a clean
wipe, but `build_lcg.ps1` still wrote into it successfully immediately
after and every downstream check (installed-package API test, all 15
examples, convergence sweep) passed clean against the result — worth a
clean re-run if it recurs, not something that needed fixing this time.
No other action items from this release; see below for prior context.

**TVP-AIDS (time-varying-parameter AIDS via a Kalman filter) is
COMPLETE** — a repo-owner-requested extension beyond the now-largely-
complete public release roadmap, staged as Stage 0–6 (see
`dev/GOLD_STANDARD_TODO.md`'s "TVP-AIDS initiative" section, now marked
complete, for the full plan and every stage's own detailed writeup).
**All six stages are committed and pushed to `origin/master`** (Stages
0–2 as `086c97a`/doc-sync `8ed8ea3`; Stage 3 as `a36c7a6`; Stage 4 as
`e3ef801`/doc-sync `3d9f3d7`; Stage 5 as `c58c511`; **Stage 6 (the final
stage) as `4129201`**). CI confirmed `success` for every one of these via
`gh run list` (Stage 6's own run id `35878499421`).

Stage 6 published real public API (`quaidsTVPFit()`/`printQuaidsTVP()`/
`quaidsTVPControlCreate()`, plus `quaidsTVPStateToFullB()`/
`quaidsTVPElasFit()` promoted from private), a real automated `sslib`
commit-pinning mechanism, docs, a new example, and a version bump to
`0.3.0`. The full release-verification pipeline
(`scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact
-InstallArtifact`) is green end to end. See Completed Work's Stage 6
entry for the full account, including six real bugs found and fixed only
once real `sslib` access became available (an array-typed struct-field
sentinel, a `string`-vs-`matrix` struct field type mismatch, a
`$|`-vs-`$+` character-matrix type mismatch found in TWO separate places,
a `diag()` vs `diagrv()` mistake, an MLE-hangs-at-`n1=4` scale limit, and
a `$+`-broadcast printer bug) — none of these were catchable by
`#include`-based testing alone, only by actually building, installing,
and running against a real `sslib` install and the installed package.
`sslib` itself is now pinned to the tagged, stable `v1.0.0` (`ad15626`)
release of `gauss-state-space`, not a commit hash on a moving branch —
see Decisions for why. Stage 3 is still Q-only MLE (H stays caller-fixed)
— see Decisions for why, a repo-owner-approved scope call made explicitly
before Stage 3 was written, not assumed.

**The latest completed work is the AIDS/QUAIDS convergence-correctness
research** in `dev/experiments/convergence-anderson/`. Its latest
result is the **completed E2 confirmation (2026-10-08/09)**, in which all
binding gates pass. See the 2026-10-09 update below and
`dev/experiments/convergence-anderson/E2_RESULTS.md`.

Commit history (all pushed; `origin/master` was `158313a` on
2026-10-09):
- **`e90a568`:** the matched 30-seed pilot.
- **`ad986fc`:** the protocol freeze.
- **`158313a`:** the E2 evidence and documentation commit.

The research itself remains experiment-only and is not shipped.

**Follow-up (2026-10-09): private AIDS direct-LM prototype.** The owner
accepted decisions Q1-Q5 in `dev/AIDS_DIRECT_LM_DESIGN.md` for a
**private** prototype, and it is now implemented:
- **Code:** `dev/prototypes/aids_direct_lm/quaidsdirect.src`, entry point
  `_quaidsDirectFit()`. It lives under `dev/` so it is never packaged
  (`build_package.ps1` copies `src/` and `tests/`, not `dev/`) and never
  loaded by `library quaids;`; a caller must `#include` it explicitly.
  No public API, version bump or release claim (`CHANGELOG.md` has only
  an `Unreleased` maintainer note).
- **Model:** AIDS only, homogeneity-only, point estimates only. No
  QUAIDS, no symmetry, no inference, no downstream integration.
- **Arithmetic:** an exact copy of the frozen E2 code.
  `quaidsdirect_evidence_test.e` reproduces every stored E2
  cold-LM cell and the `L16` selection on 7 held-out seeds with zero
  difference.
- **Other tests:** `quaidsdirect_test.e` (recovery, selection and
  no-result, determinism and RNG isolation, default compatibility) and
  five `guard_error_cases/directlm_*.e` scripts, all in `dev/prototypes/aids_direct_lm/`
  next to the code. All are wired into
  `run_source_tests.ps1`; skip the evidence test with
  `-SkipDirectEvidence`.
- **Verification (2026-10-10):** frozen manifest 34/34; the frozen/E2
  directory is unchanged since `158313a`; `quaidsdirect_test.e` 28/28;
  `quaidsdirect_evidence_test.e` 80/80 at zero difference; the five
  guards pass; `run_source_tests.ps1 -SkipBootstrap -SkipCurvature
  -SkipPubtable -SkipTVPKalman` PASS (23 tests, 23 guards); manifest,
  API and docs checks PASS. An earlier full run also passed the pubtable,
  curvature, curvature-bootstrap and TVP tests up to `quaidstvp_smooth`.
- **Known limitation: `tests/quaidstvpfit_test.e` is UNVERIFIED.** It
  did not complete in two runs (one in the full suite, one isolated with
  the `GAUSS26_CFG` override). The isolated run sat 59 minutes inside its
  first sslib/CMLmt fit before it was stopped. The test uses no file the
  prototype touched. `verify_sslib_pin.ps1` reports the shared sslib
  install has drifted from `sslib.pin.json` (`sstvp`, `sskalman`,
  `ssmain`, `ssstructural`). That drift was already there and is the
  likely cause. It was left alone (shared install) and needs a separate,
  authorized resync before the TVP fit test can be trusted again.
  **Superseded 2026-10-15:** drift was confirmed NOT to be the cause.
  The hang is the unbounded CMLmt fit on an ill-posed fixture. See
  Priority 1 above and `dev/TVP_PROMOTION_ROADMAP.md`.
- **Open owner items:** the sslib drift above; any public move needs a
  new decision (design doc §12).
The original Anderson research and benchmark were committed and pushed
earlier as `7e98720`. The methodology audit corrected the original
multi-start selector direction, added a numerically safer experimental
Anderson path, and completed the recommended 200-seed x 16-start
benchmark for both iterated AIDS and QUAIDS. The push-triggered CI run
completed successfully (`36855314457`). This work is deliberately NOT
wired into the installed package or CI test surface; no shipped `src/`
file changed. See the Completed Work entry and the experiment's
`SUMMARY.md`/`CLAUDE_HANDOFF.md` for the evidence and handoff.

**Phase 1 of that direct-nonlinear-optimization follow-up (define and
validate the objective) completed, then a Phase 1b audit corrected its
own headline conclusion.** Both are in the same experiment directory,
committed with the Phase 3 checkpoint `a776717`. Phase 1
(`direct_objective.src`'s initial objective/pack/
scale machinery, `phase1_objective_validation.e`) validated the
objective cleanly but reported a NEGATIVE recovery-gate result for
QUAIDS using a generic optimizer and a flat infeasibility penalty.
**Phase 1b found that result was confounded** by (1) a fixture whose
first-stage IV residual was itself made weakly identified by the same
noise-shrinking used to keep the GLS criterion well-conditioned, and (2)
an optimizer with no real safeguards (a flat `1e8` penalty that a
generic solver could mistake for a converged plateau). With a corrected
fixture and a real Levenberg-Marquardt solver (`_directLM`, Sinv-weighted
Gauss-Newton gradient/Hessian validated against finite differences to
~1e-9, rejects infeasible steps, requires an actual/predicted-reduction
test), direct optimization recovered the true structural parameters to
~0.001-0.002 accuracy (excluding a still-weakly-identified but
non-rank-deficient nuisance correction row) in EVERY ONE of 8 trials
(both models, 1%/5% perturbations, 2 independent STARTING-POINT
perturbations each -- one stationary-fixture dataset per model, not an
independent data redraw per trial) -- the Jacobian is full numerical
rank at truth in both models; no weak identification was found outside
the nuisance row, on this one fixture draw per model and this file's
documented parameter scaling. The fixture's truth point is NEAR- rather
than exactly stationary (observed ||FD gradient||_inf 0.084 AIDS / 0.183
QUAIDS, not zero).

**Phase 3's continuation/homotopy smoke gate then ran and PASSED on both
models**, answering Phase 1b's own open question (does a validated
direct optimizer work from far/structured starts on real noisy data, not
just Phase 1b's near-stationary fixture). On known-difficult seeds (AIDS
seed=6, QUAIDS seed=9 -- both seeds where ALL 16 structured starts x
safeguarded Anderson already failed in the completed 200-seed benchmark),
from both a Stone start and an aggressively perturbed/sign-flipped far
start: AIDS sees continuation tie cold direct LM (both beat safe
Anderson, which did not converge); QUAIDS (the harder case) sees
continuation reach a dramatically better, remarkably consistent solution
(structural recovery error 0.081, vs. cold LM's 4.02/12.29 and safe
Anderson's 2.16/0.19 depending on start) from BOTH starting points,
while cold LM lands on a different, worse, converged-but-wrong point
depending on the start. A real robustness bug (an unguarded `invpd()`
that could crash, not just mis-score, on a wild far-start step) was
found and fixed along the way, confirmed via two independent byte-for-
byte regression checks not to have changed Phase 1b's own results.

**Phase 3 hardening (2026-10-05) is done, and the strengthened smoke gate
passes cleanly in all 4 scenarios.**

What changed:
- The continuation driver now solves the real homotopy=0 problem before
  advancing.
- Failure modes have separate counters, and a new LM retcode 5 marks
  covariance-inversion failure.
- `invpd()` is fail-closed and was validated by fault injection
  (`phase3_safeguard_diag.e`, 14 of 14 checks pass).
- Endpoint differences are now reported explicitly.
- The verdict now requires an absolute correctness condition,
  structErr ≤ 1 (excluding u), with sensitivity reported at .5/1/2/5.

Corrected numbers:
- **AIDS structural error is 0.058, not 0.724.** The earlier figure was
  the u-row error, caused by a wrong row slice in the gate.
- **QUAIDS continuation endpoints are numerically equivalent, not
  identical.** The two starts differ by 5.8e-9 in scaled θ.
- **Stage zero changed no endpoint** (≤ 1.5e-8), only iteration counts.

Phase 1b's audit was rerun, and every numerical line is identical.

The hardened Phase 3 checkpoint was committed as `a776717`. The matched
pilot below is committed and pushed as `e90a568` (`HEAD`, `master`,
`origin/master`). It is experiment-only and not shipped, and it does not
authorize a 200-seed scale-up.

**Update (2026-10-05): the matched 30-seed pilot is complete.**

The pilot covers seeds 1-30, both models, and the benchmark's 16 starts.
It compares safe Anderson with its full portfolio against cold direct LM
and against hardened continuation + LM. Every one of the 2,880 cells
finished.

The Stage A harness and metric gate passed before the run:
- 27/27 Part 1 checks;
- 58/58 analysis checks;
- 6/6 resume checks;
- 377/377 checks on the final 30-seed output.

The pilot's Anderson arm reproduces the 200-seed benchmark's per-seed
results on all 60 seed/model cases. Two validation-only tolerance fixes
are recorded in `PILOT_PROTOCOL.md`; no pilot definition changed.

Selected correct counts out of 30 (benchmark-compatible `bS` ≤ 1 /
raw structural ≤ 1, excluding u):

| | Anderson | cold LM | continuation |
|---|---|---|---|
| AIDS | 19 / 22 | 27 / 30 | 26 / 30 |
| QUAIDS | 10 / 13 | 23 / 29 | 23 / 29 |

What the results show:
- **Direct LM is a genuine solver improvement over Anderson.** The
  selected and oracle results both improve.
- **Continuation ties cold LM at 16 starts.** It is more robust per start
  (QUAIDS Stone-only raw 27 vs 23; 476/480 vs 379/480 converged) at about
  12× the wall time (QUAIDS median 98 s vs 8 s per seed).
- **Solution selection is not the bottleneck.** Selected and oracle agree
  on 178 of 180 cases.
- **Most remaining `bS` failures are right-basin optima.** They fail on
  the nuisance u row, or because `quaidsFit()`'s symmetry projection is
  applied at non-fixed-point endpoints.
- **The genuine wrong-basin cases are complementary.** On QUAIDS seed 5
  continuation fails; on seed 17 cold LM fails.

**Recommendation:**
- Use a model-specific hybrid: cold LM for AIDS, and continuation (or a
  cold LM plus continuation pool) for QUAIDS.
- Decide the direct-endpoint correctness metric before any 200-seed run.
  The preferred option is to impose symmetry inside the direct
  objective.
- Do not run 200 seeds yet.

Full account: `dev/experiments/convergence-anderson/SUMMARY.md`, section
"Matched 30-seed pilot". Next entry point:
`dev/experiments/convergence-anderson/CLAUDE_HANDOFF.md`.

**Update (2026-10-09): E2 confirmation COMPLETE. All binding gates pass
for both models.**

The run used the frozen protocol from freeze commit `ad986fc`, on
held-out seeds 31-200 (n = 170 per model). `FREEZE_MANIFEST.json` verified
34/34 before and after.

| Model | Primary configuration | A2: correct / 170 (CP lower) | A4: wrong answers returned |
|---|---|---|---|
| AIDS | cold LM, 16 starts (`L16`) | **170** (0.9785 ≥ 0.90): PASS | **0**: PASS |
| QUAIDS | cold LM + continuation pool (`L16+C16`) | **160** (0.8945 ≥ 0.85): PASS | **0**: PASS |

- **A1:** all four parts pass. These are validation (340/340 benchmark
  reproductions), analysis integrity, determinism (0 differences) and
  resume (6/6).
- **A5:** no regressions against Anderson.
- **A3 (descriptive only):** improved/regressed/tied 44/0/126 for AIDS
  and 86/0/84 for QUAIDS.
- **QUAIDS has 10 visible no-result seeds:** 43, 112, 125, 134, 151,
  160, 167, 182, 187 and 190. On all of them every cold-LM start was
  infeasible or inversion-failed, and every continuation run failed.
  Anderson also failed all 10, with no result or a wrong answer.

**Limitations:**
- synthetic DGP only;
- no truth-free reliability signal;
- the 10 QUAIDS no-result seeds;
- A3 is non-confirmatory;
- timing is descriptive, not a speed conclusion (wall time per seed:
  AIDS median 4.2 s, QUAIDS median 79.5 s);
- the production symmetry design is unresolved.

The evidence package (`e2_out/`, `e2_det_out/`, `confirm_e2_out/`,
`e2_logs/`) and `E2_RESULTS.md` are committed with this update. The
exception is the two per-cell GAUSS console logs (about 1.9 MB each,
redundant with the cell CSVs), which are kept locally and git-ignored by
the owner's size decision; their SHA-256 values are recorded in
`E2_RESULTS.md`. **No
shipped `src/` file changed.**

**Update (2026-10-08; historical, superseded by the 2026-10-09 update
above): the confirmation protocol is frozen; the confirmation run has not
been run.**

`dev/experiments/convergence-anderson/CORRECTNESS_VALIDATION_PROTOCOL.md`
(FROZEN v4) was approved by the user on 2026-10-08 under a limited scope:
- it measures synthetic returned-answer accuracy for this data generator;
- it establishes **no** truth-free production reliability signal;
- the A3 comparison with Anderson is descriptive only.

What it fixes:
- **Primary metric:** projection-free structural error, which E0
  validated on seeds 1-30 with a no-solver re-score (21,594 of 21,594
  truth values exact).
- **Selection rule.**
- **Gates on held-out seeds 31-200:** A1 integrity; A2, the two-sided
  95% Clopper-Pearson lower limit (AIDS ≥ 161, QUAIDS ≥ 154 of 170);
  A4, at most 1 returned-wrong answer.

Supporting implementation:
- The analysis script `confirm_analyze.py` passes 24/24 tests with
  `ResourceWarning` as an error, and reproduces E0 exactly.
- `FREEZE_MANIFEST.json` hashes 34 files. Its SHA-256 is in
  `FREEZE_MANIFEST.sha256` and in the freeze commit message.
- The freeze is commit `ad986fc` (since pushed, with the E2 evidence
  commit `158313a`).
- At the time of the freeze, E2 had not been run. It was then approved
  and run on 2026-10-08/09 (see the update above). No shipped `src/`
  file changed.
- **Line-ending caveat (fails closed).** The current working tree
  verifies all 34 manifest files. A fresh-checkout simulation (a
  `git checkout-index` export) verifies 27 of 34, plus the manifest and
  its sidecar.
  - The 7 exceptions are `src/quaids.sdf`, `src/quaids.src`,
    `src/quaidselas.src`, `src/quaidsiv.src`, `src/quaidsslutzky.src`,
    `src/quaidsutil.src` and `tests/quaidsfixtures.src`.
  - Their blobs and current working copies are LF, but a fresh Windows
    checkout under the system's `core.autocrlf=true` rewrites them as
    CRLF.
  - The manifest then fails closed: verification and
    `confirm_analyze.py --mode confirm` refuse to run.
  - E2 ran from the verified LF working tree. Any re-run or re-analysis
    must use that tree, or a `core.autocrlf=false`/`input` checkout with
    the manifest re-verified first.
  - Root-level checkout behavior was deliberately left unchanged.

Limitations:
- 30 seeds; no population-wide claim.
- The continuation vs cold-LM differences are 1-4 seeds.
- The post-hoc pooling result is exploratory.

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
- **TVP-AIDS Stage 6 (complete, committed, and released)**: the
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
    completely missed. Committed as `4129201` and pushed to
    `origin/master`; the self-hosted push-triggered CI run completed
    `success` (confirmed via `gh run list`, run id `35878499421`, 1m28s).

- **AIDS/QUAIDS convergence-correctness research** (experiment-only,
  committed/pushed as `7e98720`): critically audited the prior Anderson
  acceleration and multi-start investigation, corrected its objective
  direction, implemented the recommended safeguarded structured
  multi-start benchmark, and retained complete scripts/output under
  `dev/experiments/convergence-anderson/`.
  - `homogCrit=-ln(det(Sigma))` must be maximized. The earlier pilot had
    minimized it; after correction, criterion selection matched the
    synthetic oracle's loose correct/wrong classification in 399 of 400
    seed/model cases in the full benchmark.
  - The safe experiment path uses true fixed-point-residual stopping,
    consecutive-difference Type-II history, a rank-truncated SVD solve,
    residual-growth restart, an accelerated-step safeguard, and a plain
    step/history reset when SVD is unusable. The SVD fallback was added
    after a reproducible QUAIDS seed-17 failure and then validated through
    that seed before the full rerun.
  - Full 200-seed x 16-start results: iterated AIDS improved from
    98 correct / 96 wrong / 6 never-converged with Stone-only safe
    Anderson to 121 / 76 / 3 with max-criterion multi-start. QUAIDS
    improved from 56 / 111 / 33 to 69 / 122 / 9. Thus multi-start
    improves reachability, but QUAIDS still converts many failures into
    converged-but-wrong answers. At the stricter `recErr <= .5` threshold,
    only 38/200 AIDS and 16/200 QUAIDS selections qualify.
  - The criterion selector is not the demonstrated bottleneck; reaching a
    good endpoint and restoring a real-data trust signal are. Substantial
    endpoint diversity was observed (71 AIDS and 119 QUAIDS seeds), so
    Anderson-induced basin collapse is not supported by this evidence.
  - Nothing here is production-ready or shipped. The next research step
    is direct nonlinear optimization plus held-out and local-rank/
    conditioning diagnostics; see `CLAUDE_HANDOFF.md` in the experiment
    directory for exact reproduction commands and sequencing.

## Decisions

- **TVP-AIDS D0 (2026-10-15): option A, sslib v2.1.0**
  (`7d0e0421a008db81fb4fda0b789bdbef21ffe5f7`), the first tag with the
  diffuse-phase termination fix. Install method: the shared install was
  already at v2.1.0 when found, installed outside `sync_sslib.ps1`. It
  was content-verified (all 9 `src/` files and `package.json` match the
  tag after line-ending normalisation) and was not reinstalled.
- **TVP optimizer caps default to CMLMT's own values (1e5 iterations,
  1e8 minutes), not finite ones.** The coordinator brief required
  backward-compatible defaults. Finite caps are opt-in and are shown in
  the docs and examples. Revisit if promotion wants bounded-by-default
  behavior (that would be a behavior change, needing a version bump and
  a changelog entry).
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
- **Historical pre-Stage-6 state:** `sslib` version pinning had been
  decided in principle ("pin to a specific commit"), but no mechanism
  existed because `package.json`'s `deps` array is a bare string list with
  no room for a commit hash. The previously documented `9132c35` pin was
  already known stale, not just
  theoretically at risk: this session found the actual installed copy at
  `C:\gauss26\pkgs\sslib` includes at least upstream commit `ae921ce`
  (confirmed by the `init_diffTVP` arity break/fix — see Completed Work's
  Stage 3 entry), and `gauss-state-space-ea`'s own concurrent session
  reported their repo (clean, matching origin) is at `7d5ed72`, LIKELY
  (not independently confirmed) close to what was then installed. This
  blocked a clean Stage 6 until the pinning mechanism and stable v1 pin in
  the next decision were implemented.
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

- Convergence research validation (2026-10-01):
  `anderson_safe_solver_test.e` passed its deliberately rank-deficient
  SVD solve (rank 2, fit difference from pseudoinverse `4.44e-16`) and a
  real QUAIDS seed-7 safe fit; `anderson_faithfulness_check.e` confirmed
  the experiment prototype at `mDepth=0` exactly reproduces shipped
  `quaidsFit()` in both `bS` and `vS`; the full benchmark output contains
  both 200-seed summaries and its explicit `run complete` marker.
- GitHub Actions for experiment commit `7e98720` completed successfully
  (`Test suite`, run `36855314457`, 2026-10-01). The experiment itself is
  not part of CI; this confirms the existing shipped suite was not
  regressed by adding the research artifacts.
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
- Full Stage 6 release verification
  (`scripts/run_release_verification.ps1 -BuildArtifact -ForceArtifact
  -InstallArtifact`) passed before the Stage 6 commit. The subsequent
  `v0.3.0` release gate also reported `GO`; artifact digest and release
  record are in Current Objective. No later work changed the installed
  package surface.
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

- **Resolved 2026-10-17 (sslib v2.1.1): sslib v2.1.0's exact diffuse
  filter was wrong for a non-diagonal H when the diffuse phase ended
  mid-period.** Found 2026-10-16 by P5 and fixed upstream in
  gauss-state-space `004c1be`. The filter now runs its diffuse test in
  LDL-decorrelation units and counts directions from `rank(P_inf_0)`.
  P5 case 3 is now a counted full-sample check (1.7e-13 per period).
  sslib v2.1.0 and earlier remain affected; see
  `tests/fixtures/tvp_loglik/README.md` ("Case 3 history").
- **Likely resolved by sslib v2.1.1, pending the P8 re-check: example 14's
  "covariance failed to invert" and sslib's fragile untransformed
  second CMLMT pass.**
  - On v2.1.0, `ssFitTVP()`'s second, untransformed CMLMT pass printed
    `CMLmt: function evaluated as complex` near a zero Q element, and
    example 14 printed "The covariance of the parameters failed to
    invert" (it still converged). In one deliberately misspecified scratch
    configuration it aborted with `G0518 Not implemented for complex` in
    `cmlmthess`, which `trap` cannot catch.
  - v2.1.1 (`717abb5`) removed that pass and converts the covariance
    by the delta method. On 2026-10-17, example 14 no longer prints either
    warning (CONVERGED, 23 iterations). The misspecified G0518
    configuration was not re-run, so P8 should confirm it.
  - The well-posed `quaidstvpfit_test.e` fixture does not trigger it.
- **RESOLVED 2026-10-15: the stale local `quaids 0.3.0.zip` no longer
  blocks the release gate.**
  - It was gitignored, in the repo root, and tripped the stale-artifact
    guard in `verify_release_artifact.ps1`.
  - Its SHA-256, `c5ce5e5ded68c64c6e5bdeaad9fdaf01bd82fca384d1aa7cef202f926d52b791`,
    matches `release_records/0.3.0-20260924-093017.json`.
  - It was moved, not deleted, to
    `C:\Users\eclow\Documents\quaids-artifact-backups\quaids 0.3.0__from-gauss-aids-model-root__moved-20261015-220303.zip`,
    with a `.README.txt` beside it.
- **The installed `quaids` package (`C:\gauss26\pkgs\quaids`) is still
  0.3.0.**
  - `examples/14_tvp_aids_estimation.e` now sets `tvpCtl.maxIters`/
    `maxTime`, which 0.3.0 lacks. A standalone
    `run_examples_smoke.ps1` run fails example 14 until 0.4.0 is
    installed.
  - The release gate (P8) installs 0.4.0 before its smoke stage, so
    this resolves itself there.
  - The example was verified against the 0.4.0 source through a scratch
    copy with absolute `#include`s.
- **`tests/quaidstvpfit_test.e` does not terminate in bounded time
  (2026-10-15: TIMEOUT at 600 s).**
  - Cause: the unbounded CMLmt defaults in `_quaidsTVPMLEFit()`, plus an
    ill-posed fixture.
  - The sslib drift below is verified to be non-causal: v1.0.0 and v2.0.0
    give bit-identical trajectories.
  - The shared install's `src/` is exactly `v2.0.0` (`0326d40`), with a
    stale 1.0.0 `package.json`. The "unverified drift" wording in the next
    item is now resolved to that.
  - Fix plan: `dev/TVP_PROMOTION_ROADMAP.md` (Priority 1).
- **`sslib` v2.0.0 released upstream (2026-09-29); the shared install at
  `C:\gauss26\pkgs\sslib` is currently a drifted/mixed tree, NOT this
  repo's pinned `v1.0.0`/`ad15626`.** Rechecked on 2026-10-01:
  `package.json` still reports `1.0.0`, but
  `scripts\verify_sslib_pin.ps1` flags `src/sstvp.src`/`sskalman.src`/
  `ssmain.src`/`ssstructural.src` as differing from the recorded v1
  hashes. A 2026-09-29 concurrent inspection identified the replacement
  source as coming from the v2 line, but the current mixed metadata means
  callers should describe it simply as unverified drift. `v1.0.0` itself
  is untouched upstream
  (`gauss-state-space`'s `v1.0.0` tag still peels to `ad15626`) — this
  repo's own pin is still valid, the INSTALLED copy is just currently
  something else. No action taken since nothing `sslib`-dependent is
  in flight; before any future TVP-AIDS/`sslib` work, re-run
  `scripts\sync_sslib.ps1 -SourceRepoPath <checkout> -Commit ad15626`
  (or deliberately move the pin to `v2.0.0`/`67bc506` — v2 requires GAUSS
  26.1.4+, which this repo's own dev environment already has, and per
  `gauss-state-space-af` every proc this repo calls by field name
  (`ssOut.final_params`/`.kfResults`/`.tvpFinal`/`.mleResults.retcode`)
  is unaffected by v2's only structural change, a new nested
  `ssOut.structural` field — worth a deliberate decision, not a default,
  since `sync_sslib.ps1`'s own hand-built `.lcg` generator does not emit
  the `typed_returns`/`keywords` annotations v2's own newer API surface
  needs; install v2 via `lib`/the GAUSS Package Manager instead if this
  repo ever moves up). Also good to know: `sslib`'s own test runner now
  defaults to an isolated install mode that never touches the shared
  directory, so the specific contention pattern that caused two earlier
  drift incidents this initiative should not recur.
- **`tsmt` package shadowing + the `library`/cross-file-global lazy-load
  gotcha** — both now documented as durable environment/language facts
  in `CLAUDE.md` (Development environment / Known GAUSS-26 language
  gotchas) rather than here, since they'll recur for any future session
  touching `sslib`. In short: any `tgauss` invocation using `tsmt` or
  `sslib` needs the `GAUSS26_CFG` override described there, and Stage 2
  code should call `ssControlCreate()` before any `sstvp.src` proc.
- `gauss-state-space`'s documented collision risk was observed repeatedly
  during TVP-AIDS development. Do not rely on an old statement about that
  repository's current working-tree state: re-check its `git status`
  before reading from it, and never install `sslib` from its live working
  tree — use `git archive` of a specific commit, as the pinning workflow
  does.
- **Historical context: the collision risk is not limited to that repo's
  own working tree — the SHARED INSTALLED COPY at
  `C:\gauss26\pkgs\sslib` itself changed mid-session** during the TVP
  work: `sstvp.src`/
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
  behaves unexpectedly in a future session. This chronology is retained
  to explain why the pin/verification tooling exists; the first Known
  Issues item is authoritative for the installed copy's current state.
- Iterated AIDS and QUAIDS (the estimator's own nonlinear iteration, not
  the TVP work) remain **experimental** per README's support-tier table
  — a known, documented estimator property, not a bug. The committed
  200-seed safeguarded multi-start benchmark improves reachability but
  does not resolve correctness: using the deliberately loose
  `recErr <= 1` synthetic cutoff, only 121/200 AIDS and 69/200 QUAIDS
  selections are correct; at `recErr <= .5`, only 38/200 and 16/200
  qualify. On real data there is no `recErr`, so a faster or more
  frequent `converged==1` result is not itself a trustworthy success
  signal. No Anderson or multi-start behavior has been shipped.

## Next Steps

**Priority 1 (2026-10-15): TVP-AIDS promotion.** P1-P7 are done, along
with the review follow-ups (roadmap section 10). P5 was closed on
2026-10-17 by Stage 1 (sslib v2.1.1), and P6 followed later that day
(Stage 2). Next comes P8. The stale-zip
blocker is resolved. Both upstream sslib risks that P5 tracked are
addressed in v2.1.1: the non-diagonal-H defect is fixed (fixed-Q likelihood checked),
and the second-pass complex/G0518 path is removed. P8 still owes the
support-tier judgment and a re-check of the G0518 configuration. The
original sequenced plan, P1-P8, in `dev/TVP_PROMOTION_ROADMAP.md`:
- P1: decision gate D0, choosing a tagged sslib baseline before any repin;
- P2: repin and resync;
- P3: bound the MLE (`maxIters`/`maxTime` on `quaidsTVPControl`);
- P4: a well-posed `quaidstvpfit_test.e`;
- P5: revalidate the numerics, plus an independent log-likelihood check;
- P6: re-measure the `n1` limit;
- P7: a per-test timeout in `run_source_tests.ps1`;
- P8: the local TVP release gate.

Item 5 below is superseded by that plan.

**Deferred and preserved, not abandoned:** the convergence-correctness
direction below. It covers AIDS direct-LM (including the frozen symmetry
study and its open owner decisions), QUAIDS direct LM / no-result seeds,
and speed benchmarking. Resume after TVP Priority 1, or by explicit owner
decision. The section below was originally written as "TVP-AIDS is
complete; the next recommended research direction is convergence
correctness for iterated AIDS/QUAIDS". That is a recommendation, not
active implementation or approved production scope.

1. Done: the direct objective and its LM solver were defined and
   validated (Phases 1, 1b, 3), and continuation was added (Phase 3 and
   its hardening).
2. Done (2026-10-05): the matched 30-seed pilot comparing safe Anderson,
   cold LM and continuation. See the Current Objective update above.

   The pilot's open items were then resolved:
   - The metric was fixed: projection-free scoring was validated by E0,
     and in-objective symmetry was waived for the protocol.
   - The model-specific hybrid was prespecified: AIDS `L16`; QUAIDS
     `L16+C16`.
   - Both are frozen in `CORRECTNESS_VALIDATION_PROTOCOL.md` (FROZEN v4,
     2026-10-08).

   Done (2026-10-08/09): E2, the 170-seed confirmation run. All binding
   gates passed (see the Current Objective update).

   **Recommended next work** (none of it started; each item needs its own
   approval):
   - Preserve the E2 evidence.
   - Prototype AIDS direct LM as an **opt-in** experimental mode, with
     the shipped default unchanged. **Done (2026-10-09) as a private
     prototype** (`dev/prototypes/aids_direct_lm/`); decisions and validation are
     in `dev/AIDS_DIRECT_LM_DESIGN.md`.
   - Investigate the 10 QUAIDS no-result seeds under a **new, separately
     prespecified exploratory protocol**.
   - Resolve the symmetry design before shipping QUAIDS (or AIDS) direct
     LM: in-objective symmetry, or a valid projection covariance.
     **AIDS symmetry study: COMPLETE (frozen 2026-10-13; stage C
     2026-10-13; G6 2026-10-15). C0-C5 and G6 all PASS; final verdict
     PROCEED** (scope and limitations at the end of this item). Protocol
     v1 (approved 2026-10-12 for the private implementation stage):
     `dev/prototypes/aids_direct_lm/SYMMETRY_PROTOCOL.md`. It covers
     in-objective symmetry against the homogeneity-only endpoint and an
     endpoint-weighted projection comparator, with fresh seed blocks
     1001-1030 (development) and 2001-2170 (confirmation). Thresholds
     approved (U6); the small-sample panel (U7) and real-data holdout
     (U8) are not run in v1. Still deferred: any shipped `quaidsFit()`
     change, the projection endpoint as a start, public exposure and new
     returned fields (U2-U5). **Implemented 2026-10-12 (private, under
     `dev/prototypes/aids_direct_lm/` only):** `quaidsdirectsym.src`
     (arm S plus the P comparator), the guarded harness
     `sym_harness.src`/`run_sym_stage.ps1`, `sym_analyze.py` and
     `sym_config.json`, and the checks (`quaidsdirectsym_test.e` 38/38,
     `test_sym_analyze.py` 11/11 including I6, five `directlm_sym_*.e`
     guards; `run_sym_checks.ps1` runs them all). A fixture-only pipeline
     smoke run (seed 602) was deterministic and agreed with the analyzer.
     **Independent pre-stage-D review (2026-10-13):** no estimator,
     math, dimension or selection defect found; guard/tooling fixes only
     (recorded in the protocol's Status): `-Full` implies `-WithEvidence`
     and clears any old stamp first; the analyzer never overwrites its
     outputs or the freeze manifest, needs a clean `dev_out` before
     freezing, rejects an incomplete determinism re-solve and maps
     asym/determinism failures to the stage D stop rule; `conf_rest` runs
     once; `truth.csv` is in key order. `test_sym_analyze.py` is now
     13/13. **I9 passed:** `run_sym_checks.ps1 -Full -WithEvidence` (all
     targeted checks, 10 guards, evidence test 80/80, and
     `run_source_tests.ps1 -SkipTVPKalman -SkipBootstrap`, which excludes
     every SSLib TVP test including `quaidstvpfit_test.e`) wrote
     `dev/prototypes/aids_direct_lm/symmetry/IMPLEMENTATION_CHECKS_PASSED.txt`
     (untracked). Any edit to a stamped file invalidates it.
     **Stage D run 2026-10-13 (seeds 1001-1030 plus the 1001/1030
     determinism re-solve; outputs kept, untracked, under
     `dev/prototypes/aids_direct_lm/symmetry/`: `dev/`, `dev_det/`,
     `dev_out/`):** the stamp re-verified first. The default invocation
     `run_sym_stage.ps1 -Stage dev` failed before any solve: under Windows
     PowerShell 5.1, `$PSScriptRoot` is empty inside the param-default of
     an advanced script (this one has `[Parameter(Mandatory)]`), so the
     `-OutRoot` default threw in `Join-Path`. No `dev`/`dev_det` directory
     or seed output existed afterwards. Rerun with the same directory
     passed explicitly (`-OutRoot <repo>\dev\prototypes\aids_direct_lm\symmetry`):
     `dev` 06:06:49-06:14:44 (7 min 55 s), `dev_det` 06:14:44-06:14:59.
     `sym_analyze.py analyze --stage dev --run symmetry/dev --det
     symmetry/dev_det --out symmetry/dev_out`: 0 integrity problems, 0
     determinism differences. **Stop rule (protocol section 8 step 3) not
     triggered:** S returned-wrong 0, S no-result 0, N2 > 0.01 on 0 seeds
     (max 3.9e-5), asym 0 exactly on every S/S1/P result.
     Descriptive (informational only; the dev-mode C1 "REJECT" is
     mechanical, 30 seeds can never reach 161): H, S, S1 and P each 30/30
     correct at cutoff 1 (max pfErr H 0.26, S 0.21, S1 0.21, P 0.21);
     N1 in [-0.0039, -0.0005]; LR in [1.6, 11.8], none above the chi2(6)
     5% point; every S start eligible; one S near-tie basin flag (seed
     1027, margin 2.8e-4); median seed time 13.7 s, projected stage C 39
     min. **Post-Stage-D wrapper repair (2026-10-13, numerically
     irrelevant; Stage D not rerun):** the `-OutRoot` default of
     `run_sym_stage.ps1` (and, as a precaution, `run_sym_checks.ps1`) now
     resolves after parameter binding; `run_sym_stage.ps1 -Preflight`
     prints the resolved paths and exits before any guard or solve, and
     `run_sym_checks.ps1` checks it from two working directories and with
     an explicit `-OutRoot` (no directory created, `symmetry/` tree
     unchanged). The old stamp was deleted, and
     `run_sym_checks.ps1 -Full -WithEvidence` (06:18:49-06:23:31) passed
     everything (I1-I9; evidence 80/80; I9 runner PASS without any SSLib
     TVP test) and wrote a new 25-file stamp. All 15 `dev/`, `dev_det/`,
     `dev_out/` files are byte-identical to before the repair; the
     execution logs of both Stage D attempts are in
     `symmetry/stageD_execution_logs/`. **Frozen 2026-10-13 06:27
     (protocol section 8 step 4):** `python -B sym_analyze.py manifest
     --gauss-version "26.1.4"` wrote `symmetry/SYM_FREEZE_MANIFEST.json`
     (25 files, identical to the stamp's hashes and `sym_config.json`'s
     `freeze_files`; git HEAD `d8b4318`; GAUSS 26.1.4, Python 3.12.0,
     NumPy 2.5.1, SciPy 1.18.0) and its sidecar
     `SYM_FREEZE_MANIFEST.sha256` (digest
     `b1036a8b04dd9a53baf8e3a20020908d6b1053413bd16bae52313ce646a0790a`).
     `verify-manifest` passes, the hashes were rechecked independently,
     and a regenerated copy is byte-identical. Do not edit any of the 25
     covered files from here on (a revision means a new protocol version
     and block 2201-2370). Only Stage D seeds (1001-1030) have been solved.
     **Stage C interim run 2026-10-13 (protocol section 8 step 5; integrity
     and runtime only):** after `verify-manifest` passed and with no
     `conf/` directory present, `run_sym_stage.ps1 -Stage conf_interim`
     (default `-OutRoot`) ran seeds 2001-2020 once, 06:29:04-06:31:46,
     into `symmetry/conf/` (complete grid: 640 cells, 21,120 raw, 80 seed
     and 900 truth rows). `sym_analyze.py analyze --stage interim --run
     symmetry/conf --out symmetry/conf_interim_out` reported 0 integrity
     problems and no accuracy or gate statistic (truth not loaded). Runtime:
     median 7.9 s per seed (max 12.8 s), projected stage C 1,339 s (about
     22 min) against the 6 h cap, so no feasibility stop. The manifest
     still verifies. **Stage C completed 2026-10-13:** `run_sym_stage.ps1
     -Stage conf_rest` (seeds 2021-2170 appended once, 06:33:06-06:51:57;
     the interim rows are byte-unchanged) and `-Stage conf_det` (seeds
     2001, 2085, 2170; 06:51:57-06:52:16). `conf/` is complete (5,440
     cells, 179,520 raw, 680 seed, 7,650 truth rows; 170 seeds, no
     duplicate keys). Total stage C wall time about 22 min against the 6 h
     cap. `sym_analyze.py analyze --stage confirm --run symmetry/conf --det
     symmetry/conf_det --manifest symmetry/SYM_FREEZE_MANIFEST.json --out
     symmetry/conf_out` (no `--real`): manifest OK, 0 integrity problems,
     0 determinism differences. **Gates:** C0 pass; C1 pass (k_S = 170/170,
     Clopper-Pearson lower 0.9785, threshold 161); C2 pass (W_S = 0); C3
     pass (0 seeds with N2 > 0.01, max N2 4.2e-4, none missing); C4 pass
     (c = 0, b = 0); C5 pass (0 regressions, so no mechanism is needed).
     G6 is pending, so the recorded verdict is `PENDING_G6`, not
     PROCEED. Descriptive: k_H = k_P = k_S1 = 170, W = 0 for every arm,
     no "none" result; median pfErr S 0.073 vs H 0.098 (paired median
     S - H -0.017); max pfErr below 0.24 for every arm; no borderline
     (0.25, 2] seed; no N1 flag; LR above the chi2(6) 5% point on 7.6% of
     seeds (preview only); 2 S near-tie flags; asym exactly 0 and
     adding-up/homogeneity errors <= 4.4e-16 on every S/S1/P result.
     **G6 completed 2026-10-15:** `verify-manifest` OK beforehand;
     `run_sym_stage.ps1 -Stage real` run once (exit 0, 06:04:00-06:04:05,
     output `symmetry/real/`); `sym_analyze.py real --run symmetry/real
     --out symmetry/real_out` wrote `symmetry/real_out/realdata.json`
     (H, S, P ok = 1; R1 = -0.667 and R2 = -0.021, both <= 0.01; R3 LR
     21.33, chi2(3) p 9.0e-5, non-inferential; R4 support 16, 0 exact or
     near ties, margin null; R5 asym, addUpErr, homogErr all 0 for S and
     P; R6 polish retcode 0 at iteration 1, critDiff 0; R7 condGN 7.0e5;
     R8 tables only). The final confirm analysis (`--real
     symmetry/real_out/realdata.json --out symmetry/conf_g6_out`; `conf_out`
     left unchanged) reported manifest OK, 0 integrity problems, 0
     determinism differences, C0-C5 pass with the same statistics as
     above, G6 pass, and **verdict PROCEED** (experimental, AIDS-only,
     point-estimate MVP evidence, subject to the protocol's section 12
     owner decisions). The optional small-sample panel was not run.
     **Artifact locations (all untracked; keep them):** the real-data
     report is `symmetry/real_out/realdata.json` (not `symmetry/real/`,
     which holds the harness CSVs). The final results are
     `symmetry/conf_g6_out/results.json` and `report.md`. In that
     `results.json`, `gates.G6.pass` is `true` but `gates.G6.statistic` is
     `null`. G6's conditions are in `realdata.json`. `symmetry/conf_out/`
     is the pre-G6 analysis (verdict `PENDING_G6`). It differs from
     `conf_g6_out` only in `gates.G6.pass` and `verdict`, and its
     `seeds.csv` is byte-identical. **Real-data descriptives (truth-free,
     never a correctness claim):** the Blanciforti86 LR is 21.33, with
     chi2(3) p about 9.0e-5. It is explicitly **non-inferential**: `u` is a
     generated regressor and no test size was validated. It is not
     evidence for or against symmetry, though it is large. S's
     Gauss-Newton condition number (condGN, R7) is about 7.0e5. G6 is
     executability only. An independent read-only audit (2026-10-15)
     re-verified the 25 manifest hashes and recomputed C0-C5 and G6 from
     the raw artifacts. It found no integrity or protocol discrepancy.
     **What PROCEED licenses:** S (in-objective symmetry) is the preferred
     symmetric path for an **experimental, private, AIDS-only,
     point-estimate** MVP. The evidence base is PROCEED, E2 and the
     I1-I9 checks kept as permanent tests. **Section 11 limitations any
     such MVP must state:** one synthetic DGP family (n = 5, T = 3000); no
     per-fit reliability signal; truth-free checks are falsification-only;
     AIDS only, unweighted, `homogenous = 1`, no user `b0`; point estimates
     only. **PROCEED does NOT license:**
     - public exposure (U2);
     - new returned fields (U3);
     - any change to shipped `quaidsFit()`, or its use as an arm (U4);
     - theta_P as a start (U5);
     - covariance, standard errors, inference or a symmetry test. Each
       needs its own prespecified study (protocol section 11 items 1-6:
       a two-step covariance, coverage, LR/Wald size and power,
       multi-start effects, heteroskedasticity/clustering/weights,
       downstream integration);
     - QUAIDS, weights, small-T or real-data correctness claims;
     - any speed claim.

     **Next owner decisions:**
     - whether to commit the untracked protocol, code and `symmetry/`
       evidence;
     - then either a private AIDS-only point-estimate MVP wrapper, or
       first prespecifying the section 11 covariance/coverage study.
   - Benchmark speed in a separate, controlled study.

   Full rationale: `E2_RESULTS.md`, "Recommended next work".
3. Add diagnostics usable without synthetic truth: fixed-point residual,
   multi-start endpoint agreement, objective gaps, Jacobian/Hessian rank
   and conditioning, theory restrictions, and training-vs-held-out fit.
   Productionization should require a defensible real-data trust signal.
4. Done in Phases 3 and 3-hardening, and confirmed in E2: continuation
   as a targeted follow-up, solving a homotopy from a linear-index stage
   to the full nonlinear objective, with parameter-block scaling. It is
   part of the confirmed QUAIDS primary configuration.
5. (Superseded 2026-10-15 by `dev/TVP_PROMOTION_ROADMAP.md` gate D0
   and step P2. Kept for history.) Before any future TVP-AIDS work, resolve the `sslib` installation
   deliberately: resync to pinned `v1.0.0`/`ad15626`, or explicitly test
   and adopt v2 using its supported installation/catalog path. The
   current shared install fails the pin verifier and must not be assumed
   reproducible.

`README.md` and `package.json` both correctly report package version
`0.3.0`; no stale package-version reference is currently known.

## Handoff Notes

- **Current state (2026-10-10, private prototype):** the freeze
  `ad986fc` and the E2 evidence commit `158313a` are both pushed
  (`origin/master` = `158313a`).

  On top of them, one local, unpushed commit adds the private AIDS
  direct-LM prototype: `dev/prototypes/aids_direct_lm/` (code and tests), and
  `dev/AIDS_DIRECT_LM_DESIGN.md`. No shipped estimator, public API,
  package manifest or version changed.

  E2 is complete, and every binding gate passes (`E2_RESULTS.md`). Next
  clean-session entry point: `dev/AIDS_DIRECT_LM_DESIGN.md` for the
  prototype; `dev/experiments/convergence-anderson/CLAUDE_HANDOFF.md`,
  then `E2_RESULTS.md`, for the research.
- **Previous state (2026-10-09, before the prototype; superseded by the
  entry above).** The text then said `ad986fc` and `158313a` were
  unpushed, which was already stale: `origin/master` was `158313a`.
- **Previous state (2026-10-08; superseded by the entry above).** Local
  `master` had the freeze commit `ad986fc` on top of `e90a568`, not
  pushed, and E2 had not yet been run.
- **Previous state (2026-10-06; superseded by the entry above).** `HEAD`,
  `master` and `origin/master` were at `e90a568`.
- **Historical (2026-10-01; superseded by the entry above).**
  Repository state before that day's status-only sync: `master` and
  `origin/master` both pointed to `7e98720`, the complete
  convergence-Anderson research commit; the only working-tree change was
  this pre-existing `PROJECT_STATUS.md` update, now reviewed and
  reconciled. The experiment commit contains 21 research code/output/
  handoff files, changes no shipped source, and has successful CI
  (`36855314457`).
- `gh run list` confirmed (as of 2026-10-01) that CI for `7e98720`, the
  `v0.3.0` release-status commit, both Stage 6 sync commits, and the
  earlier TVP commits completed with `success`. Earlier confirmed runs
  include `086c97a`, `8ed8ea3`, `a36c7a6`, `e3ef801`/`3d9f3d7`
  (Stage 4 + doc-sync), and `c58c511` (Stage 5). The push workflow runs
  with `-SkipBootstrap -SkipTVPKalman`, so CI
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
- `sslib`'s presence was reconfirmed on 2026-10-01, but its four tracked
  source files fail verification against the current `ad15626` pin even
  though the installed `package.json` still reports `1.0.0`. Presence and
  package metadata are therefore insufficient reassurance; run the pin
  verifier before relying on it and deliberately resync or repin first.
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
