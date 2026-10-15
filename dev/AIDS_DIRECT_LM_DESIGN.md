# Design decision: opt-in direct Levenberg-Marquardt estimation for iterated AIDS

- **Status:** ACCEPTED as a **private experimental prototype**, and
  implemented as one (2026-10-09). Owner decisions Q1-Q5 are recorded in
  §12. Nothing here is public API.
- **Date:** proposed and accepted 2026-10-09.
- **Evidence base:** `dev/experiments/convergence-anderson/E2_RESULTS.md`
  (E2 run under `CORRECTNESS_VALIDATION_PROTOCOL.md`, FROZEN v4, freeze
  commit `ad986fc`, evidence commit `158313a`).
- **Scope:** iterated **AIDS** only (`aCtl.linear = 1`,
  `aCtl.homogenous = 1`). QUAIDS is excluded (§11).
- **Implementation:** `dev/prototypes/aids_direct_lm/quaidsdirect.src` (deliberately outside
  `src/`; see §2). Tests, in the same directory:
  `quaidsdirect_test.e`, `quaidsdirect_evidence_test.e` and
  `guard_error_cases/directlm_*.e`.

This document does not reinterpret E2. Where a design choice goes beyond
what E2 tested, it says so and names the extra validation it needs.

## 1. What E2 supports, and what it does not

E2 supports this design work, and nothing more:

- **Primary result:** AIDS configuration `L16` (cold LM from the 16
  benchmark start recipes) was correct on 170/170 held-out seeds, with
  W = 0 (A2 and A4 PASS). `L1` (the Stone start alone) was also 0/0/170,
  but only descriptively; it was not a gated configuration.
- **Scoring:** "correct" is `pfErr ≤ 1` against the
  **homogeneity-constrained** optimum. Symmetry was not imposed or scored
  (E1 waived, protocol §12, risk R12).

E2 does **not** support the following. Code, output and documentation
must not imply any of them:

1. **Real data, other DGPs, other sizes.** Results come from one
   synthetic DGP (N = 5, tobs = 3000).
2. **A per-fit reliability signal.** None was validated (protocol §6,
   §11.4). A user cannot tell a returned wrong answer from a correct one.
3. **Better than the shipped default.** The E2 reference arm `A16` is a
   16-start experimental *Anderson* path. It is not the shipped
   single-start fixed point (`aCtl.relax = 1`). E2 made no comparison
   with shipped `quaidsFit()`.
4. **Symmetric estimates or any inference.** E2 scored point estimates
   only.
5. **Speed.** E2 timings are confounded (E2_RESULTS limitation 5).
6. **QUAIDS readiness.** 10/170 no-result seeds; see §11.

## 2. Entry point and visibility (Q1, Q3)

**A separate sibling estimator, not a mode inside `quaidsFit()`, and
private.**

```gauss
// working directory tests/; after quaids.sdf and quaids.src
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
struct _quaidsDirectOut dOut;
dOut = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
if dOut.ok; /* dOut.bFull, dOut.homogCrit, ... */ endif;
```

- **File:** `dev/prototypes/aids_direct_lm/quaidsdirect.src` holds the struct
  `_quaidsDirectOut`, the entry point `_quaidsDirectFit()`, the
  test-only seam `_quaidsDirectFitStarts()`, and the private solver
  procs.
- **Not shipped and not loaded by `library quaids;`.** The file lives
  under `dev/`, not `src/`. It is therefore not in `package.json`'s `src`
  array, not in `tests/verify_package_manifest.ps1`'s allowlist (which
  covers `src/` only), and not copied into a release artifact:
  `scripts/build_package.ps1` packages `src/`, `docs/`, `examples/`,
  `scripts/` and `tests/`, never `dev/`. The prototype's tests sit next
  to it for the same reason (`tests/` is packaged).
- **Not registered anywhere public.** Every name is underscore-prefixed.
  The struct is defined in `quaidsdirect.src`, not `src/quaids.sdf`.
  Nothing is added to `docs/public-api.json`, `docs/COMMAND_REFERENCE.md`
  or the command-reference pages. There is no version bump and no
  changelog release claim.
- **No printer.** The prototype returns a struct only. A printer is
  deferred to any future public design.
- **`quaidsControl` gains no field.** `quaidsControlCreate()` defaults
  are unchanged.
- **No `weight` argument.** E2 was unweighted.
- **No tuning knobs.** No start count, tolerance, iteration limit or step
  size is exposed. The constants are E2's (§4).

**Why not a `quaidsFit()` control flag** (the example wording in
E2_RESULTS's next-work item 2)?

1. **Downstream misuse.** A `quaidsOut` from a direct fit would flow into
   every procedure that accepts one: `quaidsElasFit()`,
   `quaidsWelfareFit()`, `quaidsCurvatureFit()`, the workflow, robust,
   survey and replicate layers. Those read `bestB`/`bestV`, `bS`/`vS`
   and `symValid`, which the prototype cannot fill validly (§7, §8). A
   separate output type makes that misuse a type error rather than a
   silent wrong standard error.
2. **Estimation-core risk.** `CLAUDE.md` prefers a sibling `.src` over
   editing `quaidsFit()`. Its variance block assumes the fixed point (§8),
   so a direct mode inside it would need branching through the variance,
   symmetry and recovery phases.

**Packaging (resolved 2026-10-10).** The first draft of the prototype
lived in `src/` and `tests/`. Because `build_package.ps1` copies those
directories whole, a release artifact would have shipped it as a
dormant file. It was moved to `dev/prototypes/aids_direct_lm/` instead
of adding an exclusion to the release tooling, which stays unchanged.

## 3. Default compatibility

- **`quaidsFit()`, `quaids()`, `quaidsControlCreate()`, `quaidsOut`,
  `src/quaids.sdf`, `package.json` and every existing procedure are
  untouched.** The shipped default behaviour and output are unchanged by
  construction. `quaidsdirect_test.e` section 5 also checks this
  directly.
- **The new code runs only when a caller `#include`s the file and calls
  `_quaidsDirectFit()`.**
- **`aCtl` field use:**

  | Field | Use |
  |---|---|
  | `linear` | Must be `1`. Otherwise fail (§5). |
  | `homogenous` | Must be `1`. Otherwise fail (§5). |
  | `alpha0` | Used, as in `quaidsFit()` (translog price index). |
  | `othnam` | Passed through to the internal Stone fit only. |
  | `b0` | Must be the default scalar `0`. A user start was never validated. Otherwise fail (§5). |
  | `maxiter`, `err`, `relax` | Not used: they control the fixed-point iteration. Ignored and not checked, so a control built for `quaidsFit()` can be reused as is. |

## 4. Algorithm and code reuse

**Objective.** Maximise `critExact = −ln det(Σ̂)` over the free
coefficients `vec(b[., 1:n1])`, where Σ̂ = R′R/N is built from the first
n−1 equations' residuals. The price index is built from the same `b`. This
is `_directObjCoreHomotopy` at homotopy 1, the objective E2 selected on,
and the same criterion as shipped `qOut.homogCrit`.

**Procedure.**

1. **Stone fit.** Call `quaidsFit()` with `aCtl.maxiter = 1` on the
   caller's data. Take `homogB` (the Stone start), `intcptFull`, `u` and
   `ng` from it. This is what the E2 harness does
   (`pilot_harness_core.src:124-140`).
2. **Scaling.** Build the scale vector (floor 0.05) from the Stone
   `homogB`.
3. **Starts (Q4).** Build the 16 benchmark recipes in recipe order, with
   recipe 1 = Stone. Noise comes from `rndKMn()` with a fixed internal
   seed (`20261009`) and an explicit state. The caller's global
   `rndseed`/`rndn` stream is never read or changed. The draws therefore
   differ from E2's, which were seeded from the synthetic DGP seed: E2
   supports the recipe *family*, not these exact draws. The test-only seam
   `_quaidsDirectFitStarts(..., startsIn)` accepts injected starts, which
   is how the evidence test feeds E2's own starts.
4. **Solve.** Run cold LM from each start with E2's constants:
   - `maxIter` 500;
   - `gradTol` 1e-3 (inf-norm of the scaled Gauss-Newton gradient);
   - finite-difference step 1e-6;
   - initial λ 1e-3;
   - acceptance ρ > 0.1, λ ×0.5 when ρ > 0.75, λ ×4 on rejection, stall at
     λ ≥ 1e12.
5. **Re-evaluate** the exact objective at each endpoint (as the harness
   did for every cell), **select** (§5) and **assemble** the output (§6).

**Code reuse.** Production cannot `#include` from `dev/`, and the frozen
experiment files stay byte-identical. So the AIDS path of these procs is
**copied** into `dev/prototypes/aids_direct_lm/quaidsdirect.src` under private names:

| Frozen proc | Prototype copy |
|---|---|
| `_directObjResidHomotopy` | `_quaidsDirectResid` |
| `_directObjCoreHomotopy` | `_quaidsDirectCore` |
| `_directObjScale`, `_directObjUnpack` | `_quaidsDirectScale`, `_quaidsDirectUnpack` |
| `_directResidJacobianFDHomotopy` | `_quaidsDirectJacobian` |
| `_directGNGradient`, `_directGNHessian` | `_quaidsDirectGNGradient`, `_quaidsDirectGNHessian` |
| `_directLMHomotopyDiag` | `_quaidsDirectLM` |
| `_pilotAddScaledNoise`, `_pilotBenchmarkStart` | `_quaidsDirectAddScaledNoise`, `_quaidsDirectStart` (RNG replaced, Q4) |

The only removals are:

- the homotopy argument, fixed at 1 (`1*.5 == .5` exactly);
- the QUAIDS branch;
- fault injection;
- the optmt / flat-penalty path;
- `verbose` printing.

Operation order, constants, the trapped fail-closed inversions and the
retcode meanings are unchanged. The evidence test confirms exact
equality (§10).

**Not reused:**
- continuation (QUAIDS only);
- analytic Jacobians;
- Anderson;
- `quaidsFit()`'s variance, symmetry and recovery blocks (§7, §8).

## 5. Selection and failure semantics

**Eligibility and selection** are protocol §4, verbatim:

- **Eligible:** retcode 0; the objective re-evaluated at the endpoint is
  feasible and finite; every endpoint coefficient is finite.
- **Selected:** the eligible start with the strictly largest `critExact`.
  There is no tolerance; exact ties go to the lower start index.
- **No silent fallback.** Ineligible starts are never selected. Nothing
  falls back to the fixed-point estimator or to any other method.

**No eligible start ("no result")** is a returned outcome, not a crash:

- `dOut.ok = 0` and `dOut.status = 1`;
- `selStart`, `bRaw`, `bFull`, `homogCrit`, `gnorm`, `iterations` and
  `marginOtherBasin` are missing scalars (`miss(0,0)`);
- the per-start table and endpoints are still filled.

This keeps E2's distinction between a visible failure (none) and a
returned answer.

**Bad input** is an immediate `errorlog` and `end`. Each case has a
`dev/prototypes/aids_direct_lm/guard_error_cases/directlm_*.e` script:

| Rejected input | Guard script |
|---|---|
| `linear ≠ 1` (the explicit QUAIDS rejection) | `directlm_quaids_rejected.e` |
| `homogenous ≠ 1` | `directlm_homogenous_zero.e` |
| `b0 ≠ 0` | `directlm_user_b0.e` |
| non-conformable shapes | `directlm_bad_shape.e` |
| missing values in the data | `directlm_missing_data.e` |

**Internal errors.** An inversion failure inside the solver is trapped
and classified (retcodes 4 and 5), exactly as in E2. It never aborts the
fit. A non-finite start ends with retcode 5 and is ineligible; that is
how the no-result path is tested.

## 6. Output and diagnostics

`_quaidsDirectOut` fields (documented in the struct definition):

- **Metadata:** `model = "AIDS direct LM (private prototype)"`,
  `experimental = 1`, `nobs`, `n`, `n1`, `ng`, `nIntFull`.
- **Outcome:** `ok`, `status` (0 = selected, 1 = no eligible start),
  `selStart`, `nEligible`, `nStarts` (16 unless starts are injected).
- **Estimates** (when `ok = 1`):
  - `bRaw`: the selected endpoint, ng x n1, in the reduced raw layout
    alpha | gamma | beta | u;
  - `bFull`: the (ng+1) x n recovery that `pfErr` uses. The n-th gamma
    row comes from homogeneity and the n-th column from adding-up. This is
    the layout of the synthetic truth, not `quaidsOut.b`'s
    absolute-price layout;
  - `homogCrit`: the selected `critExact`;
  - `gnorm` and `iterations`.
- **Per start:** `startTable` (nStarts x 13: start, retcode, eligible,
  `critExact`, `gnorm`, iterations, rejected steps, covInvFail,
  stepInvFail, infeasTrial, rhoReject, lmAccepted, final λ), and
  `startRaw` (every endpoint).
- **Context:** `bStone`, `scaleVec`, `intcptFull`, `u`. For the IV first
  stage, run `quaidsFit()` with `maxiter = 1`; the prototype does not copy
  it.
- **Multi-start diagnostics,** computed exactly as `confirm_analyze.py`
  computes them:
  - `exactTies`;
  - `marginOtherBasin`: the selected `critExact` minus the best
    `critExact` in another basin. Two endpoints are in different basins
    when the max |Δ raw| over rows 1..ng−1 (all rows except u) is above
    0.05;
  - `nearTie`: set when that margin is below 0.01.

**These diagnostics describe the multi-start run. None has been
validated as a signal that the returned estimate is correct** (protocol
§6: on development seeds every wrong direct selection was the only basin
found). The source header says so.

## 7. Symmetry treatment (Q2)

**Homogeneity-constrained output only.**

- `_quaidsDirectOut` has no `bS`, `vS` or `symValid` field.
- **Why:** E2 validated the homogeneity-constrained optimum and nothing
  else.
- **Shipped projection not applied.** The minimum-distance symmetry
  projection (`src/quaids.src:556-607`) moved alpha by 1-2 at direct
  endpoints in the pilot. Its weight matrix is also the fixed-point
  covariance, which the prototype lacks (§8).
- **Future options:**
  - **A, in-objective symmetry:** use the `R` selection-matrix idiom
    from `CLAUDE.md`. It needs its own prespecified validation (the E1
    smoke-gate design, then a held-out confirmation).
  - **B, a projection valid at direct endpoints:** this needs a valid
    covariance first.
- **Prespecified evaluation (v1, approved 2026-10-12; frozen
  2026-10-13; completed 2026-10-15):**
  `dev/prototypes/aids_direct_lm/SYMMETRY_PROTOCOL.md` evaluates A as the
  preferred MVP path against this homogeneity-only endpoint. Its only
  comparator for B is a minimum-distance projection weighted at the
  direct endpoint, used for point estimates only.
- **Outcome:** C0-C5 and G6 all PASS; **verdict PROCEED**. On seeds
  2001-2170, S was correct on 170/170 seeds (CP lower bound 0.9785), with
  0 returned-wrong, 0 N2 flags, 0 regressions and asym exactly 0.
  - **Real data (Blanciforti86, truth-free):** LR 21.33 (chi2(3) p about
    9.0e-5, non-inferential) and condGN about 7.0e5.
  - **Artifacts** (untracked, under `dev/prototypes/aids_direct_lm/symmetry/`):
    - final results are in `conf_g6_out/` (`results.json`: `gates.G6.pass`
      true, `gates.G6.statistic` null);
    - G6 details are in `real_out/realdata.json`.

  Details: `PROJECT_STATUS.md`.
- **What PROCEED licenses:** A, as the preferred path for an experimental,
  private, AIDS-only, point-estimate MVP. Such an MVP must state:
  - the evidence comes from one synthetic DGP family (n = 5, T = 3000);
  - there is no per-fit reliability signal;
  - the truth-free checks are falsification-only;
  - it is unweighted, with `homogenous = 1` and no user `b0`;
  - it gives point estimates only.

  It does **not** license public exposure, new returned fields, a change
  to shipped `quaidsFit()`, theta_P as a start (U2-U5), any covariance or
  inference (§8; protocol §11), QUAIDS, or a speed claim.

## 8. Covariance and inference boundary (Q5)

**Point estimates only.** No covariance, standard errors, t-statistics,
p-values, symmetry test or overidentification test.

Reasons:

1. **The shipped formula assumes the fixed point.** It uses
   `sse = ww − b′gw` and the matrix `J` linearised at the fixed point
   (`src/quaids.src:372-413`). The `sse` line is a normal-equation
   identity, which holds only when `b = solpd(gw, gg)`. A direct-LM
   optimum is generally not a fixed point; that is the gap the method
   closes.
2. **E2 validated no inference.**
3. **A naive alternative ignores first-stage uncertainty.** The inverse
   Gauss-Newton information treats the IV control-function residual `u`
   as known.

Any future covariance needs its own coverage simulation first.

## 9. Interaction with existing procedures

None. Every downstream procedure still accepts only `quaidsOut`, and the
prototype's result is a different type.

## 10. Validation (as implemented, 2026-10-09)

Every test uses the repository harness: one PASS/FAIL line per check,
then an `ALL N CHECKS PASSED` line. `tests/run_source_tests.ps1` runs
them all:
- `quaidsdirect_test.e` and the five guard scripts always run;
- `quaidsdirect_evidence_test.e` needs optmt and can be skipped with
  `-SkipDirectEvidence`;
- `dev/` is not shipped, so in an extracted release artifact (no
  `dev/prototypes/aids_direct_lm/`, no `.git`) the runner prints a `SKIP`
  line and skips all seven prototype checks; in a git checkout a missing
  prototype directory fails the run instead.

| ID | Check | Where | Result |
|---|---|---|---|
| T1 | **Exact arithmetic equivalence.** Frozen procs vs prototype copies on E2 seed 31 data, at 4 E2 starts: scale, unpack, residual, objective, FD Jacobian, GN gradient and Hessian, and full LM solves (all 7 outputs). Deterministic recipes 1, 11, 13, 16 are identical; the noise recipes perturb the same rows. | evidence test §1 | Zero difference. |
| T2 | **Read-only E2 reproduction.** E2's own 16 starts are injected for seeds 31, 115, 200 (determinism seeds) and 47, 90, 104, 123 (near-tie seeds). Every stored cold-LM cell is compared: retcode, eligibility, `critExact`, iterations, rejected steps, `gnorm`, the 6 LM counters, and the raw endpoint. So is the `L16` selection: start, eligible count, criterion, ties, margin, near-tie flag and `pfErr`. | evidence test §2 | Zero difference on all 7 seeds. 80/80 checks. |
| T3 | **Exact recovery.** A production-only fixture whose truth is an *exactly* stationary point: equation noise projected off every equation block of the Jacobian at the truth. The frozen `_quaidsStationaryFixtureDGP` was not used: it is only near-stationary (gradient 0.08 in `phase1b_audit_output.txt`). | core test §1 | LM started at the truth stops at iteration 1 and returns the truth exactly. The production starts recover it to 7.4e-11 (check < 1e-6). |
| T4 | **Noisy synthetic recovery.** `_quaidsSyntheticDGP` seeds 501 and 502, outside the research seeds. | core test §2 | `pfErr` 0.060 and 0.199, both ≤ 1. A descriptive `quaidsFit()` comparison is printed, not gated. |
| T5 | **Internal consistency.** The criterion is recomputed independently in absolute-price form from `bFull`; adding-up and homogeneity of `bFull` are checked; the selected start is eligible with the largest criterion. | core test §2 | Within 1e-9 and 1e-12. |
| T6 | **Selection and failure.** `_quaidsDirectSelect` unit cases: ineligible never selected, tie-break and count, none, the basin rule ignoring u, the near-tie flag. End to end: all-non-finite injected starts give `ok = 0` with missing estimates; [failed, Stone, Stone] selects start 2 with one exact tie. | core test §3 | Pass. |
| T7 | **Determinism and RNG isolation.** The caller's `rndseed` stream is unchanged, repeated fits are identical, and the internal starts are reproducible. | core test §4 | Pass. |
| T8 | **Guards.** The five `directlm_*.e` scripts. | runner | Each emits its specific message. |
| T9 | **Default compatibility.** `quaidsControlCreate()` defaults and `quaidsFit()` output before and after a direct fit, plus the full existing suite and the manifest, API and docs checks. | core test §5, runner | Pass, except that the unrelated sslib-dependent `quaidstvpfit_test.e` is unverified: it hangs because the shared sslib install has drifted from its pin (see `PROJECT_STATUS.md`). |

**Not done, by decision or scope:**
- **Installed-package test (old T10):** the prototype is not installed
  or exported.
- **Real-data smoke on Blanciforti86 (old T11):** descriptive only, and
  outside the validated envelope; deferred to any public design.

## 11. Explicit exclusions

- **QUAIDS is excluded.** `linear = 0` fails with an error. Reasons:
  - 10/170 no-result seeds;
  - it needs continuation as well as cold LM;
  - the primary configuration took a median of 79.5 s and a maximum of
    1,134 s per seed in E2.

  Any QUAIDS path waits for the separately prespecified exploratory study
  (E2_RESULTS next-work item 3) and its own confirmation. E2 must not be
  re-analysed to justify it.
- **No performance claims.** No documentation, changelog, example or
  output text may say direct LM is faster, slower or cheaper than any
  other method. E2's wall times stay context only.
- **Unsupported settings:** sampling weights, `homogenous = 0`, LA-AIDS,
  and a user `b0`.
- **Excluded features:**
  - downstream integration (§9);
  - symmetric estimates and any inference (§7, §8);
  - public API exposure;
  - a version bump;
  - a release claim.
- **No change to the default estimator.** Docs issue D1 (protocol §12) is
  also left to a separate docs-only task.

## 12. Owner decisions (accepted 2026-10-09 for the private prototype)

| ID | Decision | Accepted |
|---|---|---|
| **Q1** | Entry point | A separate sibling entry point, not a `quaidsFit()` control flag. |
| **Q2** | Symmetry | Homogeneity-only. |
| **Q3** | Exposure | Private experimental prototype only: no public API registration, package export, version bump, changelog release claim or public-release claim. |
| **Q4** | Starts | An isolated, fixed internal RNG stream for the starts, plus a test-only injected-start seam for E2 equivalence. |
| **Q5** | Inference | No covariance, standard errors or inference. |

**Still open:**
- any move to public API, which would need a new decision covering the
  printer, `docs/public-api.json`, command-reference pages, version bump,
  changelog and an installed-package test;
- symmetry option A or B as a shipped or public feature (the private
  evaluation protocol `dev/prototypes/aids_direct_lm/SYMMETRY_PROTOCOL.md`
  v1 is complete with verdict PROCEED for a private, AIDS-only,
  point-estimate MVP only; its U2-U5 stay open);
- a covariance design.

## 13. Documentation

For the private prototype:
- this document;
- the `dev/prototypes/aids_direct_lm/quaidsdirect.src` header and struct comments;
- the `tests/run_source_tests.ps1` comments;
- `PROJECT_STATUS.md`.

User-facing documentation (`docs/`, `README.md`) is deliberately
unchanged, because the prototype is not public. `CHANGELOG.md` has an
`Unreleased` maintainer note (the repository records every change there)
that names the prototype as private and unshipped and makes no release
claim.
