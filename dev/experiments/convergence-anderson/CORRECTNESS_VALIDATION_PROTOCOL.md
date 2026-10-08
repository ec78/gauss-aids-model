# Correctness-validation protocol: direct LM vs Anderson (FROZEN v4, 2026-10-08)

## Status

- **Version history:**
  - v2 drafted 2026-10-07;
  - E0, the no-solver re-score, run the same day;
  - v3 written after an adversarial pre-freeze review;
  - **v4**, 2026-10-07, after the user's scope decisions and the
    analysis implementation;
  - **FROZEN** on 2026-10-08 (America/Chicago) by the user's explicit
    approval, under the limited scope stated below. No criterion was
    changed at freeze; only the status lines, R15, §14.5 and the §15
    record were updated.
- **Inputs:**
  - the committed pilot (`e90a568`);
  - `CODEX_PILOT_REVIEW.md` (technical PASS; finding F1 since
    corrected);
  - `ORCA_CHECKPOINT.md`;
  - the E0 outputs in `e0_out/`.
- **Scope of this work:** no solver call or production-code change was
  made for v2, E0, v3, v4 or the freeze. E2 (the 170-seed confirmation
  run) has **not** been run.
- **Frozen:** this protocol is binding. The §15 checklist is ticked to
  record the user's 2026-10-08 approval and the evidence verified at
  freeze. **The frozen commit is the enclosing git commit that adds
  this text.** Its hash is reported after creation (coordinator report
  and `worker_done`). It is deliberately not written into this file:
  doing so would change the frozen content. The frozen manifest's
  SHA-256 is recorded in the sidecar `FREEZE_MANIFEST.sha256` and in
  that commit's message (§14.5).

**User decisions recorded on 2026-10-07 (reproduced, not ticked):**
- **Scope accepted.** The protocol evaluates synthetic returned-answer
  accuracy for this data generator. It does not establish a truth-free
  production reliability signal.
- **A3 is no longer a confirmatory gate.** The held-out Anderson results
  were known when it was designed, so it moves to a descriptive,
  secondary comparison.
- **Confirmation gates:** A2, A4 and the other genuinely prospective
  criteria.

What changed from v3:

| v4 change | Section |
|---|---|
| A3 demoted to descriptive and secondary. The gates are now A1, A2 and A4. A5 becomes required reporting, because its regressions are defined against the non-blind reference arm. | §7, §8 |
| Reconciliation of the AIDS development comparison: 8/0 earlier, 9/0 now. | §11.6 |
| Analysis script, configuration, tests, development-mode reproduction, E2 entry scripts and the freeze manifest are added. | §14 |
| Risks updated: R2 and R10 resolved, new items added. | §13 |
| Freeze checklist and recommendation updated. | §15 |

## 1. Decision this protocol supports

**Question.** On the synthetic DGP used throughout this research, does
a prespecified, model-specific direct-LM configuration return the
correct structural optimum often enough, and return *wrong* answers
rarely enough, to justify production-design work?

Passing this protocol does **not** authorize a change to `src/`. It also
does not show that a user can tell, from a real-data fit, whether that
fit is correct (§6).

## 2. What is validated and what is not

| Claim | Status | Evidence |
|---|---|---|
| The pilot harness, metrics, start equivalence and selector direction are correct on seeds 1-30. | **Validated** | Stage A gate 27 + 58 + 6 checks; 377/377 on the 2,880-cell grid; `CODEX_PILOT_REVIEW.md` PASS |
| The regenerated truth for seeds 1-30 is the truth the pilot used. | **Validated** | 21,594 of 21,594 stored values reproduced exactly (§11.2) |
| Projection-free scoring (`pfErr`) is computed correctly. | **Validated** | Equal to the pilot's unprojected `bsErrExclU` on all 195 symmetry-fallback cells (≤ 1.1e-13); all 6 prespecified predictions held (§11.3) |
| `pfErr` removes the projection artifact without reclassifying any converged direct-method cell. | **Validated on seeds 1-30** | 0 class changes among 1,815 converged direct cells; QUAIDS `L16+C16` 30/30 on `pfErr` vs 24/30 on projected `bS` |
| The selection rule (§4) is executed correctly by the confirmation script. | **Validated on seeds 1-30** | `confirm_analyze.py` reproduces all 2,880 cell-level `pfErr` values exactly, every E0 class table and W count, and every selected (method, start); 24/24 tests pass (§14) |
| Direct LM returns correct answers more often than Anderson **on held-out data**. | **Not validated** | Development only. Configurations were chosen on these seeds (winner's-curse risk). |
| The A2 and A4 gates pass on held-out seeds. | **Not validated** | E2 has not run. |
| Any truth-free signal identifies wrong answers. | **Not validated; evidence against** | §6, §11.4 |
| Results transfer to other DGPs, sample sizes or real data. | **Not validated** | Single DGP family, N = 5 goods, tobs = 3000 |

## 3. Arms and configurations

E2 uses the existing harness and settings unchanged (`PILOT_PROTOCOL.md`).
For every seed it runs Anderson, cold LM and continuation, each from all
16 benchmark starts. Every configuration is a deterministic subset of
those saved cells. The configurations are defined in
`confirm_config.json`.

| Configuration | Cells pooled | Role |
|---|---|---|
| `A16` | Anderson, starts 1-16 | Reference (descriptive comparisons only) |
| **AIDS primary** `L16` | cold LM, starts 1-16 | Gated (A1, A2, A4) |
| **QUAIDS primary** `L16+C16` | cold LM 1-16 ∪ continuation 1-16 | Gated (A1, A2, A4) |
| QUAIDS advisory `L16+C1` | cold LM 1-16 ∪ continuation start 1 | Same statistics reported; advisory only |
| Others (`L1`, `C1`, `L1+C1`, `C16`; QUAIDS `L16`; AIDS `L16+C16`) | as named | Reported, not gated |

Anderson cells are never pooled into a direct-method configuration.

## 4. Per-seed selection algorithm (complete specification)

The algorithm is implemented in `confirm_analyze.py` as `is_eligible`
and `select`. Apply it to one (model, seed, configuration):

1. **Candidate pool.** Take every cell of the configuration's methods ×
   starts.
2. **Eligibility.** A candidate is eligible only if all of the following
   hold:
   1. *Optimizer convergence* as recorded by the harness (`conv`):
      - **cold LM:** retcode 0 (feasible, inf-norm Gauss-Newton gradient
        < 1e-3);
      - **continuation:** homotopy reached exactly 1, the stage-zero
        solve had retcode 0, and every accepted stage had retcode 0 and
        was feasible on re-evaluation;
      - **Anderson:** `qOut.converged = 1`.
   2. A finite `critExact`.
   3. Every raw endpoint coefficient is finite.

   Ineligible candidates are never selected, whatever their objective.
   This covers continuation runs below homotopy 1, cold-LM infeasible
   starts (retcode 3), and non-finite (diverged) endpoints.
3. **Objective.** `critExact` = −ln det(Σ̂), where Σ̂ = R′R/N is built
   from the first n−1 equations' residuals at **homotopy 1**. It is
   computed by `_directObjCoreHomotopy` at the returned endpoint with
   the same data, `u` and scale for every method. **Larger is better.**
4. **Selection.** The selected endpoint is the eligible candidate with
   the strictly largest `critExact`. There is **no tolerance**. Exact
   ties go to the earliest candidate in (method ascending: cold LM
   before continuation; then start ascending) order.
5. **No eligible candidate.** The result is **"no result"** (class
   `none`).
6. **Diagnostics only** (never used to select or to decide pass/fail):
   - **Basin label:** two endpoints share a basin if the max |Δ raw| over
     all rows except u is ≤ 0.05.
   - **Near-tie flag:** set when another basin's best `critExact` is
     within ε = 0.01 of the selected value. The value was chosen because
     same-basin spread on seeds 1-30 reached 5.9e-3.

## 5. Metrics and per-seed classes

- **Primary metric, `pfErr`.** It is computed on the selected endpoint:
  1. Recover the full n-good matrix from raw columns `1:n1`:
     - by homogeneity, the n-th price row is minus the sum of the gamma
       rows;
     - by adding-up, the n-th column is minus the row sums, plus 1 on
       the constant.
  2. Take the max |error| over every row except u. Lambda is included
     for QUAIDS.

  The intercept-row count is derived from the stored layout
  (`ng − n1 − 2 − model`) and asserted (`n_int_full`).
- **Secondary metric:** `rawStructErr`.
- **Context only:** the projected `bsErr`.
- **Classes** at cutoff c (primary c = 1):
  - `correct`: an endpoint was returned and `pfErr ≤ c`;
  - `wrong`: an endpoint was returned and `pfErr > c`;
  - `none`: no result.
- **Sensitivity cutoffs:** 0.25, 0.5, 2 and 5.
- **Borderline listing:** every `pfErr` in (0.25, 2].

## 6. Optimizer convergence is not a reliability signal

- **Convergence is an eligibility rule** (§4 step 2). It means the
  method's stopping rule was met at a feasible point, and says nothing
  about which basin was reached.
- **On seeds 1-30, all wrong selections were converged**: 23 Anderson
  and every wrong direct-method selection.
- **No reliability signal worked.** E0 tested four truth-free candidates
  and none met the rule (§11.4). Every wrong direct-method selection was
  the only basin its configuration found.
- **Consequence.** The protocol measures how often the selector
  *returns* correct or wrong answers on the synthetic DGP. It validates
  no per-fit diagnostic. A real-data user has no validated way to tell a
  returned wrong answer from a correct one. The user has accepted this
  scope.

## 7. Criteria (per model; held-out seeds 31-200, n = 170)

### 7.1 Confirmation gates (all must pass for a model's primary configuration)

The model verdict is A1 ∧ A2 ∧ A4 (an intersection-union decision), so
no multiplicity adjustment is applied. Each model is decided separately.
"Correct" means `pfErr ≤ 1`, and `none` counts as not correct.

**A1 Integrity.** Four components must all pass:
1. `python pilot_analyze.py --validate e2_out`:
   - unique, complete keys;
   - counters present;
   - `critExact` consistency;
   - Anderson reproduces the benchmark's per-seed line on all 340
     held-out cases.
2. `confirm_analyze.py --mode confirm`, analysis-level integrity:
   - the manifest verifies;
   - the complete grid is present;
   - regenerated truth reproduces stored `rawStructErr`, `bsErr` and
     `critTrue` exactly;
   - the layout assertion holds;
   - no selected endpoint has a non-finite `pfErr`.
3. **Determinism:** `e2_determinism.e` re-solves seeds 31, 115 and 200
   for both models, and `pilot_analyze.py --compare e2_det_out e2_out`
   reports 0 differences apart from wall time.
4. **Resume:** `pilot_resume_test.ps1` passes.

**A2 Accuracy.**
- **Statistic:** k = the number of correct seeds, out of n = 170.
- **Interval method:** the **exact Clopper-Pearson two-sided 95%
  interval**. The criterion applies to its lower limit, which equals a
  one-sided 97.5% bound.
- **Pass thresholds:**

  | Model | Lower limit must be | Required count |
  |---|---|---|
  | AIDS | ≥ 0.90 | **k ≥ 161** (161 gives 0.9019; 160 gives 0.8945) |
  | QUAIDS | ≥ 0.85 | **k ≥ 154** (154 gives 0.8517; 153 gives 0.8447) |

**A4 Returned-wrong answers.**
- **Statistic:** W = the number of seeds where the selector returned an
  endpoint and it is wrong. "No result" seeds are not counted.
- **Pass threshold:** **W ≤ 1.** The exact two-sided 95% Clopper-Pearson
  upper limit is then ≤ 3.23%.
- **How it differs from A2.** A2 counts every failure: k = 170 − W −
  (number of `none` seeds). A4 caps only the failure a user cannot see.
  A2 alone would accept up to 9 (AIDS) or 16 (QUAIDS) failures of any
  kind, including all of them being returned-wrong. Together, A2 and A4
  force nearly all failures to be visible "no result" outcomes.
- **What it does not give.** A4 bounds a population rate under this DGP.
  It gives no per-fit protection (§6).

**Why these are prospective.** No direct-method result for seeds 31-200
exists or has been seen (§9). The gates refer only to the direct
configurations' own outcomes on those seeds.

### 7.2 Required reporting (not a gate)

**A5 Regression audit.** Every seed where `A16` is correct and the
primary configuration is not gets a written mechanism, from stored
cells and without new solves. It is not a gate because "regression" is
defined against the reference arm, whose held-out outcomes were known
when this protocol was written.

### 7.3 Descriptive and secondary comparisons (non-confirmatory; reported, never pass/fail)

**A3 Paired comparison with `A16` (DESCRIPTIVE).**
- **Definitions:** per seed, X = primary correct and Y = `A16` correct.
  - b = #(X = 1, Y = 0), c = #(X = 0, Y = 1);
  - ties = the concordant seeds.
- **Reported:** b, c, ties, the improved and regressed seed lists, and
  the exact one-sided conditional McNemar p-value P(B ≥ b | b + c, ½).
- **Not used for any decision.** The v3 thresholds (p < 0.01,
  b − c ≥ 17, c ≤ 3) are withdrawn as gates. They may be quoted only as
  "what v3 would have required".
- **Why it is non-confirmatory.** Anderson's held-out aggregate rates
  were known when A3 was designed (§9).

**A6 Robustness.** Report the primary configuration's correct counts
under `rawStructErr`, and at `pfErr` cutoffs 0.5 and 2. A divergence
from the gate verdict blocks any robustness *claim*.

**A7 Cost.** Report the median, p90 and maximum wall time per seed, and
list any cell over 900 s. Targets: median ≤ 30 s (AIDS) and ≤ 180 s
(QUAIDS).

**Advisory configuration** `L16+C1` (QUAIDS): the same statistics are
reported. It may be recommended over the primary only if it would pass
A2 and A4.

## 8. Decision mapping

| Outcome | Decision |
|---|---|
| A1 fails | The run is invalid; fix and rerun the affected cells. No conclusion is drawn. |
| A2 fails | Stop the direct-LM line for that model's configuration. |
| A2 passes, A4 fails | Returned-wrong answers are too frequent. Do search-breadth or diagnostic work before any production design. |
| A1, A2 and A4 pass | Recommend production-design work for that model's configuration, with separate approval. State §2's limits. Report A5 (required), A3, A6 and A7 (descriptive). |

There is no efficacy peeking. One interim check after the first 20
confirmation seeds may stop the run for A1 integrity reasons only.

## 9. Held-out information audit (seeds 31-200)

This audit was rechecked for v4 by read-only inspection.

**Direct methods: no information exists.**
- Every direct-method and E0 file holds seeds 1-30 only: `pilot_out/`
  1-30, `pilot_validation/` 1, 6 and 9, `e0_out/` 1-30.
- `e2_out/` and `e2_det_out/` do not exist.
- `confirm_analyze.py` was developed and tested in **dev mode**. That
  mode aborts on any row whose seed lies outside 1-30, *before* parsing
  the row's values (a test covers this).
- No LM or continuation solve has run on seeds 31-200.

**Anderson (reference arm): information exists and has been seen.**
- The 200-seed benchmark output (commit `7e98720`, which predates all
  direct-LM work) contains Anderson per-seed lines for all 340 held-out
  cases. Its aggregate rates (AIDS 60.5%, QUAIDS 34.5% on projected
  `bS`) are quoted in `SUMMARY.md`.
- Anderson's held-out `pfErr` is not known.
- Its list of held-out seeds with no converged start (AIDS 120; QUAIDS
  31, 43, 125, 134, 142, 151, 173, 182) is on disk.

**What was influenced, and what was not:**
- **Not influenced by seeds 31-200:**
  - the configurations;
  - the A2 targets and interval method;
  - A4;
  - the §4 rules, `pfErr`, the script and its tests;
  - the determinism seeds (31, 115, 200, chosen as endpoints and
    midpoint).
- **Influenced and therefore demoted:** A3's design, because Anderson's
  aggregate was known. A5's definition also depends on the reference
  arm, so it is reporting-only.

## 10. Experiments

- **E0 (done):** the no-solver re-score; results in §11.
- **E1 (recommended waived):** the symmetry-constrained objective smoke
  gate (§12). Its design is unchanged from v2.
- **E2:** the confirmation run. The run order is in §14.4, and none of
  it has been executed.

## 11. E0 results (development seeds 1-30 only)

### 11.1 Commands

```powershell
C:\gauss26\tgauss.exe -b -x e0_truth_dump.e        # > e0_out\e0_truth_dump_console.txt
$env:PYTHONIOENCODING='utf-8'; python e0_rescore.py  # > e0_out\e0_rescore_console.txt
```

No LM, continuation or Anderson routine was called.

### 11.2 Truth gate: PASS, no mismatch

All **21,594** stored values were reproduced exactly:

| Stored column | Cells |
|---|---|
| `rawStructErr`, `rawUErr`, `rawOverallErr`, `critTrue` | all 2,880 |
| `rawLambdaErr` | the 1,440 QUAIDS cells |
| `bsErr`, `bsErrExclU` | the 2,877 cells that have a `bS` |

The `pfErr` formula check passed (≤ 1.1e-13 on 195 fallback cells).

Two report-code defects (a console-encoding crash and a Stone-infeasible
summary line) were fixed before results were recorded; neither touched
any stored value. The Stone-infeasible seeds are QUAIDS 17, 18, 21 and
24, identical to the cold-LM retcode-3 Stone cells.

### 11.3 Re-scored classes (selected; none / wrong / correct at cutoff 1)

| Model | Config | `pfErr` | raw | `bS` (pilot primary) |
|---|---|---|---|---|
| AIDS | A16 | 2/7/21 | 2/6/22 | 2/9/19 |
| AIDS | L1 | 0/0/30 | 0/0/30 | 0/4/26 |
| AIDS | L16 | 0/0/30 | 0/0/30 | 0/3/27 |
| AIDS | L16+C16 | 0/0/30 | 0/0/30 | 0/3/27 |
| QUAIDS | A16 | 1/16/13 | 1/16/13 | 1/19/10 |
| QUAIDS | L1 | 4/3/23 | 4/3/23 | 4/9/17 |
| QUAIDS | C1 | 0/3/27 | 0/3/27 | 0/9/21 |
| QUAIDS | L1+C1 | 0/1/29 | 0/1/29 | 0/7/23 |
| QUAIDS | L16 | 0/1/29 | 0/1/29 | 0/7/23 |
| QUAIDS | C16 | 0/1/29 | 0/1/29 | 0/7/23 |
| QUAIDS | L16+C1 | 0/1/29 | 0/1/29 | 0/7/23 |
| QUAIDS | L16+C16 | **0/0/30** | 0/0/30 | 0/6/24 |

All 6 prespecified predictions held. There is no reclassification among
the 1,815 converged direct-method cells, and every correct direct
selection has `pfErr` ≤ 0.25.

### 11.4 Reliability-signal calibration: no candidate meets the rule

Worst cases across configurations:

| Candidate | Max wrong answers trusted | Min coverage |
|---|---:|---:|
| Support ≥ 2 | 13 | 85% |
| Basin reached by ≥ 2 methods | 1 | 0% |
| Objective ≥ Stone objective − δ | 12 | 86% |
| Objective margin Δ, with the Stone check | 7-12 | 86% |

Returned-wrong answers W on seeds 1-30:

| Config | AIDS | QUAIDS |
|---|---:|---:|
| A16 | 7 | 16 |
| L16 | 0 | 1 (seed 17) |
| C16 | 0 | 1 (seed 5) |
| L16+C1 | 0 | 1 (seed 17) |
| **L16+C16** | **0** | **0** |

### 11.5 Executability

- An independent scratch implementation reproduced every §11.3 class
  count.
- `confirm_analyze.py` then reproduced them, plus every cell-level
  `pfErr` and every selected (method, start) (§14.2).
- Measured on seeds 1-30:
  - exact ties: up to 3 per seed, always within one basin;
  - same-basin spread: ≤ 5.9e-3;
  - smallest cross-basin margin: 2.52 in the direct configurations and
    0.0056 in `A16` (AIDS seed 10, flagged as a near-tie).

### 11.6 Reconciliation: AIDS development comparison, 8/0 earlier vs 9/0 now

Both counts are correct, for different metrics. Each pairs AIDS `L16`
against `A16`, uses selected endpoints, and counts "no result" as not
correct.

| Metric | b / c | Seeds where `A16` is not correct and `L16` is |
|---|---|---|
| `rawStructErr ≤ 1` (the earlier statement, from the v2 draft) | **8 / 0** | 5, 6, 10, 13, 15, 17, 24, 29 |
| `pfErr ≤ 1` (now; the frozen metric) | **9 / 0** | the same 8, **plus 16** |

- **Why seed 16 differs.** `A16` selects Anderson start 11. Its raw
  structural error is 0.823, which is correct, but its `pfErr` is 3.294,
  which is wrong. Its projected `bsErr` is 3.266.
  - The error sits in the n-th good's price row and column, which
    `pfErr` recovers from the other goods by homogeneity and adding-up.
  - Those entries are sums of the estimated entries, so individually
    small errors accumulate.
  - `L16` is correct on seed 16 under both metrics. No direct-method
    classification changed.
- **A third, superseded count exists.** The pilot report's bS-based list
  ("Anderson wrong/never → continuation correct: seeds 5, 6, 10, 13, 15,
  16, 17, 24", also 8) used the projected `bS ≤ 1` metric, with
  continuation as the direct arm. It includes seed 16 (Anderson's `bS`
  3.27) and excludes seed 29 (continuation's `bS` 1.19, a
  projection/u-row case). It is a different metric and arm, not an
  error.
- **Verified by:** `test_confirm_analyze.TestReproducesE0.test_dev_a3_descriptive_and_reconciliation`
  asserts b = 9, c = 0, ties = 21 and the improved list
  [5, 6, 10, 13, 15, 16, 17, 24, 29].

## 12. Symmetry in the shipped estimator, and the documentation-issue log

**Finding (inspection only; nothing changed or assumed).** Shipped
`quaidsFit()` does not impose symmetry inside its objective.
- **Iteration:** homogeneity-only (`src/quaids.src:297-354`). The price
  index is computed only at line 306, from the homogeneity-only gamma.
- **Symmetry:** imposed afterwards by a one-step minimum-distance
  projection (lines 556-607). That projection feeds `bS`/`bestB` (lines
  719 and 726) and the elasticity and Slutzky outputs (lines 1090-1147).
- **Documentation agreeing with the code:**
  - the proc header (lines 74, 98-101);
  - `docs/METHODOLOGY_NOTES.md` step 5;
  - `docs/command-reference/quaidsFit.md:53-65`.
- **Plans:** none for in-objective symmetry.
- **Precedent:** `quaidsCurvatureFit()` and TVP-AIDS do impose symmetry
  inside estimation.

**Recommendation.** Waive E1 for this protocol. `pfErr` scores the
homogeneity-constrained optimum that both methods target, and it is
valid under either later symmetry design. A shipped direct-LM path would
need a separate owner decision: either in-objective symmetry (with E1 as
its smoke gate) or a valid covariance for the projection step.

### Documentation-issue log (recorded only; no file outside this directory edited)

| ID | Location | Issue | Proposed fix (not applied) |
|---|---|---|---|
| D1 | `README.md:8`, `README.md:166`, `docs/METHODOLOGY_NOTES.md:26-32` | "via iterated FGLS with cross-equation restrictions applied through a minimum-distance reparametrization" can be read as symmetry inside the iteration. The code iterates under homogeneity only and imposes symmetry by a post-iteration minimum-distance step. | Reword in a separate, approved docs-only task. |

## 13. Adversarial review: findings and unresolved risks

Severities: **H** could invalidate the decision; **M** limits the
claim; **L** is hygiene.

| ID | Sev. | Finding or risk | Status in v4 |
|---|---|---|---|
| R1 | M | **Winner's curse.** Primary configurations were chosen on seeds 1-30 (30/30). | Open. The held-out rate is expected to be lower. A2 thresholds sit below the development rate. |
| R2 | M | **A3 was not blind for the reference arm.** | **Resolved by demotion.** A3 is descriptive; A5 is reporting-only. |
| R3 | H (for production claims) | **No truth-free reliability signal.** | Accepted scope. It blocks any per-fit trust claim. |
| R4 | M | **The cutoff of 1 is lenient** compared with correct-direct `pfErr` ≤ 0.25. | Sensitivity at 0.25/0.5, borderline listing and A6. |
| R5 | M | **External validity:** one DGP family. | Accepted scope. |
| R6 | L | **AIDS and QUAIDS datasets share draws for a seed.** | Models decided separately; not independent replications. |
| R7 | M | **Seed-17-type fragility:** the correct basin was found by 1 of 32 candidates. | A4 caps W. Open. |
| R8 | L | **Hard-coded layout** in `e0_rescore.py`. | **Resolved:** `confirm_analyze.py` derives and asserts it (tested). |
| R9 | L | **Stopping-rule noise.** | Selection has no tolerance; near-tie flag is reporting only. |
| R10 | M | **No analysis program existed.** | **Resolved:** script, configuration and 24 tests; exact reproduction of E0 (§14). |
| R11 | L | **Non-finite endpoints.** | Ineligible (tested). |
| R12 | M | **E1 waived,** so E2 validates the unconstrained objective. | Owner decision (§12). |
| R13 | L | **Interim stop** is integrity-only. | §8. |
| R14 | L | **Development-mode gate verdicts are not meaningful.** With n = 30, even 30/30 gives a lower limit of 0.884, below the AIDS 0.90 target. | The script labels dev-mode verdicts "informational" (tested); only confirm mode (n = 170) is binding. |
| R15 | M | **The manifest goes stale once the approval ticks are added**, because the protocol's hash changes. A manifest digest written into this file could never verify, because the manifest hashes this file. | **Resolved at freeze (2026-10-08):** the manifest was regenerated and verified after the final ticks and wording. Its SHA-256 is kept in the sidecar `FREEZE_MANIFEST.sha256`, which the manifest does not hash, and in the freeze commit message (§14.5). |
| R16 | L | **The harness writes into `e2_out/`, which must exist first.** `e2_truth_dump.e` refuses to overwrite. | Run order §14.4 creates the directories first. |

## 14. Analysis implementation and freeze manifest

### 14.1 Files (all in this directory unless noted)

| File | Role |
|---|---|
| `confirm_analyze.py` | Implements §4, §5 and §7. `--mode dev` is restricted to seeds 1-30 and aborts on any other seed before parsing it. `--mode confirm` is restricted to seeds 31-200 and refuses to run unless `--manifest` verifies. Output is deterministic: sorted JSON keys, `repr` floats, no timestamps. Exit codes: 2 for an integrity failure, 3 for a guard refusal. |
| `confirm_config.json` | Every rule parameter: configurations, primary/advisory/reference, cutoffs, tolerances, A2 targets and interval method, A4 maximum, A7 targets. |
| `test_confirm_analyze.py` | 24 `unittest` tests (no `pytest` installed; nothing was installed). |
| `confirm_dev_out/` | Canonical development-mode outputs: `results.json`, `selected.csv`, `report.md`, `console.txt`. |
| `confirm_test_output.txt` | Saved test run. |
| `e2_run_m0.e`, `e2_run_m1.e` | E2 entry scripts, derived from `pilot_run_m*.e`. The only differences are seeds `seqa(31,1,170)`, output `e2_out`, and labels. **Not run.** |
| `e2_truth_dump.e` | No-solver truth for seeds 31-200, derived from `e0_truth_dump.e`. The only differences are the seed range, the output path and labels. **Not run.** |
| `e2_determinism.e` | A1 re-solve of seeds 31, 115 and 200 into `e2_det_out/`. **Not run.** |
| `make_freeze_manifest.py`, `FREEZE_MANIFEST.json` | SHA-256 manifest and verifier (§14.5). |

### 14.2 Tests (24/24 pass; `python -m unittest -v test_confirm_analyze`)

| Area | Tests |
|---|---|
| Eligibility | non-converged, NaN or infinite objective, and non-finite endpoint are all ineligible |
| Objective selection | the largest `critExact` wins (direction check); a non-converged candidate with a better objective is excluded |
| Exact ties | cold LM before continuation, then the lowest start |
| Non-finite endpoint | never selected, even with the best objective |
| No result | `None` selection gives class `none` |
| Error classification | `pfErr` = 1 is correct; the next float above 1 is wrong |
| `pfErr` formula | exact inversion on a hand-built consistent truth; the u row excluded; hand-derived values for a gamma perturbation (0.25) and a uniform +3 shift (12) |
| Layout assertion | derived intercept count; a bad layout raises |
| A2 intervals | k = 161 or 160, and 154 or 153, against the 0.90/0.85 targets; edge cases k = 0 and k = n; the A4 upper limit 0.0323 |
| A3 p-value | exact binomial values |
| A4 counting | synthetic seeds with a returned-wrong answer, a no-result seed and a pooled rescue. W counts only the returned-wrong seeds, never "no result". |
| Deterministic output | two runs give byte-identical `results.json`, `selected.csv` and `report.md` |
| Guards | dev mode refuses a seed-31 row; confirm mode refuses a missing or invalid manifest; a corrupted stored error fails the truth gate; an incomplete grid is detected; an empty input directory is an integrity failure (confirm mode also requires both models) |
| **E0 reproduction** (real seeds 1-30) | all 2,880 cell-level `pfErr` values equal `e0_cells_pferr.csv` **exactly** (the 2 non-finite cells match as NaN); all 48 class-table entries (§11.3: 16 configurations × 3 metrics) and all 10 W values (§11.4) equal; every selected (method, start) equals `e0_rescore.py`'s; development A3 descriptive b/c/ties (§11.6); dev gates marked non-binding |

Two of my own test mistakes were caught and fixed while writing the
suite; neither involved the script:
1. A synthetic truth fixture violated homogeneity and adding-up. It was
   replaced with a consistent construction.
2. One hand-derived expected value was wrong. With a uniform +3 shift,
   the n-th row and column intersection moves by +12, not +6. The
   derivation is in the test comment.

### 14.3 Development-mode run (`confirm_dev_out/`)

`python confirm_analyze.py --mode dev --cells pilot_out --truth e0_out --out confirm_dev_out`
- A1 analysis-level integrity: PASS, with 8,637 truth-gate values
  reproduced exactly.
- AIDS `L16`: k = 30/30, W = 0. QUAIDS `L16+C16`: k = 30/30, W = 0.
- The verdicts are marked **informational** because n = 30 (R14).
- Running twice gave byte-identical outputs.

### 14.4 E2 run order (for after the freeze; nothing here has been executed)

```powershell
python make_freeze_manifest.py --verify                     # must print MANIFEST OK
New-Item -ItemType Directory e2_out, e2_det_out
C:\gauss26\tgauss.exe -b -x e2_run_m0.e                     # concurrently with:
C:\gauss26\tgauss.exe -b -x e2_run_m1.e
C:\gauss26\tgauss.exe -b -x e2_truth_dump.e                 # no solver
python pilot_analyze.py --validate e2_out                   # A1 part 1
C:\gauss26\tgauss.exe -b -x e2_determinism.e
python pilot_analyze.py --compare e2_det_out e2_out         # A1 part 3
powershell -ExecutionPolicy Bypass -File pilot_resume_test.ps1   # A1 part 4
python confirm_analyze.py --mode confirm --cells e2_out --truth e2_out --out confirm_e2_out --manifest FREEZE_MANIFEST.json
```

### 14.5 Freeze manifest

- `FREEZE_MANIFEST.json`, written by `make_freeze_manifest.py`, records:
  - the SHA-256 of this protocol;
  - the analysis script, configuration and tests;
  - `pilot_analyze.py`;
  - every harness source E2 executes, including the shipped `src/`
    files and `tests/quaidsfixtures.src` it `#include`s;
  - the E2 entry scripts;
  - the benchmark reference output;
  - the development analysis inputs and the canonical development
    outputs;
  - `git HEAD` and the Python/NumPy/SciPy/pandas versions. The recorded
    HEAD is the **parent** of the freeze commit (`e90a568`), because the
    manifest is written before the commit exists.
- It has no timestamp, so it is reproducible.
- **Freeze procedure, executed 2026-10-08:**
  1. the §15 boxes were ticked and the wording finalized;
  2. the warning-strict test suite was rerun;
  3. `python make_freeze_manifest.py --write` was run;
  4. `python make_freeze_manifest.py --verify` was run;
  5. the sidecar `FREEZE_MANIFEST.sha256` was written in the
     conventional `<sha256>  FREEZE_MANIFEST.json` form and checked
     independently (`sha256sum -c FREEZE_MANIFEST.sha256`);
  6. everything was committed in one local commit.
- **Why the digest lives in the sidecar.** The manifest hashes this
  protocol, so its own digest cannot be written here without
  invalidating it. The manifest does not hash the sidecar.
- That manifest is the frozen record. `confirm_analyze.py --mode
  confirm` refuses to run against any file that no longer matches it.
- Before E2, run both `python make_freeze_manifest.py --verify` and
  `sha256sum -c FREEZE_MANIFEST.sha256`.

## 15. Freeze record (approved by the user 2026-10-08)

Each box is ticked to record the user's explicit 2026-10-08 approval
(America/Chicago) to freeze protocol v4 under its stated limited scope,
together with the evidence verified at freeze. No criterion was changed
when the boxes were ticked.

- [x] **Scope accepted:** synthetic returned-answer accuracy only, no
  truth-free reliability signal (§2, §6). The user stated this on
  2026-10-07 and confirmed it with the 2026-10-08 freeze approval.
- [x] **Configurations final:** AIDS `L16`; QUAIDS `L16+C16` primary,
  `L16+C1` advisory; reference `A16` descriptive only
  (§3, `confirm_config.json`).
- [x] **Selection algorithm (§4) and metric (§5) accepted** as
  implemented in `confirm_analyze.py`.
- [x] **Confirmation gates accepted:**
  - A1 as defined in §7.1;
  - A2 (two-sided 95% Clopper-Pearson lower limit; AIDS k ≥ 161,
    QUAIDS k ≥ 154 of 170);
  - A4 (W ≤ 1).
- [x] **A3 descriptive, A5 required reporting, A6 and A7 descriptive**
  (§7.2, §7.3). A3 is not a confirmatory gate.
- [x] **Reconciliation (§11.6) accepted.**
- [x] **E1 waived** for this protocol (§12).
- [x] **Information audit (§9) accepted.** It was rechecked at freeze:
  `e2_out/`, `e2_det_out/` and `confirm_e2_out/` do not exist.
- [x] **Tests 24/24 and the exact E0 reproduction reviewed** (§14.2,
  §14.3). Rerun at freeze with `ResourceWarning` as an error; the result
  is in `confirm_test_output.txt`.
- [x] **After ticking:** the manifest was regenerated and verified
  (§14.5). Its SHA-256 is recorded in the sidecar
  `FREEZE_MANIFEST.sha256`, which is checked with `sha256sum -c`, and in
  the freeze commit's message. It is not written here because the
  manifest hashes this file (R15).
- [x] **User approval recorded:**
  - date 2026-10-08 (America/Chicago);
  - the frozen commit is the enclosing local commit that adds this
    text; its hash is reported after creation (coordinator report and
    `worker_done`), not inserted here, so the frozen content is never
    amended;
  - E2 is not authorized by this approval and has not been run;
  - nothing has been pushed.
