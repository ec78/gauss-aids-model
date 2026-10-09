# Claude handoff: E2 confirmation COMPLETE (2026-10-08/09); all binding gates pass

Last updated: 2026-10-09

## Start here

This directory holds experiment-only research for `gauss-aids-model`.
Nothing in it is part of the installed package or the CI surface. No
shipped file under `src/` was changed, and `direct_objective.src` was not
changed by the pilot either.

1. Read the repository-root `CLAUDE.md`.
2. **Read `E2_RESULTS.md`.** It is the latest result: the completed E2
   confirmation, its limitations and the recommended next work.
3. Read `CORRECTNESS_VALIDATION_PROTOCOL.md` (FROZEN v4, 2026-10-08) for
   the rules E2 was run under. Its text still says E2 "has not been
   run", because it is the frozen pre-run record and is deliberately
   unchanged.
4. For background, read `SUMMARY.md`: "E2 confirmation" first, then
   "Matched 30-seed pilot (2026-10-05)" and "Phase 3 hardening". Also
   read `PILOT_PROTOCOL.md`.

Never commit or push without the user's explicit approval.

Commit history:
- **`a776717`:** the Phase 3 hardened checkpoint.
- **`e90a568`:** the matched 30-seed pilot; **pushed**.
- **`ad986fc`:** the protocol freeze; local.
- **The E2 evidence and documentation commit,** made on top of
  `ad986fc`; local.

None of the commits after `e90a568` is pushed. All of this is
experiment-only and not shipped.

## Current bottom line (2026-10-09): E2 confirmation complete

Full detail and the evidence index are in `E2_RESULTS.md`.

Under the frozen protocol, on held-out seeds 31-200 (n = 170 per model),
**A1, A2 and A4 all pass for both models**:

| Model | Primary configuration | Correct (A2) | Wrong answers returned (A4) |
|---|---|---|---|
| AIDS | cold LM, 16 starts (`L16`) | 170/170 | 0 |
| QUAIDS | cold LM + continuation pool (`L16+C16`) | 160/170 | 0 |

QUAIDS has **10 visible no-result seeds**: 43, 112, 125, 134, 151, 160,
167, 182, 187 and 190. Anderson also failed all 10, with no result or a
wildly wrong answer.

A1 (all four parts) passes. A5: no regressions against Anderson.
A3 is descriptive only: improved/regressed/tied 44/0/126 for AIDS and
86/0/84 for QUAIDS.

**Limitations:**
- synthetic DGP only;
- no truth-free reliability signal;
- QUAIDS has 10 no-result seeds;
- A3 is non-confirmatory;
- timing is descriptive, not a speed conclusion;
- the production symmetry design is unresolved.

**Next work (not started):**
1. preserve the evidence;
2. prototype AIDS direct LM as an opt-in;
3. investigate the 10 QUAIDS no-result seeds under a new exploratory
   protocol;
4. resolve symmetry before shipping QUAIDS direct LM;
5. benchmark speed separately.

The sections below are the pilot-era record, kept for history.

## Bottom line

**Stage A gate.** It passed with no change to the comparison:
- `pilot_validate.e` Part 1: 27 of 27.
- `pilot_analyze.py --validate pilot_validation`: 58 of 58.
- Resume test: 6 of 6.
- `--validate pilot_out`: 377 of 377.

Two of those reference-comparison tolerances were loosened to match the
8-significant-digit printed references; this is recorded in
`PILOT_PROTOCOL.md`.

The pilot's Anderson arm reproduces the 200-seed benchmark's per-seed
results on all 60 seed/model cases. This proves start and RNG
equivalence.

**Stage B.** All 2,880 cells completed: seeds 1-30, both models, 3
methods × 16 starts.

Results are given as selected correct counts out of 30. "bS" is primary
(benchmark-compatible, ≤ 1) and "raw" is secondary (all rows except u):

| | Anderson | cold LM | continuation |
|---|---|---|---|
| AIDS bS / raw | 19 / 22 | 27 / 30 | 26 / 30 |
| QUAIDS bS / raw | 10 / 13 | 23 / 29 | 23 / 29 |
| QUAIDS Stone-only bS / raw | 8 / 10 | 17 / 23 | 21 / 27 |

What the results show:
- **Direct LM (cold or continuation) is a genuine solver improvement over
  Anderson.** The selected and oracle results both improve. In QUAIDS,
  14 seeds improved and 1 regressed; the regression is a u-row cutoff
  effect.
- **Continuation ties cold LM at 16 starts.** Its value is per-start
  robustness:
  - QUAIDS Stone-only 27 vs 23 raw correct;
  - 476/480 vs 379/480 starts genuinely converged;
  - a single basin on 28 of 30 seeds.
- **Continuation costs much more.** The median per seed with 16 starts is
  98 s for QUAIDS, against 8 s for cold LM.
- **Solution selection is not the bottleneck.** Selected and oracle
  agree in class on 178 of 180 seed/model/method cases.
- **Most remaining bS "wrong" seeds for the direct methods are
  right-basin optima.** Raw error is ≤ 0.09 and the objective beats the
  truth's. They fail the benchmark's `bS` cutoff because of the u row,
  or because `quaidsFit()`'s symmetry projection moves alpha or gamma by
  1-2. `pilot_bs_diag.e` shows this. The projection uses a delta-method
  covariance built for the fixed-point estimator.
- **Genuine wrong basins:**
  - QUAIDS seed 5, continuation: all 16 starts go to the same wrong
    basin, while cold LM finds the right one.
  - QUAIDS seed 17, cold LM: only 1 start converges, and it is wrong;
    continuation finds the right basin from 1 start.

  The failures are complementary. Post-hoc pooling of both methods'
  candidates (exploratory, not prespecified) reaches raw 30/30 on
  QUAIDS.

## Recommendation (2026-10-05, pilot era)

This was carried out: the metric was fixed (`pfErr`), the hybrid was
prespecified and frozen, and E2 confirmed it. See `E2_RESULTS.md`.

1. **Model-specific hybrid.** Cold LM for AIDS; continuation, or a small
   pool of cold LM plus continuation, for QUAIDS.
2. **No more solution-selection work.**
3. **Before any 200-seed run, decide the direct-endpoint correctness
   metric.** The preferred option is to impose symmetry inside the
   direct objective, so `bS` needs no post-hoc projection. The
   alternatives are to prespecify an unprojected `bS` or `bsErrExclU`,
   or to derive a sandwich covariance.
4. **Then prespecify the hybrid and scale up.**

30 seeds support no population-wide claim.

**Update (2026-10-07; historical, superseded by the 2026-10-08 freeze
below):** `CORRECTNESS_VALIDATION_PROTOCOL.md` (then a draft) designs the pre-200-seed validation. It sets a projection-free
primary metric, acceptance criteria A1-A7 on held-out seeds 31-200, and a
trust-flag calibration. E0, the no-solver re-score, has run (protocol
section 9). The truth gate passed with 21,594 of 21,594 values exact, and
all 6 predictions held. No truth-free trust flag met the selection rule.
Section 10 records that the shipped quaidsFit imposes symmetry by a
post-iteration minimum-distance step, not inside the objective. Protocol
v3 adds an adversarial pre-freeze review, which recommends "do not freeze
yet". The E2 analysis script must be written and must reproduce the
seeds 1-30 results first, and the owner must decide scope items R2 and R3.
Nothing is frozen.

**Update (2026-10-07, v4):** the owner accepted the limited scope, and A3
is now descriptive. The gates are A1, A2 and A4. The analysis script
(`confirm_analyze.py`, `confirm_config.json`) passes 24/24 tests and
reproduces E0 exactly. The E2 entry scripts are written but not run.
`FREEZE_MANIFEST.json` was then a pre-freeze candidate (34 files).

**Update (2026-10-08): protocol v4 is FROZEN.** The user explicitly
approved the freeze under the limited scope: synthetic returned-answer
accuracy for this DGP, no truth-free reliability signal, and A3
descriptive only.
- Protocol section 15 is ticked.
- The final manifest (34 files) was regenerated and verified. Its SHA-256
  is in the sidecar `FREEZE_MANIFEST.sha256` and in the freeze commit
  message.
- The freeze is one local commit on top of `e90a568`. It is **not
  pushed**.
- At the time of the freeze, E2 had not been run. **Superseded:** the
  user approved E2 on 2026-10-08. It ran from the verified LF working
  tree, and its A1 resume step completed 2026-10-09. See
  `E2_RESULTS.md`.

Before any re-analysis, run `python make_freeze_manifest.py --verify` and
`sha256sum -c FREEZE_MANIFEST.sha256`, then follow protocol section 14.4.

**Line-ending caveat (known at freeze; fails closed).** This machine has
`core.autocrlf=true` at the system level.
- **The current working tree verifies all 34 manifest files.**
- **A fresh-checkout simulation verifies 27 of 34**, plus the manifest
  and the sidecar. It was a `git checkout-index` export taken after
  `dev/experiments/convergence-anderson/.gitattributes` marked the
  experiment package `-text`.
- **The 7 that do not verify** are production and reference files
  outside the experiment package:
  - `src/quaids.sdf`
  - `src/quaids.src`
  - `src/quaidselas.src`
  - `src/quaidsiv.src`
  - `src/quaidsslutzky.src`
  - `src/quaidsutil.src`
  - `tests/quaidsfixtures.src`

  Their git blobs and the current working copies are LF, which is what
  the manifest hashed. A fresh Windows checkout under
  `core.autocrlf=true` rewrites them as CRLF, so their hashes change.
- **This fails closed.** `make_freeze_manifest.py --verify` and
  `confirm_analyze.py --mode confirm` then refuse to run; they never
  pass silently.
- **E2 ran from the verified LF working tree.** Any re-run or
  re-analysis must use that tree, or a checkout made with
  `core.autocrlf=false` or `input`, and must re-verify the manifest
  first.
- No root-level `.gitattributes` or production-file checkout behavior
  was changed at freeze. That was a deliberate scope decision.

## Artifacts (all added in `e90a568`)

- `PILOT_PROTOCOL.md`
- `pilot_lib.src`, `pilot_harness_core.src`
- `pilot_validate.e`
- `pilot_run_m0.e`, `pilot_run_m1.e`
- `pilot_resume_test.e`, `pilot_resume_test.ps1`, `pilot_resume_test_output.txt`
- `pilot_analyze.py`
- `pilot_bs_diag.e`, `pilot_bs_diag_output.txt`
- `pilot_posthoc_pooled.py`
- `pilot_validation/` (Stage A cells and logs)
- `pilot_resume_test/`
- `pilot_out/`:
  - Stage B cells (6.4 MB of CSV) and console logs;
  - `pilot_report.md` and `pilot_seed_table.csv`;
  - validate, determinism and post-hoc outputs.

Exact commands are in SUMMARY.md, "Matched 30-seed pilot", "Commands".

The harness resumes per cell: rerun the same entry script and finished
`(seed, method, start)` keys are skipped. To extend seeds, copy
`pilot_run_m*.e` and change `pilotSeeds`.

## Independent review of the pilot (Claude, 2026-10-05)

The review found no correctness problems.

**What was checked:**
- **Harness:** the column indices, the `bS` and raw row ranges, and the
  convergence eligibility.
- **Analysis:** the selection direction, the oracle, the dedup and the
  classification rules all match `PILOT_PROTOCOL.md`.
- **Validators rerun:** selftest, 58 of 58 and 377 of 377.
- **Report regenerated from the raw cells:** the seed table is
  byte-identical, and the report differs only in its directory header.
- **Handoff numbers:** every headline number matches.
- **Silent misclassification:** none of the 101 symmetry-projection
  fallbacks or the 3 `bS` chain failures falls on a selected (converged)
  cell.

The pilot was not re-solved. The review relies on the saved
determinism check: 288 cells, 0 differences.

**Interpretation notes:**
- **Cold LM's QUAIDS per-start failures are mostly start-point
  failures:** 53 infeasible at homotopy 1, and 44 covariance-inversion
  failures. Continuation's per-start robustness largely comes from
  starting at homotopy 0.
- **A promising cheap hybrid (post hoc).** Running cold LM and
  continuation once each from the Stone start gives QUAIDS raw 29/30 at
  a median 7.5 s per seed. That matches cold LM with 16 starts. Treat
  it as a candidate to prespecify, not as a result.

**Raw outputs are committed as plain CSV,** not zipped. git already
compresses them to about 2.3 MB, against 2.28 MB as a zip. Plain CSV
keeps GAUSS resume and `pilot_bs_diag.e` working. `.gitattributes` marks
the cells `-text -diff`.

## Historical corrections still in force

- **AIDS structural error in the 2026-10-02 Phase 3 write-up was
  misreported as 0.724.** That figure was the u-row error. The correct
  value is 0.058, and the pilot reproduces it (V3c).
- **QUAIDS continuation endpoints from the two smoke-gate starts are
  numerically equivalent (5.8e-9), not identical.**

## Environment / gotchas

See the repo-root `CLAUDE.md`. Gotchas specific to this work:
- **Reserved words.** `line`, `key` and `le` are reserved in GAUSS and
  fail as variable names (`G0276`, `G0136`). Watch for other operator
  words too (`ge`, `eq`, `ne`, `lt`, `gt`).
- **CSV I/O.** GAUSS 26 has no `getenv` and no `fclose`; use `close`.
  - Write a numeric row with `strjoin(sprintf("%.17g", vec(x))', ",")`.
    A column vector is required; a row vector concatenates the values
    with no separator.
  - `csvReadM(file, 1)` round-trips this exactly, with `.` read back as
    missing.
- **Passing job parameters.** Because there is no `getenv`, the entry
  scripts set variables and then `#include pilot_harness_core.src`
  (top-level code). Do not `#include` that file twice in one script.
- **Editing `/* ... */` doc comments.** A partial edit can leave a
  dangling `/*` or `*/`, producing a cascading `G0008` far from the edit.
- **Row slices of `bRaw`.** Use `_pilotErrs` or `_phase3Errs`, not
  hand-written ranges.
- **Capturing output.** `tgauss ... > file` in bash captures the banner
  too; PowerShell `Out-File` captures only the script's prints.
