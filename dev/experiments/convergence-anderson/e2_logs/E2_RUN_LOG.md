# E2 run log (frozen protocol v4, CORRECTNESS_VALIDATION_PROTOCOL.md section 14.4)

Task task_2422e3bd111c / dispatch ctx_7dd5f6170ffd. Directory: dev/experiments/convergence-anderson. Times are local (America/Chicago).

## Pre-flight (2026-10-08 06:09:03)

```
$ git rev-parse HEAD
ad986fc3214bce3c86397620f51e3e7292b29430
exit=0
$ git status --porcelain (repo)
exit=1 (output above, excluding this log dir)
$ python make_freeze_manifest.py --verify
MANIFEST OK: 34 files verified; manifest sha256 ee4e359692254d46f97c595c2b7215510a335f9870a3313b44b3fe70a6d58a7b
exit=0
$ sha256sum -c FREEZE_MANIFEST.sha256
FREEZE_MANIFEST.json: OK
exit=0
$ ls e2_out e2_det_out confirm_e2_out
ls: cannot access 'e2_out': No such file or directory
ls: cannot access 'e2_det_out': No such file or directory
ls: cannot access 'confirm_e2_out': No such file or directory
```

## Step 2: create output directories (2026-10-08 06:09:11)

```
$ mkdir e2_out e2_det_out   (equivalent of New-Item -ItemType Directory e2_out, e2_det_out)
exit=0
```

## Step 3: e2_run_m0.e and e2_run_m1.e, concurrently

Launched 2026-10-08 06:09:11:
```
$ C:\gauss26\tgauss.exe -b -x e2_run_m0.e > e2_logs/e2_run_m0_console.txt 2>&1
$ C:\gauss26\tgauss.exe -b -x e2_run_m1.e > e2_logs/e2_run_m1_console.txt 2>&1
```

Completion:
```
e2_run_m0.e exit=0 finished 2026-10-08 07:34:06
e2_run_m0_console.txt: 8160 cell lines; 0 GAUSS errors; final line: e2_run_m0.e: run complete.
e2_out/cells_m0.csv rows: 8160
```
```
e2_run_m1.e exit=0 finished 2026-10-08 12:18:48
e2_run_m1_console.txt: 8160 cell lines; 0 GAUSS errors; final line: e2_run_m1.e: run complete.
e2_out/cells_m1.csv rows: 8160
```

## Step 4: no-solver truth dump (2026-10-08 12:19:08)

```
$ C:\gauss26\tgauss.exe -b -x e2_truth_dump.e > e2_logs/e2_truth_dump_console.txt 2>&1
exit=0
truth rows written: m0=170 m1=170; GAUSS errors: 0; final line: e2_truth_dump.e: run complete.
e2_out/truth_m0.csv rows=170, truth_m1.csv rows=170
```

## Step 5: A1 part 1 -- pilot_analyze.py --validate e2_out (2026-10-08 12:19:20)

```
$ python pilot_analyze.py --validate e2_out > e2_logs/pilot_analyze_validate_e2_out.txt
exit=0
PASS lines: 2023; FAIL lines: 0; SKIP lines: 4
VALIDATION: ALL 2023 CHECKS PASSED
```
V3b benchmark-reproduction checks: 1995 PASS covering 340 seed/model cases; the 4 SKIPs are the pilot smoke-gate references for dev seeds 6/9 (absent from E2 by design).

## Step 6: A1 part 3 -- determinism re-solve e2_determinism.e (seeds 31, 115, 200) (2026-10-08 12:19:35)

```
$ C:\gauss26\tgauss.exe -b -x e2_determinism.e > e2_logs/e2_determinism_console.txt 2>&1
e2_determinism.e exit=0 finished 2026-10-08 12:24:39
e2_determinism_console.txt: 288 cell lines (expected 288 = 2 models x 3 seeds x 3 methods x 16 starts); 0 GAUSS errors; final line: e2_determinism.e: run complete.
```

## Step 7: A1 part 3 -- pilot_analyze.py --compare e2_det_out e2_out (2026-10-08 12:24:54)

```
$ python pilot_analyze.py --compare e2_det_out e2_out > e2_logs/pilot_analyze_compare_det.txt
exit=0
DIFF lines: 0
DETERMINISM: 288 common cells, 0 differences (wallTime excluded)
```

## Step 8: A1 part 4 -- pilot_resume_test.ps1: DEFERRED (2026-10-08 12:25:12)

Not run yet. The frozen script re-solves dev seed 1 and rewrites the committed pilot_resume_test/ directory, which conflicts with the task's held-out-seeds-only and no-tracked-change constraints. Coordinator question msg_c4ebfbdfd084 has been pending since about 06:10, and escalation msg_645ddabb7754 was sent at 12:01. Step 9 is run first because it does not depend on step 8. A1 cannot be called PASS until part 4 is resolved.

## Step 9: confirm_analyze.py --mode confirm (2026-10-08 12:25:12)

```
$ python confirm_analyze.py --mode confirm --cells e2_out --truth e2_out --out confirm_e2_out --manifest FREEZE_MANIFEST.json > e2_logs/confirm_analyze_console.txt
exit=0
AIDS: A2 PASS (k=170/170), A4 PASS (W=0)
QUAIDS: A2 PASS (k=160/170), A4 PASS (W=0)
wrote confirm_e2_out
```

## Post-run integrity (2026-10-08 12:25:37)

```
$ python make_freeze_manifest.py --verify
MANIFEST OK: 34 files verified; manifest sha256 ee4e359692254d46f97c595c2b7215510a335f9870a3313b44b3fe70a6d58a7b
exit=0
$ sha256sum -c FREEZE_MANIFEST.sha256
FREEZE_MANIFEST.json: OK
exit=0
$ git rev-parse HEAD
ad986fc3214bce3c86397620f51e3e7292b29430
$ git status --porcelain (tracked changes only)
(end; empty = no tracked file changed)
$ git status --porcelain (untracked)
?? dev/experiments/convergence-anderson/confirm_e2_out/
?? dev/experiments/convergence-anderson/e2_det_out/
?? dev/experiments/convergence-anderson/e2_logs/
?? dev/experiments/convergence-anderson/e2_out/
```

## Summary (written at close-out)

### A1 integrity components

| Component | Outcome |
|---|---|
| Part 1: `pilot_analyze.py --validate e2_out` | **PASS**: 2,023 PASS, 0 FAIL. The 4 SKIPs are the dev-seed smoke references for seeds 6/9, absent by design. 340/340 held-out seed/model benchmark reproductions passed. |
| Part 2: `confirm_analyze.py --mode confirm` analysis-level integrity | **PASS**: manifest verified, complete grid, 48,731 truth-gate values reproduced exactly, layout asserted, no non-finite selected endpoint. |
| Part 3: determinism | **PASS**: `e2_determinism.e` re-solved 288 cells (seeds 31, 115, 200, both models); `--compare` found 0 differences apart from wall time. |
| Part 4: `pilot_resume_test.ps1` | **NOT RUN, PENDING a coordinator decision** (`msg_c4ebfbdfd084`, unanswered for about 7 h; escalation `msg_645ddabb7754`). The frozen script re-solves dev seed 1 and rewrites the committed `pilot_resume_test/`, which conflicts with this task's constraints. |

**A1 overall: not yet PASS.** Parts 1-3 pass; part 4 is unresolved.

### Gates (binding, n = 170)

| Model | Primary | A2 | A4 |
|---|---|---|---|
| AIDS | `L16` | k = 170/170, CP lower 0.9785 ≥ 0.90: **PASS** | W = 0 ≤ 1: **PASS** |
| QUAIDS | `L16+C16` | k = 160/170, CP lower 0.8945 ≥ 0.85: **PASS** | W = 0 ≤ 1: **PASS** |

QUAIDS has 10 visible "no result" seeds and no returned-wrong answers: 43, 112, 125, 134, 151, 160, 167, 182, 187, 190.

### Required reporting and descriptive results

- **A5 (required):** regressions vs `A16` = none in either model.
- **A3 (descriptive only; the reference arm is not blind):**
  - AIDS b/c/ties = 44/0/126, one-sided exact p = 5.7e-14;
  - QUAIDS 86/0/84, p = 1.3e-26.
- **A6 (descriptive):** primary correct counts are the same under raw error and at `pfErr` cutoffs 0.5 and 2 (AIDS 170, QUAIDS 160).
- **Borderline:** one primary-configuration seed per model, seed 189 (AIDS `pfErr` 0.312, QUAIDS 0.256); both are correct.
- **Near-tie flag (diagnostic):** AIDS primary seeds 47, 90, 104 and 123; QUAIDS none.
- **Advisory QUAIDS `L16+C1`:** k = 156, W = 1.
- **A7 (descriptive; timing is not a speed conclusion):** primary per-seed wall time summed over cells:
  - AIDS median 4.2 s, p90 6.0 s, max 8.2 s;
  - QUAIDS median 79.5 s, p90 238 s, max 1,134 s;
  - 0 cells over 900 s.
- **Elapsed time:**
  - AIDS grid 06:09-07:34 (about 1 h 25 min);
  - QUAIDS grid 06:09-12:18 (about 6 h 09 min; the two ran concurrently);
  - truth dump about 1 min; determinism re-solve 12:19-12:24.

### Frozen-input check

- `FREEZE_MANIFEST.json` verified 34/34 and the sidecar was OK both before and after the run.
- HEAD is unchanged at `ad986fc3214bce3c86397620f51e3e7292b29430`.
- No tracked file changed.
- New untracked directories only: `e2_out/`, `e2_det_out/`, `confirm_e2_out/`, `e2_logs/`.
- Nothing was committed or pushed.

## Step 8 (completed 2026-10-09): A1 part 4, the resume test (coordinator decision: option A)

Task task_076c4acaedeb / dispatch ctx_a6406ff85ab8. E2 was not rerun.
`e2_run_m0.e` and `e2_run_m1.e` were not touched, and neither were
`e2_out/`, `e2_det_out/` or `confirm_e2_out/`.

```
$ python make_freeze_manifest.py --verify     -> MANIFEST OK: 34 files verified (before)
$ sha256sum -c FREEZE_MANIFEST.sha256          -> FREEZE_MANIFEST.json: OK (before)
$ powershell -ExecutionPolicy Bypass -File pilot_resume_test.ps1 *> e2_logs/resume_test/console.txt
start=2026-10-09 05:43:15 end=2026-10-09 05:43:28 exit=0
PASS V6 first run writes 6 cells  rows=6; PILOT CORE DONE: cells run       6.0000000 cells skipped (already on disk)       0.0000000
PASS V6 rerun adds no rows  rows=6; PILOT CORE DONE: cells run       0.0000000 cells skipped (already on disk)       6.0000000
PASS V6 rerun leaves file byte-identical  
PASS V6 resume restores 6 rows  rows=6; PILOT CORE DONE: cells run       2.0000000 cells skipped (already on disk)       4.0000000
PASS V6 keys unique after resume  
PASS V6 recomputed cells identical except wallTime  
RESUME TEST: ALL 6 CHECKS PASSED
$ cp -p pilot_resume_test/* e2_logs/resume_test/generated/      (evidence copy, 5 files)
$ git checkout -- dev/experiments/convergence-anderson/pilot_resume_test   -> exit 0
$ git diff --quiet HEAD -- dev/experiments/convergence-anderson/pilot_resume_test  -> matches HEAD
$ git status --porcelain --untracked-files=no  -> (empty: 0 tracked modifications)
$ python make_freeze_manifest.py --verify     -> MANIFEST OK: 34 files verified (after)
$ sha256sum -c FREEZE_MANIFEST.sha256          -> FREEZE_MANIFEST.json: OK (after)
```

- The test re-solved dev seed 1 (AIDS, starts 1-2, 6 cells). These cells
  are mechanics-only, and they are excluded from every E2 gate and
  statistic.
- Evidence:
  - `e2_logs/resume_test/console.txt`;
  - `e2_logs/resume_test/exit.txt`;
  - `e2_logs/resume_test/generated/` (the regenerated fixture before
    restoration: cells_m0.csv, columns.txt, run_first.txt, run_rerun.txt,
    run_resume.txt).
- Before restoration, the run had modified 3 tracked fixture files
  (cells_m0.csv, run_first.txt, run_resume.txt). All were restored from
  HEAD.

## Final A1 verdict: **PASS**

| Part | Outcome |
|---|---|
| Part 1 | `--validate e2_out` PASS |
| Part 2 | confirm-mode analysis integrity PASS |
| Part 3 | determinism, 0 differences, PASS |
| Part 4 | resume test 6/6, exit 0, PASS |

A2 and A4 are unchanged from the 2026-10-08 run (`confirm_e2_out/`,
not recomputed):

| Model | A2 | A4 |
|---|---|---|
| AIDS `L16` | k = 170/170, CP lower 0.9785, PASS | W = 0, PASS |
| QUAIDS `L16+C16` | k = 160/170, CP lower 0.8945, PASS | W = 0, PASS |

**Under the frozen decision mapping (§8), A1, A2 and A4 all pass for
both models.**
