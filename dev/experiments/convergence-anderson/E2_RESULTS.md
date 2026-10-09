# E2 confirmation results: direct LM vs Anderson for iterated AIDS/QUAIDS

**Status: complete.** All binding gates pass for both models.
- **When:** run 2026-10-08, with the A1 resume step completed 2026-10-09.
- **Rules:** executed under the frozen protocol
  `CORRECTNESS_VALIDATION_PROTOCOL.md` (FROZEN v4), from freeze commit
  `ad986fc3214bce3c86397620f51e3e7292b29430`.
- **Manifest:** `FREEZE_MANIFEST.json` verified 34/34 before and after,
  sidecar digest
  `ee4e359692254d46f97c595c2b7215510a335f9870a3313b44b3fe70a6d58a7b`.
- **Code status:** experiment-only research. Nothing here changes the
  shipped `quaids` package or its documented behavior.

## Read this first: limitations

1. **Synthetic data only.** Every result comes from one synthetic data
   generator (`_quaidsSyntheticDGP`: N = 5 goods, tobs = 3000, the
   benchmark's 16 start recipes). Nothing here is evidence about real
   data, other DGPs, other sample sizes or other good counts.
2. **No truth-free reliability signal.** The gates measure how often the
   selector *returns* a correct or wrong answer when the truth is known.
   E0 found no truth-free signal that identifies a wrong answer, and none
   was validated. A user fitting real data still cannot tell a returned
   wrong answer from a correct one.
3. **QUAIDS has 10 visible no-result seeds out of 170.** On those seeds
   the primary configuration returns no estimate at all. That is an
   honest, visible failure, but it is still a failure: about 6% of the
   QUAIDS seeds.
4. **A3 is non-confirmatory.** The comparison with Anderson is
   descriptive only. Anderson's held-out results were known when the
   protocol was written (protocol §9).
5. **Timing is not a speed conclusion.** Wall times come from two
   concurrent processes on one 12-thread machine, with the experiment
   harness's finite-difference Jacobians. They are reported for context
   only.
6. **Production symmetry design is unresolved.** E2 scores the
   homogeneity-constrained optimum with the projection-free metric
   (`pfErr`). Shipped `quaidsFit()` imposes symmetry by a post-iteration
   minimum-distance projection. Applied at direct-LM endpoints, that
   projection was shown in the pilot to move alpha by 1-2. A shipped
   direct-LM path needs either in-objective symmetry or a valid
   projection covariance first (protocol §12).

## Binding confirmation gates (held-out seeds 31-200, n = 170 per model)

"Correct" means the selected endpoint has `pfErr ≤ 1`. "W" counts
returned answers that are wrong. A2 uses the exact Clopper-Pearson
two-sided 95% interval's lower limit.

| Model | Primary configuration | A1 | A2: correct / 170 (CP lower vs target) | A4: W (max 1) | Verdict |
|---|---|---|---|---|---|
| AIDS | `L16` (cold LM, 16 starts) | PASS | **170** (0.9785 vs 0.90): PASS | **0**: PASS | **PASS** |
| QUAIDS | `L16+C16` (cold LM + continuation, 16 starts each) | PASS | **160** (0.8945 vs 0.85): PASS | **0**: PASS | **PASS** |

The QUAIDS no-result seeds are 43, 112, 125, 134, 151, 160, 167, 182, 187
and 190.

**A1 integrity components (all PASS):**

| Component | Result |
|---|---|
| 1. `pilot_analyze.py --validate e2_out` | 2,023 PASS, 0 FAIL. The 4 SKIPs are the dev-seed smoke references, absent by design. All 340 seed/model Anderson lines reproduce the 200-seed benchmark. |
| 2. `confirm_analyze.py --mode confirm` | Manifest verified, complete grid, 48,731 truth values reproduced exactly. |
| 3. Determinism (`e2_determinism.e`, seeds 31, 115, 200) | 288 re-solved cells, 0 differences apart from wall time. |
| 4. Resume (`pilot_resume_test.ps1`) | 6/6 checks, exit 0. The dev-seed-1 cells are mechanics-only, excluded from every gate, and the committed fixture was restored afterwards. |

Under the frozen decision mapping (protocol §8), A1, A2 and A4 all pass
for both models. That supports a recommendation for production-*design*
work, subject to §2's limits. It is not a production change.

## Required and descriptive results

- **A5 (required):** regressions against `A16` are **none** in either
  model.
- **A3 (descriptive, non-confirmatory):**

  | Model | Primary improved | Regressed | Tied | Exact one-sided McNemar p |
  |---|---:|---:|---:|---:|
  | AIDS | 44 | 0 | 126 | 5.7e-14 |
  | QUAIDS | 86 | 0 | 84 | 1.3e-26 |

- **A6 (descriptive):** the primary correct counts are unchanged under
  `rawStructErr` and at `pfErr` cutoffs 0.5 and 2 (AIDS 170, QUAIDS 160).
  Correct primary selections have `pfErr` ≤ 0.25 except seed 189 in each
  model (0.31 AIDS, 0.26 QUAIDS, both correct).
- **Near-tie diagnostic:** flagged on AIDS primary seeds 47, 90, 104 and
  123; QUAIDS none. It is never used for selection.

### All configurations (selected endpoint; none / wrong / correct)

| Model | Config | `pfErr` (primary metric) | raw | projected `bS` (context only) | W |
|---|---|---|---|---|---:|
| AIDS | `A16` Anderson | 1/43/126 | 1/36/133 | 1/67/102 | 43 |
| AIDS | `L1` cold LM, Stone only | 0/0/170 | 0/0/170 | 0/38/132 | 0 |
| AIDS | `C1` continuation, Stone only | 0/0/170 | 0/0/170 | 0/33/137 | 0 |
| AIDS | **`L16`** (primary) | **0/0/170** | 0/0/170 | 0/35/135 | **0** |
| AIDS | `L16+C16` | 0/0/170 | 0/0/170 | 0/35/135 | 0 |
| QUAIDS | `A16` Anderson | 8/88/74 | 8/88/74 | 8/103/59 | 88 |
| QUAIDS | `L1` | 26/13/131 | 26/13/131 | 26/49/95 | 13 |
| QUAIDS | `C1` | 14/3/153 | 14/3/153 | 14/36/120 | 3 |
| QUAIDS | `L1+C1` | 13/2/155 | 13/2/155 | 13/38/119 | 2 |
| QUAIDS | `L16` | 20/2/148 | 20/2/148 | 20/37/113 | 2 |
| QUAIDS | `C16` | 10/1/159 | 10/1/159 | 10/37/123 | 1 |
| QUAIDS | `L16+C1` (advisory) | 13/1/156 | 13/1/156 | 13/37/120 | 1 |
| QUAIDS | **`L16+C16`** (primary) | **10/0/160** | 10/0/160 | 10/38/122 | **0** |

Every other configuration in `confirm_e2_out/report.md` (AIDS `C16`,
`L1+C1`, `L16+C1`) is 0/0/170.

Observations (descriptive; seeds 31-200 confirm what the pilot saw on
seeds 1-30):
- **QUAIDS needs both methods.** No single-method QUAIDS configuration
  reaches W = 0. Pooling cold LM and continuation does.
- **AIDS does not.** Every direct configuration is perfect, including
  cold LM from the Stone start alone.
- **The projected-`bS` metric still disagrees with `pfErr`** for direct
  endpoints: 35 AIDS and 38 QUAIDS primary seeds are "wrong" on `bS` but
  correct on `pfErr`. This is the pilot's known projection/u-row
  artifact, and it is why `pfErr` is the frozen primary metric.

### The 10 QUAIDS no-result seeds (stored outputs only; no new solves)

On every one of these seeds:
- all 16 cold-LM starts ended with retcode 3 (infeasible start) or 5
  (covariance inversion failure);
- all 16 continuation runs failed at stage zero (status 1) or exhausted
  the minimum step (status 2).

**Anderson also failed every one of them.** It returned no converged
endpoint on seeds 43, 125, 134, 151 and 182, and a wildly wrong one on
112, 160, 167, 187 and 190 (`pfErr` 380-3,855). These seeds are hard for
every method tested. The direct methods fail *visibly* (no result),
where Anderson more often fails *silently* (a wrong answer).

## Runtime (descriptive only; not a speed benchmark)

**Elapsed time:**
- AIDS grid: 06:09-07:34 on 2026-10-08, about 1 h 25 min.
- QUAIDS grid: 06:09-12:18, about 6 h 09 min. The two grids ran as
  concurrent processes.
- Truth dump: about 1 min.
- Determinism re-solve: about 5 min.
- Resume test: 13 s.

**Summed per-cell wall time by method:**

| Model | Anderson | Cold LM | Continuation |
|---|---:|---:|---:|
| AIDS | 0.01 h | 0.21 h | 1.18 h |
| QUAIDS | 0.02 h | 0.41 h | 5.73 h |

**Primary configuration per-seed wall time** (median / p90 / max):
- AIDS `L16`: 4.2 / 6.0 / 8.2 s.
- QUAIDS `L16+C16`: 79.5 / 238 / 1,134 s.

No single cell exceeded 900 s.

## Evidence package (all in this directory)

| Path | Contents |
|---|---|
| `e2_out/` | Raw harness cells `cells_m0.csv`/`cells_m1.csv` (8,160 rows each: 170 seeds × 3 methods × 16 starts), `columns.txt`, and the no-solver truth `truth_m0.csv`/`truth_m1.csv` (170 rows each). |
| `e2_det_out/` | Determinism re-solve cells (seeds 31, 115, 200; 288 cells). |
| `confirm_e2_out/` | `confirm_analyze.py --mode confirm` outputs: `report.md`, `results.json`, `selected.csv`. These are the binding gate results. |
| `e2_logs/E2_RUN_LOG.md` | Every §14.4 command with its exit status and time, plus pre- and post-run integrity checks. |
| `e2_logs/*_exit.txt`, `e2_logs/e2_truth_dump_console.txt`, `e2_logs/e2_determinism_console.txt` | Exit stamps for both grids, and GAUSS console captures for the truth dump and the determinism run. |
| `e2_logs/e2_run_m0_console.txt`, `e2_logs/e2_run_m1_console.txt` | **Deliberately not committed**, by the owner's decision of 2026-10-09 on repository size. They are kept in the local working tree only, git-ignored. Each line repeats per-cell data already in `e2_out/cells_m*.csv`. Their completion facts are recorded below and in `E2_RUN_LOG.md`. |
| `e2_logs/pilot_analyze_*.txt`, `e2_logs/confirm_analyze_console.txt` | A1 part 1 and 3 outputs and the confirm-mode console. |
| `e2_logs/resume_test/` | A1 part 4 console, exit stamp, and the regenerated fixture (`generated/`) captured before `pilot_resume_test/` was restored from HEAD. |

**Uncommitted console logs (integrity record).** Archive these files
outside git if they must be kept long term.

| File | Bytes | Lines | Cell lines | GAUSS errors | Final line | SHA-256 |
|---|---:|---:|---:|---:|---|---|
| `e2_run_m0_console.txt` | 1,877,183 | 8,172 | 8,160 | 0 | `e2_run_m0.e: run complete.` | `85edc9438db279712f18214ef40f7b73ee029e107e0ee7f183357b15adb13e57` |
| `e2_run_m1_console.txt` | 1,881,183 | 8,172 | 8,160 | 0 | `e2_run_m1.e: run complete.` | `e53ff36be80abd72263e795703a8a9894b04b32cfdc628e4c45b5bead33f8543` |

To re-derive the gate results without solving:

```
python confirm_analyze.py --mode confirm --cells e2_out --truth e2_out --out <dir> --manifest FREEZE_MANIFEST.json
```

This requires the manifest to verify. See the line-ending caveat in
`CLAUDE_HANDOFF.md` for the 7 LF production files.

**Documents intentionally unchanged:**
- `CORRECTNESS_VALIDATION_PROTOCOL.md`, `FREEZE_MANIFEST.json`/`.sha256`
  and all manifest-covered scripts and configuration are frozen. The
  protocol's own text still says E2 "has not been run", because it is
  the frozen record of the rules. This file and the run log are the
  record of the results.
- `PILOT_PROTOCOL.md`, `CODEX_PILOT_REVIEW.md`, `ORCA_CHECKPOINT.md` and
  `CODEX_HANDOFF.md` are historical records.

## Recommended next work (none of it started)

1. **Preserve the evidence.** Keep this package with the frozen commit.
   Do not overwrite `e2_out/` or `confirm_e2_out/`.
2. **Prototype AIDS direct LM as an opt-in.** For example, an
   experimental `quaidsFit()` mode behind an explicit control flag, with
   the shipped default unchanged.
   - AIDS is the strongest case: 170/170, W = 0, even from a single
     Stone start.
   - It needs its own design review, a test plan, and a decision on how
     to report symmetry-constrained output (item 4).
3. **Investigate the 10 QUAIDS no-result seeds under a new exploratory
   protocol.** Every start is infeasible, inversion-failed or
   continuation-failed on these seeds. Candidate directions:
   - feasibility-restoring starts;
   - a different continuation path or schedule;
   - scale or conditioning of the covariance.

   This must be a new, separately prespecified study. Its results must
   not be used to reinterpret E2.
4. **Resolve symmetry before shipping QUAIDS (or AIDS) direct LM.**
   Choose in-objective symmetry (the protocol's E1 smoke gate is the
   ready-made design) or a projection covariance valid at direct
   endpoints.
5. **Benchmark speed separately.** Use a dedicated, controlled timing
   study, with analytic Jacobians if pursued. E2's wall times are not
   evidence of speed.
