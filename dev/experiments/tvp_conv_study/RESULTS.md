# TVP-AIDS convergence-correctness study (Phase C): results

_Run 2026-10-20, task `task_7acc920462f5`. Protocol FROZEN v1
(`PROTOCOL.md`), freeze hash
`50f99a53048b9d29b28a8e002b7b1ede8e9e62bd240e9d66cf91faeb2c03f33e`
(`FREEZE.sha256`), committed as `e61b0b9` before any grid fit. Every
stage ran from that clean HEAD. sslib v2.1.2 (`a32885b`), strict pin
passed at all four preflights. GAUSS 26.1.4, Windows 11, one GAUSS process
at a time._

## Outcome in one paragraph

All **120 attempted grid fits** completed without a GAUSS error, outer
timeout or stop rule. All 120 returned `mleRetcode = 0` and are
**CONV-GLOBAL**: within tolerance of the best log-likelihood found across
the three starts for their dataset. On every known-truth fit (114), LL at
the fitted Q is above LL at the true Q by at least 1.31 log-likelihood
units, so the true-Q necessary check passes everywhere. No CONV-LOCAL,
CONV-FAILC3, CAPPED, OTHER-RC, ERROR or TIMEOUT cell occurred. The P6
`n1 = 3` anchor, which hit the iteration and time caps in P6, converges
from all three starts with the 2000-iteration budget, after 726-1050
iterations. The study elapsed 66.2 minutes of the 360-minute budget. These
results cover only the prespecified DGPs and anchor. They are not a
convergence guarantee and make no speed claim.

## What ran

| stage | design | cells | result |
|---|---|---|---|
| A | `n1=2`; seeds 1001-1008; tobs 200 (oracle + static H) and 500 (oracle H); 3 starts | 72 | 72 CONV-GLOBAL |
| B | `n1=3`; seeds 1001-1006; tobs 200; oracle + static H; 3 starts | 36 | 36 CONV-GLOBAL |
| C | `n1=4`; seeds **1001-1002** (reduced, see below); tobs 200; oracle H; 3 starts | 6 | 6 CONV-GLOBAL |
| anchor | P6 example data, static H; `n1=3` 2000 iters, `n1=4` 500 iters; 3 starts each | 6 | 6 CONV-GLOBAL |

- **Seed reduction (frozen rule, applied mechanically).** After Stage A,
  elapsed time was 6.2 min and the worst-case projection for B+C+anchor
  was 367.6 min, above the 360-min budget. The rule cut `n1=4` to two
  seeds, which brought the projection to 312.5 min, so `n1=3` kept all
  six seeds. The worst case assumes every fit runs to its iteration cap.
  Actual fits stopped far earlier, so the reduction proved unnecessary in
  hindsight. It is still the frozen outcome: grid 120 of the planned 126.
  `n1=4` seeds 1003-1004 were never simulated.
- **Admissibility.** Every study seed was admissible (minimum share
  across all DGP cells 0.107). No reserve substitution occurred.
- **Determinism replicate (question 6): not triggered.** The protocol
  reruns the first CAPPED cell of each stage, and no cell capped. Other
  evidence, which is not the protocol test:
  - the three Stage C seed-1001 fits rerun after the interruption below
    gave byte-identical endpoints;
  - the anchor `n1=4` shipped-start fit reproduces P6's 93 iterations
    and `Qfit` sum 3.1185578095761257 exactly.

## Non-protocol interruption (recorded)

At about 06:19 Claude Code killed the background driver process because
the machine was critically low on memory. The in-flight job,
`C-n14-t200-s1002-Horacle-q1`, died with it, after Stage C had completed
3 of its 6 cells. This is not a protocol stop rule: no GAUSS error, outer
timeout or FAILC3 occurred.

The coordinator approved one resume:
- The three completed rows moved to `results/interrupted_C_attempt1.csv`.
  Their logs and the partial output were kept as `logs/INTERRUPTED_*`
  (local; logs are gitignored).
- `run_conv.ps1 -Stage ALL` was relaunched with the unchanged frozen
  harness and the original budget clock. Stages A and B were skipped as
  done, and Stage C reran from a fresh preflight.
- The rerun seed-1001 rows match the interrupted attempt exactly:
  retcode, iterations, Fct, every LL, gradient, boundary count and
  `Qfit` sum, with the log endpoint text byte-identical. Only wall times
  differ (25-58 s before, 33-53 s after).

## Results by question

1. **`mleRetcode == 0` share:** 120/120 in every combination of `n1`,
   tobs, start and H rule. Iterations:
   - `n1=2`: 7-48;
   - `n1=3`: 6-27;
   - `n1=4`: 11-48;
   - anchor `n1=3`: 726-1050;
   - anchor `n1=4`: 61-93.
2. **Agreement with the best-found LL:** every dataset's three starts
   agree. Largest gap below the best:
   - known-truth stages: 8.1e-11 total-LL units, about 4e-7 of the
     2e-4 tolerance;
   - anchors: 9.0e-5, 0.18 of the 5e-4 tolerance (anchor `n1=3`, shipped
     start).
3. **True-Q necessary check:** passes on all 114 known-truth fits. The
   smallest LL(Qfit) - LL(Qtrue) margins are 1.31 (A), 3.18 (B) and 8.59
   (C). This is necessary, not sufficient: it rules out an endpoint below
   the truth, not a better maximum elsewhere.
4. **Boundary (descriptive), counting elements below 1e-3 x reference:**
   - known-truth fits: 0-4 of `k_states` elements per fit. True Q here is
     1e-7/1e-6, so near-zero estimates for some elements are expected
     from one realization.
   - anchor `n1=3`: 12 of 12 elements are below 1e-3 x the shipped start
     (`Qfit` sums 1.32e-3 to 1.34e-3).
   - anchor `n1=4`: 14 of 18 are (sum 3.1185).
   - The anchor starts reach the same LL within tolerance but slightly
     different Q vectors. This is a flat likelihood near the boundary, not
     a distinct maximum. Q values are not validated here.
5. **P6 anchor:** `n1=3` converges from all three starts with a
   2000-iteration budget and agrees within 9.0e-5 LL. CMLMT still prints
   "The covariance of the parameters failed to invert" at those
   boundary endpoints (3/3). `n1=4` converges from all starts as in P6,
   with the same warning in 1 of 3 fits.
6. **Determinism of capped fits:** not tested; nothing capped (see above).

Full tables: `results/analysis.md` (generated) and `results/categories.csv`.

## Timing (indicative only, no speed claim)

- Fit wall times sum to 59.6 min.
- Largest CMLMT seconds per iteration:
  - `n1=2`: 0.23 at tobs 200 and 0.55 at tobs 500;
  - `n1=3`: 0.90 at tobs 200 and 0.81 on the anchor at tobs 500;
  - `n1=4`: 2.80 at tobs 200 and 2.49 on the anchor.
- Short fits inflate the per-iteration figure, because fixed per-fit
  costs inside CMLMT's elapsed time are spread over few iterations.
- The machine was under memory pressure during the run.

## Evidence files

- `results/summary.csv`: one row per attempted fit job (120).
- `results/categories.csv`, `results/analysis.md`: analyzer output,
  reproducible with `analyze_conv.ps1`.
- `results/run_log.txt`: driver narrative, including preflights,
  projections, the reduction and the interruption note.
- `results/run_state.json`: freeze HEAD, budget start and reductions.
- `results/interrupted_C_attempt1.csv`: the three preserved interrupted
  rows.
- `results/gen/*.e`: every generated job script. The only edit is that
  the extra blank line at the end of each file was removed so that
  `git diff --check` passes; nothing else in the scripts changed.
- `results/determinism.csv`: header only, since no replicate ran.
- Per-job and preflight logs (`results/logs/*.log`) are gitignored and
  kept locally, as in P6.

## Residual risks and limits

- This is one optimizer on the prespecified DGPs. True Q was scaled to
  0.01x the draft so that shares stay in (0,1), which is a low
  signal-to-noise regime close to the boundary. Larger true Q at these
  sample sizes was infeasible for valid shares and is untested.
- Three starts per dataset bound, but do not prove, globality. A better
  maximum that no start reaches cannot be ruled out.
- The `n1=4` evidence is 2 seeds (6 fits) plus the anchor.
- The determinism of iteration-capped fits remains open.
- The anchor's "covariance failed to invert" warnings at boundary
  endpoints are consistent with the existing scope: no Q standard errors
  are claimed.
- Wall-time projections were very conservative (worst case at the cap),
  which is why the frozen rule reduced `n1=4`.
