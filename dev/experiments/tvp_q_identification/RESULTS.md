# TVP-AIDS Phase Q0: Q identification and boundary geometry -- results

_Run 2026-10-20 in two attempts.
- **Attempt 1.** Orca task `task_14e632bafa54`. Protocol FROZEN v1
  (`PROTOCOL.md`), freeze hash
  `7cecd9f9194513a972265f429c6ffc1820fd69dc370a5402e8a94262aeb69d38`,
  committed as `4a34c9d` before any study measurement. Stages A and B ran
  from that clean HEAD.
- **Amendment 1.** After measurement, this changed only the analyzer's
  number formatting.
- **Attempt 2.** Orca task `task_d62382cc106b`, an authorized resume
  recorded as Amendment 2. ANCHOR was rerun from its beginning, then the
  stage C gate, from clean HEAD `8f8ff96`. Final freeze hash:
  `dd9b17506dd20a41cb6faa13227a6d87b0a1920cd50cf262f94c17b6202f4efa`.
  Every GAUSS harness file and `run_q0.ps1` keeps its v1 hash.
- **Environment.** sslib v2.1.2 (`a32885b`); the strict pin passed at all
  four preflights. GAUSS 26, Windows 11, one GAUSS process at a time._

## Status: complete as designed

- **Measured:** stages **A** (`n1=2`, seeds 1001 and 1002), **B**
  (`n1=3`, seeds 1001 and 1002) and the **P6 `n1=3` anchor**. Each has
  three start fits, fixed-point geometry (special points, slices,
  Hessians) and every true profile.
- **Stage C (`n1=4`): excluded by the frozen gate,** applied by the
  driver. At 14:35:38 free memory was 1.55 GB, below the 4 GB rule; elapsed
  budget was 72.4 of 300 min.
- **No protocol stop rule fired** in either attempt.

## Interruption and resume (audit trail)

- **The kill.** At about 11:25, Claude Code's background-shell
  low-memory reaper killed attempt 1's driver during the first anchor
  profile job (`ANCHOR-ANCHOR-n13-t500-s204-prof-g1`). There was no GAUSS
  error, timeout, memory-cap breach or mismatch. A, B and the anchor's
  fits and geometry were done.
- **The authorization.** No coordinator answer came in that task.
  Attempt 1 was committed as partial (`067619f`). The owner then
  authorized the resume, and `PROTOCOL.md` Amendment 2 (commit `8f8ff96`)
  records:
  - the authorization;
  - the moving aside of attempt 1's ANCHOR rows and jobs to
    `results/interrupted_ANCHOR_attempt1/`;
  - unchanged harness hashes, cells, tolerances and profile method;
  - the hand edit of `run_state.json`: new `freeze_head`, and a
    `budget_start` that keeps the original 300-minute budget, counting
    the 34.9 minutes used before the kill;
  - an added stop rule: a second host-memory kill would end the study.
- **The rerun.** The driver ran detached (`Start-Process -WindowStyle
  Hidden`). It skipped A and B as done, ran a fresh strict preflight, and
  reran ANCHOR from its beginning.
- **Reproducibility.** The resumed ANCHOR fits and geometry are
  **byte-identical** to attempt 1: best, starts, elements, points, slices
  (203 rows) and all Hessian rows. The four profile rows the killed job
  had printed match the resumed rows exactly, apart from timing.
- **Which rows count.** Attempt-1 rows enter no final table.

Before any interpretation, a reminder of scope. Everything below is
description of LL geometry at these endpoints. None of it is a test, an
interval, a standard error or an identification proof (`PROTOCOL.md`,
"Claims and non-claims").

## Gates

- **Preflights:** every one passed (A, B, ANCHOR in attempt 1; ANCHOR in
  attempt 2). They checked:
  - frozen hashes, and Phase C's reused `conv_lib.src` and
    `conv_chain.inc` unchanged;
  - a clean recorded HEAD (`4a34c9d`, then `8f8ff96`) and
    `verify_sslib_pin.ps1 -Strict`;
  - Stage 6 30/30 and P5 12/12;
  - `q0_equiv.e`: the reduced fit with nothing fixed equals
    `quaidsTVPFit()` exactly, a fixed element comes back exactly, and
    `Fct` equals the LL.
- **Phase C reproduction:** all **15 of 15** final start fits (12 DGP and
  3 anchor) reproduce their Phase C rows exactly: rc, iterations,
  LL(Qfit) and sum(Qfit) text. That includes the anchor's 847, 1050 and
  726 iterations, which attempt 1 also reproduced. Every Qfit save/load
  round trip is exact.
- **Jobs (final tables):** 35 completed with no timeout, memory kill or
  GAUSS error.
  - The largest tgauss working set was 61 MB.
  - Free memory at job start was at least 1.01 GB, against the 1.0 GB
    gate.
  - The summed job wall time is 40.9 min, and the budget used is 72.4 of
    300 min.
- **Internal consistency** over all **138 profile fits** (120 DGP and 18
  anchor):
  - every fit returned rc 0, with none capped;
  - profile >= slice at every point (smallest `prof a - slice` +3.9e-7,
    at the anchor);
  - no profile above LL(qhat) + tol;
  - the two profile starts agree within 1.4e-10 LL on the DGP endpoints
    and within 3.5e-5 at the anchor (tol 5e-4);
  - each start-a point equals its kind-4 slice LL exactly.
- **Hessian step check:** stable at all five endpoints. The largest
  |eig(0.01) - eig(0.03)| is 9.4e-5.

## Results (analysis tables: `results/analysis.md`)

### Known-truth endpoints (A, B)

| endpoint | k | boundary | loss at true Q | loss at shipped start | loss at all-floor Q | min interior eig of -H (ln q) | interior condition |
|---|---|---|---|---|---|---|---|
| n1=2 s1001 | 7 | 1 | 5.37 | 32.1 | 15.0 | 0.118 | 12.2 |
| n1=2 s1002 | 7 | 0 | 3.18 | 26.0 | 80.5 | 0.721 | 2.8 |
| n1=3 s1001 | 12 | 2 | 7.11 | 46.1 | 82.1 | 0.055 | 53.5 |
| n1=3 s1002 | 12 | 1 | 3.18 | 35.7 | 178.6 | 0.488 | 5.5 |

Losses are LL(qhat) - LL(point) in total-LL units.

1. **Starts.** All three starts reach the same endpoint at every A/B
   endpoint: LL gaps of 1e-11 or less, and non-boundary elements within
   8e-7 in log10.
2. **Boundary elements are flat in their own coordinate.** Every boundary
   element (4 across these endpoints, with qhat/qRef from 6e-16 to 7e-13)
   loses at most 5e-12 when set to the floor. In ln q its curvature is 0,
   and each of these endpoints' flat Hessian directions (eig < 1e-3:
   1, 0, 2, 1) loads on exactly those boundary elements.
   - One **interior** element is also flat to the floor: `n1=3` s1001,
     element 4 (gamma), qhat/qRef 0.031, floor loss 0.041.
   - Raising the boundary elements to their true values together costs
     1.72, 2.88 and 1.96 LL.
   So these single realizations carry little information that separates
   a true variance of 1e-7 or 1e-6 from zero, element by element.
3. **Individual interior elements are weakly determined; the overall
   scale is not.**
   - Single-element slices: moving one interior element to its true value
     costs at most 1.6 LL, often < 0.1. Moving one element to 100x qRef
     costs 15-35 LL.
   - Interior log-curvatures (0.055-2.9): under the quadratic
     approximation, a factor-e change along the weakest interior
     direction costs 0.03-0.36 LL.
   - Joint scale: the ray `q = c qhat` is much sharper. log10 c = +-0.25
     (a factor 1.8) costs 0.56-3.06 LL.
4. **Profiles versus slices (cross-parameter tradeoffs).**
   - True profiles are only slightly better than fixed-path slices. The
     largest recovery is 2.76 LL: `n1=3` s1002 with gamma at the floor,
     where beta's Q sum rises 8.4x to compensate. All other recoveries
     are 1.12 LL or less, and most are below 0.6.
   - The visible tradeoff is between variance groups. With gamma pinned
     low, alpha and/or beta rise (alpha 3.4x and beta 2.9x at `n1=2`
     s1001), and with gamma pinned high they fall.
   - Interior curvature correlations are small (|c| <= 0.16).

   Grouped profile loss with a whole group at the floor:

   | endpoint | alpha | gamma | beta |
   |---|---|---|---|
   | n1=2 s1001 | 0.22 | 12.0 | 1.26 |
   | n1=2 s1002 | 21.8 | 27.1 | 22.2 |
   | n1=3 s1001 | 6.60 | 67.3 | 4.05 |
   | n1=3 s1002 | 35.7 | 118.0 | 11.7 |

   So whether a group can be removed varies by realization. In `n1=2`
   s1001 both the alpha (0.22) and beta (1.26) blocks cost less than the
   1.92 marker; elsewhere every group costs more than 4 LL. This describes
   the geometry and is not a test.
5. **True Q versus qhat.** True Q sits 3.2-7.1 LL below qhat. The best
   point on the common-shape ray `c x qRef` is the truth itself (c = 1)
   except at `n1=2` s1001 (log10 c = -0.5, loss 2.03). This matches Phase
   C's LL(Qfit) >= LL(Qtrue) and is the expected overfit of one
   realization.

### P6 `n1=3` anchor (example data, static H)

- **All 12 elements are boundary.** qhat/qRef ranges from 1e-17 to
  1.6e-5, and qRef is the shipped start (79.3 per element).
- **Static-coefficient corner.** Setting all of Q to the floor (1e-12 x
  qRef, effectively the static-coefficient model) loses only **0.112 LL**.
  Along the joint ray `c x qhat`, every c <= 1 loses <= 0.112, and c = 10
  loses 1.07.
- **Individual elements.** Each element flattens to the floor with
  |loss| <= 0.052. The largest gain on a single slice, -2.0e-4 (element 3),
  is within tol (5e-4). Raising any single element to the shipped-start
  level costs 26-429 LL, and the shipped start itself costs 1221 LL.
- **Curvature.** 10 of 12 eigenvalues of -H in ln q are flat (< 1e-3),
  each loading on boundary elements. The two non-flat ones are 0.050 and
  0.070. With no interior element, there is no interior block or
  correlation table.
- **Starts.** The three starts agree within 9.0e-5 LL, consistent with
  Phase C.
- **Grouped true profiles.** Group at 10^c x qRef; loss = LL(qhat) -
  profile. Start b ran at the floor only and took 96, 239 and 289
  iterations, all rc 0.

  | group | floor | 1e-4 | 1e-3 | 1e-2 | 1e-1 |
  |---|---|---|---|---|---|
  | alpha | -2.0e-4 (inside tol) | 0.72 | 3.45 | 13.0 | 45.9 |
  | gamma | 0.050 | 2.92 | 16.9 | 62.0 | 190.7 |
  | beta | 0.041 | 7.68 | 31.0 | 90.1 | 216.6 |

  - **At the floor.** Pinning any whole group there costs at most 0.05
    LL; for alpha the profile is 2.0e-4 above LL(qhat), within tol.
  - **Away from the floor.** Even 1e-4 x the shipped-start level
    (about 7.9e-3 per element) costs 0.7-7.7 LL, and 1e-1 costs 46-217 LL.
  - **Profiles versus slices.** Profiles improve on slices by at most
    0.14 LL.
  - **Tradeoffs.** Pinning one group high drives the other groups'
    variances toward zero: gamma or beta at 1e-2/1e-1 sends the other
    to about 1e-4 of its qhat sum. With alpha pinned, beta rises 2.4-4.3x.

On this dataset the fitted Q is therefore practically indistinguishable,
by LL, from the corner where Q is near zero. The data carry essentially no
LL information for time variation in any coefficient at this H. This is
one dataset and one H rule, and is no statement about TVP-AIDS in
general.

## Conclusions (descriptive, for these endpoints only)

1. On the known-truth DGPs, which have low signal (true Q at 1e-7/1e-6),
   the shipped optimizer finds the same endpoint from every start. 0-2
   elements per endpoint sit on the boundary; each is flat to the floor
   and accounts for that endpoint's flat curvature directions. Interior
   elements are determined only to roughly a factor of a few each, while
   the overall Q scale is much sharper.
2. Grouped true profiles track fixed-path slices closely (within 2.8 LL).
   Variance groups trade off against each other, mostly gamma against
   alpha/beta, but compensation recovers little LL. Whether a group can
   be pinned at zero cheaply varies by realization.
3. The P6 anchor endpoint is a full boundary corner: within 0.11 LL of
   all-zero Q, with 10 flat directions. The true group profiles are
   highest at or near zero for every group (<= 0.05 LL from the floor),
   and they fall steeply once a group's variance is raised above about
   1e-4 x the shipped start. Its slow, "covariance failed to invert" fits
   (Phase C) are consistent with this flat geometry and do not point to
   an optimizer defect.
4. The geometry supports Phase C's scope note: no Q standard errors or
   per-element inference. Any future boundary-aware Q inference (next
   boundary item) needs non-regular, boundary-aware machinery, not Wald or
   regular chi-square statements.

## Limits and residual risks

- **`n1=4` not measured.** The frozen memory gate excluded it
  (1.55 GB < 4 GB). The anchor's profile grid is coarse (five levels, one
  start except at the floor).
- **Small, purposive sample.** Four DGP datasets (two seeds per `n1`),
  one H rule (oracle), tobs 200, and one anchor. No static-H DGP cells, no
  tobs 500, no `n1 >= 4`.
- **Local and grid-based geometry.** Slices and profiles are on fixed
  grids. Better points off the grid, or another maximum no start reaches,
  cannot be ruled out. Finite-difference curvature in ln q is local and is
  zero by construction for boundary elements. The floor (1e-12 x qRef) is
  a numerical stand-in for zero.
- **Thresholds.** The 0.1, 1.92 and 1e-3 thresholds are descriptive
  markers chosen before measurement, not calibrated tests.
- **Machine state.** Memory pressure (1.0-1.8 GB free) caused attempt
  1's interruption and the C exclusion. The resume reproduced attempt 1
  byte for byte. Timing is not a result here.

## Evidence files

- `results/fits.csv`, `best.csv`, `starts.csv`, `elements.csv`,
  `points.csv`, `slices.csv`, `heig.csv`, `hvec.csv`, `hdiag.csv`,
  `hint.csv`, `corr.csv`, `profiles.csv`, `profile_q.csv` and `jobs.csv`;
- the analyzer outputs `results/analysis.md`, `endpoint_summary.csv`,
  `element_summary.csv` and `profile_summary.csv` (reproduce with
  `analyze_q0.ps1`). The only edit to `analysis.md` is that the
  analyzer's extra blank line at the end was removed so that
  `git diff --check` passes;
- `results/interrupted_ANCHOR_attempt1/`: attempt 1's ANCHOR rows under
  the same CSV names, `profiles_partial_prof_g1.csv` (the killed job's
  rows) and `gen/` (attempt 1's ANCHOR jobs);
- `run_log.txt` (driver narrative plus two marked manual notes on the
  interruption and the resume), `run_state.json` and `gen/*.e` (every
  final generated job);
- `calibration/pilots.txt` (the pre-freeze calibration);
- gitignored and kept locally: `results/logs/*.log` (including
  `INTERRUPTED_attempt1_*`) and `results/work/*.fmt`.
