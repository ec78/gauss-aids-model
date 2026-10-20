# TVP-AIDS speed study, Phase S0: results

_Run 2026-10-20, Orca task `task_382ca289c40d`. Protocol FROZEN v1
(`PROTOCOL.md`), freeze hash
`d8c2dc999ce5169993aaa8c8c7dfd3bb2c3f6c8f0059ee2fc0ce4a1de81c8652`
(`FREEZE.sha256`), committed as `c6fa049` before any study-cell timing.
Every stage ran from that clean HEAD. sslib v2.1.2 (`a32885b`), strict pin
passed at all three preflights. GAUSS 26, Windows 11, one GAUSS process at
a time (the owner's idle IDE was open and drew no CPU). Total elapsed time
5.3 of the 150-minute budget._

## Outcome in one paragraph

The candidate is sslib's existing `ssFitTVP` gradproc hook with a constant
diagonal `dQ`. It needs no algorithm or API choice and is **correct**.
- **Derivatives.** All 30 prespecified derivative points (6 datasets x 5
  points) pass. The largest analytic-vs-central-difference error is
  1.8e-8 against the 1e-5 gate. The gradient filter's LL equals the plain
  filter's LL exactly at every point.
- **Endpoints.** On all 9 pilot cells, the analytic arm reaches an
  **equivalent endpoint**: same retcode (0) and Phase C category
  (CONV-GLOBAL), and |LL_A - LL_N| <= 9.5e-12 against the 2e-4 or 5e-4
  tolerance. Qfit agrees within 3.2e-6 relative, bFinal within 6.5e-9.
  Both arms are byte-deterministic across repetitions.
- **Reproduction.** The numerical arm reproduced every Phase C cell
  exactly: rc, iterations, LL and `sum(Qfit)`, bit for bit.
- **Speed.** The analytic arm makes 4-8x fewer objective evaluations, but
  each one costs 5x (`n1=2`) to 8x (`n1=3`) a plain evaluation. Wall time
  shows **no consistent difference** under the frozen classification for
  either `n1=2` or `n1=3`. The median ratio N/A was below 1 (analytic
  slower) in 7 of 9 cells; range 0.66-1.06.
- **`n1=4`.** It passed preflight and the derivative gate, but the frozen
  inclusion rule **excluded** its fits: free physical memory was 2 GB,
  below 4 GB. Its profile (one gradient evaluation is 21x a plain
  evaluation, with k = 18) gives no reason to expect a gain.

No speed claim is made, and none would be justified from this pilot.

## What ran

| stage | design | derivative gate | endpoint gates E0-E3 | classification |
|---|---|---|---|---|
| A | `n1=2`, seed 1001: t200 oracle x {1,10,0.1}, t200 static x 1, t500 oracle x 1 | 15/15 points | 5/5 cells | no consistent difference (0.658-1.023) |
| B | `n1=3`, seed 1001: t200 oracle x {1,10,0.1}, t200 static x 1 | 10/10 points | 4/4 cells | no consistent difference (0.741-1.060) |
| C | `n1=4`, seed 1001: t200 oracle x {1,10,0.1} | 5/5 points | not run: excluded by the frozen memory rule | -- |

No stop rule fired. Projected C fit time was 25.4 min (within budget and
per-job limits); memory was the only failed inclusion condition. Per the
protocol, settings were not changed and nothing was rerun.

## Endpoint equivalence (repetition 1 of each arm)

| cell | cat C / N / A | iters N / A | LL_A - LL_N | tol | Qfit rel diff | bFinal abs diff | E0 (exact) | E1 | E2 | E3 |
|---|---|---|---|---|---|---|---|---|---|---|
| n12-t200-Horacle-q1 | GLOBAL x3 | 20 / 20 | 1.1e-12 | 2e-4 | 5.5e-7 | 7.4e-10 | pass (exact) | pass | pass | pass |
| n12-t200-Horacle-q10 | GLOBAL x3 | 21 / 22 | 6.1e-12 | 2e-4 | 3.1e-6 | 5.1e-9 | pass (exact) | pass | pass | pass |
| n12-t200-Horacle-q0.1 | GLOBAL x3 | 21 / 20 | 0 | 2e-4 | 1.3e-6 | 2.1e-9 | pass (exact) | pass | pass | pass |
| n12-t200-Hstatic-q1 | GLOBAL x3 | 14 / 20 | 6.8e-12 | 2e-4 | 1.1e-6 | 3.9e-9 | pass (exact) | pass | pass | pass |
| n12-t500-Horacle-q1 | GLOBAL x3 | 8 / 18 | -5.5e-12 | 5e-4 | 1.6e-6 | 2.5e-9 | pass (exact) | pass | pass | pass |
| n13-t200-Horacle-q1 | GLOBAL x3 | 20 / 24 | 9.1e-13 | 2e-4 | 1.6e-6 | 2.0e-9 | pass (exact) | pass | pass | pass |
| n13-t200-Horacle-q10 | GLOBAL x3 | 21 / 24 | -9.5e-12 | 2e-4 | 4.7e-7 | 2.0e-9 | pass (exact) | pass | pass | pass |
| n13-t200-Horacle-q0.1 | GLOBAL x3 | 21 / 20 | 2.7e-12 | 2e-4 | 6.1e-7 | 2.6e-9 | pass (exact) | pass | pass | pass |
| n13-t200-Hstatic-q1 | GLOBAL x3 | 15 / 15 | -9.1e-13 | 2e-4 | 1.2e-6 | 6.5e-9 | pass (exact) | pass | pass | pass |

All Q-agree (<= 1e-3) and B-agree (<= 1e-5) flags pass.

## Timing and evaluation counts (MLE call; median of 3 repetitions)

| cell | obj evals N | obj evals A (all with gradient) | median s N | median s A | ratio N/A |
|---|---|---|---|---|---|
| n12-t200-Horacle-q1 | 793 | 163 | 1.77 | 1.97 | 0.900 |
| n12-t200-Horacle-q10 | 857 | 199 | 1.91 | 2.28 | 0.841 |
| n12-t200-Horacle-q0.1 | 855 | 161 | 1.87 | 1.83 | 1.023 |
| n12-t200-Hstatic-q1 | 617 | 160 | 1.37 | 1.84 | 0.741 |
| n12-t500-Horacle-q1 | 470 | 137 | 2.60 | 3.95 | 0.658 |
| n13-t200-Horacle-q1 | 1666 | 299 | 5.60 | 7.55 | 0.741 |
| n13-t200-Horacle-q10 | 1781 | 294 | 6.00 | 7.65 | 0.784 |
| n13-t200-Horacle-q0.1 | 1775 | 211 | 6.07 | 5.73 | 1.060 |
| n13-t200-Hstatic-q1 | 1345 | 180 | 4.63 | 4.70 | 0.985 |

Within-arm repetition spread (max/min) was 1.01-1.09. Per-evaluation
profile at the shipped start, as seconds per evaluation (plain / gradient):
- `n1=2`, t200: 0.0024 / 0.012 (5.0x);
- `n1=2`, t500: 0.0056 / 0.029 (5.1x);
- `n1=3`, t200: 0.0034 / 0.027 (7.8x);
- `n1=4`, t200: 0.0078 / 0.167 (21x).

A Q update costs at most 0.2 ms, so the filter dominates evaluation cost.

**Why there is no gain.** With numerical gradients, CMLMT's Newton start
spends roughly `k^2` plain evaluations per Hessian. With the analytic
gradient it spends about `k` gradient evaluations instead (`cmlmthessg`).
sslib's `_tvpLFNImpl` also computes the gradient on every evaluation,
line-search steps included, because it ignores CMLMT's `ind` vector. At
`n1 = 2/3` the cost per gradient evaluation (5-8x a plain evaluation)
cancels the 4-8x cut in evaluations. Some cells also take more iterations
with the analytic gradient (e.g. t500: 8 vs 18), because the Hessian route
differs. The `n1=4` profile ratio (21x for k = 18) points the same way.
These are mechanism observations from the counts, not tested hypotheses.

Full generated tables: `results/analysis.md`.

## Caveats and residual risks

- One seed per `n1`, and 9 cells. This is pilot evidence on this machine,
  under memory pressure (2 GB free at the end), with wall-clock `hsec`
  timing at 0.01 s resolution.
- `n1=4` fits were not run (frozen memory rule), and the P6 anchors were
  excluded by design, so long and boundary fits (726-1050 iterations) are
  untested for speed. Phase C's anchor endpoints sit at the Q boundary,
  where the analytic arm's behaviour has not been checked.
- The derivative points are interior (every `p_j > 0`). The gradient at
  boundary fitted Q is not separately gate-checked. Endpoint equivalence
  covers interior fits with 0-4 near-boundary elements.
- The two arms' `covPar` differ (Hessian from the analytic gradient vs
  from function differences). It was not compared, and no Q standard
  errors are claimed.

## Implications for Phase S (owner decision; nothing done here)

Supplying the existing sslib analytic gradient, unchanged, is correct but
does not speed up TVP-AIDS fits at `n1 <= 3` under the current CMLMT
configuration. Any further speed candidate is a material algorithm,
optimizer or upstream choice and needs owner sign-off. Examples:
- an `ind`-aware objective in sslib, so function-only evaluations skip the
  gradient (upstream);
- a cheaper Q score, such as the disturbance-smoother/EM score formula
  (new derivative algorithm);
- a different CMLMT algorithm or Hessian strategy (optimizer setting).

## Evidence files

- `results/deriv.csv`, `results/profile.csv`, `results/runs.csv` (54 timed
  fits) and `results/cells.csv` (all gates per cell);
- `results/analysis.md` (generated);
- `results/run_log.txt`, `results/run_state.json` and `results/gen/*.e`
  (every generated job).
- Per-job logs (`results/logs/*.log`, with full Qfit and bFinal vectors)
  are gitignored and kept locally.
- `calibration/pilots.txt` records the pre-freeze calibration.
