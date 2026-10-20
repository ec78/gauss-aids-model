# TVP-AIDS speed study, Phase S0: analytic-gradient feasibility and speed pilot

**Status: FROZEN v1 (2026-10-20, Orca task `task_382ca289c40d`).** Freeze
hash and file hashes: `FREEZE.sha256` (see "Freeze"). The frozen files are
committed locally before any timed comparison. Any later protocol or
harness change needs a documented amendment section here and a new freeze
hash.

This is the first step of Phase S in
`../tvp_conv_study/PROTOCOL.md` ("Phase S entry gate"). It is dev-only. It
changes no shipped source, public API, package metadata, version,
documentation, tag, SSLib install or Phase C evidence. The candidate lives
only in `s0_lib.src`, which this directory's jobs `#include` explicitly.

## Claims and non-claims

S0 asks two questions on a small, prespecified set of Phase C cells:

1. Is an analytic gradient of the shipped TVP-AIDS Q likelihood, supplied
   through sslib's existing hook, correct?
2. Does it reach the same endpoint as the shipped numerical-gradient path
   on fixed Phase C datasets and starts, and how do evaluation counts and
   wall time compare?

S0 makes no speed claim for the package. A per-`n1` timing classification
is reported only if every endpoint gate for that `n1` passes, and only for
these cells, this machine and this optimizer configuration. S0 recommends
nothing about changing shipped code. Q values are not validated (Phase C
scope).

## Candidate (no algorithm or API choice)

sslib v2.1.2 `ssFitTVP()` already supports `ssControl.gradproc`. With it,
each objective evaluation runs `kalmanFilterDiffuseTVPGrad` (the
sensitivity-equation filter; per sslib's docs it uses the same diffuse-phase
decisions as `kalmanFilterDiffuseTVP`, and it is validated against finite
differences in sslib's own tests) instead of the plain filter, and CMLMT
receives the gradient rather than differencing. In the shipped model only Q
depends on the parameters, and `Q_t = diag(q)` for every period, so
`dQ_t/dq_j = e_j e_j'` and every other derivative field is zero. sslib
applies the `positive_vars` chain rule (`q = p^2`) itself.

`s0_lib.src` therefore adds only: a gradproc that copies a prebuilt
constant `dQ` field; a counting wrapper around the shipped
`_quaidsTVPQUpdate`; and `_s0MLEFit`, a line-for-line copy of
`_quaidsTVPMLEFit` that sets the gradproc when `useGrad = 1`. Objective,
transform, diffuse initialization, caps, CMLMT algorithm and settings,
final filter pass and smoother are unchanged. The numerical arm
(`useGrad = 0`) must reproduce `quaidsTVPFit()` exactly (`s0_equiv.e`).

**Known side effects (by sslib design, not choices made here):**
- sslib's `_tvpLFNImpl` ignores CMLMT's `ind` request vector, so in the
  analytic arm every objective evaluation, line-search steps included,
  also computes the gradient.
- CMLMT's default algorithm switch starts with Newton (numerical Hessian)
  and moves to BFGS. With a gradient available, CMLMT differences the
  gradient for the Hessian (`cmlmthessg`) instead of differencing the
  function twice. This holds both for the iterations and for the final
  covariance, so the arms' `covPar` differ. Only endpoints (LL, Q, bFinal)
  are compared.

## Pre-freeze calibration (recorded for transparency)

This used non-study data only: the Stage 6 fixture (`n1=2`, tobs 120, seed
7, draft Q) and reserve seed 1101 at the frozen 0.01x true Q for `n1=3/4`,
tobs 200. No S0 pilot cell was fitted and no study dataset was evaluated
before the freeze. Raw logs (`calibration/*.log`) are gitignored. Tagged
lines are in `calibration/pilots.txt`.

- **Harness equivalence:** on the fixture, the numerical arm reproduces
  `quaidsTVPFit()` exactly: Qfit and bFinal differences 0, rc equal, 10
  iterations. The analytic arm also takes 10 iterations, with LL
  difference 1.6e-12, max |dQ| 6.3e-11 and max |dbFinal| 4.3e-9.
  Objective evaluations: 534 (N) vs 100 (A, 99 with gradient).
- **Derivative errors:** at the five point types below, the largest
  relative error at step `1e-4 p` was 4.1e-8, and about 100x that at step
  `1e-3 p`, as expected for central differences. The two filters' LL
  agreed exactly (relative difference 0).
- **Per-evaluation cost:** one gradient evaluation cost 5.3x (`n1=2`,
  tobs 120), 7.1x (`n1=3`) and 26x (`n1=4`) a plain evaluation. A Newton
  iteration with a numerically differenced Hessian costs on the order of
  `k^2` plain evaluations, against about `k` gradient evaluations with the
  analytic gradient. A large `n1=4` speedup is therefore not expected
  *a priori*. This is a prior, not a result.
- **Fit-template smoke run** (fixture, all six timed fits): both arms
  converged in 10 iterations, and both arms' repetitions were
  byte-identical. MLE times were 0.65-0.73 s (N) and 0.67-0.77 s (A).

## Frozen design

All data, H rules, starts, caps and true Q are Phase C's
(`../tvp_conv_study/PROTOCOL.md`). The job reuses
`../tvp_conv_study/conv_lib.src` unchanged. Study seed **1001** is
admissible in every Phase C cell (no reserve substitution there). `qscale`
is 0.01, `maxIters` 500 and `maxTime` 600 minutes.

| stage | `n1` | pilot cells (dataset, start multiplier) | derivative datasets |
|---|---|---|---|
| A | 2 | t200 oracle H x {1, 10, 0.1}; t200 static H x 1; t500 oracle H x 1 (5 cells) | t200 oracle, t200 static, t500 oracle |
| B | 3 | t200 oracle H x {1, 10, 0.1}; t200 static H x 1 (4 cells) | t200 oracle, t200 static |
| C | 4 | t200 oracle H x {1, 10, 0.1} (3 cells), **only if the inclusion rule passes** | t200 oracle |

The `n1=3/4` P6 anchor cells (726-1050 and 61-93 iterations, 2.5-12 minutes
per numerical fit) are excluded from S0 as too long for a pilot. They are
deferred to a later S phase if S0 warrants one.

### Derivative gate (`s0_deriv.template.e`)

For each derivative dataset, the gate uses five deterministic, strictly
feasible points in natural units:
- the shipped start `q_s = 0.1*mean(diag(H))`;
- `10 q_s`;
- `0.1 q_s`;
- the true Q;
- the true Q with odd-indexed elements times 100 (`ALTQ`).

At raw `p = sqrt(q)` the job compares:
- the analytic total gradient of exactly the objective CMLMT sees (sslib
  `_tvpLFNImpl`, gradient filter, chain rule), against
- central finite differences of the plain-filter objective with step
  `h_j = 1e-4 p_j`.

**Gate per point:** `max_j |gA_j - gFD_j| / max(||gFD||_inf, 1) <= 1e-5`
**and** `|LL_gradfilter - LL_plain| / max(|LL_plain|, 1) <= 1e-12`. Step
`1e-3 p_j` is reported only. The job also profiles, at `q_s`, the seconds
per plain evaluation, per gradient evaluation (5 each) and per Q update
(20). That profile is descriptive.

### Fit cells (`s0_fit.template.e`)

One job per cell. It does one untimed warm-up objective evaluation per
arm, then six timed fits on identical data and start in the fixed order
**N A | A N | N A**: N is the shipped numerical arm, A the analytic arm,
three repetitions each. For every fit it records rc, CMLMT iterations,
objective evaluations (updateproc calls in the MLE call, including
ssFitTVP's final filter pass), gradient evaluations (gradproc calls), MLE
seconds, whole-chain seconds (MLE plus smoother), Fct, LL at Qfit
(independent filter re-evaluation, `conv_lib.src: _convLL`), CMLMT's max
|gradient| and `sum(Qfit)`. Repetition 1's Qfit and bFinal are printed in
the job log. Time comes from `hsec`, which has 0.01 s resolution.

**Gates per cell (each halts the study on failure):**
- **E0 reproduction:** the numerical arm must reproduce the Phase C cell.
  That means the same rc, the same iteration count, `|LL_N - LL_C| <= tol`
  and the same Phase C category. Exact equality of LL and `sum(Qfit)` with
  `../tvp_conv_study/results/summary.csv` is reported.
- **E1:** `rc_A == rc_N`, and the Phase C category of A equals that of N.
  The category rules are `analyze_conv.ps1`'s, applied against the Phase C
  best-found LL of that dataset (all three Phase C starts) and LL(true Q).
- **E2:** `|LL_A - LL_N| <= tol`, with `tol = 1e-6 * nobs` (the frozen
  Phase C tolerance).
- **E3:** within each arm, all three repetitions give byte-identical
  Qfit, bFinal, rc and iterations.

**Reported flags (prespecified, not gates):**
- Q-agree: `max|QA - QN| / max|QN| <= 1e-3`;
- B-agree: `max|bFinalA - bFinalN| <= 1e-5`.

Q is weakly identified near the variance boundary (Phase C). A Q flag
failure is therefore disclosed, not treated as non-equivalence, as long as
E1/E2 pass.

**Timing summary:** the median and minimum MLE seconds per arm, and the
ratio `median_N / median_A` (above 1 means the analytic arm is faster).
The frozen per-`n1` classification applies only if every cell of that
`n1` passes E0-E3:
- **analytic faster** if every ratio is `>= 1.10`;
- **analytic slower** if every ratio is `<= 1/1.10`;
- otherwise **no consistent difference**.

## Execution, budget and stop rules (`run_s0.ps1`)

Strictly sequential, one GAUSS process at a time, under
`tests/gauss26_cfg_override`. There is no parallel GAUSS: before every
job, the driver stops if any `tgauss.exe` is running or if an interactive
`gauss.exe` used more than 0.3 s of CPU in a 3 s window. (The owner's idle
IDE is allowed.) Outer per-job timeouts: 900 s for derivative jobs, and
900 / 1800 / 2400 s for `n1 = 2 / 3 / 4` fit jobs. The hard budget is
**150 minutes** from Stage A's preflight.

**Preflight before every stage:**
- the frozen hashes in `FREEZE.sha256` match;
- the Phase C harness's own `FREEZE.sha256` hashes match;
- Phase C `results/summary.csv` is unchanged (SHA-256 with CR removed,
  embedded in `run_s0.ps1`);
- HEAD equals the HEAD recorded at the first preflight, and the tree is
  clean except `results/`;
- no other GAUSS computation is running;
- `scripts/verify_sslib_pin.ps1 -Strict` passes (v2.1.2);
- `quaidstvpfit_test.e` reports ALL 30 CHECKS PASSED;
- `quaidstvp_loglik_xcheck_test.e` reports ALL 12 CHECKS PASSED;
- `s0_equiv.e` passes (the numerical arm reproduces `quaidsTVPFit()`
  exactly).

Each stage then runs its derivative jobs, then its fit cells.

**`n1 = 4` inclusion rule.** After Stage C's preflight and derivative job,
and before any C fit, each C job's time is projected as
`3 x 1.5 x (Phase C wall time of the cell)` plus
`3 x 1.5 x (Phase C iterations) x (k_states + 5) x (measured seconds per
gradient evaluation)`. C runs only if all three of these hold:
- elapsed time plus all projected C jobs is at most 150 minutes;
- every projected job is at most half its 2400 s outer timeout;
- free physical memory is at least 4 GB.

Otherwise C is recorded as excluded, which is not a failure. The same
per-job projection guards every cell against the budget.

**Stop rules, each of which halts the whole study:**
- any preflight failure;
- any GAUSS error (no `S0D_DONE`/`S0F_DONE`, or a GAUSS error code);
- any outer timeout;
- any derivative-gate failure;
- any E0, E1, E2 or E3 failure;
- a missing Phase C reference row;
- another GAUSS computation running;
- a job projected past the budget.

Settings are never changed and rerun. On a gate failure the owner (via
the Orca coordinator) is asked before anything else is done.

**Outputs** (`results/`):
- `deriv.csv` (one row per point), `profile.csv`, `runs.csv` (one row per
  timed fit) and `cells.csv` (one row per cell with every gate);
- `run_log.txt`, `run_state.json` and the generated jobs `gen/*.e`;
- `analysis.md`, from `analyze_s0.ps1`.

Per-job logs `logs/*.log` are gitignored and kept locally, as in Phase C.

## Caveats fixed in advance

- These are wall-clock timings on one Windows machine with other user
  processes present, which Phase C saw under memory pressure. Three
  repetitions bound, but do not remove, timing noise. The `hsec`
  resolution is 0.01 s.
- Iteration counts can differ between arms (the Hessian route differs), so
  time per fit, not per iteration, is the compared quantity.
- One seed per `n1` is pilot evidence only.

## Freeze

`run_s0.ps1 -Freeze` writes `FREEZE.sha256`, the SHA-256 of each frozen
file's text with CR removed. The frozen files are this file,
`s0_lib.src`, `s0_includes.inc`, `s0_equiv.e`, `s0_deriv.template.e`,
`s0_fit.template.e`, `run_s0.ps1` and `analyze_s0.ps1`. The freeze hash is
the SHA-256 of `FREEZE.sha256`, computed the same way.

Before the freeze commit, a driver smoke run (`-Smoke -Stage A
-ResultsDir <outside the repository>`) exercised the preflight, one
Stage 6 fixture derivative job and one fixture fit cell end to end. Smoke
mode skips the HEAD/clean check and E0, since no Phase C row exists for
the fixture. It is not study evidence.

Commands (from this directory):

```powershell
powershell -ExecutionPolicy Bypass -File .\run_s0.ps1 -DryRun
powershell -ExecutionPolicy Bypass -File .\run_s0.ps1 -Stage ALL
powershell -ExecutionPolicy Bypass -File .\analyze_s0.ps1
```
