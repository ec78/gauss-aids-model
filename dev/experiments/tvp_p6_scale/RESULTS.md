# TVP-AIDS P6: bounded `n1` scale re-measurement

_Run 2026-10-17 (task `task_ee20ce3a0a46`, Stage 2 of the TVP promotion
plan, `dev/TVP_PROMOTION_ROADMAP.md` step P6). Package source 0.4.0 at
`b55692b`. sslib v2.1.1 (`bf827ae`). `verify_sslib_pin.ps1 -Strict`
passed before and after the runs. GAUSS 26.1.4, Windows 11._

## Question

Earlier docs said that MLE convergence "degrades sharply above `n1=2`-`3`"
and that a 5-good, `n1=4` dataset made CMLMT "pathologically slow". That
claim was measured with an unbounded optimizer on sslib v1/v2.0. P6
re-measures `n1 = 3` and `n1 = 4` with the bounded optimizer (P3) on the
D0 baseline.

## Inputs (reproducible)

- `p6_scale_fit.template.e` + `p6_fit_body.inc`: this is
  `examples/14_tvp_aids_estimation.e`'s setup with only the number of goods
  changed. The setup is:
  - `quaidsExampleData(500, 204)`, the first `nSub = n1 + 1` goods, with
    shares renormalised;
  - `H = diagrv(eye(n1), qOut.homogSse[1:n1]/qOut.nobs)` from a static
    homogeneous LA-AIDS `quaidsFit()`, which converged in every job;
  - `q0 = 0.1*meanc(diag(H))*ones(k_states,1)`;
  - `smooth = 1`, and one `quaidsTVPFit()` call per job.
- `run_p6.ps1` writes `gen_n1<k>_<maxIters>x<maxTime>.e` (committed) and
  runs each one with `tgauss -b -x`. It sets `GAUSS26_CFG` to
  `tests/gauss26_cfg_override` and applies an **outer hard timeout of 900 s**
  per job (`taskkill /T` on expiry). Jobs run strictly one at a time.
- Commands:

  ```powershell
  cd dev\experiments\tvp_p6_scale
  powershell -ExecutionPolicy Bypass -Command "& .\run_p6.ps1 -N1 3,4 -MaxIters 500 -MaxTime 5 -Repeats 3 -HardTimeoutSeconds 900"
  powershell -ExecutionPolicy Bypass -Command "& .\run_p6.ps1 -N1 3,4 -MaxIters 100 -MaxTime 60 -Repeats 2 -HardTimeoutSeconds 900"
  ```

  Pass `-N1` through `-Command` as shown: with `-File`, `3,4` arrives as
  the single integer 34.
- Cap settings:
  - `500 x 5` (`maxIters = 500`, `maxTime = 5` min) is example 14's
    documented cap;
  - `100 x 60` makes the iteration cap the binding one, which isolates
    determinism from wall-clock timing.
- Labels: `n13` means `n1 = 3` (4 goods, `k_states = 12`); `n14` means
  `n1 = 4` (5 goods, `k_states = 18`).

## Results (`results/summary.csv`, full logs `results/*.log`)

| n1 | k_states | caps (iters x min) | repeats | `mleRetcode` | CMLMT iterations | CMLMT minutes | fit wall time (s) | outer timeout hit | warnings |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 12 | 500 x 5 | 3 | **11** (maxTime) in 3/3 | 457, 438, 437 | 5.018-5.036 | 301.1-302.2 | no | "The covariance of the parameters failed to invert" (3/3) |
| 3 | 12 | 100 x 60 | 2 | **2** (maxIters) in 2/2 | 101, 101 | 1.193-1.228 | 71.6-73.7 | no | same (2/2) |
| 4 | 18 | 500 x 5 | 3 | **0** (converged) in 3/3 | 93, 93, 93 | 3.519-3.622 | 211.2-217.4 | no | none |
| 4 | 18 | 100 x 60 | 2 | **0** (converged) in 2/2 | 93, 93 | 3.700-3.890 | 222.1-233.5 | no | none |

No job printed `function evaluated as complex` or any GAUSS error.
"Covariance failed to invert" is CMLMT's message for the non-converged
`n1 = 3` endpoint. Under `500 x 5`, CMLMT printed "covariance matrix not
computed".

**Determinism.** The `results/*.log` output from `P6_QFIT` through
`P6_DONE` (Qfit, its sum, bFinal and the `printQuaidsTVP()` report) was
compared byte for byte:
- `n1 = 4`: all 5 jobs are identical, across both cap settings
  (`Qfit` sum 3.1185578095761257).
- `n1 = 3` under `100 x 60`: both jobs are identical (`Qfit` sum
  2.5812298663137905).
- `n1 = 3` under `500 x 5`: all 3 jobs differ (`Qfit` sums 3.78e-3,
  5.54e-3, 5.61e-3). The binding cap is wall-clock time, so the stopping
  iteration (437-457) and the endpoint depend on machine timing.

**Rate.**
- `n1 = 3`: about 0.66-0.71 s per CMLMT iteration.
- `n1 = 4`: about 2.3-2.5 s per iteration.
- For reference, `n1 = 2` on the same data (example 14, Stage 1) converged
  in 23 iterations, about 11 s for the whole example.
- Wall times are from one machine with other sessions possibly active.
  Read them as indicative, not as a benchmark.

**`n1 = 4` endpoint (measured, not validated).** Nine of the 18 `Qfit`
elements are at or near the positivity boundary (below 1e-9; smallest
1.5e-15). Element 12 is 2.955. Nothing here validates these estimates:
there is no synthetic truth for this data, which is example data and not a
TVP DGP.

**`n1 = 3` trajectory.** The `Qfit` sum is 2.58 at the 101-iteration stop
and 3.8e-3 to 5.6e-3 at the 437-457-iteration stops. `Q` was therefore
still moving substantially when the time cap ended the fit.

## Conclusions (measured on this one dataset only)

1. Both `n1 = 3` and `n1 = 4` terminate within the documented caps. The
   largest overrun of the 5-minute `maxTime` was 0.036 CMLMT minutes
   (about 2 s, a few iterations at `n1 = 3`). No job needed the outer
   timeout. The old unbounded "pathological" hang does not occur
   with caps set.
2. `n1 = 4` converges (`mleRetcode = 0`) deterministically in 93 iterations
   and about 3.5-3.9 CMLMT minutes. This contradicts the old blanket
   statement that `n1 = 4` was pathologically slow on this dataset, which
   was measured on sslib before v2.1.0 with no caps.
3. `n1 = 3` does **not** converge within `500 x 5`: it stops with
   `mleRetcode = 11`, with a timing-dependent endpoint. The
   iteration-capped run stops with `mleRetcode = 2` at a deterministic
   endpoint. Whether it converges with a larger budget was not measured.
4. On this data, difficulty is not monotone in `n1`/`k_states`. Time per
   iteration grows with `k_states`, but whether a fit converges depends on
   the data subset. This is one dataset, one seed and one `q0` rule, so it
   supports no general scale guarantee in either direction.

## Not measured (out of P6 scope)

- larger budgets for `n1 = 3`;
- other seeds, `tobs`, or `q0` rules;
- `n1 >= 5`;
- the accuracy of any fitted `Q`;
- non-diagonal `H` at these sizes.
