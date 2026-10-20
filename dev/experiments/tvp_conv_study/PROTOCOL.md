# TVP-AIDS convergence and speed study protocol

**Status: FROZEN v1 (2026-10-20, task `task_7acc920462f5`).** Phase C
design approved by the owner (relayed by the Orca coordinator); the
true-Q scale was set by the coordinator-approved pre-freeze calibration
below. Freeze hash and file hashes: `FREEZE.sha256` (see "Freeze").
Any later protocol or harness change needs a documented amendment
section here and a new freeze hash.

This post-0.4.0 study addresses the next TVP boundary without changing
shipped source or public API. It first tests whether the existing optimizer
reaches the best available endpoint on known-truth data. Only after those
reference endpoints are frozen may a separate speed phase compare an
implementation change.

## Claims and non-claims

Phase C asks whether the shipped optimizer, unchanged, reaches the best
log-likelihood found across prespecified starts within a fixed iteration
budget. Results apply only to the prespecified DGPs and P6 anchor. They are
not a convergence guarantee or a general scale claim.

Phase S asks whether a candidate optimization change reduces runtime while
reaching an equivalent endpoint on every Phase C reference cell. No speed
claim is allowed before endpoint equivalence passes. `mleRetcode == 0` alone
is never treated as proof of a maximum. Phase C makes no speed claim and
changes no optimizer setting.

## Phase C questions

1. What share of known-truth fits return `mleRetcode == 0`, by `n1`, sample
   length, starting value and fixed-H rule?
2. Among those fits, which agree with the best log-likelihood found across
   the prespecified starts for the same dataset?
3. Does each returned endpoint satisfy the necessary check
   `LL(Qfit) >= LL(Qtrue) - tolerance`?
4. How often are fitted Q elements near the positivity boundary? This is
   descriptive, not a per-element accuracy gate.
5. Does the P6 `n1 = 3` example-data anchor converge under a larger,
   iteration-bound budget, and do its starts agree?
6. Are iteration-capped results deterministic?

## Pre-freeze design work (recorded for transparency)

All of this used only the reserve seeds 1101-1110, the Stage 6 test
fixture (tobs 120, seed 7) or capped smoke runs; no study seed (1001-1008)
was simulated or fitted before the freeze. Evidence: `calibration/`. As
in P6, raw GAUSS logs (`*.log`) are gitignored and kept locally; the
committed records are `calibration/screen_summary.csv` (every screened
seed) and `calibration/pilots.txt` (tagged pilot lines).

**1. Draft truth was infeasible.** With the draft true Q (alpha block
1e-5, gamma and beta blocks 1e-4) the random-walk states drive shares out
of (0,1). Reserve-seed screens (`calibration/gen_screen_q1_*.log`,
admissible / 10): `n1=2` tobs 200: 7; `n1=2` tobs 500: 2; `n1=3` tobs 200:
2; `n1=4` tobs 200: 0. An equal-share `state0` (alphas 1/n) gave 2 and 0
for `n1=3/4`, so no `state0` choice fixes it, and ten reserves cannot
replace the expected failures.

**2. Calibration rule (coordinator decision, option B).** Keep every
other factor; choose the largest common true-Q scale from the prespecified
ladder {0.1x, 0.01x} of the draft Q that makes every reserve fixture
admissible for `n1=2` at tobs 200 and 500 and `n1=3/4` at tobs 200.
Results (`calibration/gen_screen_q0.1_*`, `gen_screen_q0.01_*`):

| scale | n1=2 t200 | n1=2 t500 | n1=3 t200 | n1=4 t200 |
|---|---|---|---|---|
| 0.1x | 10/10 (min share .152) | 10/10 (.083) | 10/10 (.092) | **9/10** (seed 1104 min -.099) |
| 0.01x | 10/10 (.219) | 10/10 (.217) | 10/10 (.136) | 10/10 (.087) |

**Frozen: true-Q scale 0.01x** -- alpha block 1e-7, gamma and beta blocks
1e-6. H is unchanged (2.5e-5 * I), so the truth now has a much smaller
signal-to-noise ratio and sits closer to the positivity boundary than the
Stage 6 fixture.

**3. Harness pilots (not study evidence).** `calibration/pilot_equiv.log`:
`conv_equiv.e` reproduces `quaidsTVPFit()` on the Stage 6 fixture exactly
(Qfit and bFinal differences 0, retcode equal) and CMLMT's stored `Fct`
equals the total filter log-likelihood exactly (difference 0; Fct*nobs is
not it). `calibration/smoke_fit_mode1/2.*` exercised the fit template once
each (Stage 6 fixture with static H, and the anchor data capped at one
iteration).

## Frozen design

| Factor | Levels |
|---|---|
| `n1` | 2, 3, 4 (`k_states` 7, 12, 18); `n1 >= 5` deferred |
| `tobs` | 200 for all; 500 additionally for `n1 = 2` (oracle H only) |
| study seeds | `n1=2`: 1001-1008; `n1=3`: 1001-1006; `n1=4`: 1001-1004 |
| reserve seeds | 1101-1110, per (`n1`, tobs), in order, replacing a study seed whose DGP has any share of the n goods outside (0,1) |
| true Q | 0.01x draft: alpha block `1e-7`; gamma and beta blocks `1e-6` |
| true H | `2.5e-5 * eye(n1)` |
| H supplied to fit | oracle true H (all cells); static H (secondary, `n1=2/3`, tobs 200 only) |
| Q starts | shipped rule `q0 = 0.1*mean(diag(H))*ones`; 10x; 0.1x (H = the H supplied) |
| caps | `maxIters = 500`, `maxTime = 600` minutes (iteration cap binds) |
| P6 anchor | example-14 data `quaidsExampleData(500, 204)`, first `n1+1` goods, example-14 static H via `quaidsFit()`; `n1=3` 2000 iterations, `n1=4` 500 iterations; three starts each; `maxTime = 600` |
| grid | A 72 (`n1=2`: 8 seeds x [t200 x 2 H + t500 x 1 H] x 3 starts), B 36, C 12, anchor 6 = **126 fits** |
| hard budget | 360 minutes from the start of Stage A's preflight; staged A, B, C, anchor |

`state0` (layout `[alpha | gamma row-major upper triangle | beta]`,
`conv_lib.src: _convState0`):

```text
n1=2: .40 | .30 | .05 | -.02 | .04 | -.04 | .03          (Stage 6 fixture)
n1=3: .30 | .25 | .20 | .05 | -.02 | -.01 | .04 | -.01 | .04 | -.04 | .03 | .02
n1=4: .25 | .20 | .20 | .15 | .05 | -.02 | -.01 | -.01 | .04 | -.01 | -.01
      | .04 | -.01 | .03 | -.04 | .03 | .02 | -.01
```

Implied period-1 nth-good alpha: .30, .25, .20; gamma rows have a negative
off-diagonal pattern and positive own-price terms; implied nth-good beta
-.01 and 0 for `n1=3/4`.

**Static H** (`_convStaticH`): diagonal SSE/nobs of a static homogeneous
LA-AIDS fitted equation by equation by OLS of `w_i` on [1, relative log
prices, Stone-deflated log expenditure], using exactly the regressors
`quaidsTVPFit()` builds. (The DGP has exogenous expenditure, so no IV.)

**Data generator:** `tests/quaidsfixtures.src: _quaidsTVPRawSyntheticDGP`
(unchanged), seeded by the used seed.

## Per-fit record

Each fit is one `tgauss` job (`conv_fit.template.e`) that runs the private
chain `conv_chain.inc`, a line-for-line copy of `quaidsTVPFit()`'s
orchestration. Its tagged output and `results/summary.csv` row record: cell
identifiers (study and used seed); CMLMT retcode, iterations, minutes,
`Fct`, max absolute gradient (CMLMT's stored `Gradient`, transformed
parameters); LL at the fitted Q (independent filter re-evaluation, plus
the fit's own kfResults sum), at the true Q and at the start; tolerance;
boundary count; min/max share; full Q0, true Q, Qfit and bFinal in the job
log; wall time, CMLMT time and seconds per iteration; warnings; outer
timeout and GAUSS-error flags.

**Fct scaling (established, see pilot):** CMLMT's `Fct` is the total
log-likelihood, not a mean. All LL comparisons use totals.

**Endpoint tolerance:** `tol = 1e-6 * nobs` total-LL units (2e-4 at tobs
200, 5e-4 at tobs 500).

**Boundary count** (descriptive): number of Qfit elements below `1e-3`
times the reference Q (true Q for DGP cells, the shipped start for anchor
cells).

## Categories (`analyze_conv.ps1`)

A dataset group is (stage, `n1`, tobs, used seed, H rule); its best-found
LL is the largest LL(Qfit) among its first-run cells, whatever their
retcode. In order of precedence:

- `TIMEOUT`, `ERROR`: reported literally.
- `CONV-FAILC3`: rc 0 and LL(Qfit) < LL(Qtrue) - tol (DGP cells only).
- `CONV-LOCAL`: rc 0 and LL(Qfit) < best-found LL - tol.
- `CONV-GLOBAL`: other rc 0 (anchor cells: no true-Q check).
- `CAPPED`: rc 2.
- `OTHER-RC(rc)`: any other retcode.

Q recovery summaries are descriptive. A single realization does not
identify each Q element accurately, especially near a variance boundary.

## Determinism (question 6)

In each stage, the first cell that returns rc 2 is rerun once (replicate
`r2`) immediately, if the per-job budget guard allows; its retcode,
iteration count and everything from `CONV_FCT` to `CONV_DONE` (LLs, Q0,
true Q, Qfit, bFinal) must be byte-identical. Replicates are excluded from
categories.

## Execution, budget and stop rules (`run_conv.ps1`)

Outputs go to `results/`: `summary.csv` (one row per attempted fit job),
`categories.csv`, `analysis.md`, `determinism.csv`, `run_log.txt`,
`run_state.json` and the generated job scripts `gen/*.e` are committed;
per-job logs `logs/*.log` are gitignored and kept locally, as in P6.

Strictly sequential, one GAUSS process at a time, under the
`tests/gauss26_cfg_override` config. Outer per-job timeout: 900 s (`n1=2`),
1800 s (`n1=3`), 3600 s (`n1=4` and anchor cells).

**Preflight before every stage:** frozen file hashes match
`FREEZE.sha256`; HEAD equals the HEAD recorded at the first preflight and
the working tree is clean except `results/`;
`scripts/verify_sslib_pin.ps1 -Strict` passes (v2.1.2); TVP Stage 6
`quaidstvpfit_test.e` reports ALL 30 CHECKS PASSED; P5
`quaidstvp_loglik_xcheck_test.e` reports ALL 12 CHECKS PASSED; and
`conv_equiv.e` reproduces `quaidsTVPFit()`'s Qfit and bFinal exactly with
`Fct` equal to the total LL. A DGP stage then screens study and reserve
seeds for share admissibility and fixes the seed map.

**Stop rules -- each halts the whole study (no later stage runs):** any
preflight failure; any GAUSS error (no `CONV_DONE`, or a GAUSS error
code); any outer timeout; any `CONV-FAILC3` (flagged inside the job);
exhausted reserves; or a next job whose worst-case time would take the
elapsed budget past 360 minutes. Settings are never changed and rerun.

**Worst-case time of a cell:** overhead + (maxIters + 1) x rate(`n1`) x
tobs/500. Rate(`n1`) is the larger of the prior CMLMT s/iter at tobs 500
(`n1=2` 0.30, `n1=3` 0.75, `n1=4` 2.60; P6 maxima rounded up, `n1=2`
conservative) and the largest measured study value normalised to tobs 500.
Overhead is the larger of 30 s and the largest measured (wall - CMLMT)
seconds. A stage projection adds 120 s per preflight.

**Seed reduction (the only permitted design change):** after each stage,
if elapsed + projected worst case of the remaining stages exceeds 360
minutes, reduce `n1=4` to two seeds (1001-1002) first, re-project, then
`n1=3` to four seeds (1001-1004) if Stage B has not run. No other factor
changes. If still above budget, the per-job guard applies.

## Freeze

`run_conv.ps1 -Freeze` writes `FREEZE.sha256`: SHA-256 of each frozen
file's text with CR removed (so git line-ending conversion cannot change
it) -- this file, `conv_lib.src`, `conv_chain.inc`, `conv_includes.inc`,
`conv_fit.template.e`, `conv_screen.template.e`, `conv_equiv.e`,
`run_conv.ps1`, `analyze_conv.ps1`. The freeze hash is the SHA-256 of
`FREEZE.sha256` computed the same way. The frozen files are committed
locally before any grid fit. A post-commit harness smoke run
(`-SmokeCells`, two `n1=2` tobs-120 seed-7 cells, results outside the
repository) checks the driver end to end; it is not study evidence.

Commands (from this directory):

```powershell
powershell -ExecutionPolicy Bypass -File .\run_conv.ps1 -DryRun
powershell -ExecutionPolicy Bypass -File .\run_conv.ps1 -Stage ALL
powershell -ExecutionPolicy Bypass -File .\analyze_conv.ps1
```

## Phase S entry gate

Profile the Phase C harness before altering code: count likelihood
evaluations per optimizer iteration and measure filter versus Q-update time.
Numerical gradients are a plausible cause of per-iteration cost growth, but
that is currently an inference, not evidence. An analytic-gradient candidate
must first match numerical gradients at prespecified Q points, then preserve
Phase C categories and log-likelihood endpoints before timing is compared.

## Boundary order after this study

1. Q identification and behavior at variance boundaries.
2. Boundary-aware Q inference and standard errors.
3. Fixed non-diagonal H estimation scope, then separately estimated H.
4. IV/endogeneity treatment of total expenditure after an explicit
   methodology decision and validated inference machinery.
