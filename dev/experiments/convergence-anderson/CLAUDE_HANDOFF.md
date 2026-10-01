# Claude handoff: safeguarded Anderson + structured multi-start benchmark

Last updated: 2026-10-01

## Start here

This is experiment-only research for `gauss-aids-model`. Nothing here is
part of the installed package or CI surface, and no shipped file under
`src/` was changed. Read the repository-root `CLAUDE.md` and then
`SUMMARY.md` in this directory before editing or running GAUSS code.
Never commit or push without the user's explicit approval.

The prior handoff (`CODEX_HANDOFF.md`) asked for a methodology audit. That
audit and its recommended next step are now complete. This file records
the current state so the work is not accidentally repeated.

## Bottom line

The full 200-seed, 16-start benchmark has been completed for iterated AIDS
and QUAIDS. Correctly maximizing `homogCrit=-ln(det(Sigma))` selects the
same loose correct/wrong class as the synthetic oracle in 399 of 400
seed/model cases. Selection is therefore not the demonstrated
bottleneck; reaching a good endpoint is.

Structured multi-start improves the number of correct endpoints from 98
to 121 for AIDS and from 56 to 69 for QUAIDS relative to safeguarded
Anderson from Stone alone. It reduces total nonconvergence from 6 to 3
and 33 to 9. This is useful research evidence, but not a production-ready
fix: only 69/200 QUAIDS seeds meet the loose `recErr <= 1` threshold, and
only 16/200 meet `recErr <= .5`. A higher `converged==1` rate would still
hide many bad solutions on real data.

## Current artifacts

- `SUMMARY.md`: authoritative narrative, methodology corrections, exact
  result tables, and interpretation.
- `anderson_full_multistart_benchmark.e`: final 200-seed benchmark. It is
  configured with `nSeeds=200`, `numStarts=16`, and Anderson depth 8.
- `anderson_full_multistart_benchmark_output.txt`: complete final output;
  it ends with `anderson_full_multistart_benchmark.e: run complete.`
- `quaidsfit_anderson_prototype.src`: experiment-only estimator copy. Its
  optional `safeMode=1` path uses consecutive differences, rank-truncated
  SVD, residual-growth restart, step-size safeguard, and true
  fixed-point-residual stopping.
- `anderson_safe_solver_test.e`: isolated rank-deficient-history test plus
  a real QUAIDS seed-7 safe-fit test.
- `anderson_seed17_repro_output.txt`: evidence for the SVD-failure
  regression. QUAIDS seed 17 originally caused `svd2` to fail; the
  trapped rank-zero/plain-step fallback completes it cleanly.
- `anderson_faithfulness_check.e`: confirms `mDepth=0` remains exactly
  identical to shipped `quaidsFit()` for the tested case.
- `anderson_methodology_audit.e` and its output: 30-seed audit that fixed
  the old selector-direction error and compared Anderson/plain endpoint
  diversity.

The other scripts and outputs in this directory are historical evidence.
Do not cite the original minimized-`homogCrit` multi-start conclusion or
the old Anderson-basin-collapse hypothesis without the superseding caveats
in `SUMMARY.md`.

## Final benchmark results

| Result (of 200 seeds) | Iterated AIDS | QUAIDS |
|---|---:|---:|
| Stone-only safe AA: never / wrong / correct | 6 / 96 / 98 | 33 / 111 / 56 |
| 16-start max-criterion: never / wrong / correct | 3 / 76 / 121 | 9 / 122 / 69 |
| Oracle among converged: wrong / correct | 75 / 122 | 122 / 69 |
| Selector misses | 1 | 0 |
| Seeds with multiple endpoints | 71 | 119 |
| Converged starts / 3,200 | 3,091 | 2,639 |
| Any wrong endpoint beats truth criterion | 19 | 5 |
| Selected `recErr <= .5 / 1 / 2 / 5` | 38 / 121 / 155 / 171 | 16 / 69 / 88 / 91 |

The `recErr <= 1` split is intentionally loose. Always report the
threshold sensitivity beside the binary correct/wrong totals.

## Reproduction

Run from this directory:

```powershell
C:\gauss26\tgauss.exe -b -x anderson_safe_solver_test.e
C:\gauss26\tgauss.exe -b -x anderson_faithfulness_check.e
cmd /c "C:\gauss26\tgauss.exe -b -x anderson_full_multistart_benchmark.e > anderson_full_multistart_benchmark_output.txt 2>&1"
```

The full command performs 6,400 fits and should be run as a background or
long-timeout job. Smoke-test any changed benchmark with a few seeds
before restoring `nSeeds=200`. GAUSS batch execution may return a nonzero
shell status even when the script completes; require the explicit `run
complete` marker and inspect the output for GAUSS errors.

The experiment scripts use absolute include paths rooted at
`C:/Users/eclow/Documents/GitHub/gauss-aids-model`. Update those paths if
the repository is moved. GAUSS identifiers are case-insensitive; avoid
single-letter names that differ only by case and review the other
language-specific constraints in root `CLAUDE.md`.

## Important implementation details

`_andersonSvdMix()` intentionally returns rank zero if `svd2` fails or
produces an unusable leading singular value. The caller rejects and
clears that history, then uses a plain fixed-point step. Do not replace
this with a normal-equation inverse: that was the conditioning weakness
the safe path was introduced to remove.

The 16 start recipes are printed at the beginning of the benchmark output
and implemented in `_benchmarkStart()`. The winning-start counts show no
single recipe dominates enough to discard the portfolio without a proper
ablation. Random perturbations are generated after each deterministic
fixture seed, making the recorded sweep reproducible in the tested GAUSS
environment.

The raw-criterion replica differs from `qOut.homogCrit` by at most
`6.05e-5` (AIDS) and `4.37e-4` (QUAIDS) in the full run. Do not use
`symcCrit` as an interchangeable final-parameter criterion; the audit
found material discrepancies on bad endpoints.

## Recommended next task

Do not tune Anderson depth as the next primary task. Prototype a direct
nonlinear optimizer (LM/Gauss-Newton or the repository's `optmt`/`cmlmt`
facilities) on a clearly defined concentrated residual/GMM objective.
Sequence the work as follows:

1. Build a noiseless known-truth fixture and verify parameter/objective
   recovery before any Monte Carlo sweep.
2. Compare direct optimization with safe Anderson on the same starts and
   report both convergence and threshold-sensitive recovery.
3. Evaluate the objective on held-out draws from the same DGP; a training
   criterion beating the finite-sample truth is not an identification
   result.
4. Add local Hessian/Jacobian rank and conditioning diagnostics around
   competing endpoints.
5. Only discuss productionization if there is a real-data diagnostic that
   distinguishes trustworthy from converged-but-wrong solutions.

Keep all of that work in this experiment directory unless the user
explicitly approves a production change.
