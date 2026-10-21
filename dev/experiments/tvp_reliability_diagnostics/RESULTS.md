# TVP-AIDS Phase D0: reliability-diagnostic feasibility results

This is a dev-only retrospective analysis of committed Phase C/Q0 evidence.
It adds no public diagnostic, inference claim, or shipped source change.

## Outcome

All 5 endpoints were stable across the three starts tried, while
all 5 received `GEOMETRY-CAUTION` under the frozen descriptive rules.
The two facts coexist: optimizer endpoint agreement does not remove boundary/flat-geometry risk.
The current public result can report only the fit return code and direct scale checks;
objective agreement and geometry need extra computation.

## Endpoint table

| endpoint | k | max start LL gap / tol | max dlog10 Q | boundary | flat eig | floor < 1.92 | all-floor loss | interior cond | label |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| ANCHOR-n13-t500-s204 | 12 | 0.180195 | 0 | 12 | 10 | 12 | 0.111872 | n/a | GEOMETRY-CAUTION |
| n12-t200-s1001-Horacle | 7 | 2.50111e-08 | 4.79914e-07 | 1 | 1 | 6 | 14.986 | 12.1548 | GEOMETRY-CAUTION |
| n12-t200-s1002-Horacle | 7 | 2.50111e-08 | 5.14201e-07 | 0 | 0 | 1 | 80.5264 | 2.79093 | GEOMETRY-CAUTION |
| n13-t200-s1001-Horacle | 12 | 6.82121e-09 | 5.66852e-07 | 2 | 2 | 8 | 82.0474 | 53.5476 | GEOMETRY-CAUTION |
| n13-t200-s1002-Horacle | 12 | 3.41061e-08 | 7.72515e-07 | 1 | 1 | 2 | 178.597 | 5.51794 | GEOMETRY-CAUTION |

## Feasibility matrix

| diagnostic | one public fit | extra starts | private objective work | recommendation |
|---|---|---|---|---|
| CMLMT completion | yes: `mleRetcode` | no | no | keep as the mandatory hard status |
| Q positivity/finiteness and `Qfit/q0` scale | yes | no | no | cheap descriptive scale warning only; `q0` is caller-chosen |
| Endpoint agreement | no | yes, plus comparable LL | objective exposure is currently private | opt-in multi-start reliability layer |
| Boundary proximity against a defensible reference | only start-relative | sometimes | yes for Q0's truth/reference rule | do not call this identification from one fit |
| Score/gradient or fixed-point residual | no | no | yes | not available in `quaidsTVPOut`; investigate separately |
| Hessian rank/conditioning | no | no | yes, multiple evaluations | advanced diagnostic; boundary-aware interpretation required |
| Profile/floor loss | no | no | yes, many reduced fits | research/advanced mode, not a default fit warning |

## Findings

1. A zero `mleRetcode` and agreement from 0.1x/1x/10x starts are useful but incomplete. Every endpoint passed both and every endpoint still showed at least one frozen geometry caution.
2. The example-data anchor is the clearest case: all 12 Q elements are on Q0's boundary rule, 10 Hessian directions are flat, and the all-floor point loses only about 0.112 total LL, despite start agreement.
3. A single-fit public warning can honestly cover only solver status and caller-start-relative scale collapse. It must not claim identification, standard errors, or global optimality.
4. The next public-design candidate should be a separate result/procedure with explicit cost tiers: cheap single-fit checks, opt-in multi-start agreement, and advanced geometry diagnostics. Do not overload `mleRetcode`.
5. As a broader correctness cross-check, all 120 committed Phase C rows are `CONV-GLOBAL`. Phase S0 nevertheless found the available analytic gradient slower in 7 of 9 pilot cells, so advanced score/curvature work should not be assumed cheap.

## Limits

The set is five purposive endpoints, excludes `n1=4`, and reuses Q0's descriptive thresholds. It measures no false-positive/false-negative rate for warnings. The DGP reference Q is unavailable in applications, while a caller's `q0` is arbitrary. Therefore D0 supports an API feasibility decision, not calibrated warning thresholds for a release.

Regenerate with `python dev/experiments/tvp_reliability_diagnostics/analyze_d0.py`.
