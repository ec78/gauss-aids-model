# Matched 30-seed pilot: frozen protocol

Written 2026-10-05, before any pilot or validation result was examined
(only the Part 1 metric checks of `pilot_validate.e` had run). Changes
after that point are listed under "Deviations" at the bottom, with the
reason and when each was made.

## Design

| Item | Value |
|---|---|
| DGP | `_quaidsSyntheticDGP(3000, seed, modelFlag, 1)` (`tests/quaidsfixtures.src`), as in the 200-seed benchmark |
| Seeds | Stage A validation: 1, 6, 9. Full pilot: 1-30 |
| Models | iterated AIDS (modelFlag 0) and QUAIDS (modelFlag 1), homogeneity imposed |
| Starts | The 200-seed benchmark's 16 recipes (`_benchmarkStart`), drawn after `rndseed 100000*modelFlag + 100*seed`, in recipe order, from the Stone/LA-AIDS fit (`maxiter=1`) |
| Methods | 1 = safeguarded Anderson; 2 = cold direct LM; 3 = hardened continuation + LM |

### Method settings (unchanged from earlier phases)

- **Anderson** is `_quaidsFitAnderson(..., weight 0, depth 8, residualStop 1, safeMode 1, rankTol 1e-10, residualGrowth 10, stepFactor 10)` with `maxiter=100` and `err=1e-4`. These are exactly the benchmark's settings.
- **Cold LM** is `_directLMHomotopyDiag` at homotopy 1. It uses `maxIter 500`, `gradTol 1e-3` and FD step `1e-6`, with the scale vector taken from the Stone fit.
- **Continuation** is `_pilotContinuation`, which copies `_phase3Continuation`. It runs stage zero at homotopy 0, then uses step .25, `minStep` .01 and 500 iterations per stage, with the same LM settings as cold LM.

### Start equivalence

Every method gets the same 16 matrices `bStart` (ng x n) for a seed.
- **Anderson** receives the full matrix as `aCtl.b0`, exactly as in the benchmark.
- **The direct methods** receive `vec(bStart[.,1:n1]) ./ scaleVec`.

So the direct methods' free parameters are exactly columns 1:n1 of the same matrix. Column n, the n-th share equation, is not a parameter of the direct objective. Anderson's fixed-point map does not read it either: Tb depends only on columns 1:n-1. It enters Anderson only through the first iteration's residual and mixing history. That is the only representation difference, and it is the benchmark's own convention.

## Metrics (per cell = model x seed x method x start)

- **rawStructErr** is `_pilotErrs` (a verbatim copy of `_phase3Errs`). It is the max |raw − true| over all reduced rows **except the u row**, on columns 1:n1. It **includes lambda** for QUAIDS. It is computed on `homogB` for Anderson and on the unpacked endpoint for the direct methods.
- **bsErr** is benchmark-compatible: `max|bS − trueParams|` over the full recovered (ng+1) x n matrix, **including** the u row.
  - For Anderson, `bS` = `qOut.bS`, the benchmark's own value.
  - For a direct endpoint bD, `bS` = `_pilotRecoverBS(bIdx=bD, bCoef=bD)`. That is quaidsFit()'s own post-processing chain (delta-method variance, then symmetry minimum-distance projection, then absolute-price recovery) applied as if bD were the final iterate.
  - `bsErrExclU` (bS without its u row) is reported too.
- **critExact** = −ln det(residCov) at **homotopy = 1**. It comes from `_directObjCoreHomotopy`, evaluated at every method's returned endpoint, including continuation runs that stopped before homotopy 1. **Larger is better.**

### Genuine convergence (eligibility for selection)

| Method | Required |
|---|---|
| Anderson | `qOut.converged == 1` and critExact feasible |
| Cold LM | retcode 0 and critExact feasible |
| Continuation | reached homotopy exactly 1, stage zero converged, every accepted stage had retcode 0, every accepted stage was feasible on re-evaluation, and critExact is feasible |

## Selection and classification (per model x seed x method)

- **Operational selection** takes the genuinely converged start with the **largest critExact**. Ties go to the lowest start index. A start that did not converge is never selected, whatever its objective.
- **Oracle-best** is the genuinely converged start with the smallest error. It is computed separately for bsErr and rawStructErr. **It is a diagnostic only**, never a deployable rule.
- **Stone-only** is start 1 alone, classified the same way. It shows what a single default start costs and delivers.
- **Classes:**
  - **never**: no genuinely converged start.
  - **correct**: the selected endpoint is within the cutoff.
  - **wrong**: the selected endpoint is outside the cutoff.
- **Primary classification** uses **bsErr ≤ 1**. This is the benchmark-compatible metric and cutoff, so the results can be read against the 200-seed table.
- **Secondary classification** uses rawStructErr ≤ 1. Any seed where the two disagree is listed in the anomaly table.
- **Sensitivity** is reported at cutoffs .5, 1, 2 and 5 for both metrics.

## Tolerances (stated before the run)

- **Basin deduplication** uses leader clustering in start order over genuinely converged endpoints. A new endpoint opens a new basin if its max |Δ raw structural coefficient| (all rows except u, columns 1:n1) to every existing representative exceeds **0.05**.
  - Why 0.05: the Phase 3 T4 diagnostic showed that stopping-tolerance differences at gradTol 1e-3 reach 0.018 in the structural rows and 0.077 in the u row, which is why u is excluded.
  - The benchmark's own count is also reported for comparison: vec(bS) with tolerance 1e-4.
- **Solution equivalence** across methods (for example, "continuation and cold LM selected the same solution") uses the same 0.05 rule. **Numerically identical** means ≤ 1e-6, as in the smoke gate.

## Validation gate (Stage A)

The gate must pass before seeds 1-30 run:

| Check | What it establishes |
|---|---|
| V1 | `_pilotRecoverBS` reproduces quaidsFit()'s homogB, bS, recErr and homogCrit to ≤ 1e-8 at the same (index, coefficient) pair, and `_quaidsFitAnderson(mDepth 0)` matches quaidsFit() exactly. A polished fixed point shows the (b,b) representation equals quaidsFit()'s own. |
| V2 | `_pilotErrs` layout on synthetic matrices: the u row is excluded, lambda is included, column n is ignored. |
| V3a | Start pack/unpack round trip, and the starts are deterministic. |
| V3b | The pilot's Anderson portfolio reproduces the benchmark's per-seed line for seeds 1, 6 and 9: converged starts, selected start, selected and oracle recErr, and distinct endpoints. This proves start and RNG equivalence. |
| V3c | The pilot's start-1 cold-LM and continuation cells reproduce `phase3_smoke_gate_output.txt` (AIDS seed 6, QUAIDS seed 9). |
| V4 | critExact equals the method's own returned objective for LM and continuation endpoints at homotopy 1. For Anderson the difference is small, and it is reported. |
| V5 | Selection unit tests on synthetic candidate lists: the largest critExact wins, a non-converged better objective is excluded, ties go to the lowest start, all-failed gives "never", and a sign-reversed list is detected. |
| V6 | Resume: rerunning adds no rows, and deleting the last rows recomputes exactly those cells with identical results apart from wall time. |
| V7 | Keys are unique and complete, and the counter columns are present. |

## Outputs

- `pilot_validation/` holds Stage A, and `pilot_out/` holds seeds 1-30.
- `cells_m0.csv` and `cells_m1.csv` contain one row per cell with no header. `columns.txt` names the columns.
- After the metadata, each row carries the raw endpoint `vec(bRaw[.,1:n1])` and `vec(bS)`, so selection, deduplication and the anomaly table can all be recomputed from disk.
- `pilot_analyze.py` writes `pilot_report.md` and `pilot_seed_table.csv`.

## Deviations

1. **V3c tolerance (validation only, 2026-10-05).** The first Stage A pass
   failed one V3c check: QUAIDS seed 9 cold LM structErr 4.0194095 vs the
   reference 4.0194095, with iterations and rejected steps identical. The
   absolute tolerance 5e-9 was tighter than the 8 significant digits the
   reference was printed with. It was changed to a relative 1e-7. No
   pilot definition changed.
2. **Observation, not a change.** On QUAIDS seed 1, Anderson selection
   by critExact picks start 12, while the benchmark's `homogCrit` picked
   start 7. Both are in the same basin (eight converged starts agree on
   the objective to ~1e-9, bsErr 82.43 for each), so the order is noise
   within Anderson's 1e-4 stopping rule. The selected error is unchanged.
3. **V3b tolerance (validation only, after the 30-seed run).** Rerunning
   `--validate` on `pilot_out` (all 30 seeds) failed 10 V3b checks.
   Every one was a QUAIDS selected or oracle recErr of magnitude 100-344
   that matched the benchmark in all 8 printed significant digits (for
   example 221.38468 vs 221.38468). The absolute 1e-6 tolerance was
   tighter than the rounding of the printed reference. It was changed to
   a relative 1e-7, the same fix as deviation 1. No pilot definition
   changed, and nothing was re-solved.
