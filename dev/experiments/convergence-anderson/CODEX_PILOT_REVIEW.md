# Independent Codex review: matched 30-seed convergence pilot

Date: 2026-10-06  
Target: `e90a5686c625c3c648ad43f40dea56aea577e0ce`  
Reviewer scope: committed pilot artifacts and relevant documentation; no solver sweep rerun

## Verdict

**Technical pilot verdict: PASS.** The frozen protocol, harness, saved Stage A gate, 2,880-cell Stage B grid, metric bridge, selector, aggregate calculations, resumability evidence, and canonical generated outputs are internally consistent and reproducible from the committed cell files. I found no defect that changes a pilot result or blocks using these 30-seed results for the explicitly limited exploratory purpose stated in the documentation.

**Documentation verdict: ACTION REQUIRED (one material repository-state finding).** Several committed handoff/status statements still describe the pilot as uncommitted and identify an older `master` tip. This does not alter the numerical results, but it is materially misleading for the next operator and should be corrected before treating the handoff as current.

## Material finding

### F1 - committed handoff/status text says the committed pilot is uncommitted

Severity: material documentation defect; no numerical impact.

At review time, `git rev-parse HEAD` returned `e90a5686c625c3c648ad43f40dea56aea577e0ce`, and `git show --oneline --decorate e90a568` showed the commit at `HEAD`, `master`, `origin/master`, and `origin/HEAD`. In conflict with that state:

- `CLAUDE_HANDOFF.md:19-21` says the hardened checkpoint is `a776717` and "Everything from the pilot is uncommitted"; line 95 labels all artifacts uncommitted.
- `PROJECT_STATUS.md:140-141` says the pilot work is uncommitted; `PROJECT_STATUS.md:1065-1073` says `master` is at `a776717` with the matched pilot in the working tree. Its earlier "latest completed work" text at lines 61-70 also still points to `7e98720`.
- `SUMMARY.md:1` labels the experiment "uncommitted".

Action: update those repository-state statements to identify `e90a568` as the committed/pushed matched-pilot checkpoint and remove the obsolete working-tree instructions. Preserve the accurate experiment-only/not-shipped qualification.

## Audit results by required area

### 1. Frozen protocol and harness - PASS

- `PILOT_PROTOCOL.md:9-22` prespecifies the DGP, seeds, models, 16 starts, three methods, and method controls. `pilot_run_m0.e:21-27` and `pilot_run_m1.e:21-27` instantiate exactly 30 seeds x 16 starts x 3 methods for one model each, totaling 2,880 cells.
- `pilot_harness_core.src:31-41` matches the frozen control values. Lines 114-153 reproduce the synthetic DGP, homogeneous Stone fit, dimensions, truth mapping, scale, and start construction. Lines 186-266 dispatch the three methods with those controls.
- Every returned endpoint is independently evaluated at homotopy 1 (`pilot_harness_core.src:268-277`), and the stored `conv` flag combines solver-specific genuine convergence with exact-objective feasibility.
- The continuation predicate requires reached homotopy 1, successful stage zero, all accepted stages converged, and all accepted stages feasible (`pilot_harness_core.src:243-245`; `pilot_lib.src:145-246`).
- Scope is correctly experiment-only: `pilot_lib.src:4-9`; none of these files is wired into shipped package sources.

### 2. Harness-validation gate and saved evidence - PASS

- Part 1 tests metric recovery/layout/start representation before including the harness (`pilot_validate.e:42-239`), then runs the 288-cell validation grid (`pilot_validate.e:241-244`).
- Saved GAUSS output records `PART 1 (V1/V2/V3a): ALL 27 CHECKS PASSED` at `pilot_validation/pilot_validate_console.txt:144`, followed by 288 cells and an explicit completion marker at lines 435-436.
- The committed analyzer output ends `VALIDATION: ALL 58 CHECKS PASSED` (`pilot_validation/pilot_analyze_validate_output.txt:65`). An independent read-only rerun of `pilot_analyze.py --validate <absolute pilot_validation>` also returned 58/58.
- The full committed grid's saved output ends `VALIDATION: ALL 377 CHECKS PASSED` (`pilot_out/pilot_analyze_validate_output.txt:437`). An independent read-only rerun returned the same 377/377, including all 60 Anderson benchmark reproductions and all four Phase 3 smoke references.
- The two validation-only tolerance changes are disclosed in `PILOT_PROTOCOL.md:96-112`, are relative to printed 8-significant-digit references, and do not alter any cell definition or re-solve a cell.

### 3. Metric bridge and aggregate calculations - PASS

- `_pilotErrs` uses columns `1:n1`, excludes only the final u row from `rawStructErr`, and includes lambda for QUAIDS (`pilot_lib.src:109-135`). The synthetic V2 gate tests the u, lambda, and ignored-column layout.
- `_pilotRecoverBS` documents and implements the copied `quaidsFit()` post-estimation chain and the direct-endpoint `(bIdx, bCoef) = (bD, bD)` bridge (`pilot_lib.src:249-280,282-475`). V1 compares it to `quaidsFit()`/depth-0 Anderson at matching representations. The direct generalization uses the full residual cross-product (`pilot_lib.src:267-270,373-377`).
- No operationally selected cell and no converged Stone-only cell has `bsOk != 1`; the three full-grid QUAIDS cold-LM bS-chain failures therefore do not contaminate a selected or Stone result. Symmetry fallbacks are surfaced in the report rather than hidden.
- Independent in-memory recomputation from `cells_m0.csv` and `cells_m1.csv` produced 2,880 cells and 180 model/seed/method portfolios. Primary correct counts were AIDS 19/27/26 and QUAIDS 10/23/23; raw correct counts were AIDS 22/30/30 and QUAIDS 13/29/29 (Anderson/cold LM/continuation), matching `pilot_report.md:11-16` and the project summaries.
- The same recomputation matched the supporting claims: QUAIDS converged starts 372/379/476; QUAIDS Stone-only raw-correct 10/23/27; QUAIDS continuation had 2 multi-basin seeds (one basin on 28/30); median 16-start wall time was 8.237 s for cold LM and 98.101 s for continuation. The 178/180 selected-vs-oracle class agreement follows from the two recorded selection failures.

### 4. Equivalent starts and RNG - PASS

- `_pilotBenchmarkStart` is arithmetic-equivalent to `_benchmarkStart`, including the 16 ordered recipes (`pilot_lib.src:26-88`; `anderson_full_multistart_benchmark.e:60-118`). `_pilotAllStarts` resets exactly `rndseed 100000*modelFlag + 100*seed` and draws in recipe order (`pilot_lib.src:90-107`; benchmark lines 212-225).
- Anderson receives the full `ng x n` matrix, while direct methods pack the same matrix's free columns `1:n1` (`pilot_harness_core.src:152-166`). The protocol accurately discloses the only representation difference (`PILOT_PROTOCOL.md:24-31`).
- V3a proves pack/scale/unscale/unpack round trips and deterministic redraws. V3b independently reproduces the saved 200-seed benchmark line for every one of the 60 model/seed cases in the full pilot validation, establishing the claimed portfolio/RNG equivalence on the audited grid.

### 5. Selector and objective direction - PASS

- `critExact` is `-ln(det(residCov))`, so larger is better (`PILOT_PROTOCOL.md:40`; `pilot_harness_core.src:268-276`).
- `select_operational` sorts by start, excludes non-converged or missing-objective candidates, and replaces the incumbent only for a strictly larger criterion (`pilot_analyze.py:88-99`), giving the lowest start on exact ties.
- The independent `--selftest` rerun passed all 8 V5/dedup checks: maximize direction, sign-reversal guard, eligibility exclusion, tie behavior, all-failed behavior, NaN exclusion, oracle minimization, and basin clustering.
- Oracle selection is kept diagnostic and never feeds the operational selector (`pilot_analyze.py:102-111,193-207`).

### 6. Resumability and determinism - PASS

- The harness reads existing per-model keys, skips complete cells, appends each finished cell immediately, and adds the key to the in-memory completed set (`pilot_harness_core.src:16-23,80-112,168-175,298-301`).
- Saved V6 evidence shows: first run wrote 6; rerun wrote 0 and left the file byte-identical; after dropping two tail rows, resume wrote exactly 2; keys remained unique; recomputed fields were identical except wall time (`pilot_resume_test_output.txt:1-7`).
- Saved determinism evidence reports 288 common validation cells and 0 differences excluding wall time (`pilot_out/determinism_check.txt:1`).

### 7. Canonical outputs and documentation claims - PASS for outputs/numerics; F1 for repository-state prose

- Physical line counts are 1,440 each for `pilot_out/cells_m0.csv` and `cells_m1.csv`, 144 each for the two validation cell files, and 6 for the resume fixture. The independent validator confirmed unique complete keys, row widths, required counters, and the expected 2,880-cell grid.
- Both Stage B consoles end with `cells run 1440`, `cells skipped 0`, and explicit script completion markers (`pilot_out/pilot_run_m0_console.txt:1451-1452`; `pilot_out/pilot_run_m1_console.txt:1451-1452`).
- I copied only the two cell files and `columns.txt` to an isolated system-temp `pilot_out`, invoked `write_report('pilot_out')`, and compared SHA-256 hashes. The regenerated `pilot_report.md` and `pilot_seed_table.csv` were byte-identical to the committed canonical outputs; the temp tree was then verified under the temp root and removed (`TEMP_REMOVED=True`).
- A separate DataFrame comparison of the recomputed 180 rows to the saved seed table passed at `1e-12`; an exact in-memory-vs-parsed comparison differed only by ordinary decimal CSV round-trip ulps.
- Numeric claims in `CLAUDE_HANDOFF.md:25-79`, `PROJECT_STATUS.md:143-197`, and the matched-pilot section of `SUMMARY.md` agree with the regenerated report and saved diagnostics. The documents correctly limit inference to 30 seeds and label pooled results post hoc. Repository-state prose fails only as described in F1.

### 8. Repository boundaries, shipped source, and 200-seed execution - PASS

- `git diff-tree --no-commit-id --name-only -r e90a568` returned 38 changed paths. Filtering returned `SRC_CHANGED_COUNT=0` and no changed 200-seed/full-benchmark output path.
- The changed paths are `PROJECT_STATUS.md` plus experiment protocol, harness, analysis, diagnostics, saved cells/logs, and documentation under `dev/experiments/convergence-anderson/`. No shipped `src/` file changed.
- This review ran no GAUSS solver job, no 30-seed sweep, and no 200-seed sweep. The prior 200-seed benchmark was read only as committed reference evidence. No new 200-seed execution or artifact update is present in the target commit.
- During review, only `CODEX_PILOT_REVIEW.md` and `ORCA_CHECKPOINT.md` were edited. No commit, push, publish, install, or `PATH` change was performed.

## Commands used for decisive checks

The exact command chronology, results, recoverable sandbox failures, and next command are maintained in `ORCA_CHECKPOINT.md`. Decisive read-only checks included:

```powershell
git -C 'C:\Users\eclow\Documents\GitHub\gauss-aids-model' rev-parse HEAD
git -C 'C:\Users\eclow\Documents\GitHub\gauss-aids-model' show --stat --oneline --decorate --no-renames e90a568
git -C 'C:\Users\eclow\Documents\GitHub\gauss-aids-model' diff-tree --no-commit-id --name-only -r e90a568
python '<absolute>\pilot_analyze.py' --selftest
python '<absolute>\pilot_analyze.py' --validate '<absolute>\pilot_validation'
python '<absolute>\pilot_analyze.py' --validate '<absolute>\pilot_out'
```

No cell was resolved or recomputed by a solver during this audit.
