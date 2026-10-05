# Claude handoff: Phase 3 hardened -- strengthened smoke gate passed; stopped before the 30-seed pilot

Last updated: 2026-10-05

## Start here

This directory holds experiment-only research for `gauss-aids-model`.
Nothing in it is part of the installed package or the CI surface, and
no shipped file under `src/` was changed.

1. Read the repository-root `CLAUDE.md`.
2. Read `SUMMARY.md`, in particular its "Phase 3 hardening (2026-10-05)"
   section. That section supersedes several numbers and claims in the
   "Phase 3" section above it.

Never commit or push without the user's explicit approval. None of
Phase 1, 1b, 3 or 3-hardening is committed yet. The user said they will
review this hardened checkpoint and authorize a commit before the
30-seed pilot begins.

## Bottom line

**The strengthened Phase 3 smoke gate passes cleanly in all 4 scenarios.**
The scenarios are AIDS seed 6 and QUAIDS seed 9, each from the Stone
start and a far/structured start. Each passes all of the following:
1. reaches homotopy exactly 1;
2. stage zero genuinely converged;
3. every accepted stage returned retcode 0;
4. every accepted stage is feasible on re-evaluation;
5. structErr ≤ 1 (excluding the u row);
6. no regression versus cold LM;
7. no regression versus safe Anderson (+0.05 tolerance).

Key numbers (structErr, all rows except u):
- **AIDS:** continuation 0.058, cold LM 0.058, safe Anderson 1.52 (not
  converged).
- **QUAIDS:** continuation 0.081 from both starts, cold LM 4.02 (Stone)
  and 12.29 (far), safe Anderson 2.16 (Stone) and 0.19 (far), both not
  converged.

Corrections from the 2026-10-02 write-up:
- **AIDS structural error was misreported as 0.724.** That figure was the
  u-row error, caused by a wrong row slice. The corrected value is
  0.058.
- **QUAIDS continuation endpoints are not "identical."** The two starts
  give endpoints that differ by 5.8e-9 in scaled θ and 7.6e-10 in raw
  structural coefficients. They are numerically equivalent, not
  bit-identical.
- **Stage zero changed no endpoint.** It changes endpoints by ≤ 1.5e-8
  versus the legacy schedule; only iteration counts change.

## Artifacts

- `direct_objective.src`: the LM solver body is now
  `_directLMHomotopyDiag()`.
  - It returns an extra `diagVec`: covariance-inversion failures,
    step-inversion failures, infeasible trials, rho rejections, accepted
    steps, and final lambda.
  - Inversion is fail-closed: the destination is reset before every
    `invpd()`, and the result is checked for shape and `ismiss`.
  - New retcode 5 means a covariance-inversion failure at the current
    point; `gnorm` is then returned as -1.
  - It has a diagnostic-only fault hook (`faultMode`).
  - `_directLMHomotopy()` and `_directLM()` keep their signatures and
    are bit-identical at `faultMode=0`.
- `phase3_smoke_gate.e` / `_output.txt`: the hardened gate.
  - The driver solves stage zero (homotopy 0) first.
  - The legacy schedule runs alongside it for comparison.
  - It uses corrected error rows and has a runtime layout assertion.
  - It reports explicit endpoint differences, threshold sensitivity at
    .5/1/2/5, and the strengthened verdict.
- `phase3_safeguard_diag.e` / `_output.txt`: fault-injected
  `invpd()` failures. **14 of 14 checks pass.** The T4 check was
  revised after a first-run failure, which turned out to be a
  stopping-tolerance effect; SUMMARY.md documents it.
- Phase 1b: `phase1b_audit.e` was rerun, and every numerical line
  matches the saved `phase1b_audit_output.txt`. Only the banner and
  timestamp lines differ. The saved file was not overwritten.

Reproduce from this directory:
```powershell
C:\gauss26\tgauss.exe -b -x phase3_smoke_gate.e
C:\gauss26\tgauss.exe -b -x phase3_safeguard_diag.e
C:\gauss26\tgauss.exe -b -x phase1b_audit.e
```

## What is established

On these two hard seeds, continuation plus LM reaches a correct
endpoint (structErr ≤ 0.5) from both starts. The endpoint does not
depend on the start, to within about 1e-8. Cold LM lands in wrong
QUAIDS basins. Single-start safe Anderson does not converge.

The `invpd()` safeguards fail closed, as demonstrated by fault
injection.

## What is NOT established

- Anything about success rates. The test covers one seed per model,
  chosen as hard cases.
- That the driver's step-halving, minimum-step and stage-zero-failure
  paths work on real data. Nothing in this gate triggered them.
- Comparability with the 200-seed benchmark's `recErr`. That metric is
  `max|qOut.bS - trueParams|` on the recovered full structural matrix,
  whereas this gate's structErr is computed on the raw reduced rows
  excluding u. The threshold value matches; the metric does not.
- That continuation beats Anderson's best. Anderson ran single-start
  here, not with its 16-start portfolio.

## Next entry point: the 30-seed pilot (only after the user approves a checkpoint commit)

The pilot is a matched sweep on identical seeds and starts. It compares:
- safe Anderson with its full 16-start portfolio;
- cold LM;
- hardened continuation plus LM (with stage zero).

It should report:
- never-converged, wrong and correct counts;
- structErr **and** the benchmark-comparable `bS`-based recErr at
  .5/1/2/5;
- paired seed-level transitions;
- counts of failed stages, step reductions and minimum-step exhaustion;
- wall-clock distributions.

Follow the stage-gate discipline: no 200-seed run without the pilot
first.

## Environment / gotchas

See the repo-root `CLAUDE.md`. Gotchas specific to this work:
- **Editing `/* ... */` doc comments.** Editing only part of a doc
  comment can leave a dangling `/*` or `*/`. That produces a cascading
  `G0008` error far from the edit, so check both ends of the comment.
- **Row slices of `bRaw`.** Derive them from the layout
  `alpha(nIntFull) | gamma(n1) | beta | [lambda] | u`, with
  `ng = nIntFull + n1 + 2 + (QUAIDS)`. The 2026-10-02 slice error came
  from an off-by-model-flag expression. Use `_phase3Errs` (in the gate)
  or `_directBlockErrors` rather than hand-writing ranges.
- **Capturing output.** `tgauss ... | Out-File -Encoding utf8 file`
  captures only the script's prints. The GAUSS banner goes to the
  console.
