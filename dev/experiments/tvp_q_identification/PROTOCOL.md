# TVP-AIDS Phase Q0: Q identification and boundary geometry

**Status: FROZEN v1 (2026-10-20, Orca task `task_14e632bafa54`), with
post-measurement Amendment 1 (analyzer formatting only; see the end of
this file).** Freeze hash and file hashes: `FREEZE.sha256` (see "Freeze"). The frozen files are
committed locally before any study measurement. Any later protocol or
harness change needs a documented amendment section here and a new freeze
hash.

This is item 1 of the "Boundary order after this study" list in
`../tvp_conv_study/PROTOCOL.md`. It is dev-only. It changes no shipped
source, public API, package metadata, version, documentation, tag, SSLib
install, or Phase C / Phase S0 evidence. The profile helper lives only in
`q0_lib.src`, which this directory's jobs `#include` explicitly.

## Claims and non-claims

Q0 describes the shape of the shipped TVP-AIDS Q log-likelihood around a
small, prespecified set of endpoints: how much LL is lost moving Q away
from the fitted value, along which directions it is flat, which elements
sit on the positivity boundary, and how the other elements compensate when
one variance group is pinned.

Q0 does **not** produce, and its results must not be read as:
- a test, p-value, confidence interval or standard error for any Q
  element;
- a regular chi-square likelihood-ratio statement. Q elements at or near
  zero sit on the boundary of the parameter space, where the usual
  chi-square reference does not apply. Any threshold below is a
  **descriptive scale marker only**;
- a proof that Q is or is not identified. A flat direction found here is
  evidence of weak information in these data along that direction, not a
  proof of non-identification, and a curved direction is not a proof of
  identification;
- a statement about any DGP, `n1`, sample length or dataset other than
  those listed.

Q values are compared with the true Q where it is known, as description,
not as an accuracy gate.

## Evidence reviewed before design (no study measurement)

- **Phase C** (`../tvp_conv_study/RESULTS.md`): 120/120 grid fits
  CONV-GLOBAL; LL(Qfit) - LL(Qtrue) >= 1.31 on all known-truth fits;
  0-4 boundary elements per known-truth fit; the P6 `n1=3` anchor
  converges in 726-1050 iterations with **12 of 12** elements below
  1e-3 x the shipped start, its starts agree within 9.0e-5 LL but reach
  slightly different Q vectors, and CMLMT reports "covariance of the
  parameters failed to invert" there.
- **Phase S0** (`../tvp_speed_study/RESULTS.md`): the numerical-gradient
  path is deterministic and reproduces Phase C bit for bit; `n1=4` fits
  were excluded at 2 GB free memory.
- **P6** (`../tvp_p6_scale/RESULTS.md`): the `n1=3` anchor's Q was still
  moving when time caps ended earlier fits; at `n1=4`, 9 of 18 elements
  are below 1e-9.

## Pre-freeze calibration (non-study data only)

Calibration used **only** reserve seed 1101 (`n1=2` and `n1=3`, tobs 200,
oracle H, the frozen Phase C true Q) and non-study example data
`quaidsExampleData(500, 205)` (`n1=3`, static H, the anchor's data rule
with a different seed) for the anchor-like case. No Q0 endpoint dataset
(seeds 1001/1002 or example seed 204) was simulated, fitted or evaluated
before the freeze. A driver smoke run (reserve seed 1101, results outside
the repository) exercised every job type and the analyzer. Tagged lines:
`calibration/pilots.txt`; raw logs are gitignored.

1. **Reduced-parameter fit equals the shipped fit.** With no element
   fixed, `_q0MLEFitReduced` reproduces `_quaidsTVPMLEFit` exactly
   (Qfit difference 0, same retcode and iterations) at `n1=2` and `n1=3`.
2. **Profile fits behave.** Fixed elements come back exactly; CMLMT's
   `Fct` equals the filter LL at the rebuilt Q exactly; the two profile
   starts agree to 3e-12 (`n1=2`), 2e-11 (`n1=3`) and 1.4e-7
   (anchor-like); profile >= slice in every case; each profile fit took
   0.3-7.8 s.
3. **Floor.** LL at `c x qRef` (all elements) and at `max(qhat, c x qRef)`
   on the ladder c = 1e-2 ... 1e-14. Rule: the floor is the largest
   ladder value at which flooring qhat changes LL by at most 1e-6 in all
   three calibration endpoints. Changes at 1e-10: 5.5e-10, 5.4e-10 and
   **2.0e-5** (anchor-like); at 1e-12: 5.5e-12, 5.0e-12, 1.6e-7. **Frozen
   floor: 1e-12 x qRef** (`PHI = -12`). No filter failure occurred down
   to 1e-14.
4. **Hessian step.** Central differences in ln q at steps 0.1, 0.03,
   0.01, 0.003, 0.001. Interior eigenvalues agree to about 1e-5 relative
   between 0.03 and 0.003. Boundary-direction eigenvalues are 0 up to
   noise that grows as the step shrinks (|eig| up to 2e-7 at 0.01 and
   9e-6 at 0.001, anchor-like). **Frozen: step 0.01, with 0.03 as the
   reported reliability check.**
5. **Cost.** One LL evaluation: 4.1 ms (`n1=2`), 5.0 ms (`n1=3`), 10.5 ms
   (anchor-like, tobs 500). A two-step Hessian costs under 10 s.
6. **Implementation notes found in calibration.** `H` and `h` collide in
   GAUSS (case-insensitive), so the step is `stp`; `fn` is reserved;
   GAUSS `save` truncates a file name at its first dot, so start files
   use `q0p1` for 0.1.

### Why a reduced parameter vector, not CMLMT `active`

CMLMT's `active` switch turns fixed parameters into linear equality
constraints and moves CMLMT onto its constrained (SQP) path, a different
algorithm. Q0 instead removes the fixed elements from the CMLMT parameter
vector and inserts them in the updateproc (`_q0QUpdateReduced`). The
objective, `positive_vars` transform, diffuse initialization, CMLMT
algorithm and settings and caps are unchanged; with nothing fixed it is
bit-identical to the shipped fit (gate `q0_equiv.e`). This is the standard
concentrated-likelihood construction and makes no new algorithm or public
API choice.

## Frozen design

**Endpoints** (the Phase C cells, reproduced exactly):

| stage | endpoint | data | H | fit cap |
|---|---|---|---|---|
| A | `n1=2`, tobs 200, seeds 1001 and 1002 | Phase C DGP, true Q 0.01x draft | oracle | 500 iters |
| B | `n1=3`, tobs 200, seeds 1001 and 1002 | same | oracle | 500 iters |
| ANCHOR | P6 anchor `n1=3`: `quaidsExampleData(500, 204)`, first 4 goods | example data | static (example 14 rule) | 2000 iters |
| C | `n1=4`, tobs 200, seed 1001 (gated, see below) | Phase C DGP | oracle | 500 iters |

Seeds 1001/1002 are the first two Phase C study seeds, chosen by position,
not by outcome. Static-H DGP cells and the P6 `n1=4` anchor are not
included.

**Per endpoint:**

1. **Start fits.** Three fits, from 1x, 10x and 0.1x the shipped rule
   `q0 = 0.1 mean(diag(H))`, through Phase C's own
   `conv_chain.inc`/`conv_lib.src` (reused unchanged). Each must reproduce
   its Phase C row exactly (rc, iterations, LL(Qfit) and sum(Qfit) as
   printed text) or the study stops. **qhat** is the start endpoint with
   the largest LL (ties: earliest start).
2. **Reference and boundary.** qRef is the true Q for DGP endpoints and
   the shipped start for the anchor. An element is **boundary** if
   qhat_j < 1e-3 qRef_j (Phase C's rule). The floor is 1e-12 qRef.
3. **Special points** (LL and loss = LL(qhat) - LL): (1) qhat; (2) true
   Q (DGP); (3) shipped start; (4) all elements at the floor; (5)
   max(qhat, floor); (6) boundary elements raised to qRef.
4. **Fixed-path slices** (no reoptimization; every other element stays at
   qhat). Kind 1: joint ray q = c qhat, log10 c in {-3, -2, -1, -0.5,
   -0.25, 0, 0.25, 0.5, 1, 2, 3}. Kind 2: common-shape ray q = c qRef.
   Kind 3: each element j alone at c qRef_j. Kind 4: each variance group
   (alpha, gamma, beta blocks) at a common c qRef. Kinds 2-4 use log10 c
   in {-12 (floor), -4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2}.
5. **Local curvature.** Gradient and Hessian of LL in ln q at qhat by
   central differences, steps 0.01 and 0.03. Report all eigenvalues and
   eigenvectors of -H; the interior (non-boundary) block's eigenvalue
   range, condition number and curvature correlations
   `c_ij = -H_ij / sqrt(H_ii H_jj)` (cross-parameter tradeoffs). In ln q a
   boundary element has near-zero curvature by construction; this is
   reported, not hidden.
6. **True profiles.** For each group g and level c, group g is fixed at
   c qRef_g and every other element is reoptimized by `_q0MLEFitReduced`.
   Start a = qhat with group g replaced (the kind-4 slice point, so the
   profile can be compared with the slice); start b = the shipped start
   with group g replaced. The profile value is the larger LL. Start
   elements <= 0 would be raised to the floor (counted). Levels and caps:

   | stage | levels log10 c | start b at | profile cap |
   |---|---|---|---|
   | A, B | -12, -2, -1, 0, 1 | all levels | 500 iters |
   | ANCHOR | -12, -4, -3, -2, -1 | -12 only | 500 iters |
   | C | -12, 0 | none | 300 iters |

   A profile fit that hits its cap (rc 2) is kept and labelled: its LL is
   a lower bound on the profile, so its loss is an upper bound. It is not
   a stop rule.

## Analysis (`analyze_q0.ps1`, descriptive)

Tolerance `tol = 1e-6 x tobs` total-LL units (Phase C). Reported per
endpoint and in `results/endpoint_summary.csv`,
`element_summary.csv` and `profile_summary.csv`:
- start sensitivity: rc, iterations, LL gap to best, boundary count, and
  the largest |log10 q| difference from qhat on elements that are
  non-boundary in both;
- losses at the special points;
- per element: loss at the floor and at qRef, minimum slice loss on the
  grid (a negative value below -tol means a single-coordinate move
  improves LL, i.e. qhat is not a coordinate-wise maximum on that grid),
  and a class: **flat to floor** (floor loss < 0.1), **shallow** (< 1.92),
  **floor costly** (>= 1.92);
- joint rays;
- slice versus profile loss per group and level, profile start agreement,
  `prof a - slice` (must be >= -tol), capped profiles, and any profile
  above LL(qhat) + tol (that would mean qhat is not the best point the
  study found, and is reported as such);
- tradeoffs: at each profile point, the other groups' Q sums relative to
  qhat;
- curvature: eigenvalues, flat directions (eig < 1e-3), clearly negative
  eigenvalues (< -1e-3), whether each flat direction's largest loading is
  a boundary element, step stability (max |eig(0.01) - eig(0.03)| <=
  1e-3 x max(1, max eig)), interior condition number and the largest
  curvature correlations.

**Descriptive markers.** 0.1 LL ("flat") and 1.92 LL (half of 3.84, the
chi-square(1) 0.95 point) are scale references for the reader. They are
not tests. At a boundary the regular chi-square reference is wrong, and
no inference is drawn from either marker.

## Execution, limits and stop rules (`run_q0.ps1`)

Strictly sequential, one GAUSS process at a time, under the
`tests/gauss26_cfg_override` config; stages A, B, ANCHOR, C in that order.

**Preflight before every stage:** frozen file hashes match `FREEZE.sha256`;
the reused Phase C files (`conv_lib.src`, `conv_chain.inc`) match Phase C's
own `FREEZE.sha256`; HEAD equals the HEAD recorded at the first preflight
and the working tree is clean except this study's `results/`;
`scripts/verify_sslib_pin.ps1 -Strict` passes (v2.1.2); TVP Stage 6
`quaidstvpfit_test.e` reports ALL 30 CHECKS PASSED; P5
`quaidstvp_loglik_xcheck_test.e` reports ALL 12 CHECKS PASSED; and
`q0_equiv.e` passes (reduced fit with nothing fixed equals
`quaidsTVPFit()` exactly; a fixed element is returned exactly; `Fct`
equals the LL; the fit does not lower LL).

**Before every job:** no other GAUSS computation is running and at least
1 GB of physical memory is free (wait up to 10 minutes, else stop); and
elapsed time plus the job's outer timeout must not pass the **300-minute**
budget (else stop). **During a job:** the tgauss working set is polled
every second, and the job is killed and the study stopped above **6 GB**.

**Outer timeouts (s):** fit 600 / 900 / 2400 / 1800 (A / B / ANCHOR / C);
geometry 600 / 900 / 1800 / 1800; each profile job 1800 / 2700 / 3600 /
2400.

**Stage C inclusion gate (not a stop rule):** stage C runs only if, at its
start, at least 4 GB of physical memory is free (the S0 rule) and elapsed
time plus its 100-minute worst-case projection fits the budget. The
projection is 3 fits x 1.5 x Phase C's largest `n1=4` wall time (81 s),
plus 5 min geometry, plus 3 profile jobs x 2 levels x 301 iterations x
2.8 s (Phase C's largest `n1=4` s/iteration), plus 5 min preflight,
rounded up. Otherwise C is recorded as excluded.

**Stop rules (halt the study; nothing later runs; settings are never
changed and rerun):** any preflight failure; GAUSS error; outer timeout;
memory cap; memory/other-GAUSS wait exhausted; budget; a start fit that
does not reproduce Phase C exactly; a Qfit save/load round trip that is
not exact.

**Outputs** (`results/`): committed: `fits.csv`, `best.csv`, `starts.csv`,
`elements.csv`, `points.csv`, `slices.csv`, `heig.csv`, `hvec.csv`,
`hdiag.csv`, `hint.csv`, `corr.csv`, `profiles.csv`, `profile_q.csv`,
`jobs.csv` (wall time, outer timeout, peak working set, free memory at
start), analyzer outputs, `run_log.txt`, `run_state.json`, `gen/*.e`.
Gitignored and kept locally: `logs/*.log` and `work/*.fmt` (binary Q
vectors passed between jobs).

## Freeze

`run_q0.ps1 -Freeze` writes `FREEZE.sha256`: SHA-256 of each frozen file's
text with CR removed -- this file, `q0_lib.src`, `q0_includes.inc`,
`q0_data.inc`, `q0_design.inc`, `q0_fit.template.e`, `q0_geom.template.e`,
`q0_prof.template.e`, `q0_equiv.e`, `run_q0.ps1`, `analyze_q0.ps1`. The
freeze hash is the SHA-256 of `FREEZE.sha256` computed the same way.

Commands (from this directory):

```powershell
powershell -ExecutionPolicy Bypass -File .\run_q0.ps1 -DryRun
powershell -ExecutionPolicy Bypass -File .\run_q0.ps1 -Stage ALL
powershell -ExecutionPolicy Bypass -File .\analyze_q0.ps1
```

## Amendment 1 (2026-10-20, after measurement stopped; formatting only)

**What changed.** `analyze_q0.ps1`'s number formatter `F()` now prints
`n/a` for a missing (`$null`) or infinite value, as it already did for NaN.
Nothing else changed: no constant, threshold, classification or table, and
no GAUSS file or the driver.

**Why.** The frozen analyzer crashed (`You cannot call a method on a
null-valued expression`) on the anchor endpoint, which by design has no
true-Q point (code 2). It would also have printed `Infinity` for an endpoint
with no profile rows. This was found when the analyzer was first run on
real anchor rows, after measurement. It does not touch any measured value.

**Freeze.** v1's freeze hash
(`7cecd9f9194513a972265f429c6ffc1820fd69dc370a5402e8a94262aeb69d38`,
commit `4a34c9d`) governed every measured job. `FREEZE.sha256` was then
regenerated for the amended analyzer, and the new hash is recorded in
`RESULTS.md`. All measurement ran at HEAD `4a34c9d` under the v1 hashes.

## Amendment 2 (2026-10-20): authorized resume of the interrupted ANCHOR stage

**Authorization.** The owner, through the Orca coordinator (task
`task_d62382cc106b`), authorized resuming the interrupted Phase Q0 run
recorded at `067619f`, with no change to its scientific settings.

**What happened.** At about 11:25, Claude Code's background-shell
low-memory reaper killed the driver during the first ANCHOR profile job
(`ANCHOR-ANCHOR-n13-t500-s204-prof-g1`). There was no GAUSS error, outer
timeout, memory-cap breach or Phase C mismatch, so no protocol stop rule
fired. Stages A and B were complete. ANCHOR's three fits and its geometry
job were complete; its profiles were not. Stage C was never reached.

**Preserving attempt 1.** Every ANCHOR row from attempt 1 moved out of the
main result CSVs into `results/interrupted_ANCHOR_attempt1/` under the
same file names. The killed job's printed profile rows are in
`profiles_partial_prof_g1.csv`, and attempt 1's generated jobs are in
`gen/`. Attempt 1's ANCHOR and preflight logs were renamed
`results/logs/INTERRUPTED_attempt1_*` (gitignored, kept locally). The A
and B rows are untouched. Attempt-1 rows are evidence of the interruption
only. They enter no final table, but the resumed fits and geometry are
compared with them.

**What is unchanged.** All design, cells, starts, caps, levels, floor,
Hessian steps, tolerances, thresholds, the profile method
(`_q0MLEFitReduced`), the limits and the stop rules. Every GAUSS file and
`run_q0.ps1` keeps its v1 hash from `4a34c9d`. These hashes are identical
in v1's and the current `FREEZE.sha256`:

| file | SHA-256 (CR-stripped) |
|---|---|
| q0_lib.src | 9c1207109000de2fab1d618e374ed5e7ca3d5fea8077d87321dabc751ecbcd97 |
| q0_includes.inc | 2b80902e1103184a44843d6592dbb355c9ec2be3a31d2a6deb5af5cc804b32a6 |
| q0_data.inc | c398b27e22ce42fef520173df68e6d4dcf021e44d82e2b41941137a8f481a7d1 |
| q0_design.inc | 9c4383abbe8dd8f87e82b893f04d734f7e660882a82e9b523bc5a307429ccda8 |
| q0_fit.template.e | c6e8374d9e79714f5a4e02f1453086d86564569dd0e48234658770e8c757657a |
| q0_geom.template.e | 28860553f42e4e40ab83d88aae59c78c456fb34ea87b823bbcd9f5e6fdc3191f |
| q0_prof.template.e | fe4d6e7bc51ea2634ffd7011e4669a3c8f52108fdacdecd0b6e972353f6df5c3 |
| q0_equiv.e | d2bafdfcf995c6c79b9aab5d02a377150a886aefce782f4ff5b4ec1c22224d57 |
| run_q0.ps1 | 74456073ae0438405bb383d08be4add6c53cd41f7a6ac28fc4e1bf9d9ad10be1 |

Only `PROTOCOL.md` (these amendments) and `analyze_q0.ps1` (Amendment 1's
formatting fix) differ from v1. `FREEZE.sha256` is regenerated for them.
The v1 freeze hash was
`7cecd9f9194513a972265f429c6ffc1820fd69dc370a5402e8a94262aeb69d38`.

**Driver state, edited by hand and recorded in `results/run_log.txt`.**
Two fields in `results/run_state.json` change; the driver code does not.
- `freeze_head` becomes the new clean HEAD: this amendment's commit. The
  driver's preflight requires HEAD to equal it and the tree to be clean
  outside `results/`.
- `budget_start` moves so that the **original 300-minute budget** still
  counts the 34.9 minutes used before the kill (10:50:49 to the last
  attempt-1 output at 11:25:43). The resumed run therefore has at most
  265.1 minutes left, and the per-job budget guard applies unchanged.

The edited state keeps `stages_done = [A, B]`, so A and B are skipped.

**Resume run.**
- A fresh strict preflight runs (all v1 gates), then the ANCHOR stage
  from its beginning: three fits, which must again reproduce Phase C
  exactly, then geometry and all three profile jobs. After that comes the
  unchanged stage C gate (>= 4 GB free and the runtime projection);
  otherwise C is recorded as excluded.
- The driver is launched detached through Windows `Start-Process
  -WindowStyle Hidden`, outside Claude Code's background shell. That is
  process management, not a study setting.
- One GAUSS process at a time.

**Stop rule added for the resume.** A second host-memory kill of the
driver ends the study as recorded, with no further resume. Every original
stop rule still applies.
