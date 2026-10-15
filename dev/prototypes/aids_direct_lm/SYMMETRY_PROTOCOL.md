# Prespecified protocol: in-objective symmetry for the AIDS direct-LM prototype (v1)

## Status

- **v1, 2026-10-12. Approved for the private implementation stage; not
  frozen; nothing implemented or run.** The owner's instruction to begin
  the recommended symmetry work approves this document (U1), its
  thresholds (U6), and omits the optional small-sample panel (U7) and the
  real-data holdout (U8) from this first bounded implementation (§12).
  The freeze happens later, at §8 step 4, after stage D. No solve, sweep
  or new seed was touched to write or review this document.
- **Scope:** an experimental plan for the private prototype
  `quaidsdirect.src` (this directory). It does not change `src/`, the
  public API, `package.json`, the version or any shipped output.
- **v1 artifact-schema clarification (2026-10-12, before any stage D
  solve; coordinator-approved; no constant, gate, seed, estimator or scope
  change):** `seeds.csv` gains `critAbsErr` = |critFromFull(bFull) − crit|
  on the H, S, S1 and P rows, the input C0 needs. The harness writes the
  truth-dependent columns (`pfErr`, `rawStructErr`, `classPf`) as `.`;
  `sym_analyze.py` writes the completed `seeds.csv`, `results.json` and
  `report.md` to its own output directory, recomputes the selection
  independently from `cells.csv`/`raw.csv`, and returns INVALID on any
  disagreement. Because H is called unchanged, per-start time is not
  observable: `cells.csv` `wallSec` is `.` (not applicable), and per-arm
  time is in `seeds.csv`. The real stage also writes an intermediate
  `polish.csv` (R6), from which the analyzer builds `realdata.json`.
- **Implementation (2026-10-12, private, nothing run beyond fixtures):**
  `quaidsdirectsym.src`, `sym_harness.src`, `run_sym_stage.ps1`,
  `run_sym_checks.ps1`, `sym_analyze.py`, `sym_config.json`,
  `quaidsdirectsym_test.e`, `test_sym_analyze.py` and five
  `guard_error_cases/directlm_sym_*.e` scripts. I1-I8 pass. I9 (the full
  runner, `run_sym_checks.ps1 -Full`, which writes the stamp
  `run_sym_stage.ps1` requires before stage D) has not been run yet.
- **Independent pre-stage-D review (2026-10-13; no constant, gate, seed,
  estimator or §2-§9 change):** §2-§4 math, dimensions and GAUSS
  semantics, selection, schemas and seed guards were checked against the
  code; no estimator defect was found. Guard and tooling fixes only:
  `-Full` now implies `-WithEvidence` (I9 names the evidence test) and
  deletes any earlier stamp before running; `sym_analyze.py` refuses to
  overwrite analysis outputs or an existing freeze manifest, requires a
  clean `dev_out/results.json` before `manifest`, rejects an incomplete
  determinism re-solve, and reports asym/determinism failures in stage D
  as stop (revise) reasons; `run_sym_stage.ps1` refuses a second
  `conf_rest`; `truth.csv` rows are written in key (row, col) order. The
  I9 outcome is recorded in `PROJECT_STATUS.md` and by the stamp itself.
- **Stage D ran 2026-10-13** (outputs and execution logs kept under
  `symmetry/`; results in `PROJECT_STATUS.md`); its stop rule was not
  triggered. **Post-Stage-D wrapper repair (2026-10-13, before any freeze;
  no constant, gate, arm, seed, estimator, artifact or §2-§9 change):**
  `run_sym_stage.ps1`'s default `-OutRoot` failed under Windows PowerShell
  5.1 before any work (`$PSScriptRoot` is empty in an advanced script's
  parameter defaults), so Stage D was run with the same directory passed
  explicitly. The default now resolves after parameter binding (also in
  `run_sym_checks.ps1`), `run_sym_stage.ps1 -Preflight` prints the
  resolved paths and exits before any guard or solve, and
  `run_sym_checks.ps1` tests it with the default and an explicit
  `-OutRoot`. The repair is numerically irrelevant: Stage D was not rerun
  and none of its outputs changed. A new `-Full` stamp replaced the old
  one.
- **Not authorized here (deferred, each still an owner decision, §12):**
  a change to shipped `quaidsFit()`; θ_P as an S start; any public
  exposure; any new returned field on a shipped or public struct.
- **Relation to the frozen study.** `dev/experiments/convergence-anderson/`
  (protocol FROZEN v4 at `ad986fc`, E2 evidence at `158313a`) is read
  only. This protocol does not reinterpret E2, reuse its held-out seeds
  for a new decision, or modify any file its manifest covers.
- **Question.** Is in-objective symmetry (arm **S**) good enough to be
  the preferred path for an experimental, AIDS-only, point-estimate
  MVP? It is compared with the current homogeneity-only endpoint
  (arm **H**) and, only as a comparator, with a minimum-distance
  projection that is defensible at direct endpoints (arm **P**).

## 1. Design audit (inputs checked for this document)

| ID | Finding | Source | Consequence |
|---|---|---|---|
| F1 | The frozen protocol names E1 ("symmetry-constrained objective smoke gate", "design unchanged from v2") but contains no specification. v2 was never committed: the protocol file has a single commit (`ad986fc`), and no other committed file specifies E1. | `CORRECTNESS_VALIDATION_PROTOCOL.md` §10, §12; `git log --follow` | This protocol is self-contained. Its development stage (§8, stage D) takes the smoke-gate role. Nothing here claims to be the v2 E1 design. |
| F2 | The synthetic truth satisfies symmetry: `ga = xpnd(vech(ga))` before homogeneity and adding-up are applied. The DGP's `w = … + prices*ga` makes `ga[k, c]` the coefficient on price k in share equation c, the same orientation as G below. | `tests/quaidsfixtures.src:29-32, 58` | S is correctly specified on the E2 DGP. `pfErr` against that truth stays valid for every arm. |
| F3 | With homogeneity and adding-up, symmetry of the full n x n Γ is equivalent to symmetry of the reduced (n−1) x (n−1) block G. Proof, with γ_ij the coefficient on price j in equation i: γ_in = −Σ_{j<n} γ_ij (homogeneity) = −Σ_{j<n} γ_ji (reduced symmetry) = γ_ni (adding-up). | algebra | Only n1(n1−1)/2 restrictions, on G alone (n1 = n−1). E2's DGP has n = 5 and one intercept shifter, so nIntFull = 2, ng = 8, ng·n1 = 32: 6 restrictions, 32 → 26 free parameters. |
| F4 | `quaidsFit()`'s symmetry block builds a 0/1 matrix `X` (a local, never returned) with vec(b[.,1:n1]) = X·[non-gamma; vech(G)], via the `design(vec(xpnd(seqa(1,1,m))))` idiom. Its weight is the fixed-point delta-method covariance, which is not valid at direct endpoints. | `src/quaids.src:561-573` (X), `:585-616` (weight and its fallback); design §7, §8 | S reuses the same mapping (§2). P keeps the minimum-distance form but uses a weight evaluated at the direct endpoint (§4.3). |
| F5 | The translog index in `_quaidsDirectResid` uses only the symmetric part of G (a quadratic form). | `quaidsdirect.src:137-159` | No residual-code change: S reconstructs `bRaw` and calls the unchanged residual. |
| F6 | `_quaidsDirectJacobian`, `_quaidsDirectGNGradient`, `_quaidsDirectGNHessian` and `_quaidsDirectLM` do not depend on the parameter dimension. But `Core` calls `_quaidsDirectUnpack` by name, and `Jacobian` and `LM` call `Core` by name. | `quaidsdirect.src` | S needs copies of `Core`, `Jacobian` and `LM`, each differing only in which proc it calls (§3). `quaidsdirect.src` itself stays byte-identical. |
| F7 | Optimizer convergence is not a reliability signal, and no truth-free trust flag was validated. | frozen protocol §6, §11.4 | The real-data diagnostics in §6 are **falsification checks**. They can expose a non-global S endpoint. They never certify a correct one. |
| F8 | Seeds already solved by the prototype's own tests: 501, 502, 503 (noisy recovery, selection) and 601 (the T3 stationary fixture, noise seed 6011), plus E2 seeds 31, 47, 90, 104, 115, 123 and 200 in the evidence test. No solve of any kind has run on a seed ≥ 1001; no `symmetry/` output directory exists. | `quaidsdirect_test.e:101,133,160,214`; `quaidsdirect_evidence_test.e:100,169` | §5 seed blocks start at 1001. Fixture seeds 601 and 602 are implementation-check data only. |

## 2. Symmetric parameterization and reconstruction (exact)

**Notation** (follows `quaidsdirect.src`; every vector is a column):
- n goods, n1 = n − 1, nIntFull = intercept rows (constant plus shifters);
- `bRaw` is the **ng x n1** reduced layout with rows `[intercepts 1..nIntFull | gamma 1..n1 | beta | u]`, so ng = nIntFull + n1 + 2. It is `_quaidsDirectOut.bRaw`. Where a source matrix is ng x n (`bStone` = `quaidsOut.homogB`, every start, and `_quaidsDirectUnpack`'s padded output), only its first n1 columns are used, and its n-th column is never read;
- G = `bRaw[nIntFull+1:nIntFull+n1, 1:n1]`, where G[k, c] is the coefficient on relative price k in equation c, which is γ_ck;
- **v** = vec(`bRaw`) (column-major), the **full vector**, length ng·n1; full index of entry (r, c) = (c−1)·ng + r;
- **θ** is the **free (symmetric) vector**, length kS = ng·n1 − n1(n1−1)/2;
- `scaleVec` (length ng·n1) = `_quaidsDirectScale(bStone, n1)`, H's own scale; `scaleSym` (length kS) is defined in §3.3. A suffix "s" means a scaled coordinate: v_s = v ./ scaleVec, θ_s = θ ./ scaleSym.

**Free vector.**
θ = [θ_ng ; θ_G]:
- θ_ng: the m0 = (ng−n1)·n1 non-gamma entries of v, in ascending full-index order;
- θ_G = vech(G) in GAUSS's convention (lower triangle, row-wise; `vech`/`xpnd` are inverses), length n1(n1+1)/2.

**Index map `idx`** (ng·n1 x 1, integers in 1..kS):
- for the j-th non-gamma full position (ascending), idx = j;
- for gamma position (k, c), idx = m0 + V[k, c], where V = `xpnd(seqa(1,1,n1*(n1+1)/2))`.

T = `design(idx)`, **ng·n1 x kS**, so v = T·θ for every symmetric point. It must equal `quaidsFit()`'s `X` for the same layout; this is check I1.

**Reconstruction.**
v = θ[idx], an index gather. It is not computed as T·θ, which would spread a NaN to every row through 0·NaN. G is then symmetric **bitwise**. `bRaw = reshape(v, n1, ng)′`.

**Projection onto θ.**
proj(v) = (T′v) ./ sumc(T), length kS, with column counts 1 or 2. It is exactly (G+G′)/2 on off-diagonal pairs and the identity elsewhere. It is used for starts (§3.3) and in I2. It is an orthogonal (unweighted) projection, distinct from P's weighted one (§4.3).

## 3. Residual, objective, Jacobian, scaling and starts

All of these go in a new private sibling file
`dev/prototypes/aids_direct_lm/quaidsdirectsym.src` (names `_quaidsDirectSym*`). It is `#include`d after `quaidsdirect.src` and calls that file's unchanged procs.

### 3.1 Residual and objective

- **Residual:** reconstruct `bRaw = reshape(v, n1, ng)′ ~ zeros(ng,1)` (the same padding as `_quaidsDirectUnpack`), then call the unchanged `_quaidsDirectResid`.
- **Objective:** unchanged. critExact = −ln det(R′R/N), from the first n1 equations' residuals, larger is better; the solver minimizes ln det.
- `_quaidsDirectSymCore(thetaScaled, w, intcptFull, prices, totexp, u, aCtl, ng, n1, scaleSym, idx)` is `_quaidsDirectCore` with one difference: `bVec = (thetaScaled.*scaleSym)[idx]`.

### 3.2 Jacobian, gradient, Hessian and LM

- **Jacobian:** `_quaidsDirectSymJacobian` is a central finite-difference Jacobian over the kS coordinates of θ_s (step h = 1e-6), size (N·n1) x kS. It is `_quaidsDirectJacobian` with `Core` replaced by `SymCore`.
- **Gradient and Hessian:** `_quaidsDirectGNGradient` and `_quaidsDirectGNHessian` are reused unchanged.
- **LM:** `_quaidsDirectSymLM` is `_quaidsDirectLM` with `Core` and `Jacobian` replaced by their Sym versions, and nothing else. The E2 constants apply: maxIter 500, gradTol 1e-3 (inf-norm, now over the kS scaled coordinates), λ0 1e-3, ρ > 0.1, ×0.5 when ρ > 0.75, ×4 on rejection, stall at 1e12. The retcodes are unchanged.
- **Exactness.** Every residual entry is multilinear in the coefficients: no coefficient multiplies itself (α and γ enter through the price index, which is multiplied by β of the same equation, a different coordinate). A symmetric pair moves both G entries together and still enters linearly. Central finite differences are therefore free of truncation error. Writing J_v for the Jacobian with respect to the unscaled full v (that is, `_quaidsDirectJacobian` called with `scaleVec = ones(ng*n1,1)`), J_sym = J_v·T·diag(scaleSym) up to rounding, and H's scaled Jacobian is J_v·diag(scaleVec). This is check I4.

### 3.3 Scaling and starts

- **Scale:** scaleSym_j = max{ scaleVec[r] : idx[r] = j }. The floor of 0.05 is inherited through scaleVec.
- **Stone inputs:** S computes `bStone`, `u`, `intcptFull` and `scaleVec` exactly as `_quaidsDirectFitStarts` does (a `quaidsFit()` call with `maxiter = 1`), after the same `_quaidsDirectCheckInputs` guards. I5 checks that they equal H's `dOut.bStone`, `dOut.u`, `dOut.intcptFull` and `dOut.scaleVec` bitwise.
- **Start transformation:**
  1. Build the 16 full starts exactly as `_quaidsDirectAllStarts(bStone, nIntFull, n1)` does: the same recipes, internal seed 20261009 and draw order. Start s is the ng x n block at rows (s−1)·ng+1 .. s·ng.
  2. Take v = vec(start[., 1:n1]) and map it to θ = proj(v) (§2).
  3. θ_s = θ ./ scaleSym.

  S and H therefore start from **paired** starts drawn from the same draws, and recipe 1 is the symmetrized Stone start.

### 3.4 Output of S

`_quaidsDirectSymFit(w, intcpt, prices, totexp, instr, aCtl)` returns the existing private `struct _quaidsDirectOut`, whose definition is not changed, filled as `_quaidsDirectFitStarts` fills it, except:
- `model` = "AIDS direct LM, in-objective symmetry (private prototype)";
- `bRaw` and every block of `startRaw` are the **reconstructed** ng x n1 layouts;
- `scaleVec` holds scaleSym (kS x 1), documented in the `quaidsdirectsym.src` header;
- `bFull` = `_quaidsDirectFull(bRaw, nIntFull)`; `homogCrit` holds S's selected critExact (the field name is inherited).

No new struct and no new field is added anywhere.

### 3.5 Exclusions

No new tuning constant is introduced. No P endpoint, truth or H endpoint is ever used as an S start (U5).

## 4. Arms and selection semantics

### 4.1 H: homogeneity-only (reference)

The unchanged `_quaidsDirectFit` (L16), called as is. Its output is the current endpoint.

### 4.2 S: in-objective symmetry (candidate MVP path)

S runs the 16 paired starts through `_quaidsDirectSymLM`. Selection is design §5, verbatim:
- **Eligible:** retcode 0; the objective re-evaluated by `SymCore` at the returned θ_s is feasible and finite; every reconstructed coefficient is finite.
- **Selected:** the strictly largest critExact, with no tolerance; exact ties go to the lower start.
- **No eligible start:** "no result".
- **No fallback** to H, P, `quaidsFit()` or anything else.

The diagnostics (basin rule: max |Δ| over reconstructed `bRaw` rows 1..ng−1, that is every row except u, > 0.05; near-tie flag at margin < 0.01; exact ties) are computed with `_quaidsDirectSelect` on the reconstructed raws. They are diagnostics only.

`S1` is S restricted to start 1: it has a result when start 1 is eligible, and then its endpoint is start 1's. It is a descriptive subset of S's cells, not a separate run.

### 4.3 P: minimum-distance projection of H (comparator only)

P is evaluated only when H returned a result; otherwise P is "no result". Dimensions: v_H, v_P and scaleVec are ng·n1 x 1; J is (N·n1) x (ng·n1); B_s and B are (ng·n1) x (ng·n1); T is (ng·n1) x kS; φ_P is kS x 1.

1. v_H = vec(H.bRaw), unscaled (H.bRaw is already ng x n1).
2. Evaluate H's residual R at v_H (`_quaidsDirectCore(v_H./scaleVec, …, scaleVec)`, as H's own re-evaluation does). Sinv = `invpd(R[.,1:n1]′R[.,1:n1]/N)`, trapped; a failure gives "no result".
3. J = `_quaidsDirectJacobian(v_H./scaleVec, 1e-6, …, scaleVec)`, the Jacobian with respect to H's scaled coordinates. B_s = `_quaidsDirectGNHessian(J, Sinv, N, n1)`.
4. Unscale: since v = scaleVec .* v_s, ∂R/∂v = (∂R/∂v_s)·diag(1/scaleVec), so B = diag(1/scaleVec)·B_s·diag(1/scaleVec).
5. φ_P = `solpd(T′BT, T′B·v_H)`, trapped; a trap, a missing or wrongly sized result, or a non-finite entry gives "no result". Then v_P = φ_P[idx] and P.bRaw = reshape(v_P, n1, ng)′. P therefore minimizes (v_H − Tφ)′B(v_H − Tφ).
6. critP = −fVal from `_quaidsDirectCore(v_P, …, ones(ng*n1,1))` (unit scale, so no rounding round trip). P has a result only if critP is feasible and finite and v_P is finite.

These steps are two private procs in `quaidsdirectsym.src`: `_quaidsDirectSymPWeight` (steps 2-4, returns B or a missing scalar) and `_quaidsDirectSymPSolve(v_H, B, idx)` (step 5, returns an ok flag and v_P), so I7 can call the solve with a constructed B.

**Why it is defensible.** B is proportional to the Gauss-Newton approximation of the Hessian of the concentrated Gaussian criterion at the direct endpoint itself; at a consistent estimate it is consistent for that Hessian, because the omitted terms average residuals that are orthogonal to the regressors. Under true symmetry, a minimum-distance estimator with any consistent positive-definite weight is consistent. With a Hessian weight it is first-order equivalent to minimizing the same criterion under the restriction (that is, to S at its global optimum). The equivalence concerns point estimates and holds conditional on the first-stage `u`; it makes no variance claim.

**What it is not.** It is not shipped `quaidsFit()`'s projection, whose weight is the fixed-point covariance. P makes no covariance claim. Its scale constant is irrelevant to a point estimate.

### 4.4 Out of scope

Shipped `quaidsFit()` is **not** an arm. E2 made no shipped comparison, and adding one is owner decision U4 (deferred).

## 5. Datasets, seeds and separation

Every synthetic and real fit uses `aCtl = quaidsControlCreate()` with `linear = 1`, `homogenous = 1` and every other field at its default (as `quaidsdirect_test.e`).

| Block | Data | Seeds | Use |
|---|---|---|---|
| Fixtures | 601: the T3 stationary construction (I3). 602: `_quaidsSyntheticDGP(3000, 602, 0, 1)` | 601, 602 | Implementation checks only (§7.1) |
| Development (D) | `_quaidsSyntheticDGP(3000, s, 0, 1)`, the E2 harness DGP (`pilot_harness_core.src:122`) | **1001-1030** | Debugging and the smoke gate; gates computed as **informational** |
| Confirmation (C) | same | **2001-2170** (n = 170) | The only binding gates (§7.2) |
| Determinism | same | 2001, 2085, 2170 | C0 re-solve |
| Real data | Blanciforti86 (`tests/fixtures/published/blanciforti86_food32.csv`): n = 4, T = 32; w = the four `wFood` shares; prices = ln `pFood1..4`, totexp = ln `xFood`, instr = ln `xAgg`, intcpt = 0 (as `tests/quaids_published_validation_test.e:83-86`), so nIntFull = 1, ng = 6, n1 = 3 | none | Truth-free diagnostics (§6.2) and gate G6 |
| Reserved, not run in v1 | `_quaidsSyntheticDGP(300, s, 0, 1)` | 3001-3050 | The deferred small-sample panel (§6.3, U7) |
| Reserved revision block | `_quaidsSyntheticDGP(3000, s, 0, 1)` | 2201-2370 | Only after a revision (§8 step 3) |

**Separation rules.**
- No block overlaps E2 (1-200), the prototype tests (501-503, 601; F8) or any other block.
- No solve of any kind has run on any seed ≥ 1001 (F8).
- The harness and analysis script hard-code the block ranges and refuse any seed outside the requested block **before** reading its values. This uses the same guard pattern as `confirm_analyze.py --mode dev`.

## 6. Diagnostics (descriptive unless a gate in §7 names them)

### 6.1 Per seed (synthetic)

A quantity whose inputs include a "no result" arm is missing (`.`), never zero.

- `pfErr` and `rawStructErr` for H, S, S1 and P, using the frozen definitions (frozen §5; `confirm_analyze.py` `pf_err`, `raw_struct_err`, with model = 0). `pf_err` may be imported read-only from `confirm_analyze.py`; I6 checks the definition used.
- **class@1** = the frozen class at cutoff 1 (`correct`: result and `pfErr ≤ 1`; `wrong`: result and `pfErr > 1`; `none`: no result).
- **asym** = max |G − G′| (S and P must give exactly 0; H is reported as is).
- **addUpErr** = max over rows r of |Σ_j bFull[r, j] − e_r|, with e_1 = 1 and e_r = 0 otherwise; **homogErr** = max over equations j of |Σ_{k=1..n} bFull[nIntFull+k, j]|.
- **N1 nesting gap** = critS − critH. Because S is H's problem with restrictions added, N1 > 0.01 flags an H endpoint that is not global.
- **N2 P-dominance gap** = critP − critS. P is a feasible symmetric point, so N2 > 0.01 flags an S endpoint that is not global. The tolerance is the frozen near-tie ε; E0's same-basin spread was ≤ 5.9e-3.
- **LR** = N·(critH − critS). A χ²(n1(n1−1)/2) = χ²(6) reference is printed and **labelled non-inferential**: `u` is a generated regressor, and no size was validated.
- **condGN** = the 2-norm condition number (largest over smallest singular value) of `_quaidsDirectGNHessian` of S's **scaled** Jacobian at the selected endpoint.
- Support (the number of eligible starts in the selected basin), margin, near-tie flag and exact ties.
- Iterations and wall time.

### 6.2 Real data (Blanciforti86; truth-free, one run, everything reported)

H, S and P are each fitted once, with the internal starts. The following are reported:
- R1: N1;
- R2: N2;
- R3: LR, with a χ²(3) reference marked non-inferential;
- R4: support, margin and near-tie;
- R5: asym, addUpErr and homogErr;
- R6 polish: S's LM restarted from its selected θ_s (same scaleSym) must stop with retcode 0 at iteration 1 and the same critExact (≤ 1e-12);
- R7: condGN;
- R8: S − H and P − S coefficient differences and the symmetric own-price gammas, as tables only.

None of these is a correctness claim (F7). T = 32 is far outside the synthetic envelope.

### 6.3 Small-sample panel (deferred; not run in v1)

U7 is resolved as "not run". Seeds 3001-3050 stay reserved and unsolved. If a later version adds the panel, it must be prespecified there; this version gives no metric for it, and nothing from it could enter a decision.

## 7. Checks and gates

### 7.1 Implementation checks (on fixtures; all must pass before any D seed is solved)

| ID | Check | Tolerance |
|---|---|---|
| I1 | `design(idx)` equals `X` rebuilt by a verbatim copy of `src/quaids.src:561-573` (with nint = nIntFull − 1, nendog + nu = 2, ng = nIntFull + n1 + 2), for (nIntFull, n1) ∈ {(1,2), (1,3), (2,4)}; kS = ng·n1 − n1(n1−1)/2 | exact |
| I2 | Round-trip: for symmetric G, `(proj(v))[idx] == v`; for random θ, `proj(θ[idx]) == θ` and the reconstructed G is bitwise symmetric. | exact |
| I3 | **Exact recovery** on fixture 601, the T3 construction from `quaidsdirect_test.e` §1 (noise seed 6011) reused verbatim. Its truth is symmetric (F2), and its noise is orthogonal to every equation block of J_v at the truth, so S's gradient, diag(scaleSym)·T′ times the full unscaled gradient, also vanishes there. S's LM from proj(v*) ./ scaleSym stops at iteration 1 with retcode 0 and gnorm < 1e-8, and returns that θ_s exactly; `_quaidsDirectSymFit` recovers the truth (rows 1..ng−1). | exact; ≤ 1e-6 |
| I4 | On seed 602 at the Stone point and at one random θ: `SymCore(θ_s)` residual equals `_quaidsDirectCore((θ_s.*scaleSym)[idx], …, ones(ng*n1,1))` residual; `J_sym` vs J_v·T·diag(scaleSym) (§3.2). Relative error = max\|A − B\| / max\|B\|. | exact; ≤ 1e-7 |
| I5 | H through the harness equals a direct `_quaidsDirectFit` call; S's Stone inputs equal H's (§3.3); `quaidsdirect.src` is byte-identical to `d8b4318`. | exact; `git diff --exit-code d8b4318 -- quaidsdirect.src` |
| I6 | The `pfErr` used here equals `confirm_analyze.pf_err` on stored E0 rows (seeds 1-30, read only). | exact |
| I7 | P on seed 602: `_quaidsDirectSymPSolve` given v already symmetric returns v; given B = zeros it returns "no result"; P's v_P has asym == 0. | ≤ 1e-10; pass; exact |
| I8 | Selection unit cases (ineligible never selected, ties, none) on S; S inherits the five input guards of `_quaidsDirectCheckInputs`. | pass |
| I9 | `tests/run_source_tests.ps1`, `quaidsdirect_test.e` and `quaidsdirect_evidence_test.e` pass, and `quaidsFit()` output is unchanged. | runner |

### 7.2 Confirmatory gates (seeds 2001-2170 only; the verdict is the intersection of all, so no multiplicity adjustment)

"Correct" means `pfErr ≤ 1`, and "no result" counts as not correct.

| Gate | Statistic | Pass |
|---|---|---|
| **C0 Integrity** | Freeze manifest verifies; complete grid (§9); determinism re-solve gives 0 differences apart from `wallSec`; every S and P result has asym == 0; addUpErr and homogErr ≤ 1e-12; critS recomputed independently in absolute-price form from `bFull` (a copy of `critFromFull` in `quaidsdirect_test.e`, design T5; the S row's `critAbsErr`) ≤ 1e-9 | all |
| **C1 S accuracy** | k_S correct out of 170; exact Clopper-Pearson two-sided 95% lower limit ≥ 0.90 | **k_S ≥ 161** |
| **C2 S returned-wrong** | W_S = number of seeds with a returned and wrong endpoint | **W_S ≤ 1** |
| **C3 S optimization adequacy** (truth-free) | number of seeds with N2 > 0.01 (seeds where N2 is missing are not counted; their number is reported) | **≤ 1** |
| **C4 S vs P non-inferiority** | c = #(P correct ∧ S not), b = #(S correct ∧ P not) | **c − b ≤ 3** |
| **C5 S vs H regressions** | #(H correct ∧ S not correct), with a written mechanism for each | **≤ 3** |
| **G6 Real-data executability** (truth-free; not accuracy) | Blanciforti86: S returns `ok = 1`; P has a result; R1 ≤ 0.01; R2 ≤ 0.01; R5: asym == 0, addUpErr and homogErr ≤ 1e-12; R6 holds | all |

**Why these thresholds** (approved, U6).
- C1 and C2 reuse the frozen AIDS targets (frozen §7 A2, A4: 161 gives a lower limit of 0.9019; 160 gives 0.8945), so S is held to E2's bar.
- C3 caps undetected non-global S endpoints at the same count as A4.
- C4 and C5 margins of 3/170 (about 1.8%) are prespecified counts, not hypothesis tests.
- The seeds are fresh, so unlike frozen A3 and A5, H and P outcomes are prospective here.

### 7.3 Descriptive (never pass/fail)

- The P analogues of C1 and C2: k_P and W_P.
- H's k_H and W_H.
- S1's k and W.
- Paired `pfErr` S − H over seeds where both have a result: median difference and `scipy.stats.wilcoxon` (two-sided, `zero_method="wilcox"`, `method="auto"`) p, as context only.
- Sensitivity cutoffs 0.25, 0.5, 2 and 5.
- The borderline listing for (0.25, 2].
- The N1 flags.
- The LR empirical exceedance rate at the χ²(6) 5% point (a preview only; it validates no test).
- Wall time.

## 8. Run order, stop and falsification rules

1. **Approval (owner): done for v1 (2026-10-12).** Implementation of `quaidsdirectsym.src`, the harness, `sym_analyze.py` and `sym_config.json`, plus their unit tests, may start. Any edit to §2-§7 or §9 from here on is a new version, recorded in this file's Status.
2. **I1-I9** run on fixtures. **Stop rule:** any failure blocks stage D until it is fixed. A fix may not change a constant in §2-§4.
3. **Stage D (seeds 1001-1030; smoke gate).** Run H, S and P, then the analysis in dev mode, where verdicts are informational.
   - **Stop (revise)** if S has any returned-wrong seed, ≥ 2 no-result seeds, any N2 > 0.01, any asym ≠ 0, or a determinism difference (determinism re-solve of seeds 1001 and 1030).
   - A revision means a new protocol version **and a new confirmation block** (2201-2370). Block 2001-2170 is never solved after a revision.
   - Implementation bugs found in D may be fixed. Constants, arms, gates and seed blocks may not change.
4. **Freeze.** Write `SYM_FREEZE_MANIFEST.json` covering:
   - the SHA-256 of this protocol, `quaidsdirect.src`, `quaidsdirectsym.src`, the harness, the analysis script, the config and the tests;
   - the shipped `src/` files the harness includes;
   - `tests/quaidsfixtures.src`;
   - the GAUSS, Python and NumPy/SciPy versions and the git HEAD.

   Its digest goes in a sidecar `SYM_FREEZE_MANIFEST.sha256`, as in frozen §14.5.
5. **Stage C (seeds 2001-2170), run once,** in ascending seed order, each seed running H, S, then P.
   - **Interim:** one look after seeds 2001-2020, for integrity and runtime only. The analysis script's interim mode prints no accuracy statistic.
   - **Integrity failures:** fix and rerun **all** C seeds, keeping the failed run's outputs and a written reason. After two integrity failures the study stops and needs a revision.
   - **No efficacy peeking,** and no rerun for any other reason.
6. **G6** runs after stage C. Its result cannot change C1-C5.

**Seed-tuning protections.**
- Every seed block and the starts' internal RNG seed (20261009) are fixed here.
- H and S use paired starts.
- No constant may change after stage D begins.
- Confirmation seeds are guarded in code.
- The block is single-use.
- Every run's outputs are kept.

**Runtime cap.**
- The total for stage C, including determinism, is **6 h wall** on this machine.
- If the interim projection (median per-seed time × 170) exceeds it, stop for feasibility, record the reason and revise. Accuracy stays unseen.
- Cells over 900 s are listed.
- Timing is never a decision input and never a speed claim (design §11).

## 9. Artifacts and determinism

All future artifacts go in `dev/prototypes/aids_direct_lm/symmetry/` (`dev/` is never packaged), one subdirectory per run: `dev/` (stage D), `dev_det/`, `conf/` (stage C), `conf_det/`, `real/`, and `conf_failed_<k>/` for any integrity-failed run.

**Encoding (every CSV).** UTF-8, comma-separated, one header row, LF line ends. Missing (GAUSS missing, "no result", not applicable) is written as `.`, as in the E2 harness output that `confirm_analyze.py` reads. Integers are written as integers; reals with 17 significant digits. Rows are sorted by the key columns in the order listed.

| File | Key | Columns | Complete-grid count (stage C, ng = 8, n1 = 4) |
|---|---|---|---|
| `cells.csv` | block, seed, arm ∈ {H,S}, start 1..16 | retcode, eligible, critExact, gnorm, iterations, rejected, covInvFail, stepInvFail, infeasTrial, rhoReject, lmAccepted, lambdaFinal (`startTable` columns 2-13), wallSec (`.`: per-start time is not observable for the unchanged H) | 170 × 2 × 16 = 5,440 rows |
| `raw.csv` | block, seed, arm ∈ {H,S,P}, start, row 1..ng, col 1..n1 | value: every H and S start endpoint in the reduced ng x n1 layout (S reconstructed); P uses start = 0 and is written only when P has a result | 170 × 2 × 16 × 32 = 174,080 rows, plus 32 per P result |
| `seeds.csv` | block, seed, arm ∈ {H,S,S1,P} | ok, selStart (P: 0), nEligible, crit, pfErr, rawStructErr, classPf (class@1), asym, addUpErr, homogErr, exactTies, margin, nearTie, iterations, wallSec; N1, N2, LR and condGN on the S row only (`.` elsewhere); critAbsErr. The harness leaves pfErr, rawStructErr and classPf as `.`; the analyzer's copy fills them. S1 and P have `.` for nEligible, exactTies, margin and nearTie. | 170 × 4 = 680 rows, present even when ok = 0 |
| `truth.csv` | block, seed, row, col | the regenerated `trueParams` (no solver) | 170 × 9 × 5 = 7,650 rows |
| `results.json` | — | `protocolVersion`, `block`, `mode` (dev / interim / confirm), `nSeeds`, `gates` (for each of C0-C5 and G6: `statistic`, `threshold`, `pass`), `verdict` ∈ {PROCEED, REVISE, REJECT, INVALID}, `descriptive` (§7.3), `manifestSha256` | one file |
| `report.md` | — | gate table, descriptive tables, the borderline listing, C5 mechanisms | one file |
| `realdata.json` | — | R1-R8, plus H, S and P `ok` flags | one file (in `real/`) |
| `console_*.txt` | — | raw logs | — |

The analysis refuses to compute any gate (verdict INVALID) unless every count above is met for the block and each `seeds.csv` S and H row agrees with its `cells.csv` selection.

**Determinism.**
- JSON keys are sorted.
- Output has no timestamps except `wallSec`, which every comparison excludes.
- Two analysis runs must be byte-identical.
- RNG use is confined to the internal start stream; the caller's `rndseed` stream is untouched (as in design T7).

## 10. Decision rule

Rows are evaluated top to bottom; the first matching row decides.

| Outcome | Decision |
|---|---|
| C0 fails | **INVALID** run; no conclusion (§8 step 5). |
| C1-C5 and G6 all pass | **PROCEED.** S is the preferred symmetric path. The evidence is sufficient for an **experimental, AIDS-only, point-estimate** MVP, subject to the separate owner decisions in §12. |
| C1 or C2 fails, and P passes both (descriptive k_P ≥ 161 and W_P ≤ 1) | **REVISE.** S is not preferred. Any P-based path needs its own new protocol and confirmation block. |
| C1 or C2 fails, and P does not pass both | **REJECT** in-objective and projected symmetry for this design. The prototype stays homogeneity-only and private. |
| C3, C4 or C5 fails (C1 and C2 pass) | **REVISE.** Search breadth (C3), comparator superiority (C4) or regressions (C5) are written up first. |
| G6 fails (C1-C5 pass) | **REVISE.** A real-data executability defect. No MVP until it is fixed and G6 is rerun under a new protocol version. |

## 11. What PROCEED does and does not license

**Enough for an experimental point-estimate MVP:**
- PROCEED;
- E2 (homogeneity-constrained);
- I1-I9 kept as permanent tests next to the prototype.

The MVP must state:
- one synthetic DGP family (n = 5, T = 3000);
- no per-fit reliability signal;
- truth-free checks are falsification-only;
- AIDS only, unweighted, `homogenous = 1`, no user `b0`;
- point estimates only.

PROCEED does not itself expose anything publicly, change shipped `quaidsFit()`, or add returned fields (U2-U4).

**Still required before any covariance or inference** (outside this protocol; each needs its own prespecified study):
1. A covariance derivation for the selected S estimator that accounts for the generated regressor `u`: a two-step sandwich (Murphy-Topel or stacked GMM), or a bootstrap of the whole procedure (first stage, 16 starts and selection).
2. A coverage simulation on fresh seeds: coefficients and elasticities, several T values including small ones, and a prespecified acceptance band around 95%.
3. Size and power of a symmetry LR or Wald test, the latter under an asymmetric-truth DGP.
4. The effect of multi-start selection on inference.
5. Heteroskedasticity, clustering and weights.
6. Downstream integration (elasticities, welfare).

The LR values reported here validate none of these.

## 12. Owner decisions

| ID | Decision | State (2026-10-12) |
|---|---|---|
| U1 | Approve this document as v1 and authorize implementing §3, the harness and the analysis. | **Resolved: approved** for the private implementation stage. |
| U2 | After PROCEED: public exposure (name, printer, `docs/public-api.json`, command reference, version bump, changelog, installed-package test). | **Deferred.** Not authorized; PROCEED alone does not decide it. |
| U3 | Whether an MVP also returns the H endpoint and LR as descriptive fields. | **Deferred.** No new returned fields (§3.4). |
| U4 | Whether shipped `quaidsFit()` (`bS`) is added as a descriptive arm, or `quaidsFit()` changed. | **Deferred.** Not an arm; `quaidsFit()` unchanged. |
| U5 | Whether a later revision adds θ_P as a 17th S start. Doing so removes C3's independence and needs a new block. | **Deferred.** Not in v1 (§3.5). |
| U6 | The thresholds: C1/C2 at the E2 bar; C3 ≤ 1; C4 and C5 margins of 3. | **Resolved: approved** as stated (§7.2). |
| U7 | Whether the small-sample panel (§6.3) runs. | **Resolved: not run** in v1; seeds 3001-3050 reserved. |
| U8 | Whether a real-data holdout fit (for example, fitting on 1947-1970 and scoring 1971-1978) is wanted; it would need the first-stage `u` split too. | **Resolved: not run** in v1. |

None of U2-U5 blocks implementation or stages D and C.
