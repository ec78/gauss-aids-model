# Convergence-failure research: Anderson acceleration + multi-start (uncommitted experiment)

Status: **exploratory, not shipped**. Nothing in this directory is wired
into `package.json`'s `src` array, `tests/run_source_tests.ps1`, or any
other CI/release path. This is scratch research output, kept for a
handoff, not production code. See the repo-root `CLAUDE.md`/
`PROJECT_STATUS.md` for the actual shipped state of the library.

## Second-opinion methodology audit (2026-09-29)

This section supersedes two interpretations later in this document.  The
original raw files are retained as an audit trail, but the prior claim that
`homogCrit` cannot select correct fixed points and the hypothesis that
Anderson suppresses basin diversity are **not supported** after correcting
the pilot.

New artifacts:

- `anderson_methodology_audit.e`: 30 seeds x 2 models x 8 identical
  perturbed starts, comparing Anderson with plain iteration, both
  directions of `homogCrit`, an independently recomputed residual
  criterion, and an oracle selector.
- `anderson_methodology_audit_output.txt`: complete numeric output.
- `quaidsfit_anderson_prototype.src`: gained an opt-in
  `residualStop=1` audit mode.  The default remains zero, so the earlier
  runs and the `mDepth=0` faithfulness path retain their prior behavior.

### Audit conclusion 1: the core Anderson algebra is recognizable, but the implementation is only a prototype

The `vec()`/`reshape()` round trip is correct: `vec(b0)` stacks columns,
and `reshape(bVec, cols(b0), rows(b0))'` reconstructs the original
`ng x n` matrix under GAUSS's row-major `reshape()` behavior.

The update is an algebraically valid Type-II Anderson form *without the
ridge*: differences from every retained point to the current point span
the same subspace as the usual consecutive-difference columns, with an
invertible triangular change of basis.  The fixed `1e-10*I` ridge breaks
that basis invariance, however, so this is a particular regularized
variant rather than literally the standard consecutive-difference
algorithm.  It also solves the normal equations, which square the
condition number.  Walker and Ni explicitly recommend QR and dropping
history columns to control conditioning rather than relying on normal
equations ([paper, sections 3-4](https://users.wpi.edu/~walker/Papers/Walker-Ni,SINUM,V49,1715-1735.pdf)).
Before production use, replace the fixed ridge/normal-equation solve with
QR or rank-revealing QR/SVD, use consecutive differences, scale any
regularization, and add a residual-growth safeguard/restart.

The original accelerated branch tested the relative accelerated *step*
`b_new-b_old`, not the fixed-point residual `T(b_old)-b_old`. Anderson
history terms can make the former small through cancellation without the
latter being small.  The opt-in audit mode fixes that distinction.  On
the 30 unperturbed starts the old and residual-based flags agreed in all
60 seed/model cases, so this is a real correctness defect but **not an
empirically demonstrated explanation** of the reported 200-seed result.

The isolated linear toy and `mDepth=0` faithfulness checks were rerun
after the audit change: the toy still recovers the fixed point to
`9.3e-13`, and `mDepth=0` remains exactly identical to `quaidsFit()` in
both `bS` and `vS`.  Those checks are useful but limited: faithfulness
does not exercise the accelerated branch, and one linear problem does
not validate conditioning or globalization on this nonlinear map.

### Audit conclusion 2: the 200-seed safety-signal finding is well supported

The paired transitions in `anderson_sweep_output.txt` are stronger than
the marginal table alone:

| Baseline never-converged seed | Anderson wrong | Anderson correct | Anderson still nonconverged |
|---|---:|---:|---:|
| Iterated AIDS (78 seeds) | 59 | 17 | 2 |
| QUAIDS (109 seeds) | 62 | 8 | 39 |

Anderson also preserved almost every already-converged classification
(AIDS: 37/38 wrong stayed wrong and 84/84 correct stayed correct;
QUAIDS: 42/43 wrong stayed wrong and 47/48 correct stayed correct).
Therefore the statement "most newly converged cases are wrong" is not a
story inferred from aggregate percentages; it is directly present in
the paired seeds.  Even if some accelerated stops are better described
as small-step endpoints than rigorously verified fixed points, the
user-facing safety conclusion survives: many materially wrong returned
estimates acquire `converged==1`.

The `recErr <= 1.0` correct/wrong cutoff is deliberately loose and should
be sensitivity-tested before publication.  It does not undermine the
large failures (many QUAIDS errors are tens to hundreds), but cases close
to 1 can cross buckets under tiny numerical changes.

### Audit conclusion 3: the multi-start selection result had its objective direction reversed

`quaidsFit()` defines

```text
homogCrit = -ln(det(Sigma))
```

so **larger is better**, not smaller.  The original pilot selected the
lowest value.  That exactly explains its headline disagreement cases:
for example, QUAIDS seed 3 selects `recErr=53.07` when minimizing but
`recErr=0.62` when maximizing.

Corrected 30-seed results (residual-validated Anderson endpoints):

| Selector | Iterated AIDS correct | QUAIDS correct |
|---|---:|---:|
| maximize `homogCrit` | 16/30 | 11/30 |
| minimize `homogCrit` (old pilot) | 14/30 | 8/30 |
| maximize directly recomputed final-structural criterion | 16/30 | 11/30 |
| oracle minimum `recErr` | 16/30 | 11/30 |

Thus, on this pilot, correctly oriented `homogCrit` matches the oracle
*classification on every seed*.  The old finding that the model's own
criterion fails to discriminate correct from wrong endpoints is
withdrawn.  Reachability remains poor -- the oracle itself is only
53.3% for AIDS and 36.7% for QUAIDS -- but selection among the endpoints
that were reached is not the demonstrated problem.

An independently recomputed raw homogeneity criterion matches
`qOut.homogCrit` within `5e-5` across the pilot.  A recomputation from
final `bS` can differ materially from `qOut.symcCrit` on bad endpoints
(maximum differences 1.71 for AIDS and 29.50 for QUAIDS), so
`symcCrit`/final-structural criterion equivalence should be investigated
separately before using `symcCrit` as a selector.

### Audit conclusion 4: no current evidence that Anderson suppresses basin diversity

Using the same 0.3-relative starts and counting converged endpoints as
distinct when their final `bS` differs by more than `1e-4`:

| Model | Anderson: seeds with >1 endpoint | Plain: seeds with >1 endpoint | converged starts (Anderson / plain) |
|---|---:|---:|---:|
| Iterated AIDS | 9/30 | 3/30 | 240 / 133 |
| QUAIDS | 19/30 | 2/30 | 188 / 92 |

This is evidence against the proposed "Anderson collapses basins"
mechanism, not proof of the opposite: plain iteration is heavily censored
by nonconvergence (13 AIDS and 18 QUAIDS seeds had no converged plain
endpoint), so endpoint diversity is not directly comparable without a
longer or otherwise stabilized plain run.  Still, Anderson plainly does
not force all perturbations to one endpoint in this audit.  The 0.3
perturbation is local and multiplicative, but all Stone-start entries in
these 60 seed/model cases were nonzero, so exact-zero freezing was not
the explanation here.

### Audit conclusion 5: curvature attempt #2 has an additional evaluation-point confound

`_quaidsCurvatureSyntheticDGP()` constructs curvature at the **sample
mean**.  `_slutzkyMaxEig()` instead takes the worst eigenvalue over
**every observation**.  Those are different hypotheses.  The attempt's
own reference check reported `44.43` where the fixture comment expected
about `0.17`, and `recErr=0.465` where it expected about `0.16`; that
mismatch should have invalidated the test before interpreting its 40
identical endpoints.  Curvature attempt #2 is therefore confounded as
well as diversity-limited.  A valid retry must evaluate the true and
fitted Slutzky matrices at the same reference point (preferably the
fixture's sample mean) before studying selection.

### Recommended next step after the audit — implemented 2026-10-01

Do **not** spend the next cycle tuning Anderson depth alone.  First run a
full 200-seed corrected multi-start sweep that (a) maximizes, rather than
minimizes, `homogCrit`; (b) requires a true fixed-point residual; and
(c) uses a numerically safer AA solve/safeguard.  In parallel, expand
reachability with structured, scale-standardized starts over the
price-index-driving alpha/gamma/beta/lambda blocks, across several
perturbation scales.  The present bottleneck is reaching a good endpoint,
not recognizing it once reached.

The proposed "truth criterion" check was also run in a corrected form:
on zero of 30 AIDS seeds and zero of 30 QUAIDS seeds did any wrong
endpoint beat the true parameters under either the raw homogeneity or
directly recomputed structural residual criterion.  This is evidence
*against* the current weak-identification story and makes direct
nonlinear optimization more worth trying after the corrected multi-start
benchmark.  It is not, by itself, a proof of identification: on a finite
sample an estimate can legitimately fit better than the generating
parameter, and same-sample fit comparisons cannot distinguish weak
identification from overfitting.  A stronger identification check is a
noiseless/population-moment experiment, an objective-Hessian/rank
diagnostic near competing solutions, or held-out samples from the same
DGP.

If the corrected/structured 200-seed multi-start still leaves most seeds
without a good reachable endpoint, then move to a direct nonlinear
optimizer (LM/Gauss-Newton or `optmt`/`cmlmt`) on a clearly defined
concentrated residual/GMM objective.  Validate it first on a noiseless
known-truth fixture and compare training plus held-out criteria; do not
infer identification merely because a finite-sample optimum beats the
true parameter on its training sample.

## Full safeguarded structured multi-start benchmark (2026-10-01)

The recommended benchmark above has now been implemented and run for all
200 synthetic seeds. The implementation and complete output are:

- `anderson_full_multistart_benchmark.e`
- `anderson_full_multistart_benchmark_output.txt`
- `quaidsfit_anderson_prototype.src` (experiment-only safe-mode extension)
- `anderson_safe_solver_test.e` (isolated rank-deficient SVD and real-fit
  validation)

Each seed uses 16 starts: Stone; two relative all-parameter perturbations;
three scale-standardized all-parameter perturbations; two alpha-only and
two gamma-only perturbations; gamma and beta sign flips; beta-only noise;
two joint core-block perturbations; and a QUAIDS lambda perturbation (or a
joint gamma/beta sign flip for AIDS). Every fit uses Anderson depth 8,
true fixed-point-residual stopping, consecutive-difference history, a
rank-truncated SVD solve (`rankTol=1e-10`), a 10x residual-growth restart,
and a 10x accelerated-step safeguard. Selection maximizes
`homogCrit=-ln(det(Sigma))`.

The safe SVD helper traps an unusable decomposition and returns rank zero;
the caller then takes a plain fixed-point step and clears the history.
This path was added after QUAIDS seed 17 exposed a real SVD failure, then
reproduced successfully through that seed before restarting the complete
sweep. It prevents an ill-conditioned history from aborting the fit; it
does not label the failed accelerated step as convergence.

### Full results

| Result (of 200 seeds) | Iterated AIDS | QUAIDS |
|---|---:|---:|
| Safe Anderson, Stone only: never / wrong / correct | 6 / 96 / 98 | 33 / 111 / 56 |
| Max-`homogCrit` multi-start: never / wrong / correct | 3 / 76 / 121 | 9 / 122 / 69 |
| Oracle among converged starts: wrong / correct | 75 / 122 | 122 / 69 |
| Selector misses when oracle had a correct endpoint | 1 | 0 |
| Seeds with multiple distinct converged endpoints | 71 | 119 |
| Converged starts (of 3,200) | 3,091 | 2,639 |
| Seeds where any wrong endpoint beat truth criterion | 19 | 5 |
| Selected `recErr <= .5 / 1 / 2 / 5` | 38 / 121 / 155 / 171 | 16 / 69 / 88 / 91 |
| Oracle `recErr <= .5 / 1 / 2 / 5` | 39 / 122 / 156 / 171 | 18 / 69 / 88 / 92 |

The independent raw-criterion replica agrees with `qOut.homogCrit` to a
maximum absolute difference of `6.05e-5` for AIDS and `4.37e-4` for
QUAIDS. These differences are small relative to the large criterion gaps
at the bad endpoints and do not change the reported selector result.

### Interpretation and decision

The corrected selector is not the bottleneck. It matches the oracle's
loose correct/wrong classification on 399 of 400 seed/model cases. The
multi-start portfolio materially improves reachability for AIDS (98 to
121 correct) and QUAIDS (56 to 69 correct), and reduces complete
nonconvergence from 6 to 3 and 33 to 9 respectively. It also confirms
substantial endpoint diversity rather than Anderson-induced basin
collapse.

However, this is not ready to ship. Under the deliberately loose
`recErr <= 1` definition, the final correct rate is only 60.5% for AIDS
and 34.5% for QUAIDS; at `recErr <= .5` it is only 19% and 8%. For
QUAIDS, most newly returned results remain converged-but-wrong, so the
original safety concern persists even with safeguards and structured
starts. The Stone-only comparison here is safe Anderson, not the
original plain-iteration baseline, and should not be conflated with the
earlier baseline table.

The earlier 30-seed result that no wrong endpoint beat the true
parameters did not survive at full scale: this occurred for 19 AIDS and
5 QUAIDS seeds. That observation is compatible with finite-sample
overfit and is not proof of weak identification, but it rules out using
the simple same-sample truth comparison as an identification verdict.

The next research step should therefore be a direct nonlinear optimizer
on a precisely defined concentrated residual/GMM objective, first on a
noiseless known-truth fixture and then with held-out-sample comparison and
local Hessian/rank diagnostics. The safeguarded Anderson/multi-start
implementation should remain experiment-only until such a method restores
a trustworthy real-data failure diagnostic, not merely a higher
`converged==1` rate. See `CLAUDE_HANDOFF.md` for the current handoff.

## The problem being investigated

`quaidsFit()`'s iterated AIDS/QUAIDS modes (`aCtl.maxiter > 1`) use
successive substitution (alternating: fix the translog price index, GLS-solve
coefficients, rebuild the price index, repeat) to handle the nonlinearity in
the translog price index. This has a documented, real convergence-failure
problem, quantified by the existing committed diagnostic
`tests/quaids_convergence_sweep.e` (200 seeds, `tobs=3000`, default
settings `aCtl.relax=1`, `aCtl.err=.0001`, `aCtl.maxiter=100`):

| | Iterated AIDS | QUAIDS |
|---|---|---|
| never-converged (`qOut.converged==0`) | 39% | 54.5% |
| converged-but-**wrong** (converged, but far from the known-true synthetic DGP params) | 19% | 21.5% |
| converged-correctly | 42% | 24% |

"Converged-but-wrong" means the fixed-point iteration settled into a
self-consistent but incorrect answer — `qOut.converged==1` gives no hint
anything is wrong. See `README.md`/`docs/FEATURE_SUPPORT_MATRIX.md` for
where this is documented for library users.

## What's been tried, in order, with real empirical results

### 1. Anderson acceleration (Type-II, depth `m`) — works for convergence *speed*, not correctness

`quaidsfit_anderson_prototype.src` is a byte-faithful copy of
`src/quaids.src`'s `quaidsFit()` (verified: `mDepth=0` reproduces the real
`quaidsFit()` to floating-point identity — see
`anderson_faithfulness_check.e`), renamed `_quaidsFitAnderson()`, with one
change: the iteration loop's update step can use Anderson acceleration
(history-based extrapolation across the last `mDepth` iterates, standard
Type-II formulation: `x_new = x + beta*g - (DX + beta*DG)*mixCoef` where
`g` is the residual `T(x)-x`) instead of the existing flat
`aCtl.relax`-damped update. `aCtl.relax` is read as the Anderson mixing
parameter `beta` in the accelerated branch.

**Correctness of the math validated first, in isolation**, on a toy
linear fixed point `x_{k+1}=A*x_k+c` with spectral radius 1.4 (genuinely
divergent under plain Picard iteration) — see `anderson_toy_test.e`.
Plain iteration diverged to ~9e5 within 42 steps; Anderson with
`depth>=n` (the problem's own dimension) recovered the true fixed point
to ~1e-12 in 8 steps. This confirmed the linear algebra before touching
anything AIDS-related.

**Full 200-seed head-to-head** (`anderson_sweep_compare.e`, output in
`anderson_sweep_output.txt`, `mDepth=8`, otherwise identical settings to
the committed sweep):

| | Iterated AIDS baseline | + Anderson(8) | QUAIDS baseline | + Anderson(8) |
|---|---|---|---|---|
| never-converged | 39% | **1%** | 54.5% | **20%** |
| converged-but-wrong | 19% | **48%** | 21.5% | **52.5%** |
| converged-correctly | 42% | **51%** | 24% | **27.5%** |

**Interpretation**: Anderson nearly eliminates never-converged and
improves the net correctly-converged rate in both models. But it does
this by converting most of the never-converged failures into
converged-but-wrong ones, not into converged-correctly ones — it finds
*a* fixed point much faster, not necessarily the right one. **Practical
concern**: on real data there is no `recErr`/ground truth to check —
users only ever see `qOut.converged`. Shipping Anderson as a silent
default replacement would convert a large fraction of today's
self-flagged failures (`converged==0`, "don't trust this") into
confident-looking wrong answers (`converged==1`), which is a regression
in the *safety signal* even though raw accuracy improves. If Anderson
ships, it should not be a silent default — either opt-in, or paired with
a diagnostic that restores an equivalent trust signal (see multi-start /
curvature below, neither of which currently provides one reliably).

### 2. Multi-start, selected by `qOut.homogCrit` — does not clearly beat Anderson alone

> **Superseded by the 2026-09-29 audit above.** This pilot minimized
> `homogCrit`, but `homogCrit=-ln(det(Sigma))` must be maximized. With the
> direction corrected, it matches the oracle correct/wrong classification
> on all 30 pilot seeds in both models.

Idea: since Anderson makes each attempt cheap, refit from several
perturbed starting points (Gaussian noise, scale 0.3 relative, around the
existing deterministic Stone/LA-AIDS starting point — obtained cheaply via
`aCtl.maxiter=1`) and pick the best *converged* one using a criterion
computable on real data (no ground truth): `qOut.homogCrit`, the
already-shipped GLS log-det fit criterion (`-ln(det(S[1:n-1,1:n-1]))`).

`anderson_multistart_pilot.e` (30 seeds, `numStarts=8`, `anDepth=8`;
output `multistart_pilot_output2.txt`) reports both the realistic
(`homogCrit`-selected) outcome and an oracle upper bound (best-of-8 by
`recErr`, uses ground truth, only possible because the DGP is synthetic):

| | Iterated AIDS: homogCrit-picked | Iterated AIDS: oracle | QUAIDS: homogCrit-picked | QUAIDS: oracle |
|---|---|---|---|---|
| converged-correctly | 46.7% | 53.3% | 26.7% | 36.7% |

**Two separate negative findings, not one**:

1. **Perturbation mostly isn't reaching different basins.** On many
   seeds, all 8 perturbed starts converge to the *identical* fixed point
   (confirmed directly, e.g. `anderson_multistart_diag.e` — not copied
   here, but reproducible: seed 1, Iterated AIDS, all 8 `homogCrit`
   values agree to 4 decimal places). This is why the oracle ceiling is
   only modestly above Anderson-alone (53.3%/36.7% vs. ~51%/27.5%) — most
   of the time there's no real diversity to select among.
2. **When diversity does exist, `homogCrit` doesn't reliably pick the
   better one.** In 5 of 60 seed/model runs, a genuinely correct fixed
   point WAS reachable among the 8 starts, but `homogCrit` picked a wrong
   one instead (e.g. Iterated AIDS seed 20: oracle `recErr=0.55`
   (correct) vs. picked `recErr=7.33` (wrong); QUAIDS seed 29: oracle
   `recErr=0.85` (correct) vs. picked `recErr=26.4` (wrong)). See the raw
   per-seed lines in `multistart_pilot_output2.txt` for all cases.

**Net**: as implemented, this specific combination (0.3-scale Gaussian
perturbation + `homogCrit` selection + Anderson) is not clearly better
than Anderson alone, for 8x the compute.

**A real scripting gotcha hit and fixed here, not GAUSS-estimation-
related but worth knowing before extending this work**: GAUSS identifiers
are case-insensitive (already a documented `CLAUDE.md` gotcha). The
original pilot script used `K` for "number of starts" and `k` for the
loop counter — the SAME variable in GAUSS. Setting `k=1` inside the loop
silently reset `K` too, turning `do while k <= K` into an unconditional
infinite loop. This actually ran for ~9.7 CPU-hours before being killed
and diagnosed. Fixed by renaming to `numStarts`. **Grep any new script in
this line of work for single-letter variable names that might collide
case-insensitively with another identifier before running it
unattended.**

### 3. Curvature (Slutzky negative-semidefiniteness) as a selection criterion — inconclusive, not ruled out

> **Additional audit caveat:** attempt #2 is also evaluation-point
> confounded: the fixture imposes curvature at the sample mean while the
> script tests the worst observation in the sample. See the audit above.

Hypothesis: a wrong fixed point might be more likely to violate economic
theory (curvature) than a correct one, giving a theory-grounded,
no-truth-needed selector — reusing already-shipped `quaidsSlutzky()`
machinery (a numeric, value-returning replica, `_slutzkyMaxEig()` in both
`anderson_curvature_check.e` and `_check2.e`, validated to exactly match
`quaidsSlutzky()`'s own printed Maximum-eigenvalue output before being
trusted — see the validation block in `slutzky_eig_check.e`).

**Attempt 1** (`anderson_curvature_check.e`, output
`curvature_check_output.txt`): used the general
`_quaidsSyntheticDGP()` fixture (same one the main sweep uses).
**Confounded and invalidated**: that fixture's true gamma is NOT
constructed to be curvature-consistent at all (unconstrained random
draws) — curvature was violated in 100% of both correct (194/194) and
wrong (233/233) fits, including at the recovered-truth fixed points
themselves. This says nothing about the hypothesis; it just means the
experiment couldn't test it, because the ground truth itself doesn't
respect curvature here.

**Attempt 2** (`anderson_curvature_check2.e`, output
`curvature_check2_output.txt`): switched to
`_quaidsCurvatureSyntheticDGP()` (`tests/quaidsfixtures.src`), whose true
gamma IS constructed to be curvature-consistent at its own sample mean by
a fixed-point iteration. Real constraints: AIDS/linear only (no QUAIDS
version of this fixture exists), and it takes no `seed` argument (one
fixed internal seed=500) — so diversity has to come from many starting
points on one dataset, not many datasets.

**Inconclusive, not negative**: all 40 perturbed Anderson-accelerated
starts (0.3 relative scale, same as before) converged to the identical
fixed point (`recErr` and max Slutzky eigenvalue agree to 5+ decimals
across all 40). Zero diversity means the curvature-vs-correctness
question was never actually put to the test here. This is the SAME
basin-collapse observation as finding #2 above on this one dataset and
set of starts; the broader audit above subsequently found substantial
endpoint diversity on other seeds.

**This hypothesis was subsequently tested in the audit above and was not
supported**: with identical starts and the same iteration cap, Anderson
produced more converged endpoints and more observed endpoint diversity
than plain iteration.  The comparison remains censored by plain
iteration's high nonconvergence rate, so it does not prove that Anderson
expands basin diversity; it does rule out treating suppression as the
current leading explanation.

## What remains NOT tried after the completed full benchmark

- A systematic design/ablation of the 16-start portfolio. Structured
  block perturbations, sign flips, and broader scale-standardized noise
  have now been tried, but no recipe was uniformly dominant and the
  portfolio was not optimized for cost or cross-seed robustness.
- A fully comparable long-run multi-start without Anderson.  The audit
  did run the same starts with `mDepth=0`, but `maxiter=100` left many
  plain starts nonconverged, so its lower observed endpoint diversity is
  censored and not a definitive basin-volume comparison.
- A genuine reformulation of the estimator as direct nonlinear
  optimization (Gauss-Newton/Levenberg-Marquardt on the concentrated
  SSR/GMM objective via `optmt`/`cmlmt`, both already dependencies
  elsewhere in this repo) instead of successive substitution. This was
  flagged early as the "real" structural fix but not attempted, partly
  because the multi-start/curvature findings above suggest the deeper
  problem may be weak identification (multiple genuinely-hard-to-
  distinguish fixed points) rather than a search-strategy problem — and
  if so, a fancier optimizer would hit the same wall.  The cheap
  same-sample truth comparison has now been done (with the corrected
  direction: higher `-ln(det(Sigma))` is better) and found zero wrong
  endpoints beating truth in the 30-seed pilot; see the audit above.
  What remains is a stronger noiseless/population, held-out, and Hessian/
  rank identification analysis.

## Environment / how to run this stuff

- GAUSS 26 at `C:\gauss26`, `tgauss.exe` at `C:\gauss26\tgauss.exe`. Run
  any `.e` file: `tgauss -b -x <file>.e` (working directory matters for
  relative `#include`s — the files here use full absolute paths for their
  `#include`s of the real `src/`/`tests/` files specifically so they can
  be run from anywhere; check each file's own `#include` block).
- These files `#include` real repo source directly (absolute paths into
  `src/`/`tests/`) — they do NOT modify any shipped file. Safe to run
  against the current working tree as-is.
- **Read the repo-root `CLAUDE.md` before writing any new GAUSS code
  here** — it documents ~15 real, previously-confirmed GAUSS-26 language
  gotchas (reserved words including `gamma`/`msym`, case-insensitive
  identifiers, `reshape()` fills row-major, legacy `$+` character-matrix
  8-character cell truncation, `$|` vs `$+` character-matrix type
  incompatibility, etc.) that were each independently rediscovered the
  hard way during this exact line of work and are now durably documented
  there specifically so they don't have to be rediscovered again.
- Full runs are genuinely slow: the 200-seed sweep and the 30-seed
  multi-start pilot each took real wall-clock time (several minutes to
  run in the background). Anderson-accelerated fits are much faster
  per-attempt (~15-30 iterations) than baseline fits that exhaust
  `aCtl.maxiter=100`, but a full sweep still means hundreds of fits.
  **Always run new sweep-scale scripts in the background with a
  generous timeout, and sanity-check correctness on 1 seed / a handful of
  iterations first** — this investigation hit one genuine ~24-CPU-hour
  hang (a badly-scaled real MLE problem, TVP-AIDS-related, not part of
  this convergence work) and one ~9.7-CPU-hour infinite loop (the `K`/`k`
  case-collision bug above) by not doing this carefully enough at first.
