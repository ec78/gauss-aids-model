# Convergence-failure research: Anderson acceleration + multi-start (experiment only)

Status: **exploratory, not shipped**. Nothing in this directory is wired
into `package.json`'s `src` array, `tests/run_source_tests.ps1`, or any
other CI/release path. This is scratch research output, kept for a
handoff, not production code. See the repo-root `CLAUDE.md`/
`PROJECT_STATUS.md` for the actual shipped state of the library.

**Latest result:** see "Matched 30-seed pilot (2026-10-05)" below.

## Phase 1: direct-optimization objective validation (2026-10-01)

> **Corrected by the Phase 1b audit below (2026-10-01).** This section's
> own closing language ("direct optimization does not fix the
> convergence problem", treating the QUAIDS recovery-gate failure as
> evidence bearing on weak identification) overreached what the evidence
> actually supported and is withdrawn. The Phase 1b audit found two
> confounds entangled in this section's own Check 5: the fixture's
> first-stage residual regressor was itself made weakly identified by
> the same noise-shrinking used to keep the GLS criterion well-
> conditioned (so grading recovery against it, undifferentiated from the
> real structural blocks, was comparing against a target the fixture
> itself had made hard to estimate well), and the optimizer (generic
> optmt behind a flat infeasibility penalty, or a plain gradient-descent
> warm start) had no real safeguards and was never shown NOT to work --
> only shown that THESE TWO PARTICULAR CHOICES didn't reliably work.
> With a corrected fixture and a proper Levenberg-Marquardt solver (never
> reports convergence at an infeasible point, rejects bad steps via an
> actual/predicted-reduction test), direct optimization recovered the
> true structural parameters to ~0.001-0.002 accuracy from both 1% and
> 5% perturbations, in every one of 8 trials across both models. This
> section's own numbers and artifacts are left as-is below (an accurate
> record of what that specific recipe did), but its narrative conclusion
> should be read through the Phase 1b section instead.

Continuing from the completed 200-seed safeguarded-Anderson/multi-start
benchmark above (whose own conclusion was: selection is not the
bottleneck, reaching a good endpoint is), this phase implements and
validates the recommended next step -- a reusable, experiment-only direct
optimizer objective for the iterated AIDS/QUAIDS structural
parameterization -- BEFORE any Monte Carlo sweep, per the stage-gate
discipline in `CLAUDE_HANDOFF.md`/`CODEX_HANDOFF.md`.

**New artifacts**: `direct_objective.src` (the objective/pack/scale/
gradient-descent-warm-start module), `phase1_objective_validation.e`
(10-check validation script, both models), its saved output
`phase1_objective_validation_output.txt`, and `phase1_quaids_diag.e` (a
scratch diagnostic that found the central QUAIDS instability mechanism
below -- not part of the validation record itself, kept for evidence).

### The objective and its validation

`_quaidsDirectRawCrit()` (a renamed copy of the existing benchmark's own
`_rawBenchmarkCrit()`, one deliberate deviation -- see below) scores ANY
candidate `b` (ng x n1 reduced coefficient matrix, same row layout
`quaidsFit()`'s own fixed point and `qOut.bS` use) against
`quaidsFit()`'s own criterion `homogCrit=-ln(det(residCov))`, built from
`b`'s OWN alpha/gamma/beta rows (self-referential: the same `b` plays
both the nonlinear-price-index role and the GLS-regression-coefficient
role). A fixed point of successive substitution is a point where those
two roles agree in a GLS sense; it is not generally a stationary point of
this objective in the OTHER (price-index) direction -- that gap is
exactly what a direct optimizer is meant to close. optmt minimizes, so
the function actually handed to it is `-homogCrit = ln(det(residCov))`.

Preconditioning: `_directObjScale()` builds one scale per free parameter
element from the cheap Stone/LA-AIDS start (always available, no ground
truth needed), floored at 0.05; the optimizer always searches in units of
this scale.

A near-noiseless (not literally zero-noise -- see `direct_objective.src`'s
own header for why exact zero noise makes the criterion singular/
divergent) known-truth fixture, `_quaidsNoiselessDGP()`, was built
specifically for this validation (kept out of `tests/quaidsfixtures.src`
per this task's own scope constraint). Its noise floor (`epsScale=1e-2`)
was tuned twice from an initial `1e-4`: that smaller value made
`ln(det(.))` hypersensitive (residual covariance ~1e-32, so a ~1e-6
relative computational difference could move the criterion by O(1)) and
separately made optmt's very first trial step overflow the QUAIDS
exponential deflator -- both confirmed directly, not assumed, and both
are now documented, durable comments in the source.

Confirmed on both Iterated AIDS and QUAIDS (`phase1_objective_validation.e`,
`tobs=500`, `seed=4242`):

| Check | AIDS | QUAIDS |
|---|---|---|
| 1. Objective finite at truth | PASS (36.95) | PASS (36.95, identical -- both models share the same noise draw at this seed by this fixture's own construction, matching `_quaidsSyntheticDGP`'s existing pattern of reusing one `seed` for every draw) |
| 2. Objective matches `qOut.homogCrit` at a converged fitted endpoint | PASS (diff 3e-4 / 2e-4) | PASS |
| 3. Pack/unpack round trip exact | PASS | PASS |
| 4. FD gradient stable across step sizes 1e-4/1e-6/1e-8 | PASS | PASS |
| 5. Direct optimizer recovers truth from a 5%-perturbed start | PASS (recErr 0.078, in the run saved to disk) | **FAIL** (recErr 0.225, best of two methods) |

Check 2 needed a damped (`aCtl.relax=.3`) iterated fit to converge at all
on this near-noiseless fixture under default `relax=1` -- confirmed
directly: undamped successive substitution did NOT converge within 200
iterations here either, the same convergence-failure phenomenon this
whole research line studies, now also reproduced on a near-noiseless
fixture, not only the noisy 200-seed one. This also explains why Check
2's replica agreement is only meaningful near a genuine fixed point:
`qOut.homogCrit` is computed from the iteration's own `gg`/`m`, which were
built from the PRIOR iterate's regressors, not the final `b`'s own -- they
coincide to high precision only once `b≈b0` (converged), a previously
undocumented mechanism found while debugging this check.

### Check 5's negative result is the headline finding, not a bug

This took five rounds of real debugging before concluding it is a
genuine property of the objective, not an implementation defect:

1. A hard, flat `1e8` infeasibility penalty (when `critDet<=0`) gives
   optmt zero numerical gradient on the plateau, so it either falsely
   reports convergence after a handful of iterations or burns its whole
   budget stuck there. Fixed partially with an exponent clip (`[-30,30]`
   on the QUAIDS deflator argument before `exp()`) so a wild step cannot
   literally overflow to `inf`/`NaN`.
2. A fixed box bound in scaled units (`+/-10`) was found to EXCLUDE the
   true parameters themselves for AIDS (`thetaScaled` at truth reached
   `+/-20` there, because the u-residual-correction row's natural scale
   differs sharply from the Stone-start scale `_directObjScale()` uses) --
   widened to `+/-50` after confirming this with a direct diagnostic
   print, not a guess.
3. A targeted diagnostic (`phase1_quaids_diag.e`) at a confirmed-FEASIBLE,
   5%-perturbed QUAIDS start found: optmt's own default quasi-Newton
   search still moved the objective from -4.93 to the +1e8 penalty within
   1-3 iterations; a manual steepest-descent probe at that SAME point
   found a real descent direction exists but only within an extremely
   narrow step window (good from `alpha=1e-8` to `~1e-6`, reversed and
   blew up by `alpha=1e-4`) -- genuinely severe local ill-conditioning
   along the raw gradient direction, consistent with
   `quaidscurvature.src`'s own documented finding that QUAIDS's analogous
   outer loop is "measurably less stable" than AIDS's, both traceable to
   the same translog/quadratic-deflator nonlinearity.
4. `_directObjGradDescent()`, a small hand-rolled, strictly-monotone,
   self-shrinking backtracking gradient-descent warm start, was added to
   get past this before any optmt polish. It reliably improves the
   OBJECTIVE (QUAIDS: -1.8 to -6.8 to -8.2 across runs with more steps;
   AIDS: 12.9 to -28.7 to -29.5) but does **not** reliably improve
   PARAMETER RECOVERY -- in several runs (both models, confirmed
   reproducibly across three independent step-budget settings: 80, 400,
   and 1500 steps), continuing to improve `homogCrit` moved the fitted
   `b` FARTHER from the true parameters, and `optmt`'s subsequent polish
   step was observed to blow recErr up by three to four orders of
   magnitude from a perfectly reasonable warm-started point (once for
   QUAIDS at 1e8/3671/168478 across runs, once even for AIDS at 531 in
   the 1500-step run) -- genuine, reproduced run-to-run instability in
   optmt's own search, not a one-off fluke.

**This is direct, non-fixed-point-iteration evidence for the same
weak-identification/multiple-stationary-point concern this whole research
line has flagged from the Anderson/multi-start side**: improving this
fit criterion does not monotonically improve parameter recovery, so
"optimize harder" is not by itself a fix, independent of which search
algorithm is used to do the optimizing. It does not by itself prove weak
identification (a finite-sample objective can genuinely have a better,
wrong-parameter optimum -- the same caveat the earlier 200-seed benchmark
already carries), but it closes off "the fixed-point iteration is just a
bad optimizer for a well-behaved objective" as an explanation: a properly
validated, independently-gradient-checked direct optimizer hits the same
wall.

### Recommendation after Phase 1

Do not proceed to Phase 2's matched 30-seed pilot comparing safe Anderson
against this exact direct-optimization recipe (plain optmt, or
optmt+gradient-descent warm start) with an expectation it will
outperform -- the evidence above suggests it likely will not for QUAIDS,
and is not even reliably stable for AIDS. Two more promising directions,
both already anticipated by the original work order:

- **Phase 3's continuation/homotopy approach** is now more clearly
  motivated than before this phase ran, not less: gradually introducing
  the translog/quadratic nonlinearity from an already-solved simpler
  stage, rather than directly searching the full wild nonlinear
  landscape from a cold or lightly-perturbed start, directly targets the
  specific severe-local-ill-conditioning mechanism found in step 3 above.
- A properly derived analytic (or semi-analytic) Gauss-Newton/LM
  objective -- using the actual residual Jacobian structure rather than a
  generic quasi-Newton search on a scalar `ln(det(.))` surrogate -- may
  handle the curvature far better than either optmt or hand-rolled
  gradient descent; this phase deliberately used numerical derivatives
  throughout (the original task's own explicitly allowed starting point)
  and did not yet attempt this.

Either direction should be smoke-tested on this SAME validated objective/
fixture pair before any further Monte Carlo sweep, continuing the
stage-gate discipline this phase followed.

## Phase 1b: audit of the Phase 1 recovery-gate failure (2026-10-01)

Requested after Phase 1's own QUAIDS recovery-gate failure, specifically
to settle whether that failure reflected something real about the model
(weak identification, a genuinely hard optimization landscape) or
confounds in how Phase 1 tested it, BEFORE spending further effort on
Phase 3's continuation/homotopy work. **Bottom line: confounds, not a
hard landscape.** With both confounds fixed, direct optimization recovers
the true structural parameters reliably, for both models, from
perturbations up to 5%.

**New artifacts**: `direct_objective.src` gained roughly a dozen new
procs (see its own "Phase 1b audit additions" header comment for the
full list and the math each one implements) -- a shared `_directObjCore`
returning (objective, feasibility flag, residual matrix) instead of a
flat-penalty scalar; `_directResidJacobianFD` (central-FD residual
Jacobian); `_directGNGradient`/`_directGNHessian` (the CORRECT
Sinv-weighted gradient/Gauss-Newton-Hessian structure for the profiled
log-determinant objective, not a naive unweighted-NLS gradient);
`_directLM` (a real Levenberg-Marquardt solver: rejects infeasible
trial steps, requires an actual/predicted-reduction test to accept a
step, never reports convergence without both a feasible point and a
small gradient); `_directBlockErrors` (per-block recovery reporting);
`_directJacobianSvdReport`/`_directWeakDirectionReport` (rank/condition-
number/weak-direction diagnostics); `_directLineProfile`; and two new
fixtures, `_quaidsStatFixtureDGP` (corrected finite-sample statistical
fixture) and `_quaidsStationaryFixtureDGP` (a NEAR-stationary-point
numerical fixture -- not exactly stationary, see below -- built by
orthogonalizing equation noise against the true residual Jacobian).
`phase1b_audit.e` (the audit
script, both models, items 1-5 below) and its saved output
`phase1b_audit_output.txt`; `phase1b_smoke.e` (a scratch script that
found/fixed two real bugs before the full audit ran -- kept as evidence,
not part of the audit record itself).

### The two confounds, found and fixed

1. **The fixture's first-stage residual was weakly identified by
   construction.** Phase 1's `_quaidsNoiselessDGP()` shrank totexp's own
   noise (`epsScale`) to keep the GLS criterion well-conditioned, but
   totexp's noise IS the first-stage residual `u` becomes after
   `quaidsFit()`'s own IV stage -- shrinking it made `u` a near-
   degenerate (low-variance) regressor, which inherently gives its own
   coefficient a high-variance estimate in ANY finite sample, independent
   of whether the true endogeneity correlation is 0. `_quaidsStatFixtureDGP()`
   decouples this: `u` is drawn at a normal, informative scale (`.1*rndns`,
   matching `tests/quaidsfixtures.src`'s own convention) while the
   EQUATION-level noise stays small. `_quaidsStationaryFixtureDGP()` goes
   further: it orthogonalizes equation noise against the TRUE residual
   Jacobian before adding it, making truth a NEAR-stationary point of the
   objective, not an exactly stationary one -- the observed
   ||FD gradient||_inf at truth is 0.084 for AIDS and 0.183 for QUAIDS,
   small relative to gradient norms of several thousand seen at
   non-stationary points in Phase 1's own diagnostic, but genuinely
   nonzero, not "zero up to FD precision" (an earlier draft of this
   section overstated this -- corrected here). The residual gap has at
   least two real sources, not just FD truncation in the Jacobian used to
   orthogonalize: `_directObjCore()`'s own feasibility/objective
   evaluation uses the CLIPPED QUAIDS deflator exponent and quaidsFit()'s
   OWN fitted `u` (from a real IV first-stage regression on the
   constructed data), not the DGP's own `uResid` the orthogonalization
   was actually computed against -- neither gap was separately quantified
   in this phase, so the 0.084/0.183 figures should be read as "small
   relative to a non-stationary point" rather than "attributable
   entirely to one named cause".
2. **The optimizer had no real safeguards.** Phase 1 used a flat `1e8`
   penalty for infeasible points (a gradient-less plateau that can look
   like "convergence" to a generic optimizer) and either plain `optmt` or
   a hand-rolled gradient-descent warm start with no curvature
   information. `_directLM()` fixes this properly: it rejects (not
   penalizes) infeasible trial steps, uses the Sinv-weighted Gauss-Newton
   Hessian (always PSD, so Levenberg-Marquardt damping `B+lambda*I` is
   always well-defined), requires an actual/predicted-reduction ratio
   `rho > 0.1` to accept a step, and only reports convergence (retcode 0)
   when BOTH the current point is feasible AND the gradient's inf-norm is
   below tolerance.

The semi-analytic Gauss-Newton gradient/Hessian were validated against
the already-trusted pure finite-difference gradient before being used
for anything else, and this cross-check is now part of
`phase1b_audit.e`'s own canonical, saved-output run (promoted from the
scratch `phase1b_smoke.e`, per the Phase 1b cleanup below) for BOTH
models: max abs diff `2.19e-9` (AIDS) / `2.57e-9` (QUAIDS) at the true
parameters, and `4.26e-9` (AIDS) / `2.37e-7` (QUAIDS) at a
0.2-scaled-unit perturbed (non-stationary) point -- the "correct
gradient structure" (`g_k=(2/nobs)*tr(Sinv*R'*J_k)`, an Sinv-weighted
inner product, NOT the plain `R'*J_k` an unweighted-NLS gradient would
use) is confirmed correct, not merely plausible.

**Phase 1b cleanup (2026-10-01), before proceeding to Phase 3:** a
second pass over this section's own claims, requested explicitly,
corrected five issues:
1. "2 independent draws each" (Items 4-5 below) could be misread as 8
   independent DATASETS. It is 2 independent STARTING-POINT
   perturbations per model/scale cell, all on ONE stationary-fixture
   dataset per model (`seed=4242`) -- corrected throughout this file,
   `CLAUDE_HANDOFF.md`, and `PROJECT_STATUS.md`.
2. `_quaidsStationaryFixtureDGP()` was described as making truth "exactly
   stationary" / the gradient "zero up to FD precision". The ACTUAL
   observed ||FD gradient||_inf at truth is 0.084 (AIDS) / 0.183 (QUAIDS)
   -- genuinely nonzero. It is now called a NEAR-stationary fixture
   throughout, with the real numbers reported rather than asserted away.
3. The rank/weak-direction conclusion ("no rank deficiency... weak
   identification confined to the nuisance row") is explicitly scoped
   below to this fixture, `tobs=500`, one fixture draw per model, and
   `_directObjScale()`'s own documented parameter scaling -- not a
   claim about any other sample size, draw, or scaling choice.
4. The semi-analytic-vs-FD gradient cross-check (previously only in the
   scratch `phase1b_smoke.e`) is now part of `phase1b_audit.e` itself and
   its saved `phase1b_audit_output.txt`, for both models.
5. `_quaidsStationaryFixtureDGP()`'s noise-orthogonalization originally
   used `invpd(Jtrue'Jtrue)` (a normal-equation projection, which squares
   the Jacobian's own condition number -- observed ~1e3-4e3 at truth, so
   ~1e6-1.6e7 squared, still safely invertible in double precision, but
   not the best available option). Replaced with an orthonormal
   projection via `svd2()`'s own left singular vectors (no matrix
   inversion at all) -- the same rank-revealing decomposition this file's
   diagnostics already use elsewhere. Re-ran `phase1b_audit.e` after the
   swap: all 8 Items-4-5 trials still PASS with materially unchanged
   numbers (see `phase1b_audit_output.txt`), confirming the original
   `invpd` form was not silently driving the result.

### Item 1: which block actually drove Phase 1's 0.225 QUAIDS error

Reproducing Phase 1's exact fixture/seeds/method (the flawed
`_quaidsNoiselessDGP()` + gradient-descent-only warm start, no LM) and
running `_directBlockErrors()` on the resulting endpoint: `gamma` block
error 0.109, `alpha` block error 0.072, `beta` row errors up to 0.039,
`u`-row error 0.049 -- **error was spread across the STRUCTURAL blocks,
not cleanly isolated to the nuisance `u` row.** This is a more precise,
and more informative, finding than "it was the weak u-row": Phase 1's
plain gradient-descent warm start (no curvature information, no trust
region) simply had not converged the structural parameters either within
its iteration budget, on top of whatever the fixture's own u-row
confound contributed. This matches Items 2-5's own finding that a PROPER
optimizer (not the fixture fix alone) is necessary for reliable recovery.

### Items 2-3: fixture validation and local identification diagnostics

Both fixtures are feasible and well-conditioned at truth for both
models. Confirmed at `_quaidsStationaryFixtureDGP()`'s own truth point,
`tobs=500`:

| | AIDS | QUAIDS |
|---|---:|---:|
| Objective at truth | -9.78 | -9.79 |
| \|\|FD gradient\|\|_inf at truth | 0.084 | 0.183 |
| Jacobian shape | 2000 x 32 | 2000 x 36 |
| Jacobian numerical rank | 32 of 32 (FULL) | 36 of 36 (FULL) |
| Jacobian condition number | 980 | 3988 |
| GN-Hessian condition number | 1.17e6 | 2.12e7 |

**No rank deficiency in either model** -- the Jacobian is always full
column rank. QUAIDS is genuinely more ill-conditioned than AIDS (an
~18x larger Hessian condition number), consistent with every other
finding in this research line that QUAIDS's quadratic/lambda deflator
is harder than AIDS's, but "harder" here is quantified, not qualitative:
no near-null (truly flat) direction was found in either model.

The three smallest-singular-value directions of the Jacobian, in BOTH
models, carry essentially all their loading (>99.999%) on the `u`-row
block and essentially none on alpha/gamma/beta/lambda -- a clean, direct
confirmation that the one genuinely weakly-informed direction in this
problem is the nuisance first-stage correction, not the structural
economics. This is the Phase 1b audit's single clearest, most
load-bearing finding: **weak identification, where it exists at all
here, is confined to a nuisance row the fixture was itself responsible
for grading, not the structural parameters the convergence-correctness
research actually cares about.** It should not be read as "QUAIDS has no
identification problems anywhere" -- only as "this specific audit found
no rank deficiency or near-null direction outside the nuisance row, at
this sample size (tobs=500), this one fixture draw per model, and this
file's own documented parameter scaling (`_directObjScale()`'s
Stone-start-magnitude-based, 0.05-floored per-element scale)."

Line profiles (both models) from a 5%-perturbed start to truth are
smooth and monotonically improving at every one of 7 sampled points --
no sign of a secondary basin or discontinuity on that path. A profile
along the weakest Jacobian direction from truth (+/-2 scaled units) is
gently curved, not flat -- the objective moves by only about 0.001-0.015
over that range (confirming the direction really is low-curvature) but
never stops responding entirely (confirming it is ill-conditioned, not
literally unidentified).

### Items 4-5: the corrected local recovery gate

`_directLM()`, run from 2 independent STARTING-POINT perturbations at
each of 1% and 5% relative scale, for both models (8 trials total).
**Correction**: this is 2 independent perturbed starts per model/scale
combination, all on the SAME single stationary-fixture dataset per
model (one dataset draw per model, `seed=4242`, not independent data
redraws per trial) -- only the perturbation applied to `bRawTrue` to
build each starting point differs between trials, not the underlying
`w`/`prices`/`totexp`/`instr`. An earlier draft of this section said "2
independent draws" without that distinction, which could be misread as
independent datasets; this is now corrected throughout `SUMMARY.md`/
`CLAUDE_HANDOFF.md`/`PROJECT_STATUS.md`.

| | Every one of 8 trials |
|---|---|
| Feasible start | Yes |
| LM retcode | 0 (converged: feasible + small gradient), every trial |
| Objective gap from truth | ~1e-4, every trial |
| Structural (excl. u) recovery error | 0.0011-0.0017, every trial |
| u-row recovery error | 0.006-0.0087, every trial (reported separately, not gated) |
| Gradient norm at convergence | 1.3e-4 to 8.5e-4 (tolerance 1e-3) |

Every trial PASSES the gate as specified: a finite feasible endpoint, the
objective within 1 of the truth-level objective (in practice ~1e-4, far
tighter), a small gradient, retcode 0, and structural blocks recovered to
better than 0.05 max-abs-diff (in practice ~0.001-0.002). Recovery is
essentially IDENTICAL across both 1% and 5% starts and across the two
independent starting perturbations per cell -- strong evidence of a
single stable basin for the structural parameters FROM THESE STARTS, on
this one fixture draw per model. It is not evidence about a different
data draw or a starting point far from truth (see the Decision section
below).

### Decision and recommendation

Per the audit's own decision rule: **local LM/trust-region recovery
succeeded, and no rank deficiency was found** -- so the recommended next
step is Phase 3's continuation/homotopy work, using `_directLM()` as the
validated solver at every stage, NOT a pivot to documenting weak
identification. The audit does not establish that continuation is
unnecessary, however: `_directLM()` was only tested from already-close
(1-5%) perturbations at one sample size (tobs=500) on a fixture
specifically constructed so truth is a NEAR-stationary point (observed
gradient norm 0.084/0.183, not zero) -- exactly the easy case
continuation would also handle easily. It has NOT
been tested from the kind of cold/far/structured starts
(`anderson_full_multistart_benchmark.e`'s own 16-recipe portfolio, the
Stone start itself, real survey-scale sample sizes with real noise) that
actually caused the original convergence-correctness problem this whole
research line exists to fix. That comparison -- `_directLM()` (or a
continuation-staged version of it) against safeguarded Anderson, on the
ORIGINAL noisy fixtures and harder starts -- is the next real test, not
yet run.

## Phase 3: continuation/homotopy smoke gate (2026-10-02)

> **Correction (Phase 3 hardening, 2026-10-05) -- read "Phase 3
> hardening" below before relying on this section.** (1) The
> "structErr (excl. u)" column in this section's table used a wrong
> row slice (`1:nIntFull+nReduced+2-modelFlag`): for AIDS it INCLUDED
> the u row (so AIDS "structErr" 0.722/0.724 was really the u-row
> error; the true structural error is 0.058), and for QUAIDS it EXCLUDED
> the lambda row (QUAIDS values are unchanged after the fix because the
> lambda errors are tiny). (2) "IDENTICAL"/"EXACT SAME solution" for the
> two QUAIDS continuation endpoints was inferred from printed summaries;
> the explicit vector difference is 5.8e-9 (scaled parameters) / 7.6e-10
> (raw structural coefficients), i.e. numerically equivalent, not
> bit-identical. (3) The original verdict did not check the "every
> accepted stage feasible" claim it printed and had no absolute
> correctness condition. The original text below is kept as recorded.

Answers exactly the open question Phase 1b left: does a validated direct
optimizer, staged via continuation, work from far/structured starts on
REAL noisy data (not Phase 1b's near-stationary fixture)? **Yes, on
this tightly bounded smoke gate -- all 4 tested scenarios pass, and the
harder (QUAIDS) case shows continuation reaching a dramatically better,
remarkably consistent solution from two very different starts, while
cold direct LM and safe Anderson both do not.**

**New artifacts**: `direct_objective.src` gained a `homotopy` parameter
threaded through the objective/Jacobian/solver stack (see its own
"Phase 3 continuation/homotopy" header comments), plus a robustness fix
(trapped `invpd()` failures, found necessary by this exact smoke gate --
see below). `phase3_smoke_gate.e` / `phase3_smoke_gate_output.txt`: the
smoke gate script and its saved canonical output.

### The continuation path

ONE continuation parameter `homotopy` (t in [0,1]) scales the two
sources of self-referential nonlinearity in this objective: the translog
price index's quadratic (gamma-driven) term, and (QUAIDS only) the
exponential deflator's exponent. At t=0 the price index is purely linear
in relative prices (gamma still enters as an ordinary linear regression
coefficient, just not inside the nonlinear index) and the QUAIDS deflator
is neutral (exp(0)=1). At t=1 this is the FULL original objective.

This is exact by CONSTRUCTION, not merely by testing: `_directObjResid()`/
`_directObjCore()`/`_directResidJacobianFD()`/`_directLM()` are now thin
wrappers around their own `*Homotopy` versions at `homotopy=1` -- there is
no second, independently-written implementation that could silently
drift out of sync. `phase3_smoke_gate.e`'s own "Step 0" still checks this
numerically (max abs diff 0, both the scalar objective and the full
residual matrix) before trusting anything else, since "by construction"
is only as good as the construction.

The continuation DRIVER (`_phase3Continuation`, defined in the smoke-gate
script itself, not `direct_objective.src` -- its bookkeeping is specific
to this gate's own reporting needs) starts at `homotopy=0`, step=0.25,
and calls `_directLMHomotopy()` once per stage, warm-starting each stage
from the previous one's endpoint. A stage is accepted only if
`_directLMHomotopy()` itself returns retcode 0 (feasible AND converged);
on any other retcode the step is HALVED and retried from the last GOOD
stage (never silently advanced), down to a minimum step of 0.01 before
giving up.

### A real robustness bug found and fixed by this gate

The QUAIDS far/structured start triggered a genuine crash on the first
attempt: `error G0121 : Matrix not positive definite` inside
`_directLMHomotopy()`. Root cause: `_directObjCoreHomotopy()`'s
feasibility test only checks the SIGN of `det(residCov)`, which a wild
far-start trial step can satisfy (a tiny positive determinant) while
`residCov` is still too ill-conditioned for `invpd()` to handle --
`invpd()` does not degrade gracefully in that case, it raises a hard
error that would have crashed the whole gate, not just rejected one
step. Fixed by wrapping both `invpd()` calls in `_directLMHomotopy()`
with the same trap-based pattern this codebase already uses elsewhere
for exactly this class of failure (`_andersonSvdMix()`'s own `svd2()`
trap) -- a trapped failure is now treated as a rejected step (increase
damping, retry), never a crash, and never silently treated as success.
Confirmed this fix does not change Phase 1b's own results at all (full
byte-for-byte `phase1b_audit.e` regression check, run twice -- once
after the homotopy refactor alone, once after this robustness fix --
both identical to the original `phase1b_audit_output.txt` except
timestamps).

### Smoke gate results

Known-difficult seeds (from `anderson_full_multistart_benchmark_output.txt`'s
own "no safeguarded start converged" lines -- seeds where ALL 16
structured starts x safeguarded Anderson already failed in the completed
200-seed benchmark): **AIDS seed=6, QUAIDS seed=9**, `tobs=3000` (matching
that benchmark's own convention). Two starts per model: the Stone start,
and a structured "far" start (0.8x-relative-scale noise plus a gamma/
price-coefficient sign flip, in the same spirit as the 200-seed
benchmark's own recipe portfolio, reimplemented minimally in the smoke
-gate script rather than `#include`d from that full sweep script).

| Model/start | Method | homogCrit | structErr (excl. u) | overallErr | Notes |
|---|---|---:|---:|---:|---|
| AIDS/Stone | cold LM | -5.499 | 0.722 | 0.722 | converged |
| AIDS/Stone | continuation | -5.499 | 0.724 | 0.724 | reached homotopy=1, ties cold LM |
| AIDS/Stone | safe Anderson | -6.386 | 1.518 | 1.518 | did NOT converge (100 iters) |
| AIDS/far | cold LM | -5.499 | 0.724 | 0.724 | same basin as Stone start |
| AIDS/far | continuation | -5.499 | 0.724 | 0.724 | same basin, ties cold LM |
| AIDS/far | safe Anderson | -6.385 | 1.516 | 1.516 | did NOT converge |
| QUAIDS/Stone | cold LM | -11.345 | 4.019 | 4.019 | converged but WRONG |
| QUAIDS/Stone | continuation | **-5.577** | **0.081** | 0.515 | reached homotopy=1, dramatically better |
| QUAIDS/Stone | safe Anderson | -12.123 | 2.163 | 2.163 | did NOT converge |
| QUAIDS/far | cold LM | -11.109 | 12.291 | 12.291 | converged but WORSE than Stone-start cold LM |
| QUAIDS/far | continuation | **-5.577** | **0.081** | 0.515 | IDENTICAL to Stone-start continuation result |
| QUAIDS/far | safe Anderson | -5.774 | 0.193 | 0.544 | did NOT converge, but closest Anderson got |

For AIDS, continuation ties cold direct LM (which already worked fine
from both starts) and both clearly beat safe Anderson, which did not
converge from either start on this known-hard seed within 100
iterations. For QUAIDS -- the harder case -- continuation is dramatically
better than both alternatives and, notably, reaches the EXACT SAME
solution from two very different starting points (Stone vs. an
aggressively perturbed, sign-flipped far start), while cold LM reaches
a DIFFERENT (and for the far start, substantially worse) converged-but-
wrong point depending on the start. This is a real, direct demonstration
of exactly the mechanism continuation is meant to provide: gradually
introducing the nonlinearity avoids the bad local basin cold direct
optimization falls into, converging to what is likely the dominant/
correct basin regardless of starting point.

Cost: continuation takes meaningfully more iterations/evaluations than
cold LM when cold LM already succeeds (AIDS: ~1950 vs. a handful of
evaluations) but remains cheap in absolute wall-clock terms throughout
(under 6 seconds per solve on this machine, every scenario) -- "approx
total evals" is a clearly-labeled derived approximation
(iterations*(2k+1)+rejected steps), not an exact instrumented count (see
`direct_objective.src`'s own `_directLMHomotopy()` header for why an
exact count was not added).

### Verdict and decision

**All 4 scenarios PASS** the gate as specified: every accepted stage
feasible, continuation reaches the exact (by-construction) target
objective at homotopy=1, and performance materially improves over (or
ties) cold direct LM while not regressing against safe Anderson's best
result in any scenario. Per the user's own instruction ("do not launch
the 30-seed pilot until both model smoke cases pass or the failure
mechanism is clearly diagnosed"), **both model smoke cases passed
cleanly** -- this document stops here and does not launch a wider sweep;
that is the next session's decision to make, not this one's.

### Limitations, honestly

- One seed per model, two starts per model. This is a smoke gate, not a
  sweep -- it demonstrates continuation CAN help dramatically on a
  genuinely hard, real case, not that it reliably will across many seeds.
- The AIDS case did not actually need continuation (cold LM already
  converged to the same point) -- the smoke gate's positive result is
  carried almost entirely by the QUAIDS case. A wider sweep should
  expect, and separately report, a mix of "continuation unnecessary" and
  "continuation essential" outcomes, not assume every seed looks like
  QUAIDS seed 9.
- Safe Anderson's own 100-iteration budget matches the original 200-seed
  benchmark's convention, but Anderson was not re-run with its own
  16-start portfolio here (only the same single start each other method
  used) -- this is a fair single-start comparison, not a repeat of the
  full multi-start benchmark, and should not be read as "continuation
  beats the BEST the Anderson approach can do," only "beats Anderson
  from this one start."
- "Approx total evals" is a derived approximation, not an exact count
  (see above).

## Phase 3 hardening (2026-10-05)

A bounded hardening pass on the Phase 3 smoke gate, done before any
30-seed pilot. **Result: the strengthened smoke gate passes cleanly in
all 4 scenarios (2 models x 2 starts).** The qualitative Phase 3
conclusion stands; one reported number was wrong (AIDS structural
error) and is corrected below, and it got better, not worse.

Still one known-difficult seed per model (AIDS 6, QUAIDS 9) and two
starts per model. Nothing here says anything about population-wide
success rates.

### Code changes

- `direct_objective.src`
  - New `_directLMHomotopyDiag(..., homotopy, faultMode, faultIter,
    faultCount)` holds the LM solver body and returns a 7th output
    `diagVec = [covInvFail, stepInvFail, infeasTrial, rhoReject,
    accepted, finalLambda]`. `_directLMHomotopy()` is now a thin wrapper
    at `faultMode=0` with the same arguments and the same 6 return
    values, so `_directLM()`, `phase1b_audit.e` and `phase1b_smoke.e` are
    unchanged.
  - Fail-closed inversion: before each `invpd()` the destination is
    reset to missing, and the result is accepted only if it has the
    right shape and `not ismiss(...)`. A failed inverse can never reuse
    the previous iteration's inverse, and a failed damped-system inverse
    never moves theta.
  - New retcode **5**: covariance-inversion failure at the current
    (last accepted) point. It used to share retcode 4 ("stalled"). On
    retcode 5, `gnorm` is returned as -1, not as a stale value. No
    recorded Phase 1b or Phase 3 run ever hit this path.
  - Fault hook (diagnostic only): `faultMode=1` hands `invpd()` the
    matrix `-eye(n1)` in place of `residCov` at outer iteration
    `faultIter`; `faultMode=2` hands it `-eye(k)` in place of
    `B+lambda*I` for the first `faultCount` damped solves. Both make the
    real trap path run. No safeguard is relaxed.
- `phase3_smoke_gate.e` (rewritten)
  - **Stage zero:** the driver now solves the actual `homotopy=0`
    problem from the supplied start first, and advances only after a
    feasible, retcode-0 stage-zero endpoint exists. If stage zero fails,
    the driver stops without advancing. After that, the schedule is
    unchanged: step 0.25, warm start from the last accepted endpoint,
    step halving on any nonzero stage retcode, no advance after a failed
    stage, an honest stop (`minStepExhausted=1`) below step 0.01, and the
    last stage snapped to exactly 1.0. The old schedule (first stage at
    0.25) still runs as "legacy (b')" so the two can be compared
    directly.
  - Bookkeeping: `anyStageInfeasible` is renamed `anyStageFailedOnce`.
    New separate counters: failed stages by retcode (3 infeasible start,
    5 covariance inversion, 2 maxIter, 4 stalled), adaptive step
    reductions, `minStepExhausted`, and summed inner events (covariance
    and step inversion failures, infeasible trials, rho rejections).
    Each accepted stage's feasibility is re-evaluated independently
    (`allAcceptedFeasible`) rather than assumed.
  - **Corrected error decomposition** (`_phase3Errs`). The layout is
    `alpha (nIntFull) | gamma (n1) | beta | [lambda] | u (1 row)`, with a
    layout assertion at runtime (AIDS ng=8 with 7 structural rows,
    QUAIDS ng=9 with 8). The script reports structural (all rows except
    u), lambda, u, overall, and the old "pre-fix slice" so it can be
    traced back to the earlier output.
  - Explicit endpoint differences between the Stone and far starts, per
    model: max |Δθ| (scaled), max |Δ raw structural coefficient|, and
    |Δ objective|. "Identical" is printed only if |Δθ| = 0 exactly;
    "numerically equivalent" if both differences are ≤ 1e-6.
  - **Strengthened verdict.** A scenario is a clean pass only if ALL of
    these hold, and every one is printed:
    1. the endpoint reached homotopy exactly 1;
    2. stage zero converged genuinely;
    3. every accepted stage returned retcode 0;
    4. every accepted stage is feasible on re-evaluation;
    5. structErr ≤ 1 (primary correctness threshold);
    6. structErr ≤ cold-LM structErr + 0.05;
    7. structErr ≤ safe-Anderson structErr + 0.05.

    Rejected stages followed by step halving are allowed. A sensitivity
    table at structErr ≤ 0.5/1/2/5 is printed for all four methods.
- `phase3_safeguard_diag.e` / `_output.txt` (new): a bounded diagnostic
  of the fault hook (see below).

### Commands

From `dev/experiments/convergence-anderson/`:

```powershell
C:\gauss26\tgauss.exe -b -x phase3_smoke_gate.e      # -> phase3_smoke_gate_output.txt
C:\gauss26\tgauss.exe -b -x phase3_safeguard_diag.e  # -> phase3_safeguard_diag_output.txt (seconds)
C:\gauss26\tgauss.exe -b -x phase1b_audit.e          # regression check vs phase1b_audit_output.txt
```

Outputs were captured with `| Out-File -Encoding utf8 <file>`. In this
capture the GAUSS banner/timestamp goes to the console, not the file, so
the new output files start at the first line the script prints.

### Smoke gate results (corrected metrics)

structErr = max |raw coefficient error| over every row except the u
row, compared on the n1 estimated columns.

| Model/start | Method | homogCrit | structErr | lambdaErr | uErr | Notes |
|---|---|---:|---:|---:|---:|---|
| AIDS/Stone | cold LM | -5.49874 | 0.0583 | - | 0.722 | retcode 0 |
| AIDS/Stone | continuation (stage zero) | -5.49872 | 0.0584 | - | 0.724 | 5 stages, 33 iter, 0 failed stages |
| AIDS/Stone | continuation (legacy) | -5.49872 | 0.0584 | - | 0.724 | 4 stages, 30 iter |
| AIDS/Stone | safe Anderson | -6.386 | 1.518 | - | 0.672 | not converged (100 iter) |
| AIDS/far | cold LM | -5.49872 | 0.0584 | - | 0.724 | retcode 0 |
| AIDS/far | continuation (stage zero) | -5.49872 | 0.0584 | - | 0.724 | 5 stages, 34 iter |
| AIDS/far | continuation (legacy) | -5.49872 | 0.0584 | - | 0.724 | 4 stages, 30 iter |
| AIDS/far | safe Anderson | -6.385 | 1.516 | - | 0.673 | not converged |
| QUAIDS/Stone | cold LM | -11.345 | 4.019 | 0.020 | 0.298 | retcode 0, wrong basin |
| QUAIDS/Stone | continuation (stage zero) | **-5.5773** | **0.0815** | 5.8e-5 | 0.515 | 5 stages, 82 iter, 8 rho rejections |
| QUAIDS/Stone | continuation (legacy) | -5.5773 | 0.0815 | 5.8e-5 | 0.515 | 4 stages, 69 iter |
| QUAIDS/Stone | safe Anderson | -12.123 | 2.163 | 0.011 | 0.394 | not converged |
| QUAIDS/far | cold LM | -11.109 | 12.291 | 0.017 | 0.493 | retcode 0, wrong basin |
| QUAIDS/far | continuation (stage zero) | **-5.5773** | **0.0815** | 5.8e-5 | 0.515 | 5 stages, 72 iter, 9 rho rejections |
| QUAIDS/far | continuation (legacy) | -5.5773 | 0.0815 | 5.8e-5 | 0.515 | 4 stages, 103 iter |
| QUAIDS/far | safe Anderson | -5.774 | 0.193 | 0.001 | 0.544 | not converged |

In every continuation run, every stage was accepted on its first try:
no failed stages, no step reductions, `minStepExhausted=0`, and zero
covariance-inversion, step-inversion or infeasible-trial events. So the
step-halving and minimum-step paths were **not exercised** by this gate.
All continuation wall times were between 1.0 and 4.8 s.

Threshold sensitivity (structErr ≤ .5 / 1 / 2 / 5):

| | AIDS Stone | AIDS far | QUAIDS Stone | QUAIDS far |
|---|---|---|---|---|
| cold LM | 1/1/1/1 | 1/1/1/1 | 0/0/0/1 | 0/0/0/0 |
| continuation (stage zero) | 1/1/1/1 | 1/1/1/1 | 1/1/1/1 | 1/1/1/1 |
| continuation (legacy) | 1/1/1/1 | 1/1/1/1 | 1/1/1/1 | 1/1/1/1 |
| safe Anderson | 0/0/1/1 | 0/0/1/1 | 0/0/0/1 | 1/1/1/1 |

**Verdict:** every predicate holds in all 4 scenarios, so the output
reports `SMOKE GATE: PASS (all scenarios clean)`, 4 of 4.

### Changed numerical conclusions

- **AIDS structural error is 0.058, not 0.724.** The 0.722/0.724 figures
  were the u-row error, picked up by the wrong slice. AIDS cold LM and
  continuation are correct even at the 0.5 threshold.
- QUAIDS numbers are unchanged after the fix: lambda-row errors are
  ≤ 0.02, so excluding lambda had not mattered. The same holds for safe
  Anderson's 1.518 / 2.163 / 0.193.
- The u row carries most of the overall error for every converged
  continuation and AIDS cold-LM endpoint (0.52–0.72). This is consistent
  with Phase 1b's finding that the u row is the weakest direction on its
  fixture. It is the reason structErr excludes u.

### Endpoint equality (explicit vectors)

| Model | Pair | max \|Δθ\| scaled | max \|Δ raw structural\| | \|Δ objective\| | Classification |
|---|---|---:|---:|---:|---|
| AIDS | continuation, Stone vs far | 3.6e-9 | 2.8e-10 | 1.5e-14 | numerically equivalent (≤ 1e-6), not identical |
| AIDS | legacy continuation, Stone vs far | 1.3e-8 | 3.0e-10 | 7.2e-14 | numerically equivalent |
| AIDS | cold LM, Stone vs far | 0.21 | 4.9e-4 | 1.3e-5 | different stopping points in the same basin |
| QUAIDS | continuation, Stone vs far | 5.8e-9 | 7.6e-10 | 6.2e-15 | numerically equivalent (≤ 1e-6), not identical |
| QUAIDS | legacy continuation, Stone vs far | 6.0e-9 | 1.5e-9 | 8.9e-15 | numerically equivalent |
| QUAIDS | cold LM, Stone vs far | 8.64 | 8.27 | 0.236 | different basins |

### Did stage zero change the results?

**No endpoint changed.** Stage-zero and legacy endpoints differ by at
most 1.5e-8 in scaled θ, 8.2e-10 in raw structural coefficients, and
8.1e-14 in the objective, in all 4 scenarios. The legacy runs reproduce
the 2026-10-02 per-stage iteration counts and objective values exactly,
which also shows that the solver refactor did not change the arithmetic.

Stage zero changed cost, not the answer, and the change was mixed:
- AIDS: 33–34 vs 30 total iterations.
- QUAIDS Stone: 82 vs 69 iterations.
- QUAIDS far: 72 vs 103 iterations, with 9 vs 13 rejected steps.

On this evidence, stage zero adds robustness by construction (no stage
starts from an unsolved point) but did not alter any outcome here.

### Safeguard diagnostic (`phase3_safeguard_diag.e`)

The diagnostic uses the AIDS fixture at tobs=500, seed 1, from the Stone
start. **14 of 14 checks pass.**
- **T1 (covariance inversion fails at iteration 2,** after a valid
  inverse existed at iteration 1). The failure is trapped and the
  counter reads 1. The solver returns retcode 5 at iteration 2, and the
  returned theta/fVal equal the iteration-1 accepted point exactly. That
  is the evidence that no step was taken with the stale inverse. `gnorm`
  is -1, and fVal is the returned point's own objective.
- **T2 (covariance inversion fails at iteration 1):** retcode 5, zero
  accepted steps, and theta/fVal unchanged from the start.
- **T3 (every damped-system inversion fails):** 25 trapped failures,
  matching the damping schedule from 1e-3 to 1e12 by ×4. The solver
  makes zero accepted steps, theta/fVal are unchanged, it never
  evaluates an objective trial, and it returns retcode 4.
- **T4 (first 3 damped inversions fail):** exactly 3 trapped failures,
  each counted as a rejected step. The solver then retries with larger
  damping and converges genuinely (retcode 0, gnorm 8.4e-4 < 1e-3,
  feasible).
  - The first version of this script also required
    |f − f_baseline| < 1e-6. That check **failed** (|Δf| = 1.2e-3): the
    changed damping path stops at a different point that also satisfies
    gradTol=1e-3, differing by 0.077 in the u row and 0.018 in the
    structural rows.
  - Re-polishing both endpoints at gradTol 1e-6 brings them together
    (|Δf| 5.8e-14, max |Δθ| 1.1e-5). So the gap was the loose stopping
    tolerance, not the safeguard.
  - The check was replaced by the genuine-convergence check above, and
    the difference is printed as information. This is recorded here
    rather than silently dropped.
- **T0:** the wrapper `_directLMHomotopy` is bit-identical to
  `_directLMHomotopyDiag(faultMode 0)`, and no inversion failures occur
  without injection.

### Phase 1b regression

`phase1b_audit.e` was rerun after the refactor and compared line by line
with the saved `phase1b_audit_output.txt`. The only differences are the
GAUSS banner, job and timestamp lines (absent from the new capture, see
"Commands"). Every numerical line is identical. The saved canonical file
was left as is.

### Limitations

- One seed per model, two starts per model. This is a smoke gate, not
  evidence about success rates.
- The step-halving, minimum-step and stage-zero-failure paths of the
  driver were never triggered by real data here. Only their bookkeeping
  is in place, and the solver-level inversion paths are covered only by
  the fault-injection diagnostic.
- **The metric is not the benchmark's metric.** The threshold value
  structErr ≤ 1 is reused from the 200-seed benchmark, but that
  benchmark's `recErr` is `max|qOut.bS - trueParams|`: the recovered
  full structural matrix, including the restriction-implied n-th good.
  This gate's structErr is on the raw reduced coefficient rows except u.
  A 30-seed pilot should report the benchmark's own `bS`-based recErr
  too, so the two can be compared.
- gradTol=1e-3 leaves stopping-point differences of order 1e-3 in the
  objective along weak directions (T4 above). This is harmless for these
  correctness thresholds, but a pilot comparing objectives across
  methods should keep it in mind.
- Safe Anderson was still run single-start, not with its 16-start
  portfolio.
- "Approx total evals" is still a derived approximation.

### Next entry point

The hardened checkpoint awaits the user's review and an explicit commit
decision. After that, the next step is the 30-seed matched pilot, which
was **not started**. Its design is in CLAUDE_HANDOFF.md.

> Superseded: the hardened checkpoint was committed (`a776717`), and the
> pilot was then run and committed and pushed as `e90a568`. See
> "Matched 30-seed pilot (2026-10-05)" below.

## Matched 30-seed pilot (2026-10-05)

**Read this section first.** It is the latest result in this directory
and supersedes the "Next entry point" of the Phase 3 hardening section
above.

The pilot compares three methods on identical DGP seeds (1-30,
tobs=3000) and on the 200-seed benchmark's 16 start recipes, for both
iterated AIDS and QUAIDS:
- safeguarded Anderson with its full 16-start portfolio;
- cold direct LM;
- hardened continuation + LM (stage zero at homotopy 0).

That is 2,880 cells (model × seed × method × start). Every cell finished.

The definitions were frozen in `PILOT_PROTOCOL.md` before any result was
examined. The primary class is selected `bsErr` ≤ 1, the benchmark's own
`max|bS − trueParams|`. The secondary class is selected `rawStructErr`
≤ 1, which covers all reduced rows except u and includes lambda.

**Selection** takes the largest exact homotopy-1 objective among
genuinely converged endpoints. **Oracle** results are diagnostic only.

### Files

| File | Role |
|---|---|
| `PILOT_PROTOCOL.md` | Frozen protocol, tolerances, and 3 recorded deviations (2 validation-tolerance fixes, 1 observation). |
| `pilot_lib.src` | Verbatim copies of the benchmark's start recipes and of `_phase3Errs`; `_pilotContinuation` (a copy of `_phase3Continuation` with 2 extra info fields); `_pilotRecoverBS` (raw → `bS` through `quaidsFit()`'s own post-processing chain); CSV writer. |
| `pilot_harness_core.src` | Resumable cell loop. One CSV row is appended per finished cell; finished keys are skipped on rerun. |
| `pilot_validate.e` | Stage A Part 1 checks plus the harness on seeds 1, 6, 9 → `pilot_validation/`. |
| `pilot_resume_test.e`, `.ps1`, `_output.txt` | V6 resume test. |
| `pilot_run_m0.e`, `pilot_run_m1.e` | Stage B (seeds 1-30) → `pilot_out/`. |
| `pilot_analyze.py` | All selection, deduplication and aggregation from disk: `--selftest`, `--validate DIR`, `--compare A B`, `--report DIR`. |
| `pilot_out/pilot_report.md`, `pilot_seed_table.csv` | Full tables, including a 180-row seed table and the anomaly table. |
| `pilot_bs_diag.e`, `_output.txt` | Post-pilot diagnostic of where `bS` error comes from. |
| `pilot_posthoc_pooled.py`, `pilot_out/pilot_posthoc_pooled_output.txt` | Post-hoc exploratory pooling. Not prespecified. |

`pilot_out/` is 6.4 MB of raw CSV plus logs.

### Commands

Run from this directory. Outputs were captured as shown.

```powershell
C:\gauss26\tgauss.exe -b -x pilot_validate.e      # > pilot_validation\pilot_validate_console.txt
powershell -ExecutionPolicy Bypass -File pilot_resume_test.ps1
python pilot_analyze.py --validate pilot_validation   # > pilot_validation\pilot_analyze_validate_output.txt
C:\gauss26\tgauss.exe -b -x pilot_run_m0.e        # > pilot_out\pilot_run_m0_console.txt   (~40 min)
C:\gauss26\tgauss.exe -b -x pilot_run_m1.e        # > pilot_out\pilot_run_m1_console.txt   (~2.5 h)
python pilot_analyze.py --validate pilot_out          # > pilot_out\pilot_analyze_validate_output.txt
python pilot_analyze.py --compare pilot_validation pilot_out   # > pilot_out\determinism_check.txt
python pilot_analyze.py --report pilot_out
C:\gauss26\tgauss.exe -b -x pilot_bs_diag.e       # -> pilot_bs_diag_output.txt
python pilot_posthoc_pooled.py                     # > pilot_out\pilot_posthoc_pooled_output.txt
```

The two Stage B scripts ran concurrently as two processes, one per
model, on a 12-thread machine. The AIDS process finished first. Wall
times include that mild contention. Validation seed 1 QUAIDS took 418 s
there and 355 s in the pilot, so contention did not visibly inflate
timings.

### Stage A gate: PASSED, with no change to the comparison

| Check | Result |
|---|---|
| V1 `bS` conversion | `_pilotRecoverBS` reproduces `quaidsFit()`'s `homogB`, `bS`, `recErr` and `homogCrit` to ≤ 5e-10 in 12 of 12 cases (3 seeds × 2 models × 2 points), at the same (index, coefficient) pair. `_quaidsFitAnderson(mDepth 0)` and `quaidsFit()` give identical `bS`. A polished fixed point (AIDS seed 9) confirms that the direct (b,b) representation equals `quaidsFit()`'s own `bS` (Δ 3.7e-14). The other five polishing attempts diverged under plain iteration, which is the known instability of plain iteration, and were not scored. Across converged Anderson pilot cells, `bS(b,b)` vs `qOut.bS` differs by median 4e-7 and at most 0.021: the one-iteration index/coefficient gap at Anderson's 1e-4 stop. |
| V2 | Error-helper layout: u excluded, lambda included, column n ignored. |
| V3a | Start pack/unpack round trip ≤ 3.6e-15; starts are deterministic. |
| V3b | The pilot's Anderson arm reproduces the 200-seed benchmark's per-seed line on **all 60 seed/model cases**: converged starts, distinct endpoints, selected start, selected and oracle `recErr`, truth criterion. This proves start/RNG equivalence. |
| V3c | Start-1 cold-LM and continuation cells reproduce `phase3_smoke_gate_output.txt` (crit, structErr, iterations, rejected steps). |
| V4 | `critExact` equals the returned objective exactly for every LM and every continuation cell that reached homotopy 1. For converged Anderson endpoints, `critExact` vs `qOut.homogCrit` differs by < 1e-2. |
| V5 | Selection direction, eligibility, tie and dedup unit tests: 8 of 8. |
| V6 | Resume: 6 of 6. |
| V7 | Keys unique and complete; counter columns present. |
| Totals | `--validate`: 58 of 58 on `pilot_validation`, 377 of 377 on `pilot_out`. |
| Determinism | The 288 cells shared by `pilot_validation` and `pilot_out` are bit-identical apart from wall time. |
| Regression | `phase3_safeguard_diag.e` is unchanged (14 of 14; output identical apart from a BOM). `direct_objective.src` was not modified. |

Deviations (`PILOT_PROTOCOL.md`):
- Two reference-comparison tolerances (V3c, then V3b) were changed from
  an absolute 1e-6/5e-9 to a relative 1e-7. The cause was the 8
  significant digits of the printed reference values; every failing
  value matched in all printed digits.
- QUAIDS seed 1 Anderson: selection by `critExact` picks start 12 instead
  of the benchmark's start 7. The two are in the same basin, ~1e-9
  apart, and give the same error.
- No pilot definition changed.

### Results (30 seeds; never / wrong / correct)

| Model | Method | Stone-only bS | 16-start selected bS | oracle bS (diag.) | Stone-only raw | 16-start selected raw | converged starts |
|---|---|---|---|---|---|---|---:|
| AIDS | Anderson | 3/12/15 | 2/9/19 | 2/8/20 | 3/10/17 | 2/6/22 | 431/480 |
| AIDS | cold LM | 0/4/26 | 0/3/27 | 0/2/28 | 0/0/30 | 0/0/30 | 480/480 |
| AIDS | continuation | 0/4/26 | 0/4/26 | 0/4/26 | 0/0/30 | 0/0/30 | 480/480 |
| QUAIDS | Anderson | 5/17/8 | 1/19/10 | 1/19/10 | 5/15/10 | 1/16/13 | 372/480 |
| QUAIDS | cold LM | 4/9/17 | 0/7/23 | 0/7/23 | 4/3/23 | 0/1/29 | 379/480 |
| QUAIDS | continuation | 0/9/21 | 0/7/23 | 0/7/23 | 0/3/27 | 0/1/29 | 476/480 |

Anderson's selected correct rate on these 30 seeds (AIDS 19/30, QUAIDS
10/30) matches its 200-seed rate (60.5%, 34.5%), so the seeds are not
unusually easy or hard for Anderson.

**Threshold sensitivity** (selected error ≤ .5 / 1 / 2 / 5):

| Model | Method | bS | raw |
|---|---|---|---|
| AIDS | Anderson | 4/19/21/25 | 19/22/23/27 |
| AIDS | cold LM | 6/27/30/30 | 30/30/30/30 |
| AIDS | continuation | 6/26/30/30 | 30/30/30/30 |
| QUAIDS | Anderson | 3/10/13/14 | 13/13/14/19 |
| QUAIDS | cold LM | 8/23/28/29 | 29/29/29/29 |
| QUAIDS | continuation | 7/23/28/30 | 29/29/30/30 |

**Paired transitions** (primary bS class):

| Model | Pair | Improved to correct | Regressed |
|---|---|---|---|
| AIDS | Anderson → continuation | 8 (seeds 5, 6, 10, 13, 15, 16, 17, 24) | 1 (26) |
| AIDS | cold LM → continuation | 0 | 1 (26) |
| QUAIDS | Anderson → continuation | 14 (1, 3, 4, 6, 9, 10, 11, 13, 17, 18, 21, 23, 24, 28) | 1 (29) |
| QUAIDS | cold LM → continuation | 1 (17) | 1 (5) |

Anderson → cold LM is the same as Anderson → continuation in both
models, except that seed 26 does not regress for AIDS.

Both "regressions" against Anderson are u-row cutoff effects, not wrong
basins:
- AIDS seed 26 has continuation bsErr 1.012 against cold LM ≤ 1, and raw
  0.085 for both.
- QUAIDS seed 29 has continuation bsErr 1.08, raw 0.036; the u row is
  1.08.

**Objective-selected vs oracle.** Solution selection is **not** the
remaining problem. Selected and oracle classes agree everywhere except
two AIDS cases:
- AIDS seed 10, Anderson: selected 4.03, while the oracle reached 0.999.
- AIDS seed 22, cold LM: selected 1.006 vs oracle 0.994, a cutoff
  margin.

In QUAIDS no method has a selection failure. The continuation
selected-minus-oracle gap is at most 5.8e-6 in QUAIDS and 0.116 in AIDS.

**Cost** (wall seconds per seed; min / median / 90% / max):

| Model | Method | Stone-only | 16-start | Total over 30 seeds |
|---|---|---|---|---:|
| AIDS | Anderson | 0.006 / 0.014 / 0.045 / 0.14 | 0.13 / 0.37 / 0.69 / 2.1 | 14.5 s |
| AIDS | cold LM | 0.19 / 0.32 / 0.54 / 4.9 | 4.3 / 6.1 / 9.0 / 15 | 205 s |
| AIDS | continuation | 1.0 / 1.8 / 3.9 / 9.9 | 18 / 33 / 53 / 112 | 1,150 s |
| QUAIDS | Anderson | 0.009 / 0.026 / 0.065 / 0.40 | 0.15 / 0.49 / 0.89 / 1.5 | 15.5 s |
| QUAIDS | cold LM | 0.00 / 0.40 / 1.6 / 2.5 | 1.1 / 8.2 / 19 / 23 | 323 s |
| QUAIDS | continuation | 1.6 / 6.2 / 25 / 26 | 28 / 98 / 338 / 448 | 4,310 s |

Median iterations per cell (QUAIDS) are 27 for Anderson, 8 for cold LM
and 106 for continuation (all stages). The derived approximate
evaluation counts are 27, 584 and 7,750.

**Failure and safeguard counters.** There were no inversion failures
in AIDS, and no failed or reduced stages.

QUAIDS cold LM status counts:

| Retcode | Meaning | Cells |
|---:|---|---:|
| 0 | converged | 379 |
| 3 | infeasible start at homotopy 1 | 53 |
| 4 | stalled | 4 |
| 5 | covariance inversion failed | 44 |

QUAIDS cold LM inner events: 44 covariance-inversion failures, 178
step-inversion failures, 28 infeasible trials and 558 rho rejections.

QUAIDS continuation:
- **Status:** 476 genuine; 2 with stage zero failing (seed 21, starts
  2 and 7); 2 with the minimum step exhausted (seed 17, starts 9 and 10,
  which stopped at homotopy .969 and .812).
- **Stages:** 10 step reductions and 12 stalled (rc 4) stages; no rc
  2/3/5 stage failures.
- **Inner events:** 85 step-inversion failures, 3 infeasible trials,
  5,144 rho rejections.
- **First real-data exercise:** the step-halving, minimum-step and
  stage-zero-failure paths had never been triggered on real data before
  this pilot. They all behaved as designed, reporting honestly and never
  advancing past a failed stage.

The `bS` chain (`_pilotRecoverBS`) failed on 3 QUAIDS cold-LM cells (J
not invertible); those cells keep their raw error, and their `bsErr` is
missing. The symmetry projection used `quaidsFit()`'s own unprojected
fallback on 101 cold-LM, 2 continuation and 92 Anderson cells (QUAIDS),
and on 3 Anderson cells (AIDS).

**Basins** (raw structural, tolerance 0.05):
- **Continuation** finds a single converged basin on 30 of 30 AIDS seeds
  and 28 of 30 QUAIDS seeds.
- **Cold LM** finds more than one basin on 3 AIDS and 10 QUAIDS seeds.
- **Anderson** finds more than one basin on 7 AIDS and 14 QUAIDS seeds.

### Why bS and raw disagree (`pilot_bs_diag.e`)

Every direct-method seed that is "wrong" on bS but correct on raw sits
in the right basin:
- raw structural error ≤ 0.09;
- the selected objective is **better than the truth's** on the same
  sample, which is finite-sample overfit at the right optimum.

The diagnostic splits the `bS` error by block. It compares the
pilot's projected `bS` with an unprojected recovery (homogeneity and
adding-up only) of the same raw endpoint. The unprojected formula
agrees with `_pilotRecoverBS` to 5.7e-14 on a symmetry-fallback cell.

| Cause | Seeds | Effect |
|---|---|---|
| Symmetry projection | QUAIDS 15, 16, 20 | Projected `bS` error is 1.72, 2.18 and 1.15, all in the alpha block. Unprojected it is 0.38, 0.33 and 0.54. `quaidsFit()`'s delta-method projection metric is evaluated at a point that is not a fixed point of its map, and it moves alpha by 1-2. |
| Projection and u row | AIDS 29 | Projected `bS` error is 1.19 (gamma block). Unprojected it is 1.10, from the u row. |
| u row ≈ 1.0-1.2 | AIDS 22, 25, 26; QUAIDS 22, 25, 29 | The nuisance u row is at the cutoff, projected or not. |

So the benchmark-compatible `bS ≤ 1` rule understates the direct methods
on 6 of 7 QUAIDS seeds and on all 3 or 4 AIDS seeds where they are
classed wrong. The `bS` conversion itself is correct: it reproduces
`quaidsFit()` exactly (V1). The issue is that the projection's weighting
matrix comes from the fixed-point estimator's delta method. That is not
the right covariance for a direct-optimizer endpoint.

### Genuine wrong-basin failures of the direct methods

Only two seeds have direct-method failures with raw error > 1:
- **QUAIDS seed 5 (continuation).** All 16 starts converge to the
  **same wrong basin** (critExact −8.04, raw 1.31). Cold LM reaches the
  right basin (−5.52, raw 0.045) from 10 of 16 starts. The homotopy path
  itself leads to the wrong basin, so the start portfolio cannot help.
- **QUAIDS seed 17 (cold LM).** Only 1 of 16 starts converges, and it
  lands in a wrong basin (−17.33, raw 49). Continuation converges from
  14 starts:
  - 13 of them land in that same wrong basin;
  - only start 5 finds the right one (−5.61, raw 0.051);
  - the objective correctly selects start 5.

Continuation's start-independence works in both directions. When the
path is right, 15 extra starts are redundant. When it is wrong (seed 5),
they all fail together. Cold LM's failures are different ones.

### Post-hoc exploratory pooling (NOT prespecified)

`pilot_posthoc_pooled.py` applies the same selection rule to the pooled
converged candidates of several configurations.

| QUAIDS configuration | bS ≤ .5/1/2/5 | raw ≤ .5/1/2/5 | median / max wall per seed |
|---|---|---|---|
| cold LM, Stone only | 6/17/22/23 | 23/23/24/25 | 0.4 / 2.5 s |
| continuation, Stone only | 7/21/26/28 | 27/27/28/28 | 6.8 / 26 s |
| cold LM + continuation, Stone only (2 solves) | 8/23/28/29 | 29/29/29/29 | 7.5 / 29 s |
| cold LM, 16 starts | 8/23/28/29 | 29/29/29/29 | 8.3 / 23 s |
| continuation, 16 starts | 7/23/28/30 | 29/29/30/30 | 108 / 448 s |
| cold LM 16 + continuation 16 | 8/24/29/30 | **30/30/30/30** | 116 / 465 s |

For AIDS:
- Cold LM Stone-only already reaches raw 30/30 at every cutoff, and bS
  26/30, at a median of 0.3 s per seed.
- No configuration exceeds bS 27/30.
- The remaining bS failures are the u-row and projection cases above.

### Interpretation (per the prespecified rules)

- **AIDS.** Continuation gains nothing over cold LM:
  - both are at raw 30/30 at every cutoff;
  - bS is 26 vs 27/30, with the one-seed difference a u-row margin;
  - continuation costs about 5.6× more.

  Both direct methods are a **genuine solver improvement** over
  Anderson: selected and oracle both improve, bS 19 → 26/27 and raw
  22 → 30. Cold LM even from the Stone start alone is the efficient
  choice here.
- **QUAIDS.** Both selected and oracle results improve materially over
  Anderson: bS 10 → 23, raw 13 → 29, 14 seeds improved and 1 regressed.
  That is a **genuine solver improvement**. It is shared by cold LM and
  continuation at 16 starts, which tie (23/23 bS, 29/29 raw, one seed
  each way).

  Continuation's distinctive advantage is **per-start robustness**:
  - Stone-only bS 21 vs 17 and raw 27 vs 23;
  - 476/480 vs 379/480 genuinely converged;
  - a single basin on 28 of 30 seeds.

  The cost is a median of 98 s per seed for 16 starts.
- **Solution selection is not the bottleneck.** It matches the oracle
  class on 178 of 180 seed/model/method cases.
- **Remaining problems:**
  1. The `bS` metric for direct endpoints: the projection convention and
     the u row.
  2. Rare path-determined wrong basins for continuation (seed 5) and
     reachability failures for cold LM (seed 17). These are
     complementary.
- **Model-specific hybrid.** Explicitly evaluated, and the evidence
  favours it:
  - AIDS: cold LM.
  - QUAIDS: continuation, or a small pool of cold LM plus continuation.
    The post-hoc 2-solve Stone pool already matches cold LM's 16-start
    result.

  One method for both models is not required. Cold LM with 16 starts is
  a reasonable single default, but it misses QUAIDS seed 17.
- **No population-wide claim.** All of the above is 30 seeds. The
  continuation-vs-cold-LM differences are 1-4 seeds and not
  distinguishable from noise. The Anderson-vs-direct differences (QUAIDS
  14 improved vs 1 regressed) are large.

### Recommendation and next entry point

Do **not** launch 200 seeds as the very next step. Two decisions should
be made first, because each changes what a 200-seed run would measure:

1. **Fix the correctness metric for direct endpoints**, before any
   scale-up. Choose one of:
   - (a) impose symmetry inside the direct objective (optimize over the
     symmetric parametrization), so no post-hoc projection is needed and
     `bS` comes from homogeneity and adding-up only;
   - (b) prespecify the unprojected recovery, or `bsErrExclU`, as the
     direct-method primary metric and keep projected `bS` as
     benchmark-compatible context;
   - (c) derive the direct estimator's own sandwich covariance for the
     projection.

   Option (a) is the cleanest and matches what a shipped estimator
   would need.
2. **Prespecify the hybrid to scale up.** For example: AIDS cold LM
   Stone + k starts; QUAIDS pooled cold LM + continuation from a small
   start set. Include a seed-5-type check (a continuation-wide wrong
   basin) as a reported failure mode.

After that, a 200-seed run of the chosen hybrid is justified. The
Anderson-vs-direct gap is large and consistent on these 30 seeds, and
Anderson's 30-seed rates match its 200-seed rates. More work on
solution selection is not indicated.

To resume or extend: the harness is resumable per cell. For example,
set `pilotSeeds = seqa(31,1,170)` in a copy of `pilot_run_m*.e` with the
same `pilotOutDir`. `pilot_analyze.py --report` recomputes everything
from disk.

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
