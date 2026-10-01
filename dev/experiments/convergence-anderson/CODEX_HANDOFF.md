# Handoff: evaluate Anderson-acceleration / multi-start research for AIDS/QUAIDS convergence failure

> Historical incoming handoff. The requested audit and recommended full
> benchmark have now been completed. For the current state and next task,
> read `CLAUDE_HANDOFF.md` and `SUMMARY.md` in this directory.

## Context

This is a GAUSS library (`gauss-aids-model`, package `quaids`) estimating
Almost Ideal Demand System models. Its iterated AIDS/QUAIDS estimator
(`quaidsFit()` in `src/quaids.src`, `aCtl.maxiter > 1`) has a documented,
real convergence-reliability problem, quantified by a committed
diagnostic (`tests/quaids_convergence_sweep.e`, 200 synthetic seeds,
`tobs=3000`): iterated AIDS never-converges 39% of the time and
converges to a self-consistent but WRONG answer another 19% (58%
combined failure); QUAIDS is worse (54.5% / 21.5%, 76% combined). This is
already documented for library users in `README.md`'s "Model & Feature
Support Tiers" and `docs/FEATURE_SUPPORT_MATRIX.md` — it is not news, it
is the starting point.

A prior session (me, Claude) spent real compute investigating whether
techniques adjacent to what stabilizes optimization in deep learning
(the user's own framing — vanishing/exploding gradients) have genuine
numerical-analysis analogs that could help here, since the actual
mechanism (`quaidsFit()`'s iteration is successive substitution /
fixed-point iteration, NOT gradient descent) is governed by the same
underlying fact that governs RNN gradient stability: the spectral radius
of a Jacobian near the fixed point. That investigation's full write-up,
code, and raw numeric output live in this same directory
(`dev/experiments/convergence-anderson/`) — **read `SUMMARY.md` in this
directory first, in full, before doing anything else.** It has exact
numbers, exact file references, what was validated vs. assumed, and what
was NOT tried. Do not re-derive what's already in there.

## What I actually want from you

**Primarily: a critical evaluation, not blind continuation.** Specifically:

1. **Sanity-check the methodology.** Read `quaidsfit_anderson_prototype.src`
   (the Anderson-accelerated copy of `quaidsFit()`) and the experiment
   scripts. Is the Anderson acceleration implementation actually correct?
   (It was validated on a toy linear fixed-point problem first —
   `anderson_toy_test.e` — and the prototype was verified byte-identical
   to the real `quaidsFit()` at `mDepth=0` — `anderson_faithfulness_check.e`
   — but a second set of eyes on the linear algebra, especially the
   `vec()`/`reshape()` round-trip and the ridge-regularized least-squares
   mixing-coefficient solve, would be genuinely useful.) Are there
   confounds in the sweep/pilot/curvature experiments beyond the two I
   already found and documented (the non-curvature-consistent DGP fixture
   in curvature attempt #1; the `K`/`k` case-collision bug)?
2. **Judge whether the reported findings are actually sound**, particularly:
   - Anderson nearly eliminates never-converged but shifts most of that
     mass into converged-but-wrong rather than converged-correctly. Is
     that conclusion well-supported by the 200-seed data in
     `anderson_sweep_output.txt`, or is there a more charitable reading
     I'm missing?
   - The multi-start pilot's basin-collapse finding (most perturbed
     starts land on the identical fixed point) — is 0.3 relative Gaussian
     perturbation just too weak, or is this telling us something more
     fundamental about the shape of this optimization landscape?
   - Whether "the model's own fit criterion doesn't reliably discriminate
     correct from incorrect fixed points" (finding #2 in SUMMARY.md)
     genuinely points toward weak identification (multiple real,
     statistically-hard-to-distinguish optima), vs. some fixable artifact
     of `homogCrit` specifically (e.g. try comparing raw log-likelihood
     or a properly-normalized GMM criterion instead — `homogCrit` is a
     log-det proxy, not necessarily the best available signal).
3. **Recommend a next step, and say why**, from (at minimum) these
   candidates, all described in more depth in SUMMARY.md's "What was NOT
   tried" section:
   - More aggressive/structured multi-start perturbation.
   - Isolating whether Anderson itself suppresses basin diversity (run
     multi-start WITHOUT Anderson, same perturbation scale, compare
     diversity).
   - The cheap, not-yet-run "genuine non-identification" test: does a
     wrong fixed point ever have a strictly BETTER `homogCrit` than the
     TRUE synthetic DGP parameters' own `homogCrit` evaluated on the same
     data? (This alone would be strong evidence for real non-identification
     regardless of search strategy — cheap to check, do this early if you
     pursue anything further.)
   - A full reformulation as direct nonlinear optimization (Gauss-Newton/
     Levenberg-Marquardt via `optmt`/`cmlmt`) instead of successive
     substitution — bigger effort, and per SUMMARY.md, possibly futile if
     the real problem is non-identification rather than search strategy,
     so sequence this AFTER the identification check above, not before.
4. If you do continue the investigation, **keep working in this
   directory** (`dev/experiments/convergence-anderson/`), extend
   `SUMMARY.md` with what you find (don't just leave results in scrollback),
   and do not touch any already-shipped file under `src/` — this whole
   line of work is deliberately kept out of the installed package/CI
   surface until (if ever) a specific fix is validated and someone
   explicitly decides to productionize it. See `CLAUDE.md`'s "Don't touch
   already-shipped, tested estimation core without a strong reason"
   convention.

## Hard constraints

- **Never commit or push.** Confirm with the user first, always. This
  applies even to files inside this experiments directory.
- **Read the repo-root `CLAUDE.md` before writing or running any GAUSS
  code.** It documents ~15 real, previously-confirmed GAUSS-26 language
  gotchas — several were independently rediscovered the hard way during
  this exact investigation (case-insensitive identifiers colliding,
  `gamma`/`msym` being reserved words, `reshape()` filling row-major not
  matching `vec()`'s column-major stacking, `$|`-built character matrices
  being a different GAUSS type from `$+`-built ones, a `$+`-concatenation-
  with-a-multi-element-`ftocv()`-result broadcast bug). Do not rediscover
  these again; they're documented specifically so you don't have to.
- **Validate numerically before trusting, always.** Every step of this
  prior investigation that mattered was checked against an isolated,
  hand-verifiable case before being trusted at scale (the toy fixed-point
  problem before the real AIDS iteration; the faithfulness check before
  the sweep; the Slutzky-eigenvalue replica validated against the real
  printer's own output before the curvature checks). Keep that discipline.
- **Run anything sweep-scale in the background with a generous timeout,
  and smoke-test on one seed / a few iterations first.** This
  investigation hit two genuine multi-CPU-hour hangs from careless first
  attempts at scale (one real numerical issue in unrelated TVP-AIDS work,
  one a trivial variable-naming bug that would have been caught instantly
  by a 1-iteration smoke test). GAUSS environment:
  `C:\gauss26\tgauss.exe -b -x <file>.e`, run from any directory since the
  experiment files use absolute `#include` paths.
- This is genuinely open-ended, partially-negative research, not a
  bug fix. If your honest conclusion is "this direction is a dead end,
  here's why, here's what I'd try instead" — that is a completely valid
  and useful answer. Do not manufacture a positive result.
