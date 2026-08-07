# quaidsCurvatureBootstrapFit

## Purpose

Nonparametric i.i.d. row (pairs) bootstrap standard error for a
[quaidsCurvatureFit](quaidsCurvatureFit.md) coefficient vector -- closes a
documented gap in that proc's classical delta-method SE, which is known to
be unreliable whenever the estimated Cholesky factor lands on the boundary
of the negative-semidefinite cone. Resamples rows of `(w, intcpt, prices,
totexp, instr)` with replacement and refits the whole pipeline
(`quaidsFit()` then `quaidsCurvatureFit()`) on each resample. Silent, no
printing -- see [printQuaidsCurvatureBootstrap](printQuaidsCurvatureBootstrap.md).

## Format

```gauss
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, B);
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, B, seed=42);   // Milestone 28
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr` - same sample and shapes as
  [quaidsFit](quaidsFit.md). Must satisfy `aCtl.homogenous = 1` (checked
  up front, before any resampling).
- `aCtl` (*`quaidsControl` structure*) - passed through unchanged to every
  `quaidsFit()`/`quaidsCurvatureFit()` call, including `aCtl.relax` for
  QUAIDS (see [quaidsCurvatureFit](quaidsCurvatureFit.md)'s Remarks).
- `B` (*scalar, required*) - number of bootstrap replications to complete.
  **Required -- no default** (see Remarks for why).
- `seed` (*OPTIONAL keyword argument, default `0`*) - if `seed > 0`,
  reseeds GAUSS's random number generator before drawing any resamples
  (reproducible runs); `seed = 0` (the default) leaves the current random
  state unchanged.

## Returns

`bootOut` is a `quaidsCurvBootOut` structure:

- `b`, `seDelta` - the point estimate and existing delta-method SE from the
  base (unresampled) fit, exposed unchanged for direct comparison.
- `seBoot` - the new bootstrap standard error, same shape as `b`.
- `bBoot` - `nCompleted x (rows(b)*cols(b))` matrix of raw `vec(b)` draws
  from each completed replication (so a caller can compute percentile
  confidence intervals later if wanted -- not computed by this proc).
- `nRequested`, `nCompleted`, `nFailed`, `nAttempts` - bootstrap run
  bookkeeping. A replication counts only if both `quaidsFit()` and
  `quaidsCurvatureFit()` converge on that resample and the resulting
  coefficient vector is finite; up to `5*B` total resamples are attempted
  before giving up on reaching `B` completed replications.
- `seed` - echoed back for reproducibility records.

## Remarks

**No default `B`**: a single AIDS curvature fit measures roughly 0.9
seconds; a single QUAIDS curvature fit roughly 7.3 seconds (QUAIDS also
typically needs `aCtl.relax < 1` to converge -- see
[quaidsCurvatureFit](quaidsCurvatureFit.md)). At conventional bootstrap
replication counts (200-1000), QUAIDS alone can take 24 minutes to 2 hours.
`B` must be chosen deliberately, not defaulted to a number that quietly
means something very different for QUAIDS than for AIDS.

**Resampling**: nonparametric i.i.d. row (pairs) bootstrap -- the correct
choice for this cross-sectional data (no time/panel dependence). Refits
the *whole* pipeline on each resample (not just re-perturbs the already-
fitted curvature struct), so first-stage IV sampling variability is
captured too.

**Inclusion rule**: a resample is included only if both
[quaidsFit](quaidsFit.md) and [quaidsCurvatureFit](quaidsCurvatureFit.md)
converge on it, and the resulting coefficient vector contains no NaN/Inf.
A non-converged first-stage fit feeding the second stage is not a valid
draw from the target sampling distribution -- there is no partial-
inclusion rule.

**Failure handling**: both per-replication calls are wrapped in this
codebase's established `trap`/`scalmiss` guard (the same idiom as
[quaidsFit](quaidsFit.md)'s own Milestone 12 hardening), and
`quaidsCurvatureFit()` itself gained additional pre-call finiteness checks
around its internal eigendecomposition calls as part of this milestone --
building this bootstrap surfaced a real gap where a sufficiently degenerate
resample could crash the whole run, not just fail one replication (see
`GOLD_STANDARD_TODO.md`'s Milestone 15 section).

**Silent during the loop**: no progress printing, even on a long-running
QUAIDS bootstrap -- matches this codebase's own silent-Fit-proc convention.

**Percentile confidence intervals** are computed by the separate
[quaidsCurvatureBootstrapCI](quaidsCurvatureBootstrapCI.md) (Milestone
18) directly from `bBoot`'s raw draws -- no new resampling needed.

**A real bug, found and fixed (Milestone 18)**: `seBoot`'s individual
cells were silently scrambled relative to `b` from this proc's original
Milestone 15 release until it was caught building
`quaidsCurvatureBootstrapCI()`'s own ground-truth cross-check (GAUSS's
`reshape()` fills row-major, not column-major like `vec()` -- a subtle,
easily-missed distinction). Invisible to shape/sign/finiteness checks,
since those are permutation-invariant. If you used `seBoot` from a
version before this fix, re-run with the current release.

## Examples

```gauss
library optmt;

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, seed=42);

call printQuaidsCurvatureBootstrap(bootOut);
print "completed:" bootOut.nCompleted "failed:" bootOut.nFailed;
```

## Source

`quaidscurvature.src`

## See Also

[printQuaidsCurvatureBootstrap](printQuaidsCurvatureBootstrap.md),
[quaidsCurvatureFit](quaidsCurvatureFit.md), [quaidsFit](quaidsFit.md),
[quaidsCurvatureBootstrapCI](quaidsCurvatureBootstrapCI.md)
