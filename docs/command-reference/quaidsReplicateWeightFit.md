# quaidsReplicateWeightFit

## Purpose

Replicate-weight (jackknife/BRR-style) standard errors for
[quaidsFit](quaidsFit.md), built from a caller-supplied set of
pre-computed replicate weight columns and scale factor(s) -- reuses
Milestone 26's `quaidsFit()` `weight` argument for both the full-sample
point estimate and every replicate refit, with no new estimation logic.
Silent, no printing -- see
[printQuaidsReplicateWeight](printQuaidsReplicateWeight.md).

## Format

```gauss
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactor);
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactor, weight=surveyWeight, method="JK1");   // Milestone 28
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr` - same sample and shapes as
  [quaidsFit](quaidsFit.md).
- `aCtl` (*`quaidsControl` structure*) - passed through unchanged to
  every `quaidsFit()` call (base and every replicate).
- `replicateWeights` (*Tx R matrix, required*) - each column is one
  complete, caller-supplied alternate Tx1 weight vector for the SAME T
  observations -- e.g. one delete-one-PSU jackknife replicate, or one
  BRR half-sample replicate. This proc does not construct, infer, or
  validate design-correctness of these columns; it only refits
  `quaidsFit()` under each one exactly as supplied. Required with no
  default (not keyword-callable), matching this project's "never
  silently guess an inference-affecting parameter" precedent.
- `scaleFactor` (*positive scalar, or positive Rx1 vector, required*) -
  the replication variance scale factor(s) prescribed by the survey
  design (e.g. JK1's `(R-1)/R`, constant across all replicates; BRR's
  `1/R`). **Required -- no default**, matching
  [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
  precedent of never silently guessing an inference-affecting parameter.
  A scalar is broadcast to every replicate; an Rx1 vector supplies one
  factor per replicate, for designs where replicates are not uniformly
  scaled. Declared before `weight`/`method` below, since GAUSS requires
  every required (non-defaulted) parameter to precede any keyword-
  defaulted one -- Milestone 28 moved `replicateWeights`/`scaleFactor`
  earlier in the signature than they sat in the original (Milestone 27)
  release for exactly this reason.
- `weight` (*OPTIONAL keyword argument, default `0`*) - base/full-sample
  sampling weight, passed straight through to `quaidsFit()`'s own
  optional `weight` argument. Omit, or pass scalar `0`, to hit that
  proc's own unweighted path exactly as if the argument were omitted
  entirely. Scalar `0` is the only scalar sentinel; scalar nonzero weights
  are rejected by the underlying fit.
- `method` (*OPTIONAL keyword argument, default `"custom"`*) - purely
  descriptive (e.g. `"JK1"`, `"BRR"`). Echoed in `rOut.method` and the
  printed report; has no effect on computation, which is why (unlike
  `replicateWeights`/`scaleFactor`) a default is appropriate here.

## Returns

`rOut` is a `quaidsReplicateOut` structure:

- `b` - the base (full-sample) point estimate -- exactly
  `qOutFull.bestB`, **already in the same basis**
  [quaidsSharesFit](quaidsSharesFit.md)/[quaidsElasFit](quaidsElasFit.md)/
  [quaidsWelfareFit](quaidsWelfareFit.md) expect. Unlike
  [quaidsRobustFit](quaidsRobustFit.md)'s own reduced regressor-aligned
  basis, no expansion helper is needed before passing `rOut.b`/`rOut.v`
  to those procs.
- `se`, `v` - the replicate-weight standard errors (same shape as `b`)
  and the full covariance of `vec(b)`.
- `scaleFactor` - the per-replicate scale factor(s) actually used
  (broadcast to `Rx1` if supplied as a scalar).
- `bReplicate` - `nCompleted x (rows(b)*cols(b))` matrix of raw
  `vec(b_r)` draws from each completed replicate.
- `nRequested`, `nCompleted`, `nFailed` - replicate run bookkeeping. See
  Remarks for how a failed replicate is handled (no retry -- replicate
  weights are fixed, not random draws).
- `weighted` - echoes `qOutFull.weighted` (`1` if the base `weight` was
  a genuine, non-scalar-zero vector).
- `method` - echoed back unchanged.

## Remarks

**This is a genuinely simpler case than the bootstrap procs it otherwise
resembles.** Replicate weights are FIXED, caller-supplied columns, not
random draws -- there is no resampling loop, no `seed`, and no
"attempt up to `5*B` times" retry logic (a failed replicate cannot be
redrawn, since it is not random). If a replicate fails to converge, it
is simply dropped from the sum (`nFailed` bookkeeping) -- a documented
simplification, since the formal JK1/BRR literature does not, in
general, define how to handle a missing replicate without an adjustment
this library does not implement.

**No design is auto-detected or assumed.** This proc does not implement
jackknife, BRR, or any other specific replicate-weight method -- it
implements the shared linear form underlying all of them,
`V = sum_r c_r * vec(b_r - b_full) * vec(b_r - b_full)'`, and takes
`replicateWeights`/`scaleFactor` as required inputs. Supplying the
correct replicate columns and scale factor for your survey's actual
design is the caller's responsibility.

**A real, non-trappable crash mode was found and guarded against while
building this** (not assumed): a replicate weight column that leaves too
few *effectively*-weighted observations relative to the number of
estimated design columns (Kish's effective sample size,
`(sum w)^2 / sum(w^2)` -- the same quantity `quaidsFit()`'s own `effN`
field reports) can drive `quaidsFit()`'s internal iteration into a
rank-deficient state and fail with a plain GAUSS indexing error
(`error G0058: Index out of range`), confirmed by direct reproduction
with a replicate weight concentrated on 5 rows. `trap 1` does **not**
catch this -- the same class of non-trappable failure already documented
for `eighv()` ([quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md))
and `glm()` ([quaidsZeroFit](quaidsZeroFit.md)). Since the cause is
cheaply checkable in advance, this proc computes each replicate's own
effective sample size before ever calling `quaidsFit()` and skips
(counts as failed) any replicate falling below `2x` the number of design
columns -- a defensive, documented heuristic margin, not a formal
statistical requirement.

**`weight=0` is exactly equivalent to omitting `quaidsFit()`'s own
`weight` argument** -- both are the same keyword default, so passing
scalar `0` explicitly or leaving `weight` off entirely hit the identical
unweighted code path in `quaidsFit()`, confirmed directly, not assumed.

**`bReplicate`'s raw draws are exposed** for a caller who wants to build
percentile confidence intervals or inspect the replicate distribution
directly, mirroring
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
`bBoot` field -- no separate CI-computation proc is added by this
milestone (unlike
[quaidsCurvatureBootstrapCI](quaidsCurvatureBootstrapCI.md)); a caller
can reuse the same `quantile()`-based approach directly against
`rOut.bReplicate` if wanted.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

// replicateWeights: TxR matrix of survey-provided jackknife replicate
// weights; scaleFactorJK1 = (R-1)/R for a standard JK1 design.
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl,
    replicateWeights, scaleFactorJK1, weight=surveyWeight, method="JK1");

call printQuaidsReplicateWeight(rOut);
print "completed:" rOut.nCompleted "failed:" rOut.nFailed;

// rOut.b/rOut.v feed directly into post-estimation -- no expansion step:
sharesOut = quaidsSharesFit(rOut.b, rOut.v, intcptPt, pricesPt, totexpPt, aCtl);
```

## Source

`quaidsreplicate.src`

## See Also

[printQuaidsReplicateWeight](printQuaidsReplicateWeight.md),
[quaidsFit](quaidsFit.md), [quaidsRobustFit](quaidsRobustFit.md) (the
closed-form weighted/clustered sandwich alternative),
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) (the
nonparametric-resampling alternative this proc's structure mirrors),
[quaidsSharesFit](quaidsSharesFit.md), [quaidsElasFit](quaidsElasFit.md),
[quaidsWelfareFit](quaidsWelfareFit.md)
