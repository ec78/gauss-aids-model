# quaidsRobustBootstrapFit

## Purpose

Cluster-aware nonparametric bootstrap alternative to
[quaidsRobustFit](quaidsRobustFit.md)'s closed-form sandwich --
resamples whole clusters with replacement when `clusterId` is supplied,
or plain i.i.d. rows when it isn't (the literal `G=nobs` special case,
sharing the same resampling code path). Silent, no printing -- see
[printQuaidsRobustBootstrap](printQuaidsRobustBootstrap.md).

## Format

```gauss
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, B);
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, B, clusterId=householdId, seed=42, weight=myWeight);   // Milestone 26/28
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr` (*same sample and shapes as
  [quaidsFit](quaidsFit.md)*).
- `aCtl` (*`quaidsControl` structure*) - passed through unchanged to
  every [quaidsFit](quaidsFit.md) call.
- `B` (*positive integer scalar, required*) - number of bootstrap
  replications to complete. No default, matching
  [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
  precedent of never silently guessing an inference-affecting parameter.
  Declared before the keyword-defaulted parameters below, since GAUSS
  requires every required (non-defaulted) parameter to precede any
  keyword-defaulted one -- Milestone 28 moved `B` earlier in the
  signature than it sat in the original (Milestone 20) release for
  exactly this reason.
- `clusterId` (*OPTIONAL keyword argument, default `0`*) - `0` for a
  plain i.i.d. row bootstrap, or a `Tx1` vector of cluster group labels
  for a cluster (block) bootstrap.
- `seed` (*OPTIONAL keyword argument, default `0`*) - if `seed > 0`,
  `rndseed` is set before drawing any resamples (reproducible runs);
  `seed = 0` leaves GAUSS's current random state unchanged.
- `weight` (*OPTIONAL keyword argument, default `0`*) - Milestone 26: a
  sampling weight, same semantics as [quaidsFit](quaidsFit.md)'s own
  `weight`. Omit, or pass scalar `0`, for the unweighted bootstrap. Each
  resample carries its own rows' `weight[idx]` subvector into that
  replication's [quaidsFit](quaidsFit.md) call; the base (unresampled)
  fit and sandwich use the full `weight` the same way
  [quaidsRobustFit](quaidsRobustFit.md) would.

## Returns

`rbOut` is a `quaidsRobustBootOut` structure:

- `b` - the base (unresampled) reduced-form point estimate, same shape
  and basis as [quaidsRobustFit](quaidsRobustFit.md)'s own `b`.
- `seRobust` - the base closed-form robust/cluster SE
  ([quaidsRobustFit](quaidsRobustFit.md)'s `se`), for comparison.
- `seBoot`, `bBoot` - the bootstrap standard error and raw draws.
- `nRequested`, `nCompleted`, `nFailed`, `nAttempts` - bootstrap run
  bookkeeping. Up to `5*B` total resamples are attempted before giving up
  on reaching `B` completed replications.
- `nClusters` - number of unique clusters resampled from; equals `nobs`
  when `clusterId = 0`.

## Remarks

**Resampling**: refits only [quaidsFit](quaidsFit.md) per replication
(not [quaidsRobustFit](quaidsRobustFit.md) too), mirroring
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
"refit the estimator, not its SE stage, each rep" pattern. When
`clusterId` is supplied, each replication draws `G` cluster labels with
replacement and gathers every row belonging to each drawn cluster (the
resampled row count varies by replication, standard for a cluster
bootstrap); when it isn't, this reduces to the identical plain i.i.d. row
bootstrap [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)
already uses.

**Typically much less conservative than the closed-form sandwich**:
because this bootstrap resamples and refits the *actual*, efficient
FGLS/IV estimator, its `seBoot` is usually much closer to
`qOut.homogSE`/`symcSE` than [quaidsRobustFit](quaidsRobustFit.md)'s own
`seRobust` is -- see that proc's Remarks for the empirically-confirmed
gap between the two. Prefer `seBoot` when that gap matters for your use
case; `seRobust` is exposed alongside it for direct comparison, not
replaced.

**Post-estimation propagation**: `rbOut.bBoot` and `rbOut.seBoot` are in
the same reduced basis as [quaidsRobustFit](quaidsRobustFit.md)'s output.
Use [quaidsRobustBootstrapCovariance](quaidsRobustBootstrapCovariance.md)
to expand the empirical bootstrap covariance into `qOut.bestB`'s full
basis before passing it to [quaidsSharesFit](quaidsSharesFit.md),
[quaidsElasFit](quaidsElasFit.md), or [quaidsWelfareFit](quaidsWelfareFit.md).

**A real bug found and fixed while building this**: an early version
tracked the bootstrap point estimate in `qOut.bestB`'s full (adding-up-
recovered, `n`-column) shape while `quaidsRobustFit()`'s own `se` is in
the `n1`-column reduced form -- a genuine shape mismatch, caught by
running [printQuaidsRobustBootstrap](printQuaidsRobustBootstrap.md)
against real data (`error G0058: Index out of range`), not by re-reading
the code. Fixed by sharing a single reduction helper
(`_quaidsRobustReduceB()`) between this proc and
[quaidsRobustFit](quaidsRobustFit.md), so both report the identically-
shaped reduced form.

**The reshape/cell-position bug class Milestone 18 found twice already**
(row-major `reshape()` vs. column-major `vec()`) was guarded against
from this proc's first version, with its own regression test in
`tests/quaids_robust_bootstrap_test.e`, not found the hard way a third
time.

## Examples

```gauss
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, clusterId=householdId, seed=42);
call printQuaidsRobustBootstrap(rbOut);

{ bB, vB } = quaidsRobustBootstrapCovariance(qOut, rbOut, aCtl);
sharesB = quaidsSharesFit(bB, vB, intcptPt, pricesPt, totexpPt, aCtl);
```

## Source

`quaidsrobust.src`

## See Also

[printQuaidsRobustBootstrap](printQuaidsRobustBootstrap.md),
[quaidsRobustBootstrapCovariance](quaidsRobustBootstrapCovariance.md),
[quaidsRobustFit](quaidsRobustFit.md), [quaidsFit](quaidsFit.md),
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) (the
sibling bootstrap this proc's design mirrors)
