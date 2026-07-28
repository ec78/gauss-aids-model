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
struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, clusterId, B, seed);
```

## Parameters

- `w`, `intcpt`, `prices`, `totexp`, `instr` (*same sample and shapes as
  [quaidsFit](quaidsFit.md)*).
- `aCtl` (*`quaidsControl` structure*) - passed through unchanged to
  every [quaidsFit](quaidsFit.md) call.
- `clusterId` (*scalar `0`, or Tx1 vector*) - `0` for a plain i.i.d. row
  bootstrap, or a vector of cluster group labels for a cluster (block)
  bootstrap.
- `B` (*positive integer scalar*) - number of bootstrap replications to
  complete. Required -- no default, matching
  [quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
  precedent of never silently guessing an inference-affecting parameter.
- `seed` (*scalar*) - if `seed > 0`, `rndseed` is set before drawing any
  resamples (reproducible runs); `seed = 0` leaves GAUSS's current random
  state unchanged.

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
struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, householdId, 200, 42);
call printQuaidsRobustBootstrap(rbOut);
```

## Source

`quaidsrobust.src`

## See Also

[printQuaidsRobustBootstrap](printQuaidsRobustBootstrap.md),
[quaidsRobustFit](quaidsRobustFit.md), [quaidsFit](quaidsFit.md),
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) (the
sibling bootstrap this proc's design mirrors)
