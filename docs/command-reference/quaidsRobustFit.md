# quaidsRobustFit

## Purpose

Heteroskedasticity-robust or cluster-robust covariance for an already-
fitted [quaidsFit](quaidsFit.md) result's `bestB`, generalizing the
classic `S.*.inv(gg)` sandwich to a per-observation (robust) or
per-cluster (cluster-robust) score aggregation. Silent, no printing --
see [printQuaidsRobust](printQuaidsRobust.md).

## Format

```gauss
struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, clusterId);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md), any
  `aCtl.homogenous`.
- `w` (*TxN matrix*) - budget shares, the same sample used to fit `qOut`.
- `prices` (*TxN matrix*) - absolute log prices, same sample.
- `totexp` (*Tx1 vector*) - log total expenditure, same sample.
- `aCtl` (*`quaidsControl` structure*) - `aCtl.linear`/`aCtl.alpha0` used
  to rebuild the model-implied fitted shares; must match what was passed
  to [quaidsFit](quaidsFit.md) to produce `qOut`.
- `clusterId` (*scalar `0`, or Tx1 vector*) - `0` for heteroskedasticity-
  robust (every observation is its own cluster); a vector of cluster
  group labels for cluster-robust.

## Returns

`rOut` is a `quaidsRobustOut` structure:

- `b` - an exact (not approximate) re-expression of `qOut.bestB`'s first
  `n1` columns in a `k`-row reduced-form basis (`k = (1+nint)+n1+nendog+nu`)
  -- see Remarks for why this differs in row count from `qOut.bestB`
  itself.
- `v`, `se` - the robust/cluster-robust covariance of `vec(b)`, and
  matching standard errors.
- `nClusters` - number of unique clusters; equals `nobs` when
  `clusterId = 0`.

## Remarks

**The math**: every existing covariance in this library
(`qOut.homogV`/`symcV`, the IV first stage, `quaidsCurvatureFit`'s
`cOut.v`, `quaidsZeroFit`'s `zOut.v`) rests on a single pooled,
homoskedastic `S = sse/nobs` combined with the shared-design-matrix
`S.*.inv(gg)` Kronecker sandwich. This proc replaces that pooled `S` with
a genuine per-observation (or per-cluster) score aggregation: build the
shared regressor block `X` and residuals `U` at the converged point,
`Infl[.,(i-1)*k+1:i*k] = X.*U[.,i]` for each of the `n1` independently-
estimated equations, aggregate `Infl`'s rows by `clusterId` (or leave
ungrouped when `clusterId=0`), and sandwich with `bread =
eye(n1).*.inv(gg)` and a CR1 finite-sample correction. Robust and
cluster-robust are the same formula -- `clusterId=0` is the literal
`G=nobs` (every row its own cluster) special case, confirmed by an exact-
identity regression test (`tests/quaids_robust_test.e`) rather than
assumed.

**Simplified bread, real consequence**: this proc uses `inv(gg)` as its
bread rather than replicating [quaidsFit](quaidsFit.md)'s own intricate
nonlinear-price-index-feedback Jacobian correction -- the same class of
honestly-documented simplification [quaidsCurvatureFit](quaidsCurvatureFit.md)/
[quaidsZeroFit](quaidsZeroFit.md) already ship. **Found empirically, not
assumed**: this makes `rOut.se` dramatically more *conservative* (larger)
than `qOut.homogSE`/`qOut.symcSE` -- often by more than an order of
magnitude, especially for the IV-residual coefficient, since
`qOut`'s own SE benefits from the full cross-equation FGLS system's
efficiency gain that this proc's equation-by-equation sandwich does not
share. This is **not** a bug: an independent hand-derivation confirmed
`rOut.se` (with `clusterId=0`) lands in the same order of magnitude as a
classical `S.*.gg`-style formula built from the *same* `X`/`U` this proc
itself uses, under genuine homoskedasticity -- the large gap is entirely
against `qOut`'s different, more efficient estimator, not a defect in
this proc's own formula. [quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md)'s
own bootstrap SE, which resamples and refits the *actual* efficient
estimator, is typically much closer to `qOut`'s own SE than this proc's
closed-form sandwich is -- prefer the bootstrap when this gap matters for
your use case.

**Does not automatically propagate**: unlike a change inside
[quaidsFit](quaidsFit.md) itself, this proc's `v` does **not** flow into
`qOut.symcV` or into [quaidsElasFit](quaidsElasFit.md)/
[quaidsSharesFit](quaidsSharesFit.md)/[quaidsWelfareFit](quaidsWelfareFit.md)'s
own delta-method SEs automatically -- pass `rOut.v` into those explicitly
if you want robust/cluster-robust elasticity/share/welfare standard
errors.

**Scope**: only the `n1` independently-estimated equations are covered
(equation `n` is recovered via adding-up, never separately estimated, and
gets no independent SE here either, matching every other diagnostic in
this library).

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

// Heteroskedasticity-robust:
struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);
call printQuaidsRobust(rOut);

// Cluster-robust (householdId is a Tx1 vector of group labels):
struct quaidsRobustOut rOutCluster;
rOutCluster = quaidsRobustFit(qOut, w, prices, totexp, aCtl, householdId);
call printQuaidsRobust(rOutCluster);
```

## Source

`quaidsrobust.src`

## See Also

[printQuaidsRobust](printQuaidsRobust.md), [quaidsFit](quaidsFit.md),
[quaidsRobustBootstrapFit](quaidsRobustBootstrapFit.md) (a bootstrap
alternative that does not share this proc's conservative-bread property)
