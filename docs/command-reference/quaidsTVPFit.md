# quaidsTVPFit

## Purpose

Genuine time-varying-parameter (TVP) AIDS estimation via a Kalman filter:
coefficients evolve as a random walk over time, rather than being
constant across the sample. One call from raw data to a fitted,
optionally smoothed, state path plus a final-period coefficient snapshot.
Silent, no printing -- see [printQuaidsTVP](printQuaidsTVP.md).

**Requires the `sslib` (`gauss-state-space`) package**, installed
separately -- not a `library quaids;` dependency, so core estimation
never requires it. The supported baseline is v2.1.2; until GPM lists
that version, obtain it from the upstream v2.1.2 tag or the package
maintainer. See
[Time-Varying-Parameter Estimation](../../README.md#time-varying-parameter-estimation-optional-sslib)
in the README for setup and `examples/14_tvp_aids_estimation.e` for a
full runnable example.

This proc adds no new estimation math of its own: it is pure
orchestration over the already-independently-validated TVP-AIDS pieces
built and tested across this initiative's Stages 1-5 (reduced
homogeneity/symmetry-respecting state-vector construction, the diffuse
Kalman filter, MLE-fitted state innovation covariance `Q` against a
caller-fixed observation covariance `H`, the RTS smoother, and
state-to-coefficients recovery).

## Format

```gauss
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
```

## Parameters

- `w` (*TxN matrix*) - budget shares.
- `prices` (*TxN matrix*) - absolute log prices.
- `totexp` (*Tx1 vector*) - log total expenditure. **Unlike
  [quaidsFit](quaidsFit.md)**, this initiative's own reduced state-space
  design has no IV/endogeneity treatment -- `totexp` is used directly,
  not first-stage-purged.
- `H` (*n1 x n1 matrix*, `n1 = N-1`) - FIXED observation covariance,
  caller-supplied, never estimated. A simple starting choice: a diagonal
  matrix built from [quaidsFit](quaidsFit.md)'s own residual variances on
  the same data (see Remarks).
- `q0` (*k_states x 1 vector*, `k_states = 2*n1 + n1*(n1+1)/2`) - STARTING
  `Q` diagonal, in natural (positive variance) units. Every element must
  be strictly positive.
- `tvpCtl` (*`quaidsTVPControl` structure*) - see
  [quaidsTVPControlCreate](quaidsTVPControlCreate.md). Its `maxIters`/
  `maxTime` fields cap the `Q` MLE's CMLMT pass (see Remarks); invalid
  values stop the call with a clear error before any estimation.

## Returns

`tvOut` is a `quaidsTVPOut` structure:

- `n`, `n1`, `nobs`, `k_states` - dimensions.
- `wnam` - legacy character matrix of good names.
- `H`, `q0` - echoed inputs.
- `Qfit` - fitted `Q` diagonal (`k_states x 1`, natural units).
- `mleRetcode` - CMLMT convergence code of the estimation pass;
  `0` means converged. `2` means the `tvpCtl.maxIters` cap was reached
  and `11` that the `tvpCtl.maxTime` cap was reached; the other outputs
  are then evaluated at the last iterate, not at an optimum. With the
  pinned `sslib` v2.1.2 there is only one CMLMT pass (see Remarks).
  Checking `mleRetcode == 0` is part of the supported scope (see
  "Support tier" below).
- `filteredState` (`k_states x nobs`), `filteredStateCov` (array,
  `nobs x k_states x k_states`) - the filtered (causal) state path.
- `smoothed` - echoes `tvpCtl.smooth`.
- `smoothedState`, `smoothedStateCov` - the RTS-smoothed state path (same
  shapes as the filtered ones), or, if `tvpCtl.smooth == 0`: `0` for
  `smoothedState`, and a placeholder `1x1x1` zero array for
  `smoothedStateCov` (an `array`-typed field can't hold a bare scalar).
- `bFinal` - `(n+2) x n` full coefficient snapshot (intercept | absolute-
  price gamma | beta) at the final period's smoothed-if-available-else-
  filtered state, via [quaidsTVPStateToFullB](quaidsTVPStateToFullB.md).
  Point values only -- no standard errors (see Remarks).

## Remarks

**Support tier: Supported (scoped)** (TVP promotion P8, 2026-10-17;
release-gated locally on `sslib` v2.1.1 and repinned to source-identical
`sslib` v2.1.2, with the adoption re-gate passed 2026-10-19). Support covers only
this scope:
- a fixed, **diagonal** `H`;
- finite `tvpCtl.maxIters`/`tvpCtl.maxTime` (the defaults are
  effectively unbounded);
- a **mandatory** `tvOut.mleRetcode == 0` check on every fit. Treat any
  other value as no result.

Outside the scope, or not claimed:
- There is **no convergence guarantee**. A non-converging fit is a
  measured outcome (see Scale below).
- There is **no standard-error claim** for coefficients, elasticities or
  `Q`.
- There is **no IV/endogeneity treatment** of total expenditure.
- A **non-diagonal `H`** is validated for the fixed-`Q` log-likelihood
  only. It is not a supported estimation scope.

**No IV treatment of `totexp`.** This initiative's Stage 1 design report
scoped the reduced state-space representation without an endogeneity
correction -- a real, documented scope limitation, not an oversight. A
future extension could add one.

**`H` stays fixed, never estimated**, matching this initiative's Stage 3
scope decision (repo-owner-approved): jointly estimating both `Q` and `H`
by unconstrained MLE hits the classic state-space variance-identification
problem (documented directly in `sslib`'s own test suite, which never
converged even on a univariate local-level model given hundreds of
iterations). Estimating `Q` alone against a fixed `H` sidesteps this
entirely. A reasonable starting `H`: fit a static
[quaidsFit](quaidsFit.md) on the same data and build a diagonal matrix
from its per-equation residual variances.

**`Q` is diagonal**, not a full covariance -- one free variance per state
element, the standard TVP-VAR/TVP-AIDS simplifying assumption, and far
cheaper/safer than a full `k_states*(k_states+1)/2`-parameter covariance
given the identification risk above.

**Diagonal and non-diagonal `H`.** Against the current `sslib` v2.1.2
pin (library source identical to v2.1.1), the
fixed-`Q` log-likelihood matches an independent statsmodels
implementation to about 1e-12 over the full sample, for two diagonal
`H` cases and one non-diagonal `H` case (see
`tests/fixtures/tvp_loglik/README.md`). `sslib` v2.1.0 and earlier
could take a spurious diffuse update with a non-diagonal `H` when the
diffuse phase ended partway through a period, so do not use a
non-diagonal `H` with an older `sslib`. The non-diagonal check is one
fixed-`Q` likelihood case, not a validation of estimated `Q`, so a
non-diagonal `H` is **outside the supported estimation scope** (see
"Support tier" below). The diagonal `H` recipe is the supported one.

**No delta-method standard errors.** Propagating the state's own
filtered/smoothed covariance through the (linear) state-to-coefficients
recovery step is real, tractable future work -- not attempted yet. Use
`filteredStateCov`/`smoothedStateCov` directly if you need the raw state
covariance.

**Bound the MLE with `tvpCtl.maxIters`/`tvpCtl.maxTime`.** Their
defaults are CMLMT's own (`1e5` iterations, `1e8` minutes), so by default
a fit is unbounded and identical to one made before 0.4.0. A
non-converging fit could previously run for close to an hour per CMLMT
pass. With finite caps it returns at the cap with `mleRetcode /= 0`.
Always check `mleRetcode`, as [printQuaidsTVP](printQuaidsTVP.md) does.
[quaidsTVPDiagnostics](quaidsTVPDiagnostics.md) packages that mandatory check
with fitted/starting `Q` validity and a cheap start-relative scale warning.

The caps apply **per CMLMT pass**. With the pinned `sslib` v2.1.2,
`ssFitTVP()` runs one CMLMT pass (the covariance of the transformed `Q`
is a delta-method conversion, not a second CMLMT run), so the caps bound
the whole MLE plus a small overrun (`maxTime` is checked once per
iteration; P6 measured at most 0.036 CMLMT minutes past a 5-minute
cap), and `mleRetcode` is that pass's code. `sslib` v2.1.0 and earlier also ran an untransformed
covariance-only pass with the same control (worst case about twice each
cap).

Measured on the deterministic, correctly specified test DGP in
`tests/quaidstvpfit_test.e` (`n1=2`, `tobs=120`; originally measured on
`sslib` v2.1.1, with the same runtime code in v2.1.2): each
fit converged (`mleRetcode = 0`) in about 1 second (0.67 s for the main
fit; about 1-2 s on v2.1.0). That is evidence
for that fixture only, not a general convergence guarantee.

**`sslib` baseline.** This repository's TVP-AIDS tests are pinned to
`sslib` v2.1.2 (`sslib.pin.json`). Its library source is byte-identical
to v2.1.1; v2.1.2 updates documentation, examples and package metadata.
v2.1.0 was the first release with the
fix that correctly ends the exact diffuse Kalman filter's diffuse phase;
v2.1.1 also fixes the spurious diffuse update with a non-diagonal `H`
and replaces `ssFitTVP()`'s second CMLMT pass with a delta-method
covariance. On releases before v2.1.0, likelihoods and `Q` estimates
for this random-walk state model are affected by spurious extra diffuse
updates.

**Scale: measured behaviour, not a guarantee** (TVP promotion P6,
2026-10-17, `sslib` v2.1.1; these conclusions apply unchanged to v2.1.2
because its library source is byte-identical; evidence in
`dev/experiments/tvp_p6_scale/RESULTS.md`). These runs used example 14's
data and setup below (`quaidsExampleData(500, 204)`, `H` from a static
homogeneous LA-AIDS fit, `q0 = 0.1*mean(diag(H))`), changing only the
number of goods. Each fit used one job under an outer hard timeout,
which no job hit.

| `n1` | `k_states` | caps (`maxIters` x `maxTime` min) | `mleRetcode` | CMLMT iterations | CMLMT minutes |
|---|---|---|---|---|---|
| 2 | 7 | 500 x 5 | 0 | 23 | -- (whole example about 11 s) |
| 3 | 12 | 500 x 5 | 11 in 3/3 runs | 437-457 | 5.02-5.04 |
| 3 | 12 | 100 x 60 | 2 in 2/2 runs | 101 | 1.19-1.23 |
| 4 | 18 | 500 x 5 | 0 in 3/3 runs | 93 | 3.52-3.62 |
| 4 | 18 | 100 x 60 | 0 in 2/2 runs | 93 | 3.70-3.89 |

- `n1=4` converged with byte-identical output in all 5 runs.
- `n1=3` did not converge. When the time cap binds, the stopping
  iteration and the returned `Q` vary from run to run, because they
  depend on machine timing. When the iteration cap binds, the result is
  deterministic. CMLMT reports "The covariance of the parameters failed
  to invert" at that endpoint, and `Q` was still moving substantially
  when the fit stopped.
- Time per CMLMT iteration grew with `k_states`: about 0.7 s at
  `n1=3` and 2.3-2.5 s at `n1=4`. Whether a fit converges did not grow
  monotonically with `n1`; on this data it depends on the subset.
- This is one dataset, one seed and one `q0` rule, so it supports no
  general scale claim. A larger `n1` is unmeasured, not just untested
  for style reasons. Whether `n1=3` converges with a larger budget was
  not measured.
- The earlier claim, measured on `sslib` before v2.1.0 with no caps,
  that `n1=4` was "pathologically slow" on this data is superseded by
  these runs.

`q0`'s own scale still matters. A flat, very small start (for example
`0.001` regardless of `H`) was observed to slow the MLE before P6.
Scaling `q0` relative to `H`, as in the example below, is safer. Always
set the caps and check `mleRetcode`.

**Minimal public surface, by design.** The individual Kalman-filter,
MLE, and smoother steps this proc orchestrates are not separately public
-- this library's convention is one clean `Fit()`-style entry point per
feature (compare [quaidsFit](quaidsFit.md),
[quaidsWorkflowFit](quaidsWorkflowFit.md)).

## Examples

```gauss
library cmlmt, tsmt, sslib, quaids;
#include quaids.sdf
#include quaidstvp.src
#include quaidstvpkalman.src
#include quaidstvpmle.src
#include quaidstvpsmooth.src
#include quaidstvpelas.src
#include quaidstvpfit.src

n1 = cols(prices) - 1;
k_states = 2*n1 + n1*(n1+1)/2;

// A simple starting H from a static fit's own residual variances,
// and a small positive starting Q.
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
qOut = quaidsFit(w, 0, prices, totexp, instr, aCtl);
H = diagrv(eye(n1), qOut.homogSse[1:n1]/qOut.nobs);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);   // scaled relative to H, not a flat tiny constant -- see Remarks

tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 500;   // cap the Q MLE (defaults are effectively unbounded) ...
tvpCtl.maxTime = 5;      // ... minutes; one CMLMT pass with sslib v2.1.2
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
call printQuaidsTVP(tvOut);   // flags mleRetcode /= 0

// Elasticities at the final period's (smoothed) state.
{ er, ep, epc } = quaidsTVPElasFit(tvOut.smoothedState[., tvOut.nobs], n1, prices[tvOut.nobs, .]', totexp[tvOut.nobs], aCtl);
```

## Source

`quaidstvpfit.src`

## See Also

[printQuaidsTVP](printQuaidsTVP.md), [quaidsTVPControlCreate](quaidsTVPControlCreate.md),
[quaidsTVPStateToFullB](quaidsTVPStateToFullB.md), [quaidsTVPElasFit](quaidsTVPElasFit.md),
[quaidsTVPDiagnostics](quaidsTVPDiagnostics.md),
[quaidsTVPMultiStartDiagnostics](quaidsTVPMultiStartDiagnostics.md),
[quaidsTrendFit](quaidsTrendFit.md) (the cheap screening diagnostic to run
before committing to a full TVP-AIDS fit), [quaidsFit](quaidsFit.md)
