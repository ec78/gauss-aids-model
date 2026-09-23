# quaidsTVPFit

## Purpose

Genuine time-varying-parameter (TVP) AIDS estimation via a Kalman filter:
coefficients evolve as a random walk over time, rather than being
constant across the sample. One call from raw data to a fitted,
optionally smoothed, state path plus a final-period coefficient snapshot.
Silent, no printing -- see [printQuaidsTVP](printQuaidsTVP.md).

**Requires the `sslib` (`gauss-state-space`) package**, installed
separately -- not a `library quaids;` dependency, so core estimation
never requires it. See
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
  [quaidsTVPControlCreate](quaidsTVPControlCreate.md).

## Returns

`tvOut` is a `quaidsTVPOut` structure:

- `n`, `n1`, `nobs`, `k_states` - dimensions.
- `wnam` - legacy character matrix of good names.
- `H`, `q0` - echoed inputs.
- `Qfit` - fitted `Q` diagonal (`k_states x 1`, natural units).
- `mleRetcode` - CMLMT convergence code; `0` means converged.
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

**No delta-method standard errors.** Propagating the state's own
filtered/smoothed covariance through the (linear) state-to-coefficients
recovery step is real, tractable future work -- not attempted yet. Use
`filteredStateCov`/`smoothedStateCov` directly if you need the raw state
covariance.

**MLE convergence speed degrades sharply with `k_states`** (and so with
`n1`, the number of goods). Confirmed directly: a 5-good (`n1=4`,
`k_states=18`) dataset made the MLE step pathologically slow (still
running well past 100 seconds at `tobs` as low as 200), while the exact
same dataset restricted to 3 goods (`n1=2`, `k_states=7`) converged in
well under a minute at `tobs=500`. `q0`'s own scale matters too -- a flat,
very small starting value (e.g. `0.001` regardless of `H`'s own
magnitude) makes this worse; scaling `q0` relative to `H` (as in the
example below) is safer. Every test and example in this initiative stays
at `n1=2` or `n1=3` for exactly this reason -- treat a larger `n1` as
unproven, not just untested for style reasons.

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
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
call printQuaidsTVP(tvOut);

// Elasticities at the final period's (smoothed) state.
{ er, ep, epc } = quaidsTVPElasFit(tvOut.smoothedState[., tvOut.nobs], n1, prices[tvOut.nobs, .]', totexp[tvOut.nobs], aCtl);
```

## Source

`quaidstvpfit.src`

## See Also

[printQuaidsTVP](printQuaidsTVP.md), [quaidsTVPControlCreate](quaidsTVPControlCreate.md),
[quaidsTVPStateToFullB](quaidsTVPStateToFullB.md), [quaidsTVPElasFit](quaidsTVPElasFit.md),
[quaidsTrendFit](quaidsTrendFit.md) (the cheap screening diagnostic to run
before committing to a full TVP-AIDS fit), [quaidsFit](quaidsFit.md)
