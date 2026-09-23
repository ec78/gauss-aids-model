/*
** 14_tvp_aids_estimation.e
**
** quaidsTVPFit() fits a genuine time-varying-parameter (TVP) AIDS model:
** coefficients evolve as a random walk over time via a Kalman filter,
** rather than being constant across the sample -- distinct from
** quaidsTrendFit()'s own cheap, one-shot linear-drift screening
** diagnostic (see docs/COMMAND_REFERENCE.md's "Time-Varying-Parameter
** Screening" section). Requires the optmt-style "optional adapter"
** pattern: the sslib (gauss-state-space) package.
**
** quaidstvpfit.src (and its five sibling TVP-AIDS files) is NOT loaded by
** `library quaids;` alone -- opt-in only, since it has a hard
** compile-time dependency on sslib's struct types and core estimation
** needs no external package at all. See docs/public-api.json's
** "optional_modules" entry ("tvp").
**
** Environment note: `library ... sslib;` also pulls in `tsmt`, and this
** machine's own gauss.cfg has a documented package-shadowing issue that
** breaks plain tsmt proc resolution (see CLAUDE.md's "tsmt package
** shadowing" note) -- if this example fails with an "Undefined symbol"
** error inside a tsmt/sslib proc, set the GAUSS26_CFG environment
** variable to point at a gauss.cfg whose extra_lib_path lists
** tsmt\lib explicitly (tests/gauss26_cfg_override/gauss.cfg is a
** ready-made one this repo already ships for its own test suite).
**
** Uses only the first 3 of quaidsExampleData()'s 5 goods (n1=2), NOT the
** full 5-good dataset -- confirmed directly (not assumed) that the full
** 5-good case (n1=4, k_states=18) makes ssFitTVP()'s CMLMT optimization
** pathologically slow/non-terminating on this dataset (observed directly:
** still running after 100s+ at tobs as low as 200, vs. ~34s to convergence
** for the n1=2 subset at tobs=500) -- consistent with every other
** TVP-AIDS test/example in this initiative also staying at n1=2 or n1=3,
** never n1=4. A real, documented scale limit of this codebase's
** diffuse-initialization/random-walk-state MLE approach, not a bug in
** this example's own setup.
**
** Run via `cd examples; tgauss -b -x 14_tvp_aids_estimation.e`.
*/

new;
library cmlmt, tsmt, sslib, quaids;
#include quaids.sdf
#include quaidstvp.src
#include quaidstvpkalman.src
#include quaidstvpmle.src
#include quaidstvpsmooth.src
#include quaidstvpelas.src
#include quaidstvpfit.src
#include example_data.src

{ wFull, intcpt, pricesFull, totexp, instr } = quaidsExampleData(500, 204);

nSub = 3;
wRaw = wFull[., 1:nSub];
w = wRaw ./ sumc(wRaw');
prices = pricesFull[., 1:nSub];
goodNames = 0$+quaidsExampleGoodNames();
goodNames = goodNames[1:nSub];

n = cols(prices);
n1 = n - 1;
k_states = 2*n1 + n1*(n1+1)/2;

/* ---------------------------------------------------------------------
** H (the observation covariance) stays FIXED, never estimated -- a real,
** documented scope decision (see quaidsTVPFit.md's own Remarks). A
** simple, principled starting choice: fit a static homogeneity-
** constrained LA-AIDS model on the same data and use its own per-
** equation residual variances. q0's own scale is likewise chosen
** relative to H's, not an arbitrary tiny constant -- confirmed directly
** that starting Q many orders of magnitude smaller than H (e.g. a flat
** 0.001) contributes to the same CMLMT slowness noted above, on top of
** the n1=4 issue.
** --------------------------------------------------------------------- */

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
print "Static LA-AIDS fit converged:" qOut.converged;

H = diagrv(eye(n1), qOut.homogSse[1:n1]/qOut.nobs);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);

/* ---------------------------------------------------------------------
** The TVP-AIDS fit itself: filtered + smoothed state paths, MLE-fitted Q.
** --------------------------------------------------------------------- */

tvpCtl = quaidsTVPControlCreate();
/* 0$+ coerces quaidsExampleGoodNames()'s $|-built return value into the
   legacy character-matrix form a `matrix`-typed struct field requires --
   see src/quaids.sdf's own quaidsTVPControl header for why (a real,
   confirmed GAUSS `error G0071 : Type mismatch` otherwise); already
   applied once above when building goodNames, reused here directly. */
tvpCtl.othnam = goodNames;

tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);

print "";
call printQuaidsTVP(tvOut);

/* ---------------------------------------------------------------------
** Elasticities at the final period's smoothed state.
** --------------------------------------------------------------------- */

finalPrices = prices[tvOut.nobs, .]';
finalTotexp = totexp[tvOut.nobs];

{ er, ep, epc } = quaidsTVPElasFit(tvOut.smoothedState[., tvOut.nobs], n1, finalPrices, finalTotexp, aCtl);

print "";
print "=== Income elasticities at the final period's smoothed state ===";
print$ goodNames~ftocv(er, 0, 4);

/* ---------------------------------------------------------------------
** Filtered-only (no smoother) is available too, e.g. for a real-time/
** causal-only use case -- skips the RTS smoother entirely.
** --------------------------------------------------------------------- */

tvpCtlFiltered = quaidsTVPControlCreate();
tvpCtlFiltered.smooth = 0;
tvOutFiltered = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlFiltered);

print "";
print "=== Filtered-only fit (tvpCtl.smooth = 0) -- bFinal comes from the";
print "    filtered, not smoothed, final-period state ===";
call printQuaidsTVP(tvOutFiltered);
