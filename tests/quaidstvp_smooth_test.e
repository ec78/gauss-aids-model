/*
** quaidstvp_smooth_test.e
**
** TVP-AIDS initiative, Stage 4: validates _quaidsTVPSmoothFit()
** (src/quaidstvpsmooth.src) -- the fixed-interval (Rauch-Tung-Striebel)
** smoother wrapping sslib's ssKalmanSmoothTVP(), turning Stage 2's
** filtered state path into a full-sample smoothed one. Uses Stage 2's
** own fixed-Q/H diffuse filter (_quaidsTVPBuildModel()/
** _quaidsTVPKalmanFit()), not Stage 3's MLE -- Stage 4 has no real
** proc-level dependency on Stage 3 (see quaidstvpsmooth.src's own
** header), so this test avoids that extra dependency/cost.
**
** Requires `library cmlmt, tsmt, sslib;` and the GAUSS26_CFG override
** documented in CLAUDE.md, same as quaidstvp_kalman_test.e/
** quaidstvp_mle_test.e -- run directly with
** `tgauss -b -x quaidstvp_smooth_test.e`, GAUSS26_CFG must already point
** at tests/gauss26_cfg_override.
**
** Check 1 (algorithm invariant, exact): the smoothed state/covariance at
** the FINAL period must equal the filtered state/covariance there
** exactly -- both ssKalmanSmoothTVP's own doc comment and its backward
** recursion (initialized directly from a_F[.,nobs]/p_F[.,nobs], the loop
** only running from nobs-1 down to 1) guarantee this structurally; worth
** checking directly rather than only trusting the doc comment.
**
** Check 2 (RTS tightening property, exact -- OUTSIDE the diffuse
** initialization burn-in): every diagonal element of the smoothed
** covariance must be <= the corresponding filtered covariance diagonal
** element, from period ceil(k_states/n1) onward (the period by which
** every state has received enough observations to fully de-diffuse,
** k_states/n1 rounded up since n1 new observations arrive per period) --
** the defining property of fixed-interval smoothing (using future data
** can only reduce uncertainty about a past state, never increase it).
** Checked directly at every remaining diagonal entry, not sampled.
**
** Empirically confirmed (not assumed): WITHIN the diffuse burn-in, this
** property can genuinely fail by a small amount (~0.01 absolute at this
** seed, one state element, one period) -- the ordinary RTS backward
** recursion ssKalmanSmoothTVP() runs (unchanged from ssKalmanSmooth(),
** per its own doc comment) is not the specialized diffuse-smoother
** algorithm (e.g. de Jong's) a still-diffuse filtered covariance would
** technically call for; sslib applies the same ordinary recursion
** throughout, including across the diffuse phase. This is a property of
** sslib's own implementation, not a bug introduced by this wrapper --
** confirmed by finding the exact period of the one violation (t=3, this
** DGP's k_states=7/n1=2 needs ceil(7/2)=4 periods to fully de-diffuse)
** and by Check 3 below reproducing the identical numbers via sslib's own
** independently-established ssKalmanSmooth() fed the same input.
**
** Check 3 (internal consistency, exact -- the primary correctness
** check): on a model whose TVP matrices are actually CONSTANT across
** every period (Stage 2's own random-walk T=I/c=0/R=I plus a
** caller-fixed constant Q, exactly what _quaidsTVPBuildModel() always
** builds), feeding the SAME kalmanResult into sslib's own
** already-established time-invariant smoother (ssKalmanSmooth(), via a
** ssModel built from page 1 of tvpm's constant array fields) must
** reproduce _quaidsTVPSmoothFit()'s own output to floating-point
** precision. This isolates exactly the backward RTS recursion logic
** shared between the TVP and non-TVP smoothers (both fed the identical
** filtered input), independent of any question about the filter itself
** -- mirrors ssKalmanSmoothTVP's own doc comment, which documents being
** validated this same way internally, and this codebase's own
** convention (see quaidstvp_mle_test.e's Check 2) of never trusting a
** single check.
**
** Check 4 (loose plausibility): on a genuinely time-varying (random-walk)
** synthetic DGP with a known true state path, the smoothed path's mean
** absolute error against the true path should be no worse than the
** filtered path's own mean absolute error -- using the full sample to
** estimate each period's state should not make the average estimate
** worse than using only data up to that period, on average over the
** whole sample (a loose, expected-on-average property, not a per-period
** guarantee, which Check 2's variance property already covers exactly).
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvp_smooth_test.e
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpsmooth.src;
#include quaidsfixtures.src;

nfail = 0;
ncheck = 0;

proc (0) = check(cond, label);
    local i;
    i = ncheck + 1;
    ncheck = i;
    if cond;
        print "PASS  " $+ label;
    else;
        print "FAIL  " $+ label;
        nfail = nfail + 1;
    endif;
endp;


/* ==========================================================================
   Setup: n1 = 2 (k_states = 7), a genuinely time-varying (random-walk)
   synthetic DGP, filtered via Stage 2's own fixed-Q/H diffuse filter.
   ========================================================================== */

n1 = 2;
tobs = 200;
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

trueQdiag = 0.002*ones(k_states, 1);
trueH = 0.01*eye(n1);

{ w, pricesRel, lx, trueStateArr } = _quaidsTVPDynamicSyntheticDGP(tobs, 909, n1, trueQdiag, trueH);

Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);

struct tvpModel tvpm;
tvpm = _quaidsTVPBuildModel(Zarr, diagrv(eye(k_states), trueQdiag), trueH, n1);

struct kalmanResult rslt;
rslt = _quaidsTVPKalmanFit(tvpm, w', 1);

struct kalmanResult rsltSmooth;
{ aTS, pTS } = _quaidsTVPSmoothFit(tvpm, rslt);


/* ==========================================================================
   Check 1: dimensions, then the exact final-period invariant.
   ========================================================================== */

call check(rows(aTS) == k_states and cols(aTS) == tobs, "aTS is k_states x tobs");
ordPTS = getorders(pTS);
call check(ordPTS[1] == tobs and ordPTS[2] == k_states and ordPTS[3] == k_states,
    "pTS is tobs x k_states x k_states");

call check(maxc(abs(aTS[., tobs] - rslt.filtered_state[., tobs])) < 1e-10,
    "smoothed state at the final period exactly equals the filtered state there");
call check(maxc(maxc(abs(getmatrix(pTS, tobs) - getmatrix(rslt.filtered_state_cov, tobs)))) < 1e-10,
    "smoothed covariance at the final period exactly equals the filtered covariance there");


/* ==========================================================================
   Check 2: RTS tightening -- every smoothed variance <= filtered variance,
   at every period from the diffuse burn-in onward (see this file's own
   header for why the burn-in itself is excluded).
   ========================================================================== */

burnin = ceil(k_states/n1);

maxViolation = -1e300;
t = burnin;
do while t <= tobs;
    diagSmooth = diag(getmatrix(pTS, t));
    diagFilt = diag(getmatrix(rslt.filtered_state_cov, t));
    maxViolation = maxc(maxViolation|maxc(diagSmooth - diagFilt));
    t = t + 1;
endo;

call check(maxViolation < 1e-8,
    "every smoothed-covariance diagonal element is <= the corresponding filtered one, from period ceil(k_states/n1) onward (RTS tightening property outside the diffuse burn-in)");


/* ==========================================================================
   Check 3: internal consistency, exact -- against sslib's own established
   time-invariant ssKalmanSmooth(), fed the SAME kalmanResult, on the
   constant matrices _quaidsTVPBuildModel() always builds.
   ========================================================================== */

struct ssModel ssm;
ssm.T = getmatrix(tvpm.T, 1);
ssm.c = getmatrix(tvpm.c, 1);
ssm.R = getmatrix(tvpm.R, 1);
ssm.Q = getmatrix(tvpm.Q, 1);

{ aTSref, pTSref } = ssKalmanSmooth(ssm, rslt);

call check(maxc(maxc(abs(aTS - aTSref))) < 1e-10,
    "_quaidsTVPSmoothFit's smoothed state matches sslib's own time-invariant ssKalmanSmooth() exactly, fed the same filtered result");

maxCovDiff = -1;
t = 1;
do while t <= tobs;
    maxCovDiff = maxc(maxCovDiff|maxc(maxc(abs(getmatrix(pTS, t) - getmatrix(pTSref, t)))));
    t = t + 1;
endo;
call check(maxCovDiff < 1e-10,
    "_quaidsTVPSmoothFit's smoothed covariance matches sslib's own time-invariant ssKalmanSmooth() exactly, at every period");


/* ==========================================================================
   Check 4: loose plausibility -- smoothed path no worse, on average, than
   the filtered path against the DGP's own known true state path.
   ========================================================================== */

maeFiltered = meanc(meanc(abs(rslt.filtered_state - trueStateArr)'));
maeSmoothed = meanc(meanc(abs(aTS - trueStateArr)'));

call check(maeSmoothed < maeFiltered,
    "smoothed path's mean absolute error against the true state path is no worse than the filtered path's own");


/* ==========================================================================
   Summary
   ========================================================================== */

if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
