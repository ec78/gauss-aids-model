/*
** quaidstvp_mle_test.e
**
** TVP-AIDS initiative, Stage 3: validates _quaidsTVPMLEFit()
** (src/quaidstvpmle.src) -- hyperparameter MLE for the diagonal state
** innovation covariance Q, via sslib's ssFitTVP(), against a
** caller-FIXED observation covariance H (see that file's own header for
** why Q-only, not joint Q/H, is this stage's deliberate scope). Separate
** file from quaidstvp_kalman_test.e (Stage 2), which only exercises
** caller-supplied fixed Q/H -- no MLE at all.
**
** Requires `library cmlmt, tsmt, sslib;` and the GAUSS26_CFG override
** documented in CLAUDE.md, same as quaidstvp_kalman_test.e -- run
** directly with `tgauss -b -x quaidstvp_mle_test.e`, GAUSS26_CFG must
** already point at tests/gauss26_cfg_override.
**
** Check 1 (synthetic recovery, the primary correctness check): data
** simulated from _quaidsTVPDynamicSyntheticDGP() (tests/quaidsfixtures.src)
** -- a genuine random-walk state with KNOWN diagonal Q and KNOWN H,
** exactly the linear-Gaussian model _quaidsTVPMLEFit() assumes, so this
** is a correctly-specified-model recovery check, not just a convergence
** smoke test. CMLMT must report convergence (retcode == 0), and the
** MEAN of the fitted Q diagonal must land within a generous relative
** tolerance of the mean of the true Q diagonal. Tolerance is on the MEAN,
** not each of the k_states=7 individual diagonal elements, deliberately:
** confirmed empirically (not guessed) that individual per-state-element
** Q variances are only loosely identified even at tobs=500 (worst single
** element off by ~54% of its own true value at the seed/scale used
** here), while the AGGREGATE magnitude across all seven is much better
** identified (~7% off at this same seed) -- a known property of
** state-space variance-component MLE (the same family of identification
** issue sslib's own test/sstvpfit.inc header documents for the Q-vs-H
** case; this file's Stage 3 scope decision to fix H, not just leave it
** free, already sidesteps THAT specific problem -- this is a separate,
** milder version of the same underlying phenomenon showing up WITHIN
** Q's own diagonal instead).
**
** Check 2 (internal consistency, exact -- not approximate): rebuilding a
** tvpModel from the FITTED Q (via Stage 2's own _quaidsTVPBuildModel())
** and refiltering directly through Stage 2's own _quaidsTVPKalmanFit()
** (diffuse) must reproduce ssFitTVP()'s own sOut.kfResults.filtered_state
** to floating-point precision -- confirms _quaidsTVPQUpdate()'s Q
** reconstruction (diagrv(eye(k_states), p) plus the nobs-page wrap) is
** exactly what a caller manually replaying the fitted Q through the
** already-independently-validated Stage 2 code path would get, not an
** approximation or a different (if similar) computation.
**
** Check 3 (loose plausibility, mirroring quaidstvp_test.e's own
** real-data-plausibility check): the final-period filtered state should
** land reasonably close to the DGP's own true final-period state (which,
** unlike Stage 1/2's noiseless fixture, has itself randomly drifted over
** 500 periods) -- a generous bound (empirically confirmed max abs diff
** ~0.41 at this seed against state magnitudes up to ~1.2), not a second
** exact-recovery claim.
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvp_mle_test.e
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;
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
   synthetic DGP with known Q/H.
   ========================================================================== */

n1 = 2;
tobs = 500;
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

trueQdiag = 0.002*ones(k_states, 1);
trueH = 0.01*eye(n1);

{ w, pricesRel, lx, trueStateArr } = _quaidsTVPDynamicSyntheticDGP(tobs, 4242, n1, trueQdiag, trueH);

Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
q0 = 0.001*ones(k_states, 1);

struct ssOut sOut;
screen off;
sOut = _quaidsTVPMLEFit(Zarr, trueH, w, q0, n1);
screen on;


/* ==========================================================================
   Check 1: convergence + synthetic recovery (mean of Q's diagonal).
   ========================================================================== */

call check(rows(sOut.final_params) == k_states, "sOut.final_params is k_states x 1");
call check(sOut.mleResults.retcode == 0, "ssFitTVP (Q-only MLE) reports CMLMT convergence (retcode == 0)");
call check(minc(sOut.final_params) > 0, "every fitted Q diagonal element is strictly positive (positive_vars transform held)");
call check(maxc(sOut.final_params) < 10*maxc(trueQdiag), "no fitted Q diagonal element has blown up (sanity bound, 10x true)");

meanFittedQ = meanc(sOut.final_params);
meanTrueQ = meanc(trueQdiag);
call check(abs(meanFittedQ - meanTrueQ)/meanTrueQ < 0.35,
    "mean fitted Q diagonal within 35% (relative) of mean true Q diagonal");


/* ==========================================================================
   Check 2: internal consistency against Stage 2's independently-validated
   _quaidsTVPBuildModel()/_quaidsTVPKalmanFit() -- exact match, not loose.
   ========================================================================== */

struct tvpModel tvpmCheck;
tvpmCheck = _quaidsTVPBuildModel(Zarr, diagrv(eye(k_states), sOut.final_params), trueH, n1);

struct kalmanResult rsltCheck;
rsltCheck = _quaidsTVPKalmanFit(tvpmCheck, w', 1);

call check(maxc(maxc(abs(rsltCheck.filtered_state - sOut.kfResults.filtered_state))) < 1e-8,
    "refiltering at the fitted Q via Stage 2's own _quaidsTVPKalmanFit() exactly reproduces ssFitTVP's own kfResults.filtered_state");


/* ==========================================================================
   Check 3: loose plausibility -- final filtered state vs. the DGP's own
   (randomly drifted) true final-period state.
   ========================================================================== */

finalStateHat = sOut.kfResults.filtered_state[., tobs];
call check(maxc(abs(finalStateHat - trueStateArr[., tobs])) < 2,
    "final filtered state stays within a generous bound of the true (drifted) final-period state");


/* ==========================================================================
   Summary
   ========================================================================== */

if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
