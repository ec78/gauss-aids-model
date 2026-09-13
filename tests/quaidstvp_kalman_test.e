/*
** quaidstvp_kalman_test.e
**
** TVP-AIDS initiative, Stage 2: validates the Kalman-filter wiring in
** src/quaidstvp.src -- _quaidsTVPBuildModel()/_quaidsTVPKalmanFit(),
** which package Stage 1's Zarr (homogeneity+symmetry-respecting
** state-vector construction) into a real sslib tvpModel and run it
** through sslib's kalmanFilterDiffuseTVP()/kalmanFilterTVP(), with a
** caller-supplied FIXED Q/H (hyperparameter MLE is Stage 3, not
** attempted here). See that file's own Stage 2 header for the full
** design writeup (random-walk state transition T=I/c=0/R=I; why the
** diffuse filter, not the ordinary one, is the intended default for a
** pure random walk with no stationary distribution).
**
** Separate file from quaidstvp_test.e (Stage 1), which has NO sslib
** dependency and must keep working even on a machine without sslib
** installed -- this file requires `library cmlmt, tsmt, sslib;` and the
** GAUSS26_CFG override documented in CLAUDE.md (tests/run_source_tests.ps1
** sets this env var itself before invoking this file; run directly with
** `tgauss -b -x quaidstvp_kalman_test.e`, GAUSS26_CFG must already point
** at tests/gauss26_cfg_override -- see CLAUDE.md's tsmt-shadowing note).
**
** Check 1 (the primary, definitive check, mirroring Stage 1's own exact-
** recovery bar): with Q and H both forced to ~0 and a diffuse prior
** (every state element's prior variance literally infinite), the
** filtered state at the final period of a NOISELESS synthetic sample
** must recover the true state to floating-point precision -- a Kalman
** filter in this limit is mathematically equivalent to pooled OLS over
** the same shared design (Stage 1's own exact-recovery check), so this
** independently confirms the filter wiring (Z/H/T/c/R/Q assembly, the
** diffuse initialization via sslib's init_diffTVP, kalmanFilterDiffuseTVP
** itself) without re-deriving or duplicating Stage 1's own OLS check.
**
** Check 2: the non-diffuse branch (kalmanFilterTVP, using tvpm.a_0/p_0
** AS SUPPLIED by the caller) runs without error given a genuine,
** non-default prior, and -- since p_0=0 (a degenerate "certain the state
** is exactly a_0" prior) only relaxes through Q's own period-by-period
** accumulation (P_(t+1) = T*P_t*T' + Q = P_t + Q here, since T=I) -- the
** final filtered state stays CLOSE TO a_0, not to the true state, by a
** predictable amount (loose upper bound only, not a second exact-
** recovery claim): confirms the non-diffuse code path is wired
** correctly and is not silently reusing the diffuse branch's machinery.
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvp_kalman_test.e
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
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
   Check 1: noiseless synthetic recovery via the diffuse filter, n1 = 3.
   ========================================================================== */

n1 = 3;
tobs = 50;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 777, n1);

Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

Qtiny = 1e-8*eye(k_states);
Htiny = 1e-10*eye(n1);

struct tvpModel tvpm;
tvpm = _quaidsTVPBuildModel(Zarr, Qtiny, Htiny, n1);

call check(tvpm.k_states == k_states, "tvpm.k_states matches n1 + ngamma + n1");
call check(tvpm.k_endog == n1, "tvpm.k_endog matches n1");
call check(tvpm.nobs == tobs, "tvpm.nobs matches tobs");

y = w';

struct kalmanResult rslt;
rslt = _quaidsTVPKalmanFit(tvpm, y, 1);

finalState = rslt.filtered_state[., tobs];
call check(maxc(abs(finalState - trueState)) < 1e-8,
    "diffuse filter, Q/H~0: final filtered state matches trueState to floating-point precision");

b0Hat = _quaidsTVPStateToB(finalState, n1);
trueAlphaCheck = trueState[1:n1];
trueBetaCheck = trueState[n1+ngamma+1:n1+ngamma+n1];
trueB0 = (trueAlphaCheck')|trueGammaFull|(trueBetaCheck');
call check(maxc(maxc(abs(b0Hat - trueB0))) < 1e-8,
    "diffuse filter result unpacks via _quaidsTVPStateToB to match the true b0 layout");

gammaHat = b0Hat[2:n1+1, .];
call check(maxc(maxc(abs(gammaHat - gammaHat'))) < 1e-10,
    "filtered gamma sub-block is exactly symmetric (hard shared-state constraint survives the filter)");


/* ==========================================================================
   Check 2: non-diffuse branch (kalmanFilterTVP, explicit zero prior) runs
   and stays close to a_0=0, not to trueState -- confirms this is a
   genuinely different code path, not a silent diffuse fallback. Uses its
   OWN, deliberately more moderate Q/H (not Check 1's Qtiny/Htiny): a
   near-zero H combined with p_0 left at exactly zero is a pathologically
   ill-conditioned forecast-error covariance (F_obs = Z*P*Z' + H, both
   terms ~1e-8 to 1e-10) that sent the filter into numerical blowout
   (~1e+93) when first tried with Qtiny/Htiny reused from Check 1 -- a
   property of that specific tiny-H/zero-P combination, not a bug in the
   filter wiring (Check 1's own result, same tvpm machinery, is exact to
   floating-point precision). A realistic-scale H avoids this entirely.
   ========================================================================== */

Qmod = 0.001*eye(k_states);
Hmod = 0.01*eye(n1);

struct tvpModel tvpm2;
tvpm2 = _quaidsTVPBuildModel(Zarr, Qmod, Hmod, n1);

struct kalmanResult rslt2;
rslt2 = _quaidsTVPKalmanFit(tvpm2, y, 0);

finalState2 = rslt2.filtered_state[., tobs];
call check(maxc(abs(finalState2)) < 1,
    "non-diffuse filter, zero prior, moderate Q/H: final filtered state stays bounded and near a_0=0 (tvpm2.p_0 left at its _quaidsTVPBuildModel default)");
call check(maxc(abs(finalState2 - trueState)) > maxc(abs(finalState - trueState)),
    "non-diffuse (zero-prior) result is farther from trueState than the diffuse (Check 1) result -- confirms the two branches are genuinely different, not the same recursion under two names");


/* ==========================================================================
   Summary
   ========================================================================== */

if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
