/*
** quaidstvpfit_test.e
**
** TVP-AIDS initiative, Stage 6: validates quaidsTVPFit()/printQuaidsTVP()
** (src/quaidstvpfit.src) -- the consolidated public entry point. Adds no
** new estimation math of its own (see that file's own header), so the
** core checks here are INTERNAL CONSISTENCY (exact, not loose):
** quaidsTVPFit()'s output must exactly match the SAME already-validated
** Stage 1-5 procs (_quaidsTVPStoneIndex/_quaidsTVPBuildZ/_quaidsTVPMLEFit/
** _quaidsTVPSmoothFit/quaidsTVPStateToFullB) called directly on the
** identical inputs. Synthetic recovery of Q is quaidstvp_mle_test.e's job;
** the Qfit band in check 1 is only a single-seed sanity check.
**
** Requires `library cmlmt, tsmt, sslib;` and the GAUSS26_CFG override
** documented in CLAUDE.md, same as quaidstvp_kalman_test.e/
** quaidstvp_mle_test.e/quaidstvp_smooth_test.e -- run directly with
** `tgauss -b -x quaidstvpfit_test.e`, GAUSS26_CFG must already point at
** tests/gauss26_cfg_override.
**
** Data (TVP-AIDS promotion P4): a deterministic, correctly specified
** random-walk TVP-AIDS DGP in quaidsTVPFit()'s raw input form, from
** tests/quaidsfixtures.src's _quaidsTVPRawSyntheticDGP (n = 3, tobs =
** 120, seed 7; every share stays inside (0,1)). It replaced an ill-posed
** fixture of random shares with a fixed H near the sample share variance
** and a flat q0, on which the Q MLE never converged: the uncapped fits
** could each sit in CMLMT for close to an hour (see
** dev/TVP_PROMOTION_ROADMAP.md section 2). On this DGP each fit converges
** in a few seconds; the checks below assert that (retcode 0, a wall-clock
** budget) instead of assuming it.
**
** MLE runs: the main fit (smooth = 1, custom names), the smooth = 0 fit,
** one direct _quaidsTVPMLEFit() reference call, and two deliberately
** capped fits that stop after about one CMLMT iteration. The main fit is
** reused for checks 1-3 rather than refit.
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvpfit_test.e
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;
#include ../src/quaidstvpsmooth.src;
#include ../src/quaidstvpelas.src;
#include ../src/quaidstvpfit.src;
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
   Setup: well-posed TVP-AIDS DGP, n1 = 2 (3 goods, k_states = 7).
   ========================================================================== */

n = 3;
n1 = n - 1;
tobs = 120;
k_states = 2*n1 + n1*(n1+1)/2;

state0 = .40|.30|.05|-.02|.04|-.04|.03;
trueQ = 1e-5*ones(n1, 1)|1e-4*ones(k_states-n1, 1);
H = 2.5e-5*eye(n1);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);

{ w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, 7, state0, trueQ, H);

/* Per-fit wall-clock budget, seconds. Measured at about 1-2 s per fit on
   this DGP (sslib v2.1.0, GAUSS 26.1.4; 0.67 s for the main fit on
   v2.1.1, which drops the second CMLMT pass); 60 s leaves wide headroom while
   still catching a regression to the old multi-minute non-converging
   path. */
budget = 60;
capIters = 500;
capMinutes = 2;

call check(minc(minc(w)) > 0 and maxc(maxc(w)) < 1 and maxc(abs(sumc(w') - 1)) < 1e-12,
    "fixture: every share is strictly inside (0,1) and each row sums to 1");
call check(maxc(abs(_quaidsTVPStoneIndex(w, prices, totexp, n) - lxTrue)) < 1e-10,
    "fixture: quaidsTVPFit()'s Stone index reproduces the simulation regressor exactly");


/* ==========================================================================
   Check 0: control defaults and cap resolution.
   ========================================================================== */

tvpCtl = quaidsTVPControlCreate();
call check(tvpCtl.smooth == 1, "quaidsTVPControlCreate() default: smooth = 1");
call check(tvpCtl.maxIters == 1e5 and tvpCtl.maxTime == 1e8,
    "quaidsTVPControlCreate() default caps equal CMLMT's own defaults (maxIters 1e5, maxTime 1e8 minutes)");

struct quaidsTVPControl tvpCtlBare;
tvpCtlBare.smooth = 1;
tvpCtlBare.othnam = 0;
miBare = 0;
mtBare = 0;
{ miBare, mtBare } = _quaidsTVPResolveCaps(tvpCtlBare);
call check(miBare == 1e5 and mtBare == 1e8,
    "a control struct built without quaidsTVPControlCreate() (empty caps) resolves to the CMLMT defaults");


/* ==========================================================================
   Check 1: the main fit -- tvpCtl.smooth = 1, a caller-supplied othnam,
   finite practical caps -- matches a direct _quaidsTVPStoneIndex +
   _quaidsTVPBuildZ + _quaidsTVPMLEFit chain run at CMLMT's DEFAULT caps,
   exactly. A converged fit must not depend on a cap it never reached.
   ========================================================================== */

/* 0$+ coerces the $|-built character matrix into the legacy form a
   `matrix`-typed struct field requires -- a bare $|-built value fails
   with a real, confirmed `error G0071 : Type mismatch` otherwise (found
   via examples/14_tvp_aids_estimation.e; see quaidsTVPControl's own
   header in src/quaids.sdf). */
customNames = 0$+("Alpha" $| "Beta" $| "Gamma");
tvpCtl.othnam = customNames;
tvpCtl.maxIters = capIters;
tvpCtl.maxTime = capMinutes;

screen off;
t0 = hsec;
struct quaidsTVPOut tvOut;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
secsMain = (hsec - t0)/100;

lxDirect = _quaidsTVPStoneIndex(w, prices, totexp, n);
pricesRelDirect = prices[., 1:n1] - prices[., n];
ZarrDirect = _quaidsTVPBuildZ(pricesRelDirect, lxDirect, n1);

struct ssOut sOutDirect;
sOutDirect = _quaidsTVPMLEFit(ZarrDirect, H, w[., 1:n1], q0, n1);
screen on;

print "" $+ "main fit: " $+ ftocv(secsMain, 1, 2) $+ " s, retcode " $+ ftocv(tvOut.mleRetcode, 1, 0);

call check(tvOut.n == n and tvOut.n1 == n1 and tvOut.nobs == tobs and tvOut.k_states == k_states,
    "tvOut dimensions match n/n1/nobs/k_states");
call check(tvOut.mleRetcode == 0, "well-posed DGP: the Q MLE converges (tvOut.mleRetcode == 0) within the finite caps");
call check(secsMain < budget, "main fit finishes inside the per-fit wall-clock budget");
call check(maxc(abs(tvOut.Qfit - sOutDirect.final_params)) < 1e-10,
    "tvOut.Qfit (finite caps) exactly matches a direct _quaidsTVPMLEFit() call at CMLMT's default caps");
call check(tvOut.mleRetcode == sOutDirect.mleResults.retcode,
    "tvOut.mleRetcode exactly matches a direct _quaidsTVPMLEFit() call's mleResults.retcode");
call check(maxc(maxc(abs(tvOut.filteredState - sOutDirect.kfResults.filtered_state))) < 1e-10,
    "tvOut.filteredState exactly matches a direct _quaidsTVPMLEFit() call's kfResults.filtered_state");
call check(minc(tvOut.Qfit) > 0 and minc(tvOut.Qfit ./ trueQ) > 0.2 and maxc(tvOut.Qfit ./ trueQ) < 5,
    "single-seed sanity band: every fitted Q element is positive and within 0.2x-5x of its true value");
call check(tvOut.wnam[1] $== "Alpha" and tvOut.wnam[2] $== "Beta" and tvOut.wnam[3] $== "Gamma",
    "caller-supplied tvpCtl.othnam (a real character matrix) is honored in tvOut.wnam");


/* ==========================================================================
   Check 2: tvOut.smoothedState matches a direct _quaidsTVPSmoothFit() call
   on (sOutDirect.tvpFinal, sOutDirect.kfResults), exactly.
   ========================================================================== */

struct tvpModel tvpFinalDirect;
tvpFinalDirect = sOutDirect.tvpFinal;

aTSDirect = 0;
pTSDirect = 0;
{ aTSDirect, pTSDirect } = _quaidsTVPSmoothFit(tvpFinalDirect, sOutDirect.kfResults);

call check(tvOut.smoothed == 1, "tvOut.smoothed echoes tvpCtl.smooth (1)");
call check(maxc(maxc(abs(tvOut.smoothedState - aTSDirect))) < 1e-10,
    "tvOut.smoothedState exactly matches a direct _quaidsTVPSmoothFit() call");


/* ==========================================================================
   Check 3: tvOut.bFinal matches a direct quaidsTVPStateToFullB() call on
   the final period's smoothed state column, exactly.
   ========================================================================== */

bFinalDirect = quaidsTVPStateToFullB(aTSDirect[., tobs], n1);
call check(rows(tvOut.bFinal) == n+2 and cols(tvOut.bFinal) == n, "tvOut.bFinal shape is (n+2) x n");
call check(maxc(maxc(abs(tvOut.bFinal - bFinalDirect))) < 1e-10,
    "tvOut.bFinal exactly matches a direct quaidsTVPStateToFullB() call at the final smoothed state");


/* ==========================================================================
   Check 4: tvpCtl.smooth = 0 with the default othnam -- smoothedState/
   smoothedStateCov are the documented placeholders, bFinal is derived
   from the FILTERED final state, and the MLE itself is unchanged.
   ========================================================================== */

tvpCtlFiltered = quaidsTVPControlCreate();
tvpCtlFiltered.smooth = 0;
tvpCtlFiltered.maxIters = capIters;
tvpCtlFiltered.maxTime = capMinutes;

screen off;
t0 = hsec;
struct quaidsTVPOut tvOutFiltered;
tvOutFiltered = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlFiltered);
secsFiltered = (hsec - t0)/100;
screen on;

print "" $+ "smooth = 0 fit: " $+ ftocv(secsFiltered, 1, 2) $+ " s";

call check(secsFiltered < budget, "smooth = 0 fit finishes inside the per-fit wall-clock budget");
call check(tvOutFiltered.wnam[1] $== "W1", "default tvpCtl.othnam (0): tvOut.wnam auto-generates W1..Wn");
call check(maxc(abs(tvOutFiltered.Qfit - tvOut.Qfit)) < 1e-10 and tvOutFiltered.mleRetcode == tvOut.mleRetcode,
    "tvpCtl.smooth = 0: the Q MLE is identical to the smooth = 1 fit's");
call check(tvOutFiltered.smoothed == 0, "tvpCtl.smooth = 0: tvOut.smoothed echoes 0");
call check(tvOutFiltered.smoothedState == 0, "tvpCtl.smooth = 0: tvOut.smoothedState is 0");
covOrd = getorders(tvOutFiltered.smoothedStateCov);
call check(covOrd[1] == 1 and covOrd[2] == 1 and covOrd[3] == 1 and getmatrix(tvOutFiltered.smoothedStateCov, 1) == 0,
    "tvpCtl.smooth = 0: tvOut.smoothedStateCov is the documented 1x1x1 zero-array placeholder");

bFinalFilteredDirect = quaidsTVPStateToFullB(tvOutFiltered.filteredState[., tobs], n1);
call check(maxc(maxc(abs(tvOutFiltered.bFinal - bFinalFilteredDirect))) < 1e-10,
    "tvpCtl.smooth = 0: tvOut.bFinal exactly matches quaidsTVPStateToFullB() at the final FILTERED state");


/* ==========================================================================
   Check 5: an exhausted cap returns promptly with mleRetcode /= 0 instead
   of running on. maxIters = 1 stops CMLMT on its iteration limit; a
   1e-6-minute maxTime stops it after its first iteration on the time
   limit (CMLMT retcode 11, "Maximum Time Exceeded").
   ========================================================================== */

tvpCtlCapIter = quaidsTVPControlCreate();
tvpCtlCapIter.smooth = 0;
tvpCtlCapIter.maxIters = 1;

screen off;
t0 = hsec;
struct quaidsTVPOut tvOutCapIter;
tvOutCapIter = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlCapIter);
secsCapIter = (hsec - t0)/100;
screen on;

print "" $+ "maxIters = 1 fit: " $+ ftocv(secsCapIter, 1, 2) $+ " s, retcode " $+ ftocv(tvOutCapIter.mleRetcode, 1, 0);

call check(tvOutCapIter.mleRetcode /= 0, "maxIters = 1: the fit returns with mleRetcode /= 0 (iteration cap reported, not hidden)");
call check(secsCapIter < budget, "maxIters = 1: the capped fit returns inside the wall-clock budget");
call check(rows(tvOutCapIter.Qfit) == k_states and minc(tvOutCapIter.Qfit) > 0,
    "maxIters = 1: the capped fit still returns a full, positive Qfit");

tvpCtlCapTime = quaidsTVPControlCreate();
tvpCtlCapTime.smooth = 0;
tvpCtlCapTime.maxTime = 1e-6;

screen off;
t0 = hsec;
struct quaidsTVPOut tvOutCapTime;
tvOutCapTime = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlCapTime);
secsCapTime = (hsec - t0)/100;
screen on;

print "" $+ "maxTime = 1e-6 fit: " $+ ftocv(secsCapTime, 1, 2) $+ " s, retcode " $+ ftocv(tvOutCapTime.mleRetcode, 1, 0);

call check(tvOutCapTime.mleRetcode == 11, "maxTime = 1e-6 minutes: the fit returns with CMLMT's time-limit code (mleRetcode == 11)");
call check(secsCapTime < budget, "maxTime = 1e-6 minutes: the capped fit returns inside the wall-clock budget");


/* ==========================================================================
   Check 6: printQuaidsTVP() runs without error (smoke), including the
   non-converged branch.
   ========================================================================== */

screen off;
call printQuaidsTVP(tvOut);
call printQuaidsTVP(tvOutFiltered);
call printQuaidsTVP(tvOutCapIter);
screen on;
call check(1, "printQuaidsTVP() runs without error on the smoothed, filtered-only and capped (non-converged) fits");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "QUAIDS TVP STAGE 6 TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "QUAIDS TVP STAGE 6 TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
