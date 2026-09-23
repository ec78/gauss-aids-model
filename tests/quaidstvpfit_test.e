/*
** quaidstvpfit_test.e
**
** TVP-AIDS initiative, Stage 6: validates quaidsTVPFit()/printQuaidsTVP()
** (src/quaidstvpfit.src) -- the consolidated public entry point. Adds no
** new estimation math of its own (see that file's own header), so every
** check here is INTERNAL CONSISTENCY (exact, not loose): quaidsTVPFit()'s
** output must exactly match the SAME already-validated Stage 1-5 procs
** (_quaidsTVPStoneIndex/_quaidsTVPBuildZ/_quaidsTVPMLEFit/
** _quaidsTVPSmoothFit/quaidsTVPStateToFullB) called directly on the
** identical inputs -- no synthetic-recovery-of-truth claim is needed or
** attempted here, since that correctness was already established by
** Stages 1-5's own test files.
**
** Requires `library cmlmt, tsmt, sslib;` and the GAUSS26_CFG override
** documented in CLAUDE.md, same as quaidstvp_kalman_test.e/
** quaidstvp_mle_test.e/quaidstvp_smooth_test.e -- run directly with
** `tgauss -b -x quaidstvpfit_test.e`, GAUSS26_CFG must already point at
** tests/gauss26_cfg_override.
**
** Synthetic data here is deliberately NOT a validated AIDS/QUAIDS DGP
** fixture (unlike tests/quaidsfixtures.src's own
** _quaidsTVPStaticSyntheticDGP/_quaidsTVPDynamicSyntheticDGP, which both
** already return pre-built relative-price/Stone-deflated pieces, not the
** raw absolute-price/totexp/full-w inputs quaidsTVPFit()'s own public
** contract actually takes) -- any plausible, positive, row-summing-to-1
** share matrix works equally well for an internal-consistency check, so a
** small self-contained generator is used instead of forcing a fixture
** mismatch.
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
   Setup: a small, plausible (not a validated true DGP -- see this file's
   own header) synthetic dataset, n1 = 2 (3 goods).
   ========================================================================== */

rndseed 909;

n = 3;
n1 = n - 1;
tobs = 60;
k_states = 2*n1 + n1*(n1+1)/2;

prices = 1 + 0.1*rndn(tobs, n);
totexp = 5 + 0.2*rndn(tobs, 1);
raw = 0.1 + abs(rndn(tobs, n));
w = raw ./ sumc(raw');

H = 0.01*eye(n1);
q0 = 0.001*ones(k_states, 1);


/* ==========================================================================
   Check 1: quaidsTVPFit() (tvpCtl.smooth = 1, the default) matches a
   direct _quaidsTVPStoneIndex + _quaidsTVPBuildZ + _quaidsTVPMLEFit call
   chain on the identical inputs, exactly.
   ========================================================================== */

tvpCtl = quaidsTVPControlCreate();
call check(tvpCtl.smooth == 1, "quaidsTVPControlCreate() default: smooth = 1");

screen off;
struct quaidsTVPOut tvOut;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);

lxDirect = _quaidsTVPStoneIndex(w, prices, totexp, n);
pricesRelDirect = prices[., 1:n1] - prices[., n];
ZarrDirect = _quaidsTVPBuildZ(pricesRelDirect, lxDirect, n1);

struct ssOut sOutDirect;
sOutDirect = _quaidsTVPMLEFit(ZarrDirect, H, w[., 1:n1], q0, n1);
screen on;

call check(tvOut.n == n and tvOut.n1 == n1 and tvOut.nobs == tobs and tvOut.k_states == k_states,
    "tvOut dimensions match n/n1/nobs/k_states");
call check(maxc(abs(tvOut.Qfit - sOutDirect.final_params)) < 1e-10,
    "tvOut.Qfit exactly matches a direct _quaidsTVPMLEFit() call's final_params");
call check(tvOut.mleRetcode == sOutDirect.mleResults.retcode,
    "tvOut.mleRetcode exactly matches a direct _quaidsTVPMLEFit() call's mleResults.retcode");
call check(maxc(maxc(abs(tvOut.filteredState - sOutDirect.kfResults.filtered_state))) < 1e-10,
    "tvOut.filteredState exactly matches a direct _quaidsTVPMLEFit() call's kfResults.filtered_state");
call check(tvOut.wnam[1] $== "W1", "default tvpCtl.othnam (0): tvOut.wnam auto-generates W1..Wn");


/* ==========================================================================
   Check 1b: a caller-supplied tvpCtl.othnam (a real character matrix, not
   the default 0) is honored -- the ONE place in this codebase's history
   that field is ever exercised with a real value (found via a real GAUSS
   `error G0071 : Type mismatch` when this field was still `string`-typed,
   like quaidsControl.othnam -- see quaidsTVPControl's own header in
   quaids.sdf for why it is `matrix`-typed instead).
   ========================================================================== */

/* 0$+ coerces the $|-built character matrix into the legacy form a
   `matrix`-typed struct field requires -- a bare $|-built value fails
   with a real, confirmed `error G0071 : Type mismatch` otherwise (found
   via examples/14_tvp_aids_estimation.e; see quaidsTVPControl's own
   header in src/quaids.sdf). */
customNames = 0$+("Alpha" $| "Beta" $| "Gamma");
tvpCtlNamed = quaidsTVPControlCreate();
tvpCtlNamed.othnam = customNames;

screen off;
struct quaidsTVPOut tvOutNamed;
tvOutNamed = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlNamed);
screen on;

call check(tvOutNamed.wnam[1] $== "Alpha" and tvOutNamed.wnam[2] $== "Beta" and tvOutNamed.wnam[3] $== "Gamma",
    "caller-supplied tvpCtl.othnam (a real character matrix) is honored in tvOut.wnam");


/* ==========================================================================
   Check 2: tvOut.smoothedState/smoothedStateCov match a direct
   _quaidsTVPSmoothFit() call on (sOutDirect.tvpFinal, sOutDirect.kfResults),
   exactly.
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
   Check 4: tvpCtl.smooth = 0 -- smoothedState/smoothedStateCov are 0,
   bFinal is derived from the FILTERED final state instead.
   ========================================================================== */

tvpCtlFiltered = quaidsTVPControlCreate();
tvpCtlFiltered.smooth = 0;

screen off;
struct quaidsTVPOut tvOutFiltered;
tvOutFiltered = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtlFiltered);
screen on;

call check(tvOutFiltered.smoothed == 0, "tvpCtl.smooth = 0: tvOut.smoothed echoes 0");
call check(tvOutFiltered.smoothedState == 0, "tvpCtl.smooth = 0: tvOut.smoothedState is 0");
covOrd = getorders(tvOutFiltered.smoothedStateCov);
call check(covOrd[1] == 1 and covOrd[2] == 1 and covOrd[3] == 1 and getmatrix(tvOutFiltered.smoothedStateCov, 1) == 0,
    "tvpCtl.smooth = 0: tvOut.smoothedStateCov is the documented 1x1x1 zero-array placeholder");

bFinalFilteredDirect = quaidsTVPStateToFullB(tvOutFiltered.filteredState[., tobs], n1);
call check(maxc(maxc(abs(tvOutFiltered.bFinal - bFinalFilteredDirect))) < 1e-10,
    "tvpCtl.smooth = 0: tvOut.bFinal exactly matches quaidsTVPStateToFullB() at the final FILTERED state");


/* ==========================================================================
   Check 5: printQuaidsTVP() runs without error (smoke).
   ========================================================================== */

screen off;
call printQuaidsTVP(tvOut);
call printQuaidsTVP(tvOutFiltered);
screen on;
call check(1, "printQuaidsTVP() runs without error on both the smoothed and filtered-only fits");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "QUAIDS TVP STAGE 6 TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "QUAIDS TVP STAGE 6 TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
