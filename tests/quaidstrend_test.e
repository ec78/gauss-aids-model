/*
** quaidstrend_test.e
**
** TVP-AIDS initiative, Stage 0: validates quaidsTrendFit() (a cheap,
** one-shot LA-AIDS-based screening diagnostic for whether a demand
** system's coefficients show a linear trend over time -- see
** src/quaidstrend.src's own header for the full design rationale, and
** the project's TVP design report for why this exists as Stage 0 ahead
** of the much larger Kalman-filter-based TVP-AIDS effort).
**
** Checks, in order: the exact adding-up/homogeneity identities on BOTH
** the level and trend-slope coefficient blocks (this proc's central
** mathematical claim -- these hold by construction of the shared-design
** GLS mechanism, not because of any separate restriction-imposition
** step, so this is a direct regression guard on that claim); SIZE (does
** the joint trend-block Wald test correctly fail to reject on a
** no-true-trend DGP?); POWER (does it correctly reject on a genuine-
** trend DGP?); the aCtl.homogenous guard; and that printQuaidsTrend()
** runs without error.
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstrend_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidstrend.src
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
   Exact-identity checks: adding-up and homogeneity on BOTH the level and
   trend-slope blocks, to floating-point precision -- the central claim
   this whole design rests on (see quaidstrend.src's Remarks). Checked on
   the no-true-trend fixture; the same identities are re-checked on the
   genuine-trend fixture below to confirm they hold regardless of whether
   there actually IS a trend in the data.
   ========================================================================== */

{ w, intcpt, prices, totexp, instr } = _quaidsTrendSyntheticDGP(400, 12345, 0);

aCtl = quaidsControlCreate;
aCtl.homogenous = 1;

tOut = quaidsTrendFit(w, intcpt, prices, totexp, instr, aCtl);

tol = 1e-8;
call check(abs(sumc(tOut.b0[1, .]') - 1) < tol, "level block: constant row sums to exactly 1");
call check(maxc(abs(sumc(tOut.b0[2:1+tOut.n1, .]'))) < tol, "level block: price rows sum to exactly 0 (homogeneity)");
call check(abs(sumc(tOut.b0[2+tOut.n1, .]')) < tol, "level block: expenditure row sums to exactly 0");
call check(abs(sumc(tOut.b1[1, .]')) < tol, "trend block: intercept-trend row sums to exactly 0");
call check(maxc(abs(sumc(tOut.b1[2:1+tOut.n1, .]'))) < tol, "trend block: price-trend rows sum to exactly 0 (per price)");
call check(abs(sumc(tOut.b1[rows(tOut.b1), .]')) < tol, "trend block: expenditure-trend row sums to exactly 0");

call check(tOut.n == 4, "metadata: n == 4");
call check(tOut.n1 == 3, "metadata: n1 == n-1 == 3");
call check(tOut.ng == 1+0+3+1+1, "metadata: ng == 1+nint+n1+1(lx)+nu");
call check(tOut.ngTrend == 1+3+1, "metadata: ngTrend == 1(t)+n1+1(t*lx)");
call check(rows(tOut.b0) == tOut.ng and cols(tOut.b0) == tOut.n, "b0 shape: ng x n");
call check(rows(tOut.b1) == tOut.ngTrend and cols(tOut.b1) == tOut.n, "b1 shape: ngTrend x n");
call check(rows(tOut.se0) == rows(tOut.b0) and cols(tOut.se0) == cols(tOut.b0), "se0 same shape as b0");
call check(rows(tOut.se1) == rows(tOut.b1) and cols(tOut.se1) == cols(tOut.b1), "se1 same shape as b1");
call check(not scalmiss(tOut.trendStat) and tOut.trendStat >= 0, "trendStat finite and non-negative");


/* ==========================================================================
   SIZE: no true coefficient drift by construction -- the joint trend test
   should fail to reject.
   ========================================================================== */

call check(tOut.trendPval > 0.05, "trend test SIZE: fails to reject on a no-true-trend DGP (pval > 0.05)");
call check(tOut.trendDf == tOut.ngTrend*(tOut.n-1), "trend test df == ngTrend*(n-1)");


/* ==========================================================================
   POWER: genuine linear-in-t drift by construction -- the joint trend
   test should reject, decisively.
   ========================================================================== */

{ wD, intcptD, pricesD, totexpD, instrD } = _quaidsTrendSyntheticDGP(400, 12345, 1);

tOutD = quaidsTrendFit(wD, intcptD, pricesD, totexpD, instrD, aCtl);

call check(tOutD.trendPval < 0.001, "trend test POWER: rejects decisively on a genuine-trend DGP (pval < 0.001)");
call check(abs(sumc(tOutD.b1[1, .]')) < tol, "trend block still sums to exactly 0 (intercept-trend) even with a true trend present");
call check(maxc(abs(sumc(tOutD.b1[2:1+tOutD.n1, .]'))) < tol, "trend block still sums to exactly 0 (price-trend) even with a true trend present");


/* ==========================================================================
   Guard: aCtl.homogenous must be 1 (see quaidstrend.src's Remarks for why).
   tgauss's own "run;" halting on an unguarded errorlog means this must be
   checked as a separate guard-error-case script, not inline here (the
   same reasoning documented for every other guard test in this project,
   the tests under tests/guard_error_cases) -- see
   tests/guard_error_cases/trend_requires_homogenous.e.
   ========================================================================== */


/* ==========================================================================
   printQuaidsTrend() runs without error on a real fit.
   ========================================================================== */

call printQuaidsTrend(tOut);
call check(1, "printQuaidsTrend() ran without error");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "QUAIDS TREND TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "QUAIDS TREND TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
