/*
** quaids_welfare_test.e
**
** Milestone 11: validates quaidsWelfareFit() (src/quaidswelfare.src) --
** exact compensating variation (CV) and equivalent variation (EV) for a
** price change, holding nominal expenditure fixed.
**
** Checked via exact algebraic identities (like Milestone 5's Engel/
** Cournot/homogeneity checks), not tolerance-based approximations, plus
** one limiting-case numerical check tying the exact measure back to a
** well-known approximation:
**
**   1. Zero price change: CV == EV == 0 exactly.
**   2. Round-trip exactness: plugging e(p1,u0) back into the ORIGINAL
**      indirect utility formula at p1 returns u0 exactly -- this is the
**      defining property of an expenditure function (it is genuinely the
**      inverse of the indirect utility function), so it should hold to
**      floating-point precision if the closed-form inversion is correct.
**   3. First-order consistency: for a small price change on one good,
**      CV/EV should approach the standard Marshallian first-order
**      approximation (share-weighted expenditure change) as the price
**      change shrinks toward zero.
**
** Exercised on both a QUAIDS fit (lambda != 0) and an AIDS/LA-AIDS fit
** (aCtl.linear=1, no lambda term) from the same _quaidsSyntheticDGP
** fixture already used by tests/quaids_elasticities_test.e and others --
** the formula is unified across both model choices (verified during
** implementation: setting lambda=0 collapses the QUAIDS expenditure
** function to the simpler, independently-verified AIDS one).
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_welfare_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidswelfare.src;
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

/* Recomputes ln V(x,p) directly from the model equations (not via
   quaidsWelfareFit()), for the round-trip identity check -- mirrors
   quaidscurvature.src's/tests' own "recompute independently, don't just
   call the function under test" validation style. */
proc (1) = modelLnV(b, nint, n, alpha0, intcptPt, pricesPt, totexpPt, linear);
    local alpha, gama, beta, lambda, a_p, b_p, lamP, lx;
    alpha = b[1:1+nint, .];
    gama = b[1+nint+1:1+nint+n, .];
    beta = b[1+nint+n+1, .];
    a_p = alpha0 + intcptPt'alpha*pricesPt + .5*pricesPt'gama*pricesPt;
    b_p = exp(beta*pricesPt);
    if linear;
        lamP = 0;
    else;
        lambda = b[1+nint+n+2, .];
        lamP = lambda*pricesPt;
    endif;
    lx = totexpPt - a_p;
    retp(1/(b_p/lx + lamP));
endp;

/* Recomputes the model-implied share at a point (same helper pattern as
   tests/quaids_elasticities_test.e's modelShareAt()), for the first-order
   approximation check. */
proc (1) = modelShareAt(b, nint, n, alpha0, intcptPt, pricesPt, totexpPt, linear);
    local alpha, gama, beta, lambda, a_p, lx, b_p, lx2, wPt;
    alpha = b[1:1+nint, .];
    gama = b[1+nint+1:1+nint+n, .];
    beta = b[1+nint+n+1, .];
    a_p = alpha0 + intcptPt'alpha*pricesPt + .5*pricesPt'gama*pricesPt;
    lx = totexpPt - a_p;
    wPt = alpha'intcptPt + gama'pricesPt + beta'lx;
    if not linear;
        lambda = b[1+nint+n+2, .];
        b_p = exp(beta*pricesPt);
        lx2 = (lx^2)./b_p;
        wPt = wPt + lambda'lx2;
    endif;
    retp(wPt);
endp;

/* Runs the full battery of checks for one fitted model (QUAIDS or AIDS). */
proc (0) = checkWelfare(label, struct quaidsOut qOut, struct quaidsControl aCtl,
        prices, totexp, intcptFullMean);
    local n, nint, pricesPt0, pricesPt1, totexpPt0, b, v, lnE_p1_u0,
        u0_roundtrip, dlnp, pricesPtTiny, wPt, approx;
    struct quaidsWelfareOut wOut, wOutZero, wOutTiny;

    n = qOut.n;
    nint = qOut.nint;
    b = qOut.bestB;
    v = qOut.bestV;

    pricesPt0 = meanc(prices);
    totexpPt0 = meanc(totexp);
    pricesPt1 = pricesPt0;
    pricesPt1[1] = pricesPt1[1] + 0.05;

    /* --- Check 1: zero price change --- */
    wOutZero = quaidsWelfareFit(b, v, intcptFullMean, pricesPt0, pricesPt0, totexpPt0, aCtl);
    call check(wOutZero.cv == 0, label $+ ": CV == 0 exactly at zero price change");
    call check(wOutZero.ev == 0, label $+ ": EV == 0 exactly at zero price change");

    /* --- Check 2: round-trip exactness --- */
    wOut = quaidsWelfareFit(b, v, intcptFullMean, pricesPt0, pricesPt1, totexpPt0, aCtl);
    lnE_p1_u0 = ln(wOut.cv + exp(totexpPt0));
    u0_roundtrip = modelLnV(b, nint, n, aCtl.alpha0, intcptFullMean, pricesPt1, lnE_p1_u0, aCtl.linear);
    call check(abs(u0_roundtrip - wOut.u0) < 1e-8,
        label $+ ": e(p1,u0) fed back into V(p1,.) returns u0 exactly (expenditure-function inverse property)");

    /* --- Check 3: first-order consistency --- */
    dlnp = 1e-5;
    pricesPtTiny = pricesPt0;
    pricesPtTiny[1] = pricesPtTiny[1] + dlnp;
    wOutTiny = quaidsWelfareFit(b, v, intcptFullMean, pricesPt0, pricesPtTiny, totexpPt0, aCtl);
    wPt = modelShareAt(b, nint, n, aCtl.alpha0, intcptFullMean, pricesPt0, totexpPt0, aCtl.linear);
    approx = exp(totexpPt0)*wPt[1]*dlnp;
    call check(abs(wOutTiny.cv - approx) < 0.01*abs(approx) + 1e-6,
        label $+ ": CV for a tiny price change matches the first-order (share-weighted) approximation");
    call check(abs(wOutTiny.ev - approx) < 0.01*abs(approx) + 1e-6,
        label $+ ": EV for a tiny price change matches the first-order (share-weighted) approximation");

    /* --- Shape/finiteness checks on standard errors --- */
    call check(wOut.seCV >= 0 and wOut.seCV == wOut.seCV, label $+ ": seCV is non-negative and not missing/NaN");
    call check(wOut.seEV >= 0 and wOut.seEV == wOut.seEV, label $+ ": seEV is non-negative and not missing/NaN");

    /* --- CV and EV agree in sign (both represent the same welfare
       direction for a given price change) --- */
    call check((wOut.cv > 0 and wOut.ev > 0) or (wOut.cv < 0 and wOut.ev < 0),
        label $+ ": CV and EV agree in sign for this price change");

    call printQuaidsWelfare(wOut);
    call check(1, label $+ ": printQuaidsWelfare() runs without error");
endp;


{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 204, 1, 1);

struct quaidsControl aCtlQ;
aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 1;
aCtlQ.err = .0001;

struct quaidsOut qOutQ;
qOutQ = quaidsFit(w, intcpt, prices, totexp, instr, aCtlQ);
call check(qOutQ.converged == 1, "QUAIDS starting fit converged");

intcptPtQ = meanc(qOutQ.intcptFull);
call checkWelfare("QUAIDS", qOutQ, aCtlQ, prices, totexp, intcptPtQ);


struct quaidsControl aCtlL;
aCtlL = quaidsControlCreate();
aCtlL.linear = 1;
aCtlL.maxiter = 100;
aCtlL.homogenous = 1;
aCtlL.err = .0001;

struct quaidsOut qOutL;
qOutL = quaidsFit(w, intcpt, prices, totexp, instr, aCtlL);
call check(qOutL.converged == 1, "AIDS starting fit converged");

intcptPtL = meanc(qOutL.intcptFull);
call checkWelfare("AIDS", qOutL, aCtlL, prices, totexp, intcptPtL);


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "WELFARE TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "WELFARE TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
