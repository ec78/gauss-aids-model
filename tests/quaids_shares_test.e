/*
** quaids_shares_test.e
**
** Milestone 16: validates quaidsSharesFit()/printQuaidsShares()
** (src/quaidsshares.src) -- the model-implied predicted budget share
** vector (and its delta-method covariance/standard errors) at an
** arbitrary evaluation point.
**
** Checks, in order: (1) the point estimate matches a fresh, independently
** hand-evaluated version of the share formula, written directly in this
** test rather than calling any src/ proc, on both an AIDS fixture
** (nint=0) and a QUAIDS fixture (nint=1) -- not just "it runs"; (2) exact
** adding-up (sum(w)==1), which follows from qOut.bestB's own adding-up
** construction and holds regardless of evaluation point; (3) shape/
** finiteness/non-negativity of se and v, and that se is exactly
** sqrt(diag(v)); (4) that different evaluation points give genuinely
** different predicted shares (non-vacuous); (5) printQuaidsShares() runs
** without error.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_shares_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
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

/* Independent, hand-evaluated share formula -- deliberately NOT calling
   _quaidsSharesAt() or quaidsElas_(), so this is a real cross-check
   against a fresh derivation, not a tautology. */
proc (1) = handComputedShare(b, intcptPt, pricesPt, totexpPt, struct quaidsControl aCtl);
    local nint, n, alpha, gama, _beta, lambda, a_p, lx, lx2, b_p, wPt;
    nint = rows(intcptPt);
    n = rows(pricesPt);
    alpha = intcptPt'b[1:nint, .];
    gama = b[nint+1:nint+n, .];
    _beta = b[nint+n+1, .];
    a_p = aCtl.alpha0 + alpha*pricesPt + .5*pricesPt'gama*pricesPt;
    lx = totexpPt - a_p;
    wPt = alpha' + gama'pricesPt + _beta'lx;
    if not aCtl.linear;
        b_p = exp(_beta*pricesPt);
        lambda = b[nint+n+2, .];
        lx2 = (lx^2)./b_p;
        wPt = wPt + lambda'lx2;
    endif;
    retp(wPt);
endp;


/* ===================================================================
   AIDS block (nint=0)
   =================================================================== */

tobs = 3000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsCurvatureSyntheticDGP(tobs);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "AIDS: starting fit converged (precondition)");

n = qOut.n;
nint = qOut.nint;
intcptMean = meanc(qOut.intcptFull);
pricesMean = meanc(prices);
totexpMean = meanc(totexp);

struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);

wHand = handComputedShare(qOut.bestB, intcptMean, pricesMean, totexpMean, aCtl);
call check(maxc(abs(sharesOut.w - wHand)) < 1e-10, "AIDS: quaidsSharesFit.w matches an independent hand-evaluated formula");

call check(abs(sumc(sharesOut.w) - 1) < 1e-8, "AIDS: predicted shares sum to 1 (adding-up)");

call check(rows(sharesOut.w) == n and cols(sharesOut.w) == 1, "AIDS: w has shape n x 1");
call check(rows(sharesOut.se) == n and cols(sharesOut.se) == 1, "AIDS: se has shape n x 1");
call check(rows(sharesOut.v) == n and cols(sharesOut.v) == n, "AIDS: v has shape n x n");
call check(minc(sharesOut.se) >= 0, "AIDS: se are all non-negative");
call check(sumc(sharesOut.se .== sharesOut.se) == n, "AIDS: se contains no NaN/missing values");
call check(maxc(abs(sharesOut.se - sqrt(diag(sharesOut.v)))) < 1e-10, "AIDS: se exactly equals sqrt(diag(v))");

/* Non-vacuous: a different evaluation point gives a genuinely different
   predicted share. */
pricesShift = pricesMean;
pricesShift[1] = pricesShift[1] + ln(1.2);
struct quaidsSharesOut sharesShift;
sharesShift = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesShift, totexpMean, aCtl);
call check(maxc(abs(sharesShift.w - sharesOut.w)) > 1e-6, "AIDS: a shifted price point gives a genuinely different predicted share");

call printQuaidsShares(sharesOut);
call check(1, "AIDS: printQuaidsShares() runs without error");


/* ===================================================================
   QUAIDS block (nint=1)
   =================================================================== */

tobsQ = 3000;
{ wQ, intcptQ, pricesQ, totexpQ, instrQ, trueParamsQ } = _quaidsSyntheticDGP(tobsQ, 204, 1, 1);

struct quaidsControl aCtlQ;
aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 1;
aCtlQ.err = .0001;

struct quaidsOut qOutQ;
qOutQ = quaidsFit(wQ, intcptQ, pricesQ, totexpQ, instrQ, aCtlQ);
call check(qOutQ.converged == 1, "QUAIDS: starting fit converged (precondition)");

nQ = qOutQ.n;
nintQ = qOutQ.nint;
intcptMeanQ = meanc(qOutQ.intcptFull);
pricesMeanQ = meanc(pricesQ);
totexpMeanQ = meanc(totexpQ);

struct quaidsSharesOut sharesOutQ;
sharesOutQ = quaidsSharesFit(qOutQ.bestB, qOutQ.bestV, intcptMeanQ, pricesMeanQ, totexpMeanQ, aCtlQ);

wHandQ = handComputedShare(qOutQ.bestB, intcptMeanQ, pricesMeanQ, totexpMeanQ, aCtlQ);
call check(maxc(abs(sharesOutQ.w - wHandQ)) < 1e-10, "QUAIDS: quaidsSharesFit.w matches an independent hand-evaluated formula");

call check(abs(sumc(sharesOutQ.w) - 1) < 1e-8, "QUAIDS: predicted shares sum to 1 (adding-up)");

call check(rows(sharesOutQ.w) == nQ and cols(sharesOutQ.w) == 1, "QUAIDS: w has shape n x 1");
call check(rows(sharesOutQ.se) == nQ and cols(sharesOutQ.se) == 1, "QUAIDS: se has shape n x 1");
call check(rows(sharesOutQ.v) == nQ and cols(sharesOutQ.v) == nQ, "QUAIDS: v has shape n x n");
call check(minc(sharesOutQ.se) >= 0, "QUAIDS: se are all non-negative");
call check(sumc(sharesOutQ.se .== sharesOutQ.se) == nQ, "QUAIDS: se contains no NaN/missing values");
call check(maxc(abs(sharesOutQ.se - sqrt(diag(sharesOutQ.v)))) < 1e-10, "QUAIDS: se exactly equals sqrt(diag(v))");

call printQuaidsShares(sharesOutQ);
call check(1, "QUAIDS: printQuaidsShares() runs without error");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "SHARES TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "SHARES TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
