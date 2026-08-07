/*
** quaids_elasticities_test.e
**
** Milestone 5: validates the quaidsElasFit()/printQuaidsElas() split
** (parity against the pre-split quaidsElas() on the same point) and the
** "arbitrary evaluation point" capability the roadmap asked to generalize
** -- quaidsElas_()/quaidsElasFit() already accepted any point as an
** argument; what was missing was a silent, struct-returning entry point
** and confidence that it's actually correct away from the four points
** quaids() has always used (mean/Q1/median/Q3).
**
** The "is it actually correct" check uses three EXACT algebraic
** identities that any valid AIDS/QUAIDS elasticity set must satisfy
** given adding-up/homogeneity hold (which this estimator always imposes
** by construction) -- not tolerance-based approximations, since these are
** consequences of the functional form, not separately estimated
** quantities:
**   Engel aggregation:     sum_i(w_i * er_i)        = 1
**   Cournot aggregation:   sum_i(w_i * ep_ij) + w_j  = 0, for each price j
**   Homogeneity (elast.):  sum_j(ep_ij) + er_i       = 0, for each good i
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_elasticities_test.e
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

/* Checks the three exact elasticity identities at one point, given the
   model-implied share vector wPt (NOT the noisy observed share -- the
   identities are properties of the estimated coefficients evaluated at a
   point, they hold regardless of how noisy any one observation's actual
   share is). */
proc (0) = checkIdentities(label, struct quaidsElasOut elasOut, wPt);
    local engel, cournot, homog;
    engel = sumc(wPt.*elasOut.er);
    cournot = (wPt'elasOut.ep)' + wPt;
    homog = sumc(elasOut.ep') + elasOut.er;
    call check(abs(engel - 1) < 1e-8, label $+ ": Engel aggregation sum(w*er)==1");
    call check(maxc(abs(cournot)) < 1e-8, label $+ ": Cournot aggregation sum(w*ep)+w==0");
    call check(maxc(abs(homog)) < 1e-8, label $+ ": elasticity homogeneity sum(ep)+er==0");
endp;

/* Milestone 16: the model-implied share at a point used to be recomputed
   here via a private modelShareAt() helper duplicating quaidsElas_()'s
   own formula, needed only so this test could check the identities below
   against the RIGHT w. Now uses the public quaidsSharesFit() instead --
   the Engel/Cournot/homogeneity identity checks below are themselves a
   real cross-check that its point estimate is correct (they are exact
   algebraic identities that only hold given the true model-implied w). */

aCtl = quaidsControlCreate;
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 204, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);


/* --- Parity: quaidsElasFit()/printQuaidsElas() vs. quaidsElas_() directly,
   at the sample mean (the same point quaids() has always used first). --- */

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

elasMean = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);

{ erDirect, epDirect, epcDirect } = quaidsElas_(qOut.bestB, intcptMean, pricesMean, totexpMean, aCtl);
call check(maxc(abs(elasMean.er - erDirect)) == 0, "quaidsElasFit.er matches quaidsElas_ exactly");
call check(maxc(maxc(abs(elasMean.ep - epDirect))) == 0, "quaidsElasFit.ep matches quaidsElas_ exactly");
call check(maxc(maxc(abs(elasMean.epc - epcDirect))) == 0, "quaidsElasFit.epc matches quaidsElas_ exactly");
call check(rows(elasMean.ser) == n, "quaidsElasFit.ser has n rows");
call check(minc(minc(elasMean.sep)) >= 0, "quaidsElasFit.sep (standard errors) are all non-negative");
call check(minc(minc(elasMean.sepc)) >= 0, "quaidsElasFit.sepc (standard errors) are all non-negative");

sharesMean = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call checkIdentities("at sample mean", elasMean, sharesMean.w);


/* --- Genuinely arbitrary evaluation points: NOT mean/Q1/median/Q3. --- */

/* A specific observation's own data. */
pt = 500;
intcptPt1 = qOut.intcptFull[pt,.]';
pricesPt1 = prices[pt,.]';
totexpPt1 = totexp[pt];
elasObs = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt1, pricesPt1, totexpPt1, aCtl);
sharesObs = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt1, pricesPt1, totexpPt1, aCtl);
call checkIdentities("at observation 500", elasObs, sharesObs.w);

/* A fully synthetic counterfactual point (e.g. a hypothetical 20% price
   increase on good 1 relative to the sample mean), demonstrating the
   point doesn't need to come from the data at all. */
intcptPt2 = intcptMean;
pricesPt2 = pricesMean;
pricesPt2[1] = pricesPt2[1] + ln(1.20);
totexpPt2 = totexpMean;
elasCounterfactual = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt2, pricesPt2, totexpPt2, aCtl);
sharesCf = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt2, pricesPt2, totexpPt2, aCtl);
call checkIdentities("at a synthetic counterfactual (20% price-1 increase)", elasCounterfactual, sharesCf.w);

call check(maxc(maxc(abs(elasCounterfactual.er - elasMean.er))) > 1e-6,
    "counterfactual point gives genuinely different elasticities than the sample mean");


/* --- printQuaidsElas() runs without error on a custom-point struct. --- */
call printQuaidsElas(elasObs);
call check(1, "printQuaidsElas() runs without error on a custom-point quaidsElasOut");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "ELASTICITIES TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "ELASTICITIES TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
