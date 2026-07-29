/*
** quaids_workflow_test.e
**
** Milestone 21 seed: validates quaidsWorkflowFit() as a thin composition
** layer over quaidsFit(), quaidsSharesFit(), quaidsElasFit(), and
** quaidsRobustFit(). The point is parity with the explicit manual workflow,
** not new econometrics.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_workflow_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
#include ../src/quaidsworkflow.src
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

seed = 204;
tobs = 1000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, 1, 1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsSharesOut sharesOut;
struct quaidsElasOut elasOut;
struct quaidsRobustOut robustOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
robustOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);

call check(wfOut.converged == 1 and wfOut.postValid == 1 and wfOut.robustValid == 1,
    "workflow fit converged and computed post-estimation outputs");
call check(wfOut.model $== qOut.model and wfOut.n == qOut.n and wfOut.nobs == qOut.nobs,
    "workflow metadata matches quaidsFit");
call check(maxc(maxc(abs(wfOut.bestB - qOut.bestB))) == 0,
    "workflow bestB matches quaidsFit exactly");
call check(maxc(abs(wfOut.evalIntcpt - intcptMean)) == 0 and maxc(abs(wfOut.evalPrices - pricesMean)) == 0 and wfOut.evalTotexp == totexpMean,
    "workflow sample-mean evaluation point matches manual construction");
call check(maxc(abs(wfOut.shares - sharesOut.w)) == 0,
    "workflow predicted shares match quaidsSharesFit exactly");
call check(maxc(abs(wfOut.sharesSE - sharesOut.se)) == 0,
    "workflow share SE match quaidsSharesFit exactly");
call check(maxc(abs(wfOut.incomeElas - elasOut.er)) == 0,
    "workflow income elasticities match quaidsElasFit exactly");
call check(maxc(maxc(abs(wfOut.priceElas - elasOut.ep))) == 0,
    "workflow uncompensated price elasticities match quaidsElasFit exactly");
call check(maxc(maxc(abs(wfOut.compPriceElas - elasOut.epc))) == 0,
    "workflow compensated price elasticities match quaidsElasFit exactly");
call check(maxc(maxc(abs(wfOut.robustSE - robustOut.se))) == 0 and wfOut.robustNClusters == robustOut.nClusters,
    "workflow robust SE match quaidsRobustFit exactly");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "WORKFLOW TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "WORKFLOW TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
