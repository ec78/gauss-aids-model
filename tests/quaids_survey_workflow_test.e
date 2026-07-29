/*
** quaids_survey_workflow_test.e
**
** Milestone 25 seed: validates quaidsSurveyWorkflowFit() as an opt-in
** survey/microdata workflow layer. Sampling weights affect the
** representative evaluation point and post-estimation summaries; the
** underlying quaidsFit() estimator remains the same as quaidsWorkflowFit().
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_survey_workflow_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstests.src
#include ../src/quaidswelfare.src
#include ../src/quaidsrobust.src
#include ../src/quaidsdiagnostics.src
#include ../src/quaidsworkflow.src
#include ../src/quaidssurvey.src
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

weight = seqa(1, 1, tobs);
weight = weight/sumc(weight)*tobs;

struct quaidsWorkflowOut wfBase;
struct quaidsWorkflowOut wfSurvey;
wfBase = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);
wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0, weight);

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
wgt = weight/sumc(weight);
mW = (wgt'(qOut.intcptFull~prices~totexp))';
intcptW = mW[1:1+nint];
pricesW = mW[1+nint+1:1+nint+n];
totexpW = mW[1+nint+n+1];

struct quaidsSharesOut sharesOut;
struct quaidsElasOut elasOut;
struct quaidsSharesOut sharesRobustOut;
struct quaidsElasOut elasRobustOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptW, pricesW, totexpW, aCtl);
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptW, pricesW, totexpW, aCtl);
sharesRobustOut = quaidsSharesFit(wfSurvey.robustBestB, wfSurvey.robustBestV, intcptW, pricesW, totexpW, aCtl);
elasRobustOut = quaidsElasFit(wfSurvey.robustBestB, wfSurvey.robustBestV, intcptW, pricesW, totexpW, aCtl);

call check(wfSurvey.converged == wfBase.converged and maxc(maxc(abs(wfSurvey.bestB - wfBase.bestB))) == 0,
    "survey workflow leaves the underlying estimator unchanged");
call check(wfSurvey.surveyWeighted == 1 and wfSurvey.surveyWeightValid == 1,
    "survey workflow marks valid weighted evaluation");
call check(wfSurvey.surveyWeightSum == sumc(weight) and wfSurvey.surveyWeightNPositive == tobs,
    "survey workflow weight sum/count diagnostics match input");
call check(wfSurvey.surveyWeightMin == minc(weight) and wfSurvey.surveyWeightMax == maxc(weight),
    "survey workflow weight range diagnostics match input");
call check(maxc(abs(wfSurvey.evalIntcpt - intcptW)) < 1e-12 and
    maxc(abs(wfSurvey.evalPrices - pricesW)) < 1e-12 and
    abs(wfSurvey.evalTotexp - totexpW) < 1e-12,
    "survey workflow evaluation point matches manual weighted mean");
call check(maxc(abs(wfSurvey.evalPrices - wfBase.evalPrices)) > 1e-8 or
    abs(wfSurvey.evalTotexp - wfBase.evalTotexp) > 1e-8,
    "nonconstant weights change the representative evaluation point");
call check(maxc(abs(wfSurvey.shares - sharesOut.w)) == 0,
    "survey workflow predicted shares match direct weighted-point call");
call check(maxc(abs(wfSurvey.incomeElas - elasOut.er)) == 0,
    "survey workflow income elasticities match direct weighted-point call");
call check(maxc(maxc(abs(wfSurvey.priceElas - elasOut.ep))) == 0,
    "survey workflow uncompensated price elasticities match direct weighted-point call");
call check(maxc(maxc(abs(wfSurvey.compPriceElas - elasOut.epc))) == 0,
    "survey workflow compensated price elasticities match direct weighted-point call");
call check(wfSurvey.postRobustValid == 1 and maxc(abs(wfSurvey.sharesRobustSE - sharesRobustOut.se)) == 0,
    "survey workflow robust share SE match weighted-point robust covariance propagation");
call check(maxc(abs(wfSurvey.incomeElasRobustSE - elasRobustOut.ser)) == 0,
    "survey workflow robust income elasticity SE match weighted-point robust covariance propagation");
call check(maxc(maxc(abs(wfSurvey.priceElasRobustSE - elasRobustOut.sep))) == 0,
    "survey workflow robust price elasticity SE match weighted-point robust covariance propagation");

weightOnes = ones(tobs, 1);
struct quaidsWorkflowOut wfOnes;
wfOnes = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0, weightOnes);
call check(maxc(abs(wfOnes.evalPrices - wfBase.evalPrices)) < 1e-12 and abs(wfOnes.evalTotexp - wfBase.evalTotexp) < 1e-12,
    "constant weights reproduce the default workflow evaluation point");
call check(maxc(abs(wfOnes.shares - wfBase.shares)) < 1e-10 and maxc(abs(wfOnes.incomeElas - wfBase.incomeElas)) < 1e-10,
    "constant weights reproduce default workflow post-estimation point estimates");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "SURVEY WORKFLOW TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "SURVEY WORKFLOW TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
