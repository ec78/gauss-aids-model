/*
** quaids_workflow_test.e
**
** Milestone 21 seed: validates quaidsWorkflowFit() as a thin composition
** layer over quaidsFit(), quaidsSharesFit(), quaidsElasFit(), and
** quaidsRobustFit(); Milestone 24 also validates the workflow's compact
** quaidsPreflight() summary fields. quaidsWorkflowScenarioFit() is the
** explicit welfare-scenario extension. The point is parity with the
** explicit manual workflow, not new econometrics.
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
#include ../src/quaidstests.src
#include ../src/quaidswelfare.src
#include ../src/quaidsrobust.src
#include ../src/quaidsdiagnostics.src
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

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0);

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, 0);

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
robustOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);
{ bRobustFull, vRobustFull } = quaidsRobustCovariance(qOut, robustOut, aCtl);
sharesRobustOut = quaidsSharesFit(bRobustFull, vRobustFull, intcptMean, pricesMean, totexpMean, aCtl);
elasRobustOut = quaidsElasFit(bRobustFull, vRobustFull, intcptMean, pricesMean, totexpMean, aCtl);

call check(wfOut.converged == 1 and wfOut.postValid == 1 and wfOut.robustValid == 1,
    "workflow fit converged and computed post-estimation outputs");
call check(wfOut.model $== qOut.model and wfOut.n == qOut.n and wfOut.nobs == qOut.nobs,
    "workflow metadata matches quaidsFit");
call check(wfOut.preflightOk == pOut.ok and wfOut.preflightErrors == pOut.nErrors and wfOut.preflightWarnings == pOut.nWarnings,
    "workflow preflight status summary matches quaidsPreflight");
call check(wfOut.preflightDesignInvOk == pOut.designInvOk and wfOut.preflightIVFstat == pOut.ivFstat and wfOut.preflightNClusters == pOut.nClusters,
    "workflow preflight design/IV/cluster summary matches quaidsPreflight");
call check(wfOut.symValid == qOut.symValid and wfOut.symStat == qOut.symStat and wfOut.symPval == qOut.symPval and wfOut.symDf == qOut.symDf,
    "workflow symmetry summary matches quaidsFit");
call check(wfOut.overidValid == 0 and scalmiss(wfOut.overidPvf),
    "workflow overidentification summary is missing when exactly identified");
call check(wfOut.quadraticValid == 0 and scalmiss(wfOut.quadraticPval),
    "workflow quadratic summary is missing when model is AIDS");
call check(wfOut.welfareValid == 0 and scalmiss(wfOut.cv),
    "base workflow leaves welfare scenario fields missing");
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
call check(wfOut.postRobustValid == 1 and maxc(maxc(abs(wfOut.robustBestV - vRobustFull))) == 0,
    "workflow expanded robust covariance matches quaidsRobustCovariance");
call check(maxc(abs(wfOut.sharesRobustSE - sharesRobustOut.se)) == 0,
    "workflow robust share SE match manual robust covariance propagation");
call check(maxc(abs(wfOut.incomeElasRobustSE - elasRobustOut.ser)) == 0,
    "workflow robust income elasticity SE match manual robust covariance propagation");
call check(maxc(maxc(abs(wfOut.priceElasRobustSE - elasRobustOut.sep))) == 0,
    "workflow robust price elasticity SE match manual robust covariance propagation");

pricesPt1 = pricesMean;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl, intcptMean, pricesMean, pricesPt1, totexpMean);
welfareOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, pricesPt1, totexpMean, aCtl);
welfareRobustOut = quaidsWelfareFit(bRobustFull, vRobustFull, intcptMean, pricesMean, pricesPt1, totexpMean, aCtl);

call check(wfScenario.welfareValid == 1,
    "workflow scenario computed welfare outputs");
call check(wfScenario.cv == welfareOut.cv and wfScenario.ev == welfareOut.ev,
    "workflow scenario CV/EV match quaidsWelfareFit exactly");
call check(wfScenario.seCV == welfareOut.seCV and wfScenario.seEV == welfareOut.seEV,
    "workflow scenario welfare SE match quaidsWelfareFit exactly");
call check(wfScenario.welfareRobustValid == 1 and wfScenario.seCVRobust == welfareRobustOut.seCV and wfScenario.seEVRobust == welfareRobustOut.seEV,
    "workflow scenario robust welfare SE match manual robust covariance propagation");
call check(maxc(abs(wfScenario.scenarioPrices1 - pricesPt1)) == 0 and wfScenario.scenarioTotexp0 == totexpMean,
    "workflow scenario echoes welfare scenario inputs");

aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 0;
aCtlQ.err = .0001;

wfQ = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtlQ, 0);
qOutQ = quaidsFit(w, intcpt, prices, totexp, instr, aCtlQ);
{ statQ, pvalQ, dfQ } = quaidsQuadraticTest(qOutQ);

call check(wfQ.quadraticValid == 1,
    "workflow quadratic summary is valid for unconstrained QUAIDS");
call check(wfQ.quadraticStat == statQ and wfQ.quadraticPval == pvalQ and wfQ.quadraticDf == dfQ,
    "workflow quadratic summary matches quaidsQuadraticTest exactly");

/* Milestone 26: quaidsWorkflowFit()'s own optional estimator-level weight
   argument -- distinct from surveyWeighted (Milestone 25's evaluation-
   point-only weighting, already covered by quaids_survey_workflow_test.e).
   Checks parity against explicit quaidsFit()/quaidsPreflight()/
   quaidsRobustFit() calls, the same standard every other field in this
   file is held to. */
rndseed 55;
wgtWf = 1 + rndu(tobs, 1)*3;

wfWeighted = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0, wgtWf);

qOutWfWeighted = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, wgtWf);

pOutWfWeighted = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, wgtWf);

robustWfWeighted = quaidsRobustFit(qOutWfWeighted, w, prices, totexp, aCtl, 0, wgtWf);

call check(wfWeighted.weighted == 1 and wfWeighted.weightSum == qOutWfWeighted.weightSum and wfWeighted.effN == qOutWfWeighted.effN,
    "workflow with an estimator weight echoes quaidsFit's own weight diagnostics");
call check(wfWeighted.preflightWeightValid == pOutWfWeighted.weightValid and wfWeighted.preflightWeightSum == pOutWfWeighted.weightSum
    and wfWeighted.preflightEffN == pOutWfWeighted.effN,
    "workflow preflight weight summary matches a direct quaidsPreflight call with the same weight");
call check(maxc(maxc(abs(wfWeighted.bestB - qOutWfWeighted.bestB))) == 0,
    "workflow with an estimator weight matches a direct weighted quaidsFit call exactly");
call check(maxc(maxc(abs(wfWeighted.robustSE - robustWfWeighted.se))) == 0,
    "workflow with an estimator weight matches a direct weighted quaidsRobustFit call exactly");
call check(maxc(maxc(abs(wfWeighted.bestB - wfOut.bestB))) > 1e-6,
    "a genuinely unequal estimator weight changes the workflow's point estimate (non-vacuous)");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "WORKFLOW TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "WORKFLOW TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
