/*
** 12_applied_workflow.e
**
** quaidsWorkflowFit() bundles the individual pieces demonstrated across
** this suite (preflight, quaidsFit, predicted shares/elasticities at
** the sample mean, robust standard errors) into one silent, struct-
** returning call -- the fastest path from raw data to a usable applied
** result. quaidsWorkflowScenarioFit() adds an explicit CV/EV price-
** change scenario on top. See docs/USAGE_GUIDE.md's "Choosing An API"
** section.
**
** Run from the examples/ directory:
**   tgauss -b -x 12_applied_workflow.e
*/

new;
library quaids;
#include example_data.src

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);
goodNames = quaidsExampleGoodNames();

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

/* ---------------------------------------------------------------------
** One call: preflight summary, fit, mean-point shares/elasticities, and
** heteroskedasticity-robust standard errors.
** --------------------------------------------------------------------- */

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);

print "Model:" wfOut.model;
print "Converged:" wfOut.converged "in" wfOut.iterations "iterations";
print "";
print "Preflight (echoed inside the workflow struct):";
print "  ok:" wfOut.preflightOk "  warnings:" wfOut.preflightWarnings;
print "  weak instrument flag:" wfOut.preflightWeakIV;
print "";
print "Model/restriction summary:";
print "  symmetry rejected at 5%?" wfOut.symPval < .05;
print "  overidentification test valid (ninst > nu)?" wfOut.overidValid;
print "  (this example uses exactly one instrument, so it is exactly";
print "   identified -- the overidentification test needs a genuinely";
print "   extra instrument, and is not applicable here.)";
print "  quadratic-term test valid (needs an unconstrained QUAIDS fit)?" wfOut.quadraticValid;
print "  (not applicable here -- this workflow fit imposes homogeneity/";
print "   symmetry; see 04_hypothesis_tests.e for the unconstrained case.)";

print "";
print "Predicted shares at the sample mean:";
print goodNames';
print wfOut.shares';
print "Income elasticities at the sample mean:";
print wfOut.incomeElas';
print "";
print "Robust (heteroskedasticity-robust) coefficient SE mean:" meanc(vec(wfOut.robustSE));
print "Robust-propagated predicted-share SE:" wfOut.sharesRobustSE';

/* ---------------------------------------------------------------------
** Add an explicit CV/EV scenario: a 10% price increase on Housing,
** evaluated at the same sample-mean point the base workflow used.
** --------------------------------------------------------------------- */

pricesPt1 = wfOut.evalPrices;
pricesPt1[2] = pricesPt1[2] + ln(1.10);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);

print "";
print "=== 10% Housing price increase scenario ===";
if wfScenario.welfareValid;
    print "CV:" wfScenario.cv "  se:" wfScenario.seCV;
    print "EV:" wfScenario.ev "  se:" wfScenario.seEV;
endif;
if wfScenario.welfareRobustValid;
    print "Robust se -- CV:" wfScenario.seCVRobust "  EV:" wfScenario.seEVRobust;
endif;

print "";
print "Next: 13_pubtable_reporting.e -- exporting these results to";
print "LaTeX/Markdown/CSV tables.";
