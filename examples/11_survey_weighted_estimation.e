/*
** 11_survey_weighted_estimation.e
**
** quaidsFit() accepts an optional sampling-weight argument for a
** genuine weighted point estimate -- useful when your sample was drawn
** with unequal selection probabilities (e.g. a household survey that
** oversamples some strata). This example uses data drawn with a real,
** deliberate sampling bias (see example_data.src's own header) and
** compares naive (unweighted) estimation against weighted estimation on
** the SAME biased sample, against the TRUE population parameters --
** the whole point of the weight argument. See docs/USAGE_GUIDE.md's
** "Choosing An API" (quaidsSurveyWorkflowFit) and "Limitations"
** sections.
**
** Run from the examples/ directory:
**   tgauss -b -x 11_survey_weighted_estimation.e
*/

new;
library quaids;
#include example_data.src

{ w, intcpt, prices, totexp, instr, sampwt } = quaidsExampleSurveyData(6000, 11);

print "Sampled observations:" rows(w);
print "Sampling weight range:" minc(sampwt) "to" maxc(sampwt);
print "(a weight of, say, 4 means that row represents about 4 population";
print "households -- rows from the under-sampled stratum get bigger weights.)";

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

/* ---------------------------------------------------------------------
** 1. quaidsFit()'s own weight= argument -- the estimator-level primitive
** --------------------------------------------------------------------- */

qOutNaive = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
qOutWeighted = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, weight=sampwt);

print "";
print "Naive (unweighted) fit converged:" qOutNaive.converged;
print "Weighted fit converged:" qOutWeighted.converged;
print "Weighted fit's effective sample size (Kish):" qOutWeighted.effN "of" rows(w) "rows";

print "";
print "Max abs difference between weighted and naive gamma coefficients:";
print maxc(maxc(abs(qOutWeighted.bestB - qOutNaive.bestB)));
print "(a real, non-trivial difference -- this sample is genuinely biased,";
print "not just re-weighted for no reason; see example_data.src's header";
print "for how and why.)";

/* ---------------------------------------------------------------------
** 2. quaidsSurveyWorkflowFit() -- weight the estimator AND evaluate
** post-estimation quantities at a sampling-weighted representative point
**
** Since Milestone 26, the same weight argument does double duty: it
** fits quaidsFit() as a genuine weighted estimator (like step 1 above)
** AND changes the evaluation point used for predicted shares/
** elasticities to the weighted mean, with SE recomputed there too.
** --------------------------------------------------------------------- */

wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, sampwt);

print "";
print "=== quaidsSurveyWorkflowFit() ===";
print "weighted:" wfSurvey.weighted "  weightSum:" wfSurvey.weightSum "  effN:" wfSurvey.effN;
print "Weighted evaluation point (total expenditure):" wfSurvey.evalTotexp;
print "Predicted shares at the weighted representative point:";
print wfSurvey.shares';
print "Income elasticities at the weighted representative point:";
print wfSurvey.incomeElas';
