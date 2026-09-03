/*
** 05_elasticities_shares_slutzky.e
**
** Post-estimation diagnostics evaluated at an arbitrary point: income and
** price elasticities (quaidsElasFit), the model-implied predicted budget
** share (quaidsSharesFit), and the Slutzky negative-semidefiniteness
** diagnostic (quaidsSlutzky) -- first at the sample mean, then at a
** synthetic counterfactual price scenario. See docs/USAGE_GUIDE.md's
** "Elasticities At Any Point" and "Predicted Budget Shares At Any Point"
** sections.
**
** Run from the examples/ directory:
**   tgauss -b -x 05_elasticities_shares_slutzky.e
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

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

if not qOut.converged;
    print "quaidsFit() did not converge.";
    end;
endif;

/* Always evaluate against qOut.bestB/qOut.bestV -- "whichever is the
   most-constrained estimate actually fit" (symmetry-constrained here,
   since aCtl.homogenous = 1). */
n = qOut.n;
nint = qOut.nint;

/* ---------------------------------------------------------------------
** 1. At the sample mean
** --------------------------------------------------------------------- */

m_ = meanc(qOut.intcptFull~prices~totexp);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

print "=== At the sample mean ===";

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call printQuaidsElas(elasOut);

/* Predicted shares below can fall outside [0,1] -- a known property of
   this synthetic dataset's noise/price scale, not a bug; see
   example_data.src's header comment and 01_basic_estimation.e. Adding-up
   (they sum to exactly 1) is what's structurally guaranteed. */
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call printQuaidsShares(sharesOut);
print "Sum of predicted shares (adding-up, exact by construction):" sumc(sharesOut.w);

call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);

/* ---------------------------------------------------------------------
** 2. At a counterfactual: a 20% price increase on Food, all else at the
** sample mean
**
** quaidsElasFit()/quaidsSharesFit() take an explicit evaluation point,
** not a fixed set of sample statistics -- this is what makes policy
** questions like "what if Food prices rose 20%?" directly answerable.
** --------------------------------------------------------------------- */

pricesCf = pricesMean;
pricesCf[1] = pricesCf[1] + ln(1.20);

print "";
print "=== At a hypothetical 20% Food price increase (mean otherwise) ===";

elasOutCf = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesCf, totexpMean, aCtl);
call printQuaidsElas(elasOutCf);

sharesOutCf = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesCf, totexpMean, aCtl);
call printQuaidsShares(sharesOutCf);

print "";
print "Change in Food's own predicted share:" sharesOutCf.w[1] - sharesOut.w[1];
print "Food's own-price uncompensated elasticity at the mean:" elasOut.ep[1,1];
print "(a more negative elasticity here would predict a bigger share drop above)";
