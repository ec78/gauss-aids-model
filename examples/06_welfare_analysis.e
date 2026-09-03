/*
** 06_welfare_analysis.e
**
** Exact compensating variation (CV) and equivalent variation (EV) for a
** hypothetical price change, holding nominal expenditure fixed --
** answers "how much would this household need to be compensated (or
** could be taxed) to leave them exactly as well off as before this
** price change?" Works for any model choice (LA-AIDS, iterated AIDS,
** QUAIDS), no extra package required. See docs/USAGE_GUIDE.md's
** "Welfare Analysis" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 06_welfare_analysis.e
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

n = qOut.n;
intcptPt = meanc(qOut.intcptFull);
pricesPt0 = meanc(prices);
totexpPt0 = meanc(totexp);

/* ---------------------------------------------------------------------
** Scenario: a 10% price increase on Housing, all else unchanged
** --------------------------------------------------------------------- */

pricesPt1 = pricesPt0;
pricesPt1[2] = pricesPt1[2] + ln(1.10);

wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
call printQuaidsWelfare(wOut);

print "";
print "CV/EV are positive when the price change REDUCES welfare, negative";
print "when it IMPROVES welfare. CV and EV should always agree with EACH";
print "OTHER in sign (checked below); they are NOT guaranteed to be positive";
print "just because this is a price INCREASE -- that additionally requires";
print "the fitted system to be well-behaved (Slutzky negative semidefinite)";
print "at this point, which an unconstrained/uncurved fit is not guaranteed";
print "to satisfy (see 10_curvature_imposition.e, which imposes exactly";
print "this). If CV/EV below come out negative despite a price increase,";
print "that is a real, honest consequence of this fit's own curvature";
print "properties at the reference point, not a sign error in the formula.";
print "CV sign matches EV sign:" (wOut.cv .> 0) == (wOut.ev .> 0);

/* ---------------------------------------------------------------------
** A sanity check built into the formula itself: no price change at all
** must give exactly zero CV and EV.
** --------------------------------------------------------------------- */

wOutNoChange = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt0, totexpPt0, aCtl);
print "";
print "CV/EV for pricesPt1 == pricesPt0 (should both be exactly 0):";
print wOutNoChange.cv wOutNoChange.ev;

/* ---------------------------------------------------------------------
** A larger price change, for comparison -- a 25% increase on Housing
** --------------------------------------------------------------------- */

pricesPt2 = pricesPt0;
pricesPt2[2] = pricesPt2[2] + ln(1.25);

wOutBig = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt2, totexpPt0, aCtl);
print "";
print "10% Housing price increase -- CV:" wOut.cv "  EV:" wOut.ev;
print "25% Housing price increase -- CV:" wOutBig.cv "  EV:" wOutBig.ev;
