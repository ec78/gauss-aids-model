/*
** 10_curvature_imposition.e
**
** quaidsSlutzky() always DIAGNOSES curvature (Slutzky negative
** semidefiniteness, the demand-theory condition for a well-behaved
** consumer) but never imposes it. quaidsCurvatureFit() can impose it
** locally, at the sample mean, via the Diewert-Wales (1987) Cholesky
** reparametrization -- for LA-AIDS/AIDS and (since Milestone 13) QUAIDS.
** Requires the optmt package. See docs/USAGE_GUIDE.md's "Imposing
** Curvature (Diewert-Wales)" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 10_curvature_imposition.e
*/

new;
library optmt, quaids;
#include example_data.src

/* ---------------------------------------------------------------------
** 1. QUAIDS, on the standard example dataset
**
** quaidsCurvatureFit() needs an already-converged, homogeneity+symmetry-
** constrained starting fit (aCtl.homogenous = 1). QUAIDS's curvature
** outer loop is measurably less stable than AIDS's own -- aCtl.relax =
** .25 is effectively required here (undamped runs on this dataset
** diverge) -- and it takes noticeably longer to converge (dozens to
** ~200 outer iterations, vs. AIDS's ~10-20).
** --------------------------------------------------------------------- */

{ wQ, intcptQ, pricesQ, totexpQ, instrQ } = quaidsExampleData(3000, 204);

aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 1;
aCtlQ.err = .0001;

qOutQ = quaidsFit(wQ, intcptQ, pricesQ, totexpQ, instrQ, aCtlQ);
print "QUAIDS fit converged:" qOutQ.converged;

aCtlQ.relax = .25;
cOutQ = quaidsCurvatureFit(qOutQ, wQ, pricesQ, totexpQ, aCtlQ);

print "";
print "=== QUAIDS, after imposing curvature ===";
call printQuaidsCurvature(cOutQ);
print "Curvature-constrained eigenvalues at the mean (all should be <= 0):" cOutQ.eigenvalues';

/* ---------------------------------------------------------------------
** 2. AIDS, on a dataset built so the TRUE gamma is already curvature-
** consistent -- lets us show a clean "before: violates curvature,
** after: fixed" comparison. (quaidsExampleData() above works fine for
** QUAIDS curvature imposition, as just shown, but its own much larger,
** uncurved Slutzky violations -- needed for that dataset's general-
** purpose role throughout the rest of this suite -- make AIDS curvature
** imposition on it numerically unreliable; confirmed directly, not
** assumed. quaidsExampleCurvatureData() below is a separate, smaller-
** scale AIDS-only dataset built specifically for this.)
** --------------------------------------------------------------------- */

{ w, intcpt, prices, totexp, instr } = quaidsExampleCurvatureData(3000);

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
print "";
print "AIDS fit converged:" qOut.converged;

intcptPt = meanc(qOut.intcptFull)';
pricesPt = meanc(prices)';
totexpPt = meanc(totexp)';

print "";
print "=== AIDS Slutzky eigenvalues at the mean, BEFORE imposing curvature ===";
call quaidsSlutzky(qOut.bestB, intcptPt, pricesPt, totexpPt, aCtl);
print "(any POSITIVE eigenvalue above is a real violation of this fit's own";
print "curvature -- even though this dataset's TRUE gamma is curvature-";
print "consistent by construction, the ESTIMATED fit is not guaranteed to be.)";

cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
print "";
print "=== AIDS, after imposing curvature ===";
call printQuaidsCurvature(cOut);
print "Curvature-constrained eigenvalues at the mean (all should be <= 0):" cOut.eigenvalues';

/* A real, documented caveat: the estimated Cholesky factor often has
   entries at exactly zero -- the constrained optimum sits on the EDGE
   of the negative-semidefinite cone, where classical delta-method
   inference is known to be unreliable. Point estimates and the exact
   curvature property are unaffected. See docs/USAGE_GUIDE.md#limitations. */

/* ---------------------------------------------------------------------
** 3. A bootstrap alternative to the delta-method standard error
**
** A small B keeps this example fast -- a single AIDS curvature fit
** takes under a second, but conventional replication counts (200-1000)
** would take minutes; see docs/USAGE_GUIDE.md#limitations before
** choosing B for real work. There is deliberately no default.
** --------------------------------------------------------------------- */

print "";
print "=== AIDS curvature bootstrap (B=15) ===";
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 15, seed=42);
call printQuaidsCurvatureBootstrap(bootOut);

{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, 0.05);
print "";
print "95% percentile CI for good 1's first coefficient:" ciLower[1,1] "to" ciUpper[1,1];
print "(point estimate:" bootOut.b[1,1] ")";
