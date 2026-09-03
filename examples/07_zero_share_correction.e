/*
** 07_zero_share_correction.e
**
** Real survey/microdata routinely has corner solutions -- some
** households report zero expenditure on some goods. Fitting quaidsFit()
** directly on such data is a known source of bias, since the linear/
** log-linear share equation has no mechanism for a censored dependent
** variable. quaidsZeroFit() corrects for this via the Shonkwiler & Yen
** (1999) two-step procedure: a per-good first-stage probit (the
** probability of a non-zero share) followed by a corrected second-stage
** GLS fit. See docs/USAGE_GUIDE.md's "Zero Budget Shares (Corner
** Solutions)" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 07_zero_share_correction.e
*/

new;
library quaids;
#include example_data.src

{ w, intcpt, prices, totexp, instr } = quaidsExampleZeroData(3000, 1);
goodNames = quaidsExampleGoodNames();

print "Budget categories:" goodNames';
print "Fraction of zero shares per good:" (meanc(w .== 0))';

/* ---------------------------------------------------------------------
** 1. Unconstrained (aCtl.homogenous = 0)
** --------------------------------------------------------------------- */

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 0;
aCtl.err = .0001;

print "";
print "=== Unconstrained Shonkwiler-Yen correction ===";
zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);

print "";
print "All 5 first-stage probits converged:" prodc(zOut.probitConverged) == 1;
print "Observed shareZeroFrac matches the raw fraction printed above:" zOut.shareZeroFrac';

/* ---------------------------------------------------------------------
** 2. Homogeneity + symmetry imposed on top of the correction
**
** zOut.bS's recovered n x n gamma block is exactly symmetric by
** construction of the combined symmetry+diagonal-delta restriction --
** see src/quaidszerocorrect.src's own header for the derivation.
** --------------------------------------------------------------------- */

aCtl.homogenous = 1;

print "";
print "=== Homogeneity+symmetry-constrained Shonkwiler-Yen correction ===";
zOutH = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOutH);

nint = zOutH.nint;
n = zOutH.n;
gammaBS = zOutH.bS[1+nint+1:1+nint+n, .];
print "";
print "Recovered gamma is exactly symmetric:" maxc(maxc(abs(gammaBS - gammaBS'))) < 1e-6;

/* ---------------------------------------------------------------------
** 3. A real, documented limitation: adding-up is only approximate here
**
** Each equation is independently rescaled by its own good-specific
** first-stage probability, so the corrected shares do not sum to
** exactly 1 the way an uncorrected quaidsFit() always does -- this is a
** genuine property of the Shonkwiler-Yen method itself, not a bug.
** --------------------------------------------------------------------- */

print "";
print "Note: unlike quaidsFit(), the Shonkwiler-Yen correction does not";
print "make adding-up hold exactly for its own corrected coefficients --";
print "a documented property of the method, not a defect -- see";
print "docs/USAGE_GUIDE.md's 'Zero Budget Shares' section.";
