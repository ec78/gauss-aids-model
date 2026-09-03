/*
** 04_hypothesis_tests.e
**
** Standalone Wald tests of the demand-theory restrictions this library
** can otherwise impose by construction: homogeneity, symmetry (jointly
** with homogeneity), and whether the QUAIDS quadratic log-expenditure
** term is actually needed. All three require an UNCONSTRAINED fit
** (aCtl.homogenous = 0) -- quaidsFit() itself already reports a
** symmetry-given-homogeneity test as part of a constrained fit; these
** are for testing the restrictions themselves, not assuming them. See
** docs/USAGE_GUIDE.md's "Homogeneity, Symmetry, and Overidentification"
** and "Is QUAIDS Needed?" sections.
**
** Run from the examples/ directory:
**   tgauss -b -x 04_hypothesis_tests.e
*/

new;
library quaids;
#include example_data.src

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);

aCtl = quaidsControlCreate();
aCtl.linear = 0;          // QUAIDS -- required for quaidsQuadraticTest below
aCtl.maxiter = 100;
aCtl.homogenous = 0;      // UNCONSTRAINED -- required by all three tests
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

if not qOut.converged;
    print "quaidsFit() did not converge; hypothesis tests need a converged fit.";
    end;
endif;

/* ---------------------------------------------------------------------
** 1. Homogeneity: does each equation's own price coefficients sum to 0?
** --------------------------------------------------------------------- */

{ statH, pvalH, dfH } = quaidsHomogeneityTest(qOut);
print "Homogeneity test:  chi2(" dfH ") =" statH "  p-value =" pvalH;
if pvalH < .05;
    print "  -> rejected at the 5% level.";
else;
    print "  -> not rejected at the 5% level.";
endif;

/* ---------------------------------------------------------------------
** 2. Joint homogeneity + symmetry
** --------------------------------------------------------------------- */

{ statJ, pvalJ, dfJ } = quaidsJointTest(qOut);
print "";
print "Joint homogeneity+symmetry test:  chi2(" dfJ ") =" statJ "  p-value =" pvalJ;
if pvalJ < .05;
    print "  -> rejected at the 5% level.";
else;
    print "  -> not rejected at the 5% level.";
endif;

/* ---------------------------------------------------------------------
** 3. Is the QUAIDS quadratic term actually needed?
**
** The synthetic data here is genuinely QUAIDS (a real, nonzero true
** lambda), so this test SHOULD reject -- a real, non-vacuous check that
** the test has power, not just that it runs.
** --------------------------------------------------------------------- */

{ statQ, pvalQ, dfQ } = quaidsQuadraticTest(qOut);
print "";
print "Quadratic-term test:  chi2(" dfQ ") =" statQ "  p-value =" pvalQ;
if pvalQ < .05;
    print "  -> rejected: the quadratic term is needed (expected here -- the";
    print "     synthetic data was generated with a real quadratic term).";
else;
    print "  -> not rejected: plain AIDS (aCtl.linear = 1) would fit as well,";
    print "     and is simpler and more numerically stable -- see";
    print "     docs/USAGE_GUIDE.md's 'Choosing A Model' section.";
endif;

print "";
print "Note: quaidsHomogeneityTest()/quaidsJointTest()/quaidsQuadraticTest()";
print "all error clearly if called on a qOut with aCtl.homogenous == 1, and";
print "quaidsQuadraticTest() additionally requires qOut.linear == 0 (an AIDS";
print "fit never estimates a quadratic term, so there is nothing to test).";
