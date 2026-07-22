/*
** quaids_published_validation_test.e
**
** Milestone 3 published-data validation: fits quaidsFit() (LA-AIDS,
** aCtl.linear=1, aCtl.maxiter=1, Stone price index) on the Blanciforti86
** food-demand dataset (see tests/fixtures/published/SOURCE.md) and checks
** the recovered coefficients against an independent reference computed in
** R via micEconAids::aidsEst(..., instNames=...) (3SLS, instrumenting
** log(xFood) with log(xAgg) -- the same identification strategy GAUSS
** uses, via a different estimation algorithm: GAUSS uses a control-
** function/residual-inclusion approach, R's 3SLS projects regressors onto
** the instrument space directly. Both are valid, asymptotically
** consistent IV estimators for the same model; they are not expected to
** match bit-for-bit, only closely.
**
** This test is also what caught a real, now-fixed bug: quaidsFit()'s
** Stone-index starting-value computation previously used a mutated
** (partly-relative, partly-absolute) price matrix, producing a materially
** wrong deflator for the aCtl.maxiter=1 (LA-AIDS) case specifically. See
** GOLD_STANDARD_TODO.md Milestone 3 for the full writeup.
**
** Milestone 9 addition: extends the same published-data cross-check to
** iterated (nonlinear translog price index) AIDS -- aCtl.linear=1,
** aCtl.maxiter>1 -- against R's aidsEst(method="IL", ...) (the Iterated
** Linear Least Squares Estimator, Blundell & Robin 1999), closing a gap
** the Milestone 9 gate found: only the LA-AIDS special case had
** independent published-data validation. micEconAids's IL method does not
** support instrumental variables (instNames is silently ignored for
** method="LA"/3SLS only -- confirmed empirically: combining method="IL"
** with instNames segfaults R's aidsEst rather than erroring cleanly, so
** this is not exercised or shipped as a working example anywhere), while
** GAUSS's iterated fit still always instruments log(xFood) -- so this
** comparison spans both a different estimation algorithm AND an IV-vs-
** no-IV difference, and is expected to (and does) show a wider gap than
** the tighter LA-AIDS-vs-3SLS-IV comparison above. QUAIDS (the quadratic
** extension) has no such reference available: micEconAids does not
** implement a quadratic log-expenditure term at all. QUAIDS's validation
** remains the known-true synthetic-DGP recovery in
** tests/quaids_synthetic_validation_test.e -- documented as a real but
** different tier of evidence (recovery of independently-generated
** known-true parameters, not cross-implementation agreement on real
** data) in docs/FEATURE_SUPPORT_MATRIX.md.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_published_validation_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsformula.src;

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

/* Tolerance: an 0.05 absolute-difference bound. The observed max
   difference against the R reference below is ~0.021 (alpha, good 4);
   0.05 gives real headroom above that without being so loose it would
   fail to catch a genuine regression. This is a cross-implementation
   agreement bound between two different (but both valid) IV algorithms,
   not a bit-identical-output bound -- see the file header. */
tol = 0.05;

data = loadd("fixtures/published/blanciforti86_food32.csv");
call check(rows(data) == 32, "loaded 32 years of Blanciforti86 food data");

w = data[., "wFood1" "wFood2" "wFood3" "wFood4"];
prices = ln(data[., "pFood1" "pFood2" "pFood3" "pFood4"]);
totexp = ln(data[., "xFood"]);
instr = ln(data[., "xAgg"]);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 1;
aCtl.maxiter = 1;
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, 0, prices, totexp, instr, aCtl);

call check(qOut.model $== "LA-AIDS", "model == LA-AIDS (aCtl.linear=1, aCtl.maxiter=1)");
call check(rows(qOut.bS) == 7 and cols(qOut.bS) == 4, "bS is 7 x 4 (alpha, 4 gamma rows, beta, u-coef)");

alphaGauss = qOut.bS[1, .];
gammaGauss = qOut.bS[2:5, .];
betaGauss = qOut.bS[6, .];

/* Adding-up / homogeneity / symmetry sanity checks, independent of the
   external reference -- these must hold by construction of the estimator
   regardless of any published-data comparison. */
call check(abs(sumc(alphaGauss') - 1) < 1e-8, "sum(alpha) == 1 (adding-up)");
call check(maxc(abs(sumc(gammaGauss'))) < 1e-8, "gamma row sums == 0 (homogeneity)");
call check(maxc(abs(sumc(gammaGauss))) < 1e-8, "gamma col sums == 0 (adding-up across equations)");
call check(abs(sumc(betaGauss')) < 1e-8, "sum(beta) == 0 (adding-up)");
call check(maxc(maxc(abs(gammaGauss - gammaGauss'))) < 1e-8, "gamma is symmetric");

/* --- Reference: R micEconAids::aidsEst(priceNames, shareNames, "xFood", --- */
/* --- data=Blanciforti86[1:32,], priceIndex="Ls", instNames=c(log      --- */
/* --- prices, "lxAgg")), estMethod dispatches to 3SLS since instNames  --- */
/* --- is supplied. Computed 2026-07-20, R 4.5.0, micEconAids 0.6-20.   --- */
let alphaR[1, 4] = -0.329388137 0.032896406 0.292524909 1.003966821;
let betaR[1, 4] = 0.373054464 0.100376128 -0.093012355 -0.380418237;
let gammaR[4, 4] =
     0.095575402 -0.140321326 -0.012528920  0.057274844
    -0.140321326  0.140292705  0.008121780 -0.008093159
    -0.012528920  0.008121780  0.013275950 -0.008868810
     0.057274844 -0.008093159 -0.008868810 -0.040312876;

call check(maxc(maxc(abs(alphaGauss - alphaR))) <= tol, "alpha matches R (3SLS, log(xAgg) instrument) within tolerance");
call check(maxc(maxc(abs(betaGauss - betaR))) <= tol, "beta matches R (3SLS, log(xAgg) instrument) within tolerance");
call check(maxc(maxc(abs(gammaGauss - gammaR))) <= tol, "gamma matches R (3SLS, log(xAgg) instrument) within tolerance");

print;
print "alpha (GAUSS):" alphaGauss;
print "alpha (R)    :" alphaR;
print "beta  (GAUSS):" betaGauss;
print "beta  (R)    :" betaR;
print "max abs diff, alpha/beta/gamma:";
print maxc(maxc(abs(alphaGauss-alphaR)))~maxc(maxc(abs(betaGauss-betaR)))~maxc(maxc(abs(gammaGauss-gammaR)));


/* --- Iterated AIDS (aCtl.linear=1, aCtl.maxiter>1, translog price index) --- */

/* Tolerance is wider than the LA-AIDS check above (0.05): this comparison
   spans both a different estimation algorithm AND an IV-vs-no-IV
   difference (see file header). Observed max abs difference is ~0.112
   (alpha); 0.15 gives real headroom above that without being so loose it
   would fail to catch a genuine regression. */
tolIL = 0.15;

struct quaidsControl aCtlIL;
aCtlIL = quaidsControlCreate;
aCtlIL.linear = 1;
aCtlIL.maxiter = 100;
aCtlIL.homogenous = 1;
aCtlIL.err = .0001;

struct quaidsOut qOutIL;
qOutIL = quaidsFit(w, 0, prices, totexp, instr, aCtlIL);

call check(qOutIL.model $== "AIDS", "model == AIDS (aCtl.linear=1, aCtl.maxiter>1)");
call check(qOutIL.converged == 1, "iterated AIDS fit converged");

alphaGaussIL = qOutIL.bS[1, .];
gammaGaussIL = qOutIL.bS[2:5, .];
betaGaussIL = qOutIL.bS[6, .];

call check(abs(sumc(alphaGaussIL') - 1) < 1e-6, "iterated AIDS: sum(alpha) == 1 (adding-up)");
call check(maxc(abs(sumc(gammaGaussIL'))) < 1e-6, "iterated AIDS: gamma row sums == 0 (homogeneity)");
call check(maxc(maxc(abs(gammaGaussIL - gammaGaussIL'))) < 1e-6, "iterated AIDS: gamma is symmetric");

/* --- Reference: R micEconAids::aidsEst(priceNames, shareNames, "xFood", --- */
/* --- data=Blanciforti86[1:32,], method="IL", priceIndex="Ls") -- SUR,   --- */
/* --- no instrument (see file header for why). Computed 2026-07-21,     --- */
/* --- R 4.5.0, micEconAids 0.6-20, 7 IL iterations to convergence.      --- */
let alphaRIL[1, 4] = -0.2591720 0.1244717 0.2726351 0.8620652;
let betaRIL[1, 4] = 0.33056489 0.04694985 -0.08110712 -0.29640762;
let gammaRIL[4, 4] =
    -0.08546936 -0.169939590  0.034689651  0.220719303
    -0.16993959  0.156128397  0.011038062  0.002773131
     0.03468965  0.011038062  0.001339222 -0.047066935
     0.22071930  0.002773131 -0.047066935 -0.176425499;

call check(maxc(maxc(abs(alphaGaussIL - alphaRIL))) <= tolIL, "iterated AIDS: alpha matches R IL within tolerance");
call check(maxc(maxc(abs(betaGaussIL - betaRIL))) <= tolIL, "iterated AIDS: beta matches R IL within tolerance");
call check(maxc(maxc(abs(gammaGaussIL - gammaRIL))) <= tolIL, "iterated AIDS: gamma matches R IL within tolerance");

print;
print "alpha (GAUSS, iterated AIDS):" alphaGaussIL;
print "alpha (R, method=IL)        :" alphaRIL;
print "max abs diff, alpha/beta/gamma (iterated AIDS vs. R IL):";
print maxc(maxc(abs(alphaGaussIL-alphaRIL)))~maxc(maxc(abs(betaGaussIL-betaRIL)))~maxc(maxc(abs(gammaGaussIL-gammaRIL)));

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "PUBLISHED VALIDATION TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "PUBLISHED VALIDATION TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
