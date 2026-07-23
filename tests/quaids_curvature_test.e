/*
** quaids_curvature_test.e
**
** Milestone 10: validates quaidsCurvatureFit() (src/quaidscurvature.src),
** the Diewert-Wales Cholesky reparametrization imposing local curvature
** (Slutzky negative semidefiniteness) at the sample mean, scoped to
** LA-AIDS/AIDS (aCtl.linear=1).
**
** Uses tests/quaidsfixtures.src's _quaidsCurvatureSyntheticDGP(): a
** synthetic dataset whose TRUE gamma is curvature-consistent at its own
** actual sample mean (built from a known Cholesky factor), verified
** (in that fixture's own header comment) to be a genuinely non-vacuous
** test case: the existing UNCONSTRAINED quaidsFit() recovers the true
** parameters reasonably but its own Slutzky matrix at the sample mean has
** a positive eigenvalue -- i.e. it really does violate curvature, even
** though the truth does not, so imposing curvature has something real to
** fix.
**
** Checks, in order: (1) adding-up/homogeneity/symmetry hold for the
** curvature-constrained gamma exactly, as they must by construction;
** (2) the Slutzky matrix at the reference point is negative semidefinite
** to the outer iteration's own convergence tolerance -- not exact
** floating-point zero, since (unlike the homogeneity/symmetry minimum-
** distance restrictions elsewhere in this codebase) this is an iterated
** nonlinear fit, but comfortably tighter than the violation size found
** below; (3) the unconstrained fit's own Slutzky matrix at the same
** point is genuinely positive somewhere (non-vacuous test); (4) the
** curvature-constrained gamma recovers the true gamma at least as well
** as the unconstrained fit (a real, calibrated improvement, not just "not
** obviously worse"); (5) shape/finiteness checks on the returned
** standard errors.
**
** Standard errors are NOT checked for "reasonable size": the Diewert-
** Wales Cholesky reparametrization frequently converges with some
** Cholesky-factor entries at exactly zero (a boundary solution -- the
** unconstrained optimum lies on the edge of the negative-semidefinite
** cone rather than its interior), where classical delta-method inference
** is known to be unreliable (the same boundary-inference complication
** that arises for variance components constrained to be non-negative
** elsewhere in econometrics). This run's estimated Cholesky factor has
** exactly this property (confirmed directly, not assumed) -- checked
** only for non-negativity/finiteness here, not magnitude. See
** src/quaidscurvature.src's header comment and
** GOLD_STANDARD_TODO.md's Milestone 10 section.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_curvature_test.e
*/

new;
library optmt;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidscurvature.src;
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

tobs = 3000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsCurvatureSyntheticDGP(tobs);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "starting homogeneity+symmetry fit converged");

n = qOut.n;
n1 = n - 1;
nint = qOut.nint;
gaTrue = trueParams[1+1:1+n, .];
gaUnc = qOut.bestB[1+nint+1:1+nint+n, .];

/* --- Non-vacuous test: the unconstrained fit really violates curvature
   at the sample mean, even though the true DGP does not. --- */

pricesPt = meanc(prices);
totexpPt = meanc(totexp);
intcptPt = meanc(qOut.intcptFull);
alphaUnc = qOut.bestB[1:1+nint, .];
betaUnc = qOut.bestB[1+nint+n+1, .];

a_pPtUnc = intcptPt'alphaUnc*pricesPt + .5*pricesPt'gaUnc*pricesPt;
lxPtUnc = totexpPt - a_pPtUnc;
wPtUnc = alphaUnc'intcptPt + gaUnc'pricesPt + betaUnc'lxPtUnc;
wepcUnc = -diagrv(eye(n), wPtUnc) + wPtUnc*wPtUnc' + gaUnc + (betaUnc'betaUnc)*lxPtUnc;
eigUnc = eigh(wepcUnc);

call check(maxc(eigUnc) > 0, "unconstrained fit genuinely violates curvature at the sample mean (non-vacuous test)");

/* --- Fit the curvature-constrained model. --- */

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);

call check(cOut.converged == 1, "curvature-constrained fit converged");
call check(cOut.n == n and cOut.n1 == n1, "quaidsCurvOut metadata matches qOut");


/* --- Adding-up/homogeneity/symmetry hold exactly, by construction. --- */

call check(maxc(abs(sumc(cOut.gama'))) < 1e-8, "curvature-constrained gamma: row sums == 0 (homogeneity)");
call check(maxc(abs(sumc(cOut.gama))) < 1e-8, "curvature-constrained gamma: column sums == 0 (adding-up)");
call check(maxc(maxc(abs(cOut.gama - cOut.gama'))) < 1e-8, "curvature-constrained gamma is symmetric");


/* --- Negative semidefinite at the reference point, to the outer
   iteration's own convergence tolerance. --- */

call check(maxc(cOut.eigenvalues) < 1e-3,
    "curvature-constrained fit: Slutzky matrix at the reference point is NSD (to iteration tolerance)");
call check(maxc(cOut.eigenvalues) < maxc(eigUnc),
    "curvature-constrained fit's largest eigenvalue is smaller than the unconstrained fit's");


/* --- Recovery: the curvature-constrained gamma is at least as close to
   the true (curvature-consistent) gamma as the unconstrained fit --
   empirically calibrated, not guessed (observed ~0.108 vs ~0.163). --- */

diffConstrained = maxc(maxc(abs(cOut.gama - gaTrue)));
diffUnconstrained = maxc(maxc(abs(gaUnc - gaTrue)));

call check(diffConstrained < 0.20, "curvature-constrained gamma recovers the true gamma within tolerance");
call check(diffConstrained <= diffUnconstrained,
    "curvature-constrained gamma is at least as close to the truth as the unconstrained fit");


/* --- Shape/finiteness checks on the returned coefficients/SEs. --- */

call check(rows(cOut.b) == 1+nint+n+1 and cols(cOut.b) == n, "cOut.b has the expected shape");
call check(rows(cOut.se) == rows(cOut.b) and cols(cOut.se) == cols(cOut.b), "cOut.se matches cOut.b's shape");
call check(minc(minc(cOut.se)) >= 0, "cOut.se are all non-negative");
call check(sumc(sumc(cOut.se .== cOut.se)) == rows(cOut.se)*cols(cOut.se), "cOut.se contains no NaN/missing values");

/* This run's estimated Cholesky factor genuinely has boundary (exactly
   zero) entries -- confirmed directly, not assumed, since it changes
   what can honestly be claimed about the accompanying standard errors
   (see file header). */
call check(minc(abs(vech(cOut.cholA))) < 1e-6,
    "estimated Cholesky factor has at least one boundary (near-zero) entry, as documented");


/* --- printQuaidsCurvature() runs without error. --- */
call printQuaidsCurvature(cOut);
call check(1, "printQuaidsCurvature() runs without error");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "CURVATURE TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "CURVATURE TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
