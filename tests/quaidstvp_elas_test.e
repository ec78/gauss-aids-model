/*
** quaidstvp_elas_test.e
**
** TVP-AIDS initiative, Stage 5: validates src/quaidstvp.src's
** _quaidsTVPStateToFullB() and src/quaidstvpelas.src's
** _quaidsTVPElasFit() -- no sslib dependency (plain #include, no
** `library` statement), since a state here is just a plain vector, not
** any sslib struct type (see quaidstvpelas.src's own header).
**
** Check 1 (the primary, definitive check -- noiseless synthetic recovery
** of a KNOWN FULL n-good system, not just the n1-equation reduced system
** Stage 1's own quaidstvp_test.e already checks): builds a full n x n
** absolute-price gamma matrix directly (symmetric, row-sum-zero by
** construction via double-centering -- an INDEPENDENT construction, not
** _quaidsTVPStateToFullB()'s own recovery formula) plus full true
** alpha/beta vectors (equation n's own values set by the adding-up
** identities directly, in the test, not by calling any library proc).
** The n1-equation REDUCED system used to generate data is then only a
** trivial submatrix/subvector extraction from these true full values, so
** an exact match between the recovered bFull and the independently-built
** true full system is a real correctness claim about
** _quaidsTVPStateToFullB()'s own recovery logic, not a tautology.
**
** Check 2 (regression guard): homogeneity (gamma row sums), symmetry,
** and adding-up (alpha/beta/gamma column sums) all hold exactly on the
** RECOVERED bFull -- confirms adding-up-on-gamma really does fall out
** automatically from homogeneity+symmetry alone, as
** _quaidsTVPStateToFullB()'s own header claims, rather than merely
** trusting that derivation.
**
** Check 3 (internal consistency, exact -- the second independent check
** required for new estimation logic per CLAUDE.md's Testing
** expectations): _quaidsTVPElasFit()'s own output must exactly match a
** direct call to _quaidsElas() (the already-correct, already-tested
** sibling proc it wraps) fed the same recovered bFull -- isolates "does
** the wrapper plumb its inputs through correctly" from "is the recovery
** math correct" (Check 1's own concern).
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvp_elas_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstvp.src;
#include ../src/quaidstvpelas.src;
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


/* ==========================================================================
   Setup: n1 = 3 (n = 4 goods). Build the TRUE FULL system independently of
   any library recovery logic.
   ========================================================================== */

n1 = 3;
n = n1 + 1;
tobs = 60;
rndseed 4242;

/* True full n x n absolute-price gamma: symmetric and row-sum-zero by
   double-centering a random symmetric matrix -- for symmetric M, row i's
   mean equals column i's mean, so centering by both row and column means
   (and adding back the grand mean) preserves symmetry while forcing every
   row (and, by symmetry, every column) to sum to exactly zero. This is an
   INDEPENDENT construction of a valid true system, not
   _quaidsTVPStateToFullB()'s own row-sum recovery formula. */
rawMat = .1*round(rndn(n, n)*10)/10;
symMat = (rawMat + rawMat')/2;
rowMeanVec = sumc(symMat')/n;
gMeanVal = meanc(rowMeanVec);
trueGammaAbs = symMat - rowMeanVec*ones(1, n) - ones(n, 1)*rowMeanVec' + gMeanVal;

call check(maxc(abs(sumc(trueGammaAbs'))) < 1e-10, "constructed trueGammaAbs has exact zero row sums (homogeneity, by construction)");
call check(maxc(maxc(abs(trueGammaAbs - trueGammaAbs'))) < 1e-10, "constructed trueGammaAbs is exactly symmetric (by construction)");

trueAlpha = .1*round(rndn(n1, 1)*10)/10 + .15;
trueAlphaN = 1 - sumc(trueAlpha);
trueAlphaFull = trueAlpha|trueAlphaN;

trueBeta = .05*round(rndn(n1, 1)*10)/10;
trueBetaN = -sumc(trueBeta);
trueBetaFull = trueBeta|trueBetaN;

/* The n1-equation REDUCED system is just a trivial submatrix/subvector
   extraction from the true full system above -- see this proc's own
   derivation in _quaidsTVPStateToFullB()'s header for why the relative-
   price coefficient on (p_j - p_n) in equation i equals gamma_abs[i,j]
   exactly whenever homogeneity holds. */
trueGammaRel = trueGammaAbs[1:n1, 1:n1];

trueBFull = (trueAlphaFull')|trueGammaAbs|(trueBetaFull');


/* ==========================================================================
   Generate noiseless n1-equation data from the reduced system, exactly as
   tests/quaidstvp_test.e's own Check 1 does, then recover stateHat via
   plain pooled OLS across the shared per-period design.
   ========================================================================== */

pricesRel = .1*rndn(tobs, n1);
lx = .2*rndn(tobs, 1);

w = zeros(tobs, n1);
t = 1;
do while t <= tobs;
    i = 1;
    do while i <= n1;
        w[t, i] = trueAlpha[i] + pricesRel[t, .]*trueGammaRel[i, .]' + trueBeta[i]*lx[t];
        i = i + 1;
    endo;
    t = t + 1;
endo;

Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

Xstack = zeros(tobs*n1, k_states);
ystack = zeros(tobs*n1, 1);
row = 1;
t = 1;
do while t <= tobs;
    Zt = getmatrix(Zarr, t);
    i = 1;
    do while i <= n1;
        Xstack[row, .] = Zt[i, .];
        ystack[row] = w[t, i];
        row = row + 1;
        i = i + 1;
    endo;
    t = t + 1;
endo;

stateHat = invpd(Xstack'Xstack)*Xstack'ystack;


/* ==========================================================================
   Check 1: noiseless recovery of the FULL n-good system, against the
   independently-constructed true full system (see Setup above).
   ========================================================================== */

bFullHat = _quaidsTVPStateToFullB(stateHat, n1);

call check(rows(bFullHat) == n+2 and cols(bFullHat) == n, "bFullHat shape is (n+2) x n");
call check(maxc(maxc(abs(bFullHat - trueBFull))) < 1e-8,
    "noiseless recovery: _quaidsTVPStateToFullB's bFull matches the independently-constructed true full system exactly");


/* ==========================================================================
   Check 2: homogeneity/symmetry/adding-up hold exactly on the RECOVERED
   bFull -- confirms adding-up-on-gamma falls out of homogeneity+symmetry
   automatically, rather than merely trusting the derivation.
   ========================================================================== */

alphaHat = bFullHat[1, .];
gammaHat = bFullHat[2:n+1, .];
betaHat = bFullHat[n+2, .];

call check(maxc(maxc(abs(gammaHat - gammaHat'))) < 1e-8, "recovered gamma is exactly symmetric");
call check(maxc(abs(sumc(gammaHat'))) < 1e-8, "recovered gamma has exact zero row sums (homogeneity)");
call check(maxc(abs(sumc(gammaHat))) < 1e-8, "recovered gamma has exact zero column sums (adding-up, NOT separately imposed -- falls out of homogeneity+symmetry)");
call check(abs(sumc(alphaHat') - 1) < 1e-8, "recovered alpha sums to exactly 1 (adding-up)");
call check(abs(sumc(betaHat')) < 1e-8, "recovered beta sums to exactly 0 (adding-up)");


/* ==========================================================================
   Check 3: internal consistency, exact -- _quaidsTVPElasFit()'s own
   output must exactly match a direct call to _quaidsElas() (the
   already-correct sibling proc it wraps) fed the same recovered bFull.
   ========================================================================== */

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;

pricesEval = .05*ones(n, 1);
totexpEval = .1;

{ erWrap, epWrap, epcWrap } = _quaidsTVPElasFit(stateHat, n1, pricesEval, totexpEval, aCtl);
{ erDirect, epDirect, epcDirect } = _quaidsElas(bFullHat, 1, pricesEval, totexpEval, aCtl);

call check(maxc(abs(erWrap - erDirect)) < 1e-12, "_quaidsTVPElasFit's income elasticities exactly match a direct _quaidsElas() call on the same recovered bFull");
call check(maxc(maxc(abs(epWrap - epDirect))) < 1e-12, "_quaidsTVPElasFit's uncompensated price elasticities exactly match a direct _quaidsElas() call");
call check(maxc(maxc(abs(epcWrap - epcDirect))) < 1e-12, "_quaidsTVPElasFit's compensated price elasticities exactly match a direct _quaidsElas() call");

/* Loose sanity check: income elasticities should be finite, non-degenerate
   numbers (not NaN/Inf from a shape mismatch slipping through silently). */
call check(sumc(erWrap .== erWrap) == n, "income elasticities are all finite (no NaN)");


/* ==========================================================================
   Summary
   ========================================================================== */

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "QUAIDS TVP STAGE 5 TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "QUAIDS TVP STAGE 5 TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
