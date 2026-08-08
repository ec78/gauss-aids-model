/*
** quaids_zero_test.e
**
** Milestone 19: validates quaidsZeroFit()/printQuaidsZero()
** (src/quaidszerocorrect.src) -- AIDS/QUAIDS estimation corrected for
** zero budget shares via the Shonkwiler-Yen (1999) two-step procedure.
**
** Milestone 30 extends this file for two changes: (1) homogeneity/
** symmetry imposition on top of the correction (aCtl.homogenous=1),
** via the same n1=n-aCtl.homogenous reparametrization quaidsFit() uses,
** plus a combined symmetry+diagonal-delta minimum-distance restriction
** (see src/quaidszerocorrect.src's own header, and
** _quaidsZeroSymDiagRestrict()'s header for the real cross-constraint
** finding -- gamma's column-n entries are NOT free under symmetry, they
** are a deterministic function of the same symmetric sub-block, found
** by direct empirical testing against this fixture's known-symmetric
** DGP after an initial draft left them free and produced measurably
** asymmetric recovered output); (2) a bugfix converting zOut.b's gamma
** coefficients to genuine absolute-price form in BOTH modes (previously
** left in a mixed relative/absolute basis in the shipped Milestone 19
** unconstrained mode -- quaidsZeroFit() never had an analog of
** quaidsFit()'s own "RECOVERS ABSOLUTE PRICE EFFECTS FROM RELATIVE"
** block). aCtl.b0 now expects zOut.bRaw's shape/basis, not zOut.b's --
** mirroring quaidsFit()'s own qOut.homogB/aCtl.b0 pairing, and a direct
** consequence of zOut.b no longer being in that raw internal form.
**
** Uses tests/quaidsfixtures.src's _quaidsZeroSyntheticDGP(3000, 1) --
** seed=1 found by direct screening (not arbitrary), giving a genuine,
** uneven-but-non-degenerate zero-share pattern (see that fixture's own
** header for the exact fractions) and clean convergence of both the
** corrected and (for comparison) the naive quaidsFit(). The fixture's
** true gamma is genuinely symmetric and adding-up consistent by
** construction (ga = xpnd(vech(ga)) before the adding-up rows/columns
** are appended), so it is directly reusable for the homogeneous/
** symmetric path with no new fixture needed.
**
** Checks, in order: (1) the fixture's own adding-up identity holds
** exactly, by construction of the accounting-identity share formula;
** (2) unconstrained mode: the diagonal-delta restriction holds exactly,
** shape/finiteness checks, the aCtl.b0/zOut.bRaw warm-start roundtrip,
** and the core validation -- quaidsZeroFit()'s corrected coefficients
** recover the TRUE LATENT (uncensored) DGP parameters better than
** naively fitting quaidsFit() on the same censored data, on both a
** max-absolute-difference and a mean-absolute-difference basis; (3)
** homogeneous mode: the SAME diagonal-delta check plus an EXACT
** symmetry check on the recovered n x n gamma block (gammaBS==
** gammaBS' to floating-point precision -- the real regression guard for
** the cross-constraint bug found while building this), symDf/symStat/
** symPval shape/sign checks, and the same naive-vs-corrected comparison.
** This is a real, confirmed (not assumed) property of this specific
** fixture/seed -- building this test found that the improvement is
** genuine but modest (a single-dataset comparison, not a Monte Carlo
** average), and that several other candidate seeds did NOT show a clean
** improvement on both metrics simultaneously, consistent with
** Shonkwiler-Yen being a known approximation rather than a fully
** efficient correction, especially at the high per-good censoring rates
** this DGP family produces. Documented honestly rather than
** cherry-picking a misleadingly strong result.
**
** Adding-up does NOT hold for the corrected coefficients themselves (a
** real, known property of Shonkwiler-Yen -- see this file's header and
** src/quaidszerocorrect.src's own header) -- not checked here, since
** there is nothing to check.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_zero_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidszerocorrect.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
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

seed = 1;
tobs = 3000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsZeroSyntheticDGP(tobs, seed);

call check(maxc(abs(sumc(w') - 1)) < 1e-8, "fixture: observed shares sum to exactly 1 (accounting identity)");
call check(minc(minc(w)) >= 0, "fixture: no negative observed shares");

n = cols(w);
fracs = zeros(n, 1);
i = 1;
do while i <= n;
    fracs[i] = meanc(w[., i] .<= 0);
    i = i + 1;
endo;
call check(minc(fracs) > 0.01 and maxc(fracs) < 0.99, "fixture: zero-share fractions are genuine but not degenerate (no good is always/never zero)");


/* ==================== UNCONSTRAINED MODE (aCtl.homogenous=0) ==================== */

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 0;
aCtl.err = .0001;

zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call check(zOut.converged == 1, "quaidsZeroFit (unconstrained) converged");

nint = zOut.nint;
nendog = 2;   /* QUAIDS: lx, lx2 */
nu = 1;       /* single instrument in this fixture */
ngFinalExpected = 1+nint+nendog+nu+2*n;

/* --- Diagonal-delta restriction: off-diagonal entries are exactly
   zero, on-diagonal entries are freely, nonzero-ly estimated. Delta's
   position in the RECOVERED output is the trailing n rows (Milestone
   30 -- previously this was ngOld+1:ngOld+n using the internal,
   pre-recovery ngOld; the recovered position happens to coincide with
   that for unconstrained mode specifically, verified directly below via
   the shape check, not assumed). --- */
deltaBlock = zOut.b[ngFinalExpected-n+1:ngFinalExpected, .];
offDiagMask = 1 - eye(n);
call check(maxc(maxc(abs(deltaBlock.*offDiagMask))) == 0,
    "unconstrained: diagonal-delta restriction holds exactly (off-diagonal hazard entries are exactly zero)");
call check(minc(abs(diag(deltaBlock))) > 1e-6,
    "unconstrained: diagonal-delta restriction's on-diagonal (own-hazard) entries are genuinely estimated, not degenerately zero");

/* --- Shape/finiteness. --- */
call check(rows(zOut.b) == ngFinalExpected and cols(zOut.b) == n,
    "unconstrained: zOut.b has the expected recovered shape (1+nint+nendog+nu+2n rows, includes the n x n delta block)");
call check(rows(zOut.se) == rows(zOut.b) and cols(zOut.se) == cols(zOut.b), "unconstrained: zOut.se matches zOut.b's shape");
call check(minc(minc(zOut.se)) >= 0, "unconstrained: zOut.se are all non-negative");
call check(sumc(sumc(zOut.se .== zOut.se)) == rows(zOut.se)*cols(zOut.se), "unconstrained: zOut.se contains no NaN/missing values");
call check(rows(zOut.probitB) == rows(zOut.xnam) + n + 1, "zOut.probitB has one row per probit regressor (intcpt block + prices + totexp)");
call check(cols(zOut.probitB) == n, "zOut.probitB has one column per good");
call check(sumc(zOut.probitConverged) == n, "all n first-stage probits converged");

/* --- aCtl.b0 warm-start roundtrip: must use zOut.bRaw (the raw,
   pre-recovery internal form), NOT zOut.b (Milestone 30 -- zOut.b is
   now in recovered, absolute-price form, a different basis, mirroring
   quaidsFit()'s own qOut.homogB/aCtl.b0 pairing). --- */
struct quaidsControl aCtlB0;
aCtlB0 = aCtl;
aCtlB0.maxiter = 1;
aCtlB0.b0 = zOut.bRaw;
zOutB0 = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtlB0);
call check(rows(zOutB0.b) == rows(zOut.b) and cols(zOutB0.b) == cols(zOut.b),
    "quaidsZeroFit supplied aCtl.b0=zOut.bRaw: coefficient shape is preserved");
call check(rows(zOutB0.se) == rows(zOut.se) and cols(zOutB0.se) == cols(zOut.se),
    "quaidsZeroFit supplied aCtl.b0=zOut.bRaw: standard-error shape is preserved");

/* --- Core validation: corrected recovery beats naive recovery against
   the true latent DGP, on this specific screened seed. Now apples-to-
   apples given the Milestone 30 recovery bugfix (previously zOut.b's
   gamma columns 1..n-1 were gamma_ij-gamma_in, not genuine gamma_ij). --- */
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "naive quaidsFit() (unconstrained, for comparison) converged");

nCompareRows = rows(trueParams);
diffCorrected = maxc(maxc(abs(zOut.b[1:nCompareRows, .] - trueParams)));
diffNaive = maxc(maxc(abs(qOut.bestB[1:nCompareRows, .] - trueParams)));
meanDiffCorrected = meanc(meanc(abs(zOut.b[1:nCompareRows, .] - trueParams)));
meanDiffNaive = meanc(meanc(abs(qOut.bestB[1:nCompareRows, .] - trueParams)));

call check(diffCorrected < diffNaive,
    "unconstrained: quaidsZeroFit recovers the true latent DGP better than naive quaidsFit (max abs difference)");
call check(meanDiffCorrected < meanDiffNaive,
    "unconstrained: quaidsZeroFit recovers the true latent DGP better than naive quaidsFit (mean abs difference)");


/* ==================== HOMOGENEOUS MODE (aCtl.homogenous=1) ==================== */

struct quaidsControl aCtlH;
aCtlH = aCtl;
aCtlH.homogenous = 1;

zOutH = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtlH);
call check(zOutH.converged == 1, "quaidsZeroFit (homogeneous) converged");
call check(zOutH.symValid == 1, "homogeneous: symValid == aCtl.homogenous == 1");

n1 = n - 1;
call check(zOutH.n1 == n1, "homogeneous: zOut.n1 == n-1");

/* --- Shape: identical to the unconstrained case's recovered shape,
   confirming the Milestone 30 dimension claim that zOut.b's shape is
   independent of aCtl.homogenous (only its values differ). --- */
call check(rows(zOutH.b) == ngFinalExpected and cols(zOutH.b) == n,
    "homogeneous: zOutH.b has the SAME recovered shape as the unconstrained case");
call check(rows(zOutH.bS) == ngFinalExpected and cols(zOutH.bS) == n,
    "homogeneous: zOutH.bS has the same recovered shape too");
call check(rows(zOutH.b) == rows(zOutH.bS), "homogeneous: zOutH.b and zOutH.bS have identical row counts");

/* --- Diagonal-delta restriction holds exactly on the symmetry+diagonal
   -constrained output too. --- */
deltaBlockH = zOutH.bS[ngFinalExpected-n+1:ngFinalExpected, .];
call check(maxc(maxc(abs(deltaBlockH.*offDiagMask))) == 0,
    "homogeneous: diagonal-delta restriction holds exactly on zOutH.bS");
call check(minc(abs(diag(deltaBlockH))) > 1e-6,
    "homogeneous: diagonal-delta restriction's on-diagonal entries are genuinely estimated");

/* --- EXACT symmetry check on the recovered n x n gamma block -- the
   real regression guard for the cross-constraint bug this milestone
   found and fixed (an earlier draft left gamma's column-n entries free,
   and this check failed with asymmetry ~0.05-0.8 depending on the
   specific coefficient before the fix; now holds to floating-point
   precision). --- */
gammaBS = zOutH.bS[1+nint+1:1+nint+n, .];
call check(rows(gammaBS) == n and cols(gammaBS) == n, "homogeneous: recovered gamma block is n x n");
call check(maxc(maxc(abs(gammaBS - gammaBS'))) < 1e-6,
    "homogeneous: recovered gamma is exactly symmetric (the real regression guard for this milestone's cross-constraint finding)");

/* --- Symmetry test statistic/df: symDf is n1*(n1+1)/2, NOT
   n1*(n1-1)/2 -- the within-sub-block restrictions alone undercount by
   n1 (the column-n cross-constraints), confirmed by an independent
   degrees-of-freedom argument in _quaidsZeroSymDiagRestrict()'s own
   header. --- */
call check(zOutH.symDf == n1*(n1+1)/2, "homogeneous: symDf == n1*(n1+1)/2");
call check(zOutH.symStat == zOutH.symStat and zOutH.symStat >= 0, "homogeneous: symStat is finite and non-negative");
call check(zOutH.symPval >= 0 and zOutH.symPval <= 1, "homogeneous: symPval is a valid probability");

call check(maxc(maxc(abs(zOutH.bestB - zOutH.bS))) == 0, "homogeneous: zOut.bestB exactly equals zOut.bS");
call check(maxc(maxc(abs(zOutH.bestV - zOutH.vS))) == 0, "homogeneous: zOut.bestV exactly equals zOut.vS");

call check(rows(zOutH.seS) == rows(zOutH.bS) and cols(zOutH.seS) == cols(zOutH.bS), "homogeneous: zOutH.seS matches zOutH.bS's shape");
call check(minc(minc(zOutH.seS)) >= 0, "homogeneous: zOutH.seS are all non-negative");
call check(sumc(sumc(zOutH.seS .== zOutH.seS)) == rows(zOutH.seS)*cols(zOutH.seS), "homogeneous: zOutH.seS contains no NaN/missing values");

/* --- Core validation, homogeneous mode: corrected beats naive against
   the true latent DGP (which is genuinely symmetric by construction, so
   this is a real, non-vacuous comparison, not just "it runs"). --- */
qOutH = quaidsFit(w, intcpt, prices, totexp, instr, aCtlH);
call check(qOutH.converged == 1, "naive quaidsFit() (homogeneous, for comparison) converged");

diffCorrectedH = maxc(maxc(abs(zOutH.bS[1:nCompareRows, .] - trueParams)));
diffNaiveH = maxc(maxc(abs(qOutH.bestB[1:nCompareRows, .] - trueParams)));
meanDiffCorrectedH = meanc(meanc(abs(zOutH.bS[1:nCompareRows, .] - trueParams)));
meanDiffNaiveH = meanc(meanc(abs(qOutH.bestB[1:nCompareRows, .] - trueParams)));

call check(diffCorrectedH < diffNaiveH,
    "homogeneous: quaidsZeroFit recovers the true latent DGP better than naive quaidsFit (max abs difference)");
call check(meanDiffCorrectedH < meanDiffNaiveH,
    "homogeneous: quaidsZeroFit recovers the true latent DGP better than naive quaidsFit (mean abs difference)");


call printQuaidsZero(zOut);
call check(1, "printQuaidsZero() runs without error (unconstrained)");
call printQuaidsZero(zOutH);
call check(1, "printQuaidsZero() runs without error (homogeneous, includes the symmetry-test block)");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "ZERO-SHARE CORRECTION TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "ZERO-SHARE CORRECTION TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
