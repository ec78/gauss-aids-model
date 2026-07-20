/*
** quaids_hypothesis_tests_test.e
**
** Milestone 4: validates the new standalone quaidsHomogeneityTest() and
** quaidsJointTest() (src/quaidstests.src), strengthens coverage of the
** existing symmetry-given-homogeneity test with a power check, and
** exercises the overidentification test for the first time in this
** repo's history (every prior fixture used exactly-identified instruments,
** ninst==nu, so qOut.overidValid was always 0 -- an untested code path
** until now).
**
** Every hypothesis test here is checked for both SIZE (does it correctly
** fail to reject when the null is true by construction?) and, where
** practical, POWER (does it correctly reject when the null is false by
** construction?) -- a test that always passes regardless of the data
** would be worthless; checking only size cannot catch that.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_hypothesis_tests_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstests.src;
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
   Homogeneity and joint tests: size and power.
   Base data: the standard homogeneity-true-by-construction fixture
   (seed=204, matching Milestone 3), fit UNCONSTRAINED (aCtl.homogenous=0)
   so quaidsHomogeneityTest()/quaidsJointTest() have something to test.
   ========================================================================== */

struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 0;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 204, 1, 1);

struct quaidsOut qOutTrue;
qOutTrue = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ stat, pval, df } = quaidsHomogeneityTest(qOutTrue);
call check(df == qOutTrue.n - 1, "homogeneity test: df == n-1 (4)");
call check(pval > 0.05, "homogeneity test SIZE: fails to reject on a true-homogeneity DGP (pval > 0.05)");

{ statJ, pvalJ, dfJ } = quaidsJointTest(qOutTrue);
call check(dfJ == (qOutTrue.n-1) + (qOutTrue.n-1)*(qOutTrue.n-2)/2, "joint test: df == (n-1)+(n-1)(n-2)/2 (10)");
call check(pvalJ > 0.05, "joint test SIZE: fails to reject on a true-homogeneity+symmetry DGP (pval > 0.05)");

/* Power: inject a clean, deliberate homogeneity violation into the
   OBSERVED shares (a share responding to one price without a fully
   offsetting response elsewhere is precisely a homogeneity violation --
   scaling all prices together no longer leaves w invariant). Preserves
   adding-up exactly (added to one share, subtracted from another). */
wViol = w;
c = 1.5;
wViol[.,1] = wViol[.,1] + c*prices[.,1];
wViol[.,2] = wViol[.,2] - c*prices[.,1];

struct quaidsOut qOutViol;
qOutViol = quaidsFit(wViol, intcpt, prices, totexp, instr, aCtl);

{ statV, pvalV, dfV } = quaidsHomogeneityTest(qOutViol);
call check(pvalV < 0.01, "homogeneity test POWER: rejects on a deliberately homogeneity-violating dataset (pval < 0.01)");

{ statJV, pvalJV, dfJV } = quaidsJointTest(qOutViol);
call check(pvalJV < 0.01, "joint test POWER: rejects on the same homogeneity-violating dataset (pval < 0.01)");

call check(qOutTrue.homogenous == 0, "quaidsHomogeneityTest/quaidsJointTest operate on the unconstrained (homogenous=0) fit");


/* ==========================================================================
   Symmetry-given-homogeneity test: strengthen with a power check.
   Build a gamma matrix that satisfies adding-up AND homogeneity (both row
   and column sums exactly zero, by construction of the fixup below) but is
   NOT symmetric -- isolates a pure symmetry violation. (A SYMMETRIC gamma
   with adding-up automatically satisfies homogeneity too: row i sum =
   sum_j gamma_ij = sum_j gamma_ji [symmetry] = column i sum = 0. So a
   homogeneity-preserving symmetry violation necessarily starts from an
   asymmetric core -- that is what this DGP does, deliberately.)
   ========================================================================== */

seedSym = 305;
tobsSym = 3000;
N = 5;

al = round(rndns(1,N-1,seedSym)*10)/10;
al = al~(1-sumc(al'));
al1 = .5*round(rndns(1,N-1,seedSym)*10)/10;
al1 = al1~(-sumc(al1'));

/* Asymmetric (N-1)x(N-1) core -- NOT run through xpnd(vech(.)). The
   row-then-column fixup below still forces both row and column sums of
   the FULL NxN matrix to zero exactly (homogeneity + adding-up both
   hold), while leaving the matrix asymmetric (symmetry violated). */
gaCore = round(rndns(N-1,N-1,seedSym)*10)/10;
ga = gaCore|(-sumc(gaCore)');
ga = ga~(-sumc(ga'));
call check(maxc(abs(sumc(ga))) < 1e-8, "symmetry-violation fixture: column sums (adding-up) are exactly 0");
call check(maxc(abs(sumc(ga'))) < 1e-8, "symmetry-violation fixture: row sums (homogeneity) are exactly 0");
call check(maxc(maxc(abs(ga - ga'))) > 0.1, "symmetry-violation fixture: gamma is NOT symmetric (max|gamma-gamma''| > 0.1)");

be = .5*round(rndns(1,N-1,seedSym)*10)/10;
be = be~(-sumc(be'));
ro = round(rndns(1,N-1,seedSym)*10)/10;

pricesSym = 1+rndns(tobsSym,N,seedSym);
instrSym = 5+5*rndns(tobsSym,1,seedSym);
intcptSym = 2+2*rndns(tobsSym,1,seedSym);
uSym =  .1*rndns(tobsSym,1,seedSym);
totexpSym = .85*instrSym + uSym;
eSym = 2*rndns(tobsSym,N-1,seedSym) + uSym*ro;
eSym = eSym~(-sumc(eSym'));

a_pSym = sumc( (pricesSym.*(al+intcptSym*al1) )') + .5*sumc(((pricesSym*ga).*pricesSym)');
lxSym = totexpSym - a_pSym;
wSym = al + pricesSym*ga + lxSym*be + eSym + intcptSym*al1;

struct quaidsControl aCtlHom;
aCtlHom = quaidsControlCreate;
aCtlHom.linear = 1;
aCtlHom.maxiter = 100;
aCtlHom.homogenous = 1;
aCtlHom.err = .0001;

struct quaidsOut qOutSym;
qOutSym = quaidsFit(wSym, intcptSym, pricesSym, totexpSym, instrSym, aCtlHom);

call check(qOutSym.symValid == 1, "symmetry test runs when aCtl.homogenous=1");
call check(qOutSym.symPval < 0.01, "symmetry-given-homogeneity test POWER: rejects on a deliberately asymmetric-but-homogeneous DGP (pval < 0.01)");

/* Size, for contrast: the ORIGINAL (symmetric-by-construction) fixture,
   already implicitly checked in quaids_synthetic_validation_test.e, but
   confirmed explicitly here alongside the power check for a complete
   picture in one place. Must fit with aCtl.linear=0 to match: `w` was
   generated with quadratic=1 (a true QUAIDS DGP); fitting it with
   aCtlHom.linear=1 would omit the true quadratic term and bias every
   coefficient (including gamma) from specification error alone, which
   would spuriously reject symmetry for a reason that has nothing to do
   with symmetry -- confirmed by hitting exactly that failure mode first
   and tracing it back to this mismatch. */
struct quaidsControl aCtlHomQ;
aCtlHomQ = quaidsControlCreate;
aCtlHomQ.linear = 0;
aCtlHomQ.maxiter = 100;
aCtlHomQ.homogenous = 1;
aCtlHomQ.err = .0001;

struct quaidsOut qOutSymTrue;
qOutSymTrue = quaidsFit(w, intcpt, prices, totexp, instr, aCtlHomQ);
call check(qOutSymTrue.symPval > 0.05, "symmetry-given-homogeneity test SIZE: fails to reject on a truly symmetric DGP (pval > 0.05)");


/* ==========================================================================
   Overidentification test: first-ever exercise of this code path.
   Every prior fixture in this repo used exactly-identified instruments
   (ninst==nu), so qOut.overidValid was always 0. Build a 2-instrument
   variant (ninst=2 > nu=1) so the ninst>nu branch actually runs.
   ========================================================================== */

seedOI = 406;
tobsOI = 3000;

alO = round(rndns(1,N-1,seedOI)*10)/10;
alO = alO~(1-sumc(alO'));
al1O = .5*round(rndns(1,N-1,seedOI)*10)/10;
al1O = al1O~(-sumc(al1O'));
gaCoreO = round(rndns(N-1,N-1,seedOI)*10)/10;
gaO = xpnd(vech(gaCoreO));
gaO = gaO|(-sumc(gaO)');
gaO = gaO~(-sumc(gaO'));
beO = .5*round(rndns(1,N-1,seedOI)*10)/10;
beO = beO~(-sumc(beO'));
roO = round(rndns(1,N-1,seedOI)*10)/10;

pricesO = 1+rndns(tobsOI,N,seedOI);
/* Two valid, independent instruments, both genuinely driving totexp. */
instrO = (5+5*rndns(tobsOI,1,seedOI))~(3+3*rndns(tobsOI,1,seedOI));
intcptO = 2+2*rndns(tobsOI,1,seedOI);
uO =  .1*rndns(tobsOI,1,seedOI);
totexpO = .60*instrO[.,1] + .25*instrO[.,2] + uO;
eO = 2*rndns(tobsOI,N-1,seedOI) + uO*roO;
eO = eO~(-sumc(eO'));

a_pO = sumc( (pricesO.*(alO+intcptO*al1O) )') + .5*sumc(((pricesO*gaO).*pricesO)');
lxO = totexpO - a_pO;
wO = alO + pricesO*gaO + lxO*beO + eO + intcptO*al1O;

struct quaidsOut qOutOI;
qOutOI = quaidsFit(wO, intcptO, pricesO, totexpO, instrO, aCtlHom);

call check(qOutOI.ninst == 2, "overID fixture: ninst == 2");
call check(qOutOI.overidValid == 1, "overID test runs (ninst=2 > nu=1) -- first time this code path has been exercised");
call check(qOutOI.overidDf == 1, "overID test: df == ninst-nu (1)");
call check(rows(qOutOI.overidGamma) == 2 and cols(qOutOI.overidGamma) == 5, "overID test: overidGamma is ninst x n (2x5)");
call check(minc(qOutOI.overidPvf) >= 0 and maxc(qOutOI.overidPvf) <= 1, "overID test: all per-equation p-values are valid probabilities");
/* SIZE check: both instruments are valid by construction, so the test
   should typically NOT reject overidentification for most/all equations
   -- not a hard guarantee for every single equation at any one seed, but
   the MEAN p-value across equations should be comfortably away from 0. */
call check(meanc(qOutOI.overidPvf) > 0.10, "overID test SIZE: mean p-value across equations is not close to 0 (both instruments valid by construction)");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "HYPOTHESIS TESTS: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "HYPOTHESIS TESTS: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
