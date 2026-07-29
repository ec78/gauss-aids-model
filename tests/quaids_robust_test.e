/*
** quaids_robust_test.e
**
** Milestone 20: validates quaidsRobustFit()/printQuaidsRobust()
** (src/quaidsrobust.src) -- heteroskedasticity-robust and cluster-robust
** standard errors for an already-fitted quaidsFit() result, generalizing
** the classic S.*.inv(gg) sandwich to a per-observation (robust) or
** per-cluster (cluster-robust) score aggregation, unified through one
** clusterId argument (0 = robust, a Tx1 group-label vector = cluster-
** robust with a CR1 correction).
**
** Checks, in order: (1) the point estimate matches a fresh,
** independently hand-evaluated version of the sandwich formula (written
** directly in this test, not calling any src/ helper), confirming the
** implementation matches its own documented derivation exactly, not just
** that it runs; (2) the exact-identity regression guard that clusterId=0
** and an explicit seqa(1,1,nobs) per-row cluster label give byte-
** identical output -- confirming "robust is the G=nobs special case of
** cluster-robust" actually holds in code, not just on paper; (3) the
** reshape/cell-position regression guard, written from day one (unlike
** quaidsCurvatureFit()/quaidsCurvatureBootstrapFit(), where this exact
** class of bug was found and fixed only at Milestone 18, after shipping
** with it silently present since Milestone 10/15); (4) shape/finiteness/
** non-negativity of se/v; and (5) the core non-vacuous check -- on
** tests/quaidsfixtures.src's new _quaidsClusterSyntheticDGP fixture
** (genuine within-cluster-correlated noise, unlike every other fixture
** in this file), the cluster-robust se is measurably larger than the
** naive (clusterId=0) se, the standard textbook consequence of ignoring
** real clustering.
**
** A real, honestly-documented finding from building this test (see
** src/quaidsrobust.src's own header and docs/METHODOLOGY_NOTES.md): this
** proc's simplified bread makes its closed-form SE dramatically more
** CONSERVATIVE (larger) than quaidsFit()'s own classical homogSE/symcSE
** (which benefits from the full cross-equation FGLS system's efficiency
** and the nonlinear-price-index Jacobian correction this proc
** deliberately does not replicate) -- confirmed directly, not assumed,
** by comparing against an independently-built classical S.*.gg formula
** computed from the SAME regressors/residuals this proc itself uses
** (check 1's own hand-evaluation), which lands in the same order of
** magnitude, unlike qOut's own SE. This test does NOT assert closeness to
** qOut.homogSE/symcSE anywhere -- that comparison is apples-to-oranges
** between two genuinely different-efficiency estimators, not a property
** this proc's own formula is meant to have.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_robust_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidswelfare.src
#include ../src/quaidsrobust.src
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

/* --- Plain i.i.d. fixture: point-estimate cross-check, exact-identity
   regression guard, reshape regression guard, shape/finiteness. --- */

seed = 204;
tobs = 1000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, 1, 1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "prerequisite quaidsFit() converged");

struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);

{ bRobustFull, vRobustFull } = quaidsRobustCovariance(qOut, rOut, aCtl);
call check(maxc(maxc(abs(bRobustFull - qOut.bestB))) == 0,
    "quaidsRobustCovariance returns qOut.bestB as the full-basis point estimate");
call check(rows(vRobustFull) == rows(qOut.bestV) and cols(vRobustFull) == cols(qOut.bestV),
    "quaidsRobustCovariance returns a covariance shaped like qOut.bestV");
call check(maxc(maxc(abs(vRobustFull - vRobustFull'))) < 1e-8,
    "quaidsRobustCovariance returns a symmetric full-basis covariance");

/* Fresh, independent hand-evaluation of the sandwich formula --
   deliberately re-derived here, not calling any src/ helper, matching
   this project's established point-estimate cross-check convention
   (e.g. tests/quaids_shares_test.e). */
n = qOut.n;
n1 = qOut.n1;
nint = qOut.nint;
intcptFull = qOut.intcptFull;
nint1 = cols(intcptFull);

alphaHat = intcptFull*qOut.bestB[1:nint1, .];
gamaHat = qOut.bestB[nint1+1:nint1+n, .];
betaHat = qOut.bestB[nint1+n+1, .];
a_pHat = aCtl.alpha0 + sumc((prices.*alphaHat)') + .5*sumc(((prices*gamaHat).*prices)');
lxHat = totexp - a_pHat;
fittedW = alphaHat + prices*gamaHat + lxHat*betaHat;

uResidHand = w[., 1:n1] - fittedW[., 1:n1];
pricesHybridHand = prices;
pricesHybridHand[., 1:n-1] = prices[., 1:n-1] - prices[., n];
Xhand = intcptFull~pricesHybridHand[., 1:n1]~lxHat~qOut.u;
khand = cols(Xhand);
gghand = moment(Xhand, 0)/tobs;

InflHand = zeros(tobs, n1*khand);
i = 1;
do while i <= n1;
    InflHand[., (i-1)*khand+1:i*khand] = Xhand.*uResidHand[., i];
    i = i + 1;
endo;
Kdofhand = n1*khand;
chand = (tobs/(tobs-1))*((tobs-1)/(tobs-Kdofhand));
MeatHand = chand*(InflHand'InflHand)/tobs;
breadHand = eye(n1).*.invpd(gghand);
vHand = breadHand*MeatHand*breadHand;
seHand = reshape(sqrt(diag(vHand)), n1, khand)';

call check(maxc(maxc(abs(rOut.se - seHand))) < 1e-8, "quaidsRobustFit point estimate matches a fresh, independent hand-evaluation");
call check(maxc(maxc(abs(rOut.v - vHand))) < 1e-6, "quaidsRobustFit covariance matches the fresh hand-evaluation");

/* Independent classical (S.*.gg-style) formula built from the SAME
   X/U this proc itself uses -- confirms the robust formula stays in the
   same order of magnitude as its own classical special case under
   genuine homoskedasticity (this fixture has plain i.i.d. noise), unlike
   qOut's own, much more efficient, full-system SE (see this file's own
   header for why that comparison would be apples-to-oranges). */
Sclassical = uResidHand'uResidHand/tobs;
vClassical = breadHand*(Sclassical.*.gghand)*breadHand;
seClassical = reshape(sqrt(diag(vClassical)), n1, khand)';
ratio = rOut.se./seClassical;
call check(minc(minc(ratio)) > 0.5 and maxc(maxc(ratio)) < 5, "robust se stays within the same order of magnitude as the classical special case (same X/U)");

/* Exact-identity regression guard: clusterId=0 is the literal G=nobs
   special case of the cluster-robust formula. */
struct quaidsRobustOut rOutExplicit;
rOutExplicit = quaidsRobustFit(qOut, w, prices, totexp, aCtl, seqa(1, 1, tobs));
call check(maxc(maxc(abs(rOut.se - rOutExplicit.se))) == 0, "clusterId=0 and an explicit per-row seqa(1,1,nobs) label give byte-identical se");
call check(maxc(maxc(abs(rOut.v - rOutExplicit.v))) == 0, "clusterId=0 and an explicit per-row seqa(1,1,nobs) label give byte-identical v");
call check(rOut.nClusters == tobs, "clusterId=0: nClusters == nobs");

/* Reshape/cell-position regression guard -- written from day one. */
seCheck = reshape(sqrt(diag(rOut.v)), cols(rOut.b), rows(rOut.b))';
call check(maxc(maxc(abs(rOut.se - seCheck))) == 0, "se's cells are correctly positioned relative to v (reshape/transpose regression guard)");

/* Shape/finiteness/non-negativity. */
call check(rows(rOut.b) == rows(rOut.se) and cols(rOut.b) == cols(rOut.se), "rOut.se matches rOut.b's shape");
call check(cols(rOut.b) == n1, "rOut.b has n1 columns (equation n excluded, recovered via adding-up)");
call check(minc(minc(rOut.se)) >= 0, "rOut.se are all non-negative");
call check(sumc(sumc(rOut.se .== rOut.se)) == rows(rOut.se)*cols(rOut.se), "rOut.se contains no NaN/missing values");

nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsSharesOut sharesClassical;
struct quaidsSharesOut sharesRobust;
sharesClassical = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
sharesRobust = quaidsSharesFit(bRobustFull, vRobustFull, intcptMean, pricesMean, totexpMean, aCtl);
call check(maxc(abs(sharesRobust.w - sharesClassical.w)) == 0,
    "robust covariance propagation leaves predicted share point estimates unchanged");
call check(maxc(abs(sharesRobust.se - sharesClassical.se)) > 1e-8,
    "robust covariance propagation changes predicted share SE");

struct quaidsElasOut elasClassical;
struct quaidsElasOut elasRobust;
elasClassical = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
elasRobust = quaidsElasFit(bRobustFull, vRobustFull, intcptMean, pricesMean, totexpMean, aCtl);
call check(maxc(abs(elasRobust.er - elasClassical.er)) == 0 and maxc(maxc(abs(elasRobust.ep - elasClassical.ep))) == 0,
    "robust covariance propagation leaves elasticity point estimates unchanged");
call check(maxc(abs(elasRobust.ser - elasClassical.ser)) > 1e-8,
    "robust covariance propagation changes income elasticity SE");

pricesPt1 = pricesMean;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

struct quaidsWelfareOut welfareClassical;
struct quaidsWelfareOut welfareRobust;
welfareClassical = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, pricesPt1, totexpMean, aCtl);
welfareRobust = quaidsWelfareFit(bRobustFull, vRobustFull, intcptMean, pricesMean, pricesPt1, totexpMean, aCtl);
call check(welfareRobust.cv == welfareClassical.cv and welfareRobust.ev == welfareClassical.ev,
    "robust covariance propagation leaves welfare point estimates unchanged");
call check(abs(welfareRobust.seCV - welfareClassical.seCV) > 1e-8 or abs(welfareRobust.seEV - welfareClassical.seEV) > 1e-8,
    "robust covariance propagation changes welfare SE");

call printQuaidsRobust(rOut);
call check(1, "printQuaidsRobust() runs without error");


/* --- Cluster fixture: the core non-vacuous check -- ignoring genuine
   within-cluster correlation understates sampling variability. --- */

nClustersTrue = 40;
{ wc, intcptc, pricesc, totexpc, instrc, clusterId } = _quaidsClusterSyntheticDGP(tobs, seed, nClustersTrue);

struct quaidsOut qOutC;
qOutC = quaidsFit(wc, intcptc, pricesc, totexpc, instrc, aCtl);
call check(qOutC.converged == 1, "cluster-fixture quaidsFit() converged");

struct quaidsRobustOut rOutNaive;
rOutNaive = quaidsRobustFit(qOutC, wc, pricesc, totexpc, aCtl, 0);

struct quaidsRobustOut rOutCluster;
rOutCluster = quaidsRobustFit(qOutC, wc, pricesc, totexpc, aCtl, clusterId);

call check(rOutCluster.nClusters == nClustersTrue, "cluster fixture: nClusters matches the true number of clusters");
call check(meanc(meanc(rOutCluster.se)) > meanc(meanc(rOutNaive.se)), "cluster-robust se is measurably larger than the naive (per-row) se on genuinely clustered data");
call check(rOutCluster.se[1, 1] > rOutNaive.se[1, 1], "cluster-robust se exceeds naive se at a specific cell too (not just on average)");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "ROBUST/CLUSTER STANDARD ERROR TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "ROBUST/CLUSTER STANDARD ERROR TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
