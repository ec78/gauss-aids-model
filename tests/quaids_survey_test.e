/*
** quaids_survey_test.e
**
** Milestone 26: validates the new estimator-level sampling-weight argument
** threaded through quaidsFit()/_quaidsIVFirstStage() (src/quaids.src,
** src/quaidsiv.src), quaidsPreflight() (src/quaidsdiagnostics.src), and
** quaidsRobustFit()/quaidsRobustBootstrapFit() (src/quaidsrobust.src).
**
** Checks, in order: (1) exact-identity regression guard -- weight omitted
** and an explicit uniform weight (ones(nobs,1)) give byte-identical
** bestB/bestV to the pre-Milestone-26 unweighted quaidsFit(), and the same
** for quaidsPreflight()/quaidsRobustFit()/quaidsRobustBootstrapFit(); (2)
** quaidsPreflight()'s weight validation correctly flags an invalid
** (negative) weight as a hard preflight error while leaving a valid one
** passing; (3) the single easiest detail in this milestone to get
** backwards, checked directly: quaidsRobustFit()'s sandwich scales the
** bread by sqrt(weight) (matching quaidsFit()'s own WLS design) but the
** per-observation score contribution by PLAIN weight (the Horvitz-Thompson
** pweight-robust convention) -- confirmed by comparing against two fresh,
** independent hand-evaluations, one using the documented (correct)
** convention and one deliberately using the wrong one; (4) the core
** non-vacuous check -- on tests/quaidsfixtures.src's new
** _quaidsSurveyWeightedDGP fixture (informative, biased sampling on an
** unobserved error term), weighting the estimator by the fixture's own
** Horvitz-Thompson weight recovers the true population parameters
** measurably better than ignoring the weight on the same biased sample,
** mirroring Milestone 19's own "corrected-beats-naive" validation pattern.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_survey_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
#include ../src/quaidsdiagnostics.src
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

/* --- Exact-identity regression guard: weight omitted / explicit uniform
   weight must reproduce the pre-Milestone-26 unweighted behavior exactly,
   across every proc that gained the new argument. --- */

seed = 204;
tobs = 1000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, 1, 1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsOut qOutBase;
qOutBase = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOutBase.converged == 1, "prerequisite unweighted quaidsFit() converged");
call check(qOutBase.weighted == 0 and qOutBase.weightSum == tobs and qOutBase.effN == tobs,
    "unweighted quaidsFit reports weighted=0, weightSum=effN=nobs");

wgtOnes = ones(tobs, 1);
struct quaidsOut qOutOnes;
qOutOnes = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, wgtOnes);
call check(qOutOnes.weighted == 1 and qOutOnes.weightSum == tobs and qOutOnes.effN == tobs,
    "explicit uniform weight reports weighted=1, weightSum=effN=nobs");
call check(maxc(maxc(abs(qOutOnes.bestB - qOutBase.bestB))) == 0,
    "explicit uniform weight reproduces unweighted bestB exactly");
call check(maxc(maxc(abs(qOutOnes.bestV - qOutBase.bestV))) == 0,
    "explicit uniform weight reproduces unweighted bestV exactly");

struct quaidsPreflightOut pOutBase;
struct quaidsPreflightOut pOutOnes;
pOutBase = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, 0);
pOutOnes = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, wgtOnes);
call check(pOutBase.weightValid == 1 and pOutBase.weightSum == tobs and pOutBase.effN == tobs,
    "quaidsPreflight with weight=0 reports weightValid=1, weightSum=effN=nobs");
call check(pOutOnes.weightValid == 1 and pOutOnes.weightSum == tobs and pOutOnes.effN == tobs,
    "quaidsPreflight with an explicit uniform weight reports the same diagnostics");
call check(pOutBase.ok == pOutOnes.ok and pOutBase.nErrors == pOutOnes.nErrors,
    "quaidsPreflight ok/nErrors unaffected by an explicit uniform weight");

struct quaidsRobustOut rOutBase;
struct quaidsRobustOut rOutOnes;
rOutBase = quaidsRobustFit(qOutBase, w, prices, totexp, aCtl, 0);
rOutOnes = quaidsRobustFit(qOutBase, w, prices, totexp, aCtl, 0, wgtOnes);
call check(maxc(maxc(abs(rOutOnes.b - rOutBase.b))) == 0 and maxc(maxc(abs(rOutOnes.se - rOutBase.se))) == 0,
    "quaidsRobustFit with an explicit uniform weight reproduces the unweighted sandwich exactly");

struct quaidsRobustBootOut rbOutBase;
struct quaidsRobustBootOut rbOutOnes;
rbOutBase = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 10, seed=42);
rbOutOnes = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 10, seed=42, weight=wgtOnes);
call check(maxc(maxc(abs(rbOutOnes.b - rbOutBase.b))) == 0 and maxc(maxc(abs(rbOutOnes.seRobust - rbOutBase.seRobust))) == 0,
    "quaidsRobustBootstrapFit with an explicit uniform weight reproduces the unweighted base fit/sandwich exactly");

/* --- quaidsPreflight() weight validation: invalid weight is a hard
   preflight error; a valid one is not. --- */

wgtBad = wgtOnes;
wgtBad[1] = -1;
struct quaidsPreflightOut pOutBad;
pOutBad = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, wgtBad);
call check(pOutBad.weightValid == 0 and pOutBad.ok == 0,
    "a negative weight entry is flagged as a hard preflight error");
call check(pOutBad.nErrors == pOutBase.nErrors + 1,
    "an invalid weight adds exactly one preflight error relative to the clean baseline");

/* --- The single easiest detail in this milestone to get backwards:
   quaidsRobustFit()'s bread uses sqrt(weight), its per-observation score
   contribution (meat) uses PLAIN weight. Checked directly against two
   fresh, independent hand-evaluations -- one matching the documented
   (correct) convention, one deliberately using the wrong one. --- */

rndseed 77;
wgtUneq = 1 + rndu(tobs, 1)*3;
wgtUneq = wgtUneq*tobs/sumc(wgtUneq);

struct quaidsOut qOutW;
qOutW = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, wgtUneq);
call check(qOutW.converged == 1, "prerequisite weighted quaidsFit() converged");

struct quaidsRobustOut rOutW;
rOutW = quaidsRobustFit(qOutW, w, prices, totexp, aCtl, 0, wgtUneq);

n = qOutW.n;
n1 = qOutW.n1;
nint = qOutW.nint;
intcptFull = qOutW.intcptFull;
nint1 = cols(intcptFull);

alphaHat = intcptFull*qOutW.bestB[1:nint1, .];
gamaHat = qOutW.bestB[nint1+1:nint1+n, .];
betaHat = qOutW.bestB[nint1+n+1, .];
a_pHat = aCtl.alpha0 + sumc((prices.*alphaHat)') + .5*sumc(((prices*gamaHat).*prices)');
lxHat = totexp - a_pHat;
fittedW = alphaHat + prices*gamaHat + lxHat*betaHat;

uResidHand = w[., 1:n1] - fittedW[., 1:n1];
pricesHybridHand = prices;
pricesHybridHand[., 1:n-1] = prices[., 1:n-1] - prices[., n];
Xhand = intcptFull~pricesHybridHand[., 1:n1]~lxHat~qOutW.u;
khand = cols(Xhand);
rootWtHand = sqrt(wgtUneq);
Kdofhand = n1*khand;
chand = (tobs/(tobs-1))*((tobs-1)/(tobs-Kdofhand));
breadHand = eye(n1).*.invpd(moment(rootWtHand.*Xhand, 0)/tobs);

/* Correct convention: meat's score contribution scaled by PLAIN weight. */
InflCorrect = zeros(tobs, n1*khand);
i = 1;
do while i <= n1;
    InflCorrect[., (i-1)*khand+1:i*khand] = wgtUneq.*(Xhand.*uResidHand[., i]);
    i = i + 1;
endo;
vCorrect = breadHand*(chand*(InflCorrect'InflCorrect)/tobs)*breadHand;
seCorrect = reshape(sqrt(diag(vCorrect)), n1, khand)';

/* Deliberately WRONG convention: meat's score contribution scaled by
   sqrt(weight) instead of plain weight -- the natural mistake this
   milestone's own header warns about. */
InflWrong = zeros(tobs, n1*khand);
i = 1;
do while i <= n1;
    InflWrong[., (i-1)*khand+1:i*khand] = rootWtHand.*(Xhand.*uResidHand[., i]);
    i = i + 1;
endo;
vWrong = breadHand*(chand*(InflWrong'InflWrong)/tobs)*breadHand;
seWrong = reshape(sqrt(diag(vWrong)), n1, khand)';

call check(maxc(maxc(abs(rOutW.se - seCorrect))) < 1e-8,
    "quaidsRobustFit matches the documented sqrt(weight)-bread / plain-weight-meat convention");
call check(maxc(maxc(abs(rOutW.se - seWrong))) > 1e-4,
    "quaidsRobustFit does NOT match the sqrt(weight)-meat (wrong) convention");

/* --- Core non-vacuous check: weighted recovery beats naive on a
   genuinely informatively-sampled fixture. --- */

nPop = 6000;
seedSurvey = 11;
{ wS, intcptS, pricesS, totexpS, instrS, weightS, trueParamsS } = _quaidsSurveyWeightedDGP(nPop, seedSurvey);

struct quaidsControl aCtlS;
aCtlS = quaidsControlCreate();
aCtlS.linear = 1;
aCtlS.maxiter = 100;
aCtlS.homogenous = 1;
aCtlS.err = .0001;

struct quaidsOut qOutNaive;
qOutNaive = quaidsFit(wS, intcptS, pricesS, totexpS, instrS, aCtlS);
struct quaidsOut qOutWeighted;
qOutWeighted = quaidsFit(wS, intcptS, pricesS, totexpS, instrS, aCtlS, weightS);

call check(qOutNaive.converged == 1 and qOutWeighted.converged == 1,
    "survey fixture: both naive and weighted quaidsFit() converge");

nRowsS = rows(trueParamsS);
nuS = qOutNaive.nu;
rowsStructural = seqa(1, 1, nRowsS - nuS);

maxDiffNaive = maxc(maxc(abs(qOutNaive.bestB[rowsStructural, .] - trueParamsS[rowsStructural, .])));
maxDiffWeighted = maxc(maxc(abs(qOutWeighted.bestB[rowsStructural, .] - trueParamsS[rowsStructural, .])));
meanDiffNaive = meanc(meanc(abs(qOutNaive.bestB[rowsStructural, .] - trueParamsS[rowsStructural, .])));
meanDiffWeighted = meanc(meanc(abs(qOutWeighted.bestB[rowsStructural, .] - trueParamsS[rowsStructural, .])));

call check(maxDiffWeighted < maxDiffNaive,
    "weighted estimation recovers the true population structural coefficients better (max abs diff) than naive on biased-sample data");
call check(meanDiffWeighted < meanDiffNaive,
    "weighted estimation recovers the true population structural coefficients better (mean abs diff) than naive on biased-sample data");
call check(qOutWeighted.weighted == 1 and qOutWeighted.effN < qOutWeighted.weightSum,
    "survey fixture: weighted fit reports weighted=1 and effN < weightSum (unequal weighting genuinely reduces effective sample size)");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "SURVEY WEIGHT TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "SURVEY WEIGHT TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
