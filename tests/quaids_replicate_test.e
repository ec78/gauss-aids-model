/*
** quaids_replicate_test.e
**
** Milestone 27: validates quaidsReplicateWeightFit()/
** printQuaidsReplicateWeight() (src/quaidsreplicate.src) -- replicate-
** weight (jackknife/BRR-style) standard errors built from a caller-
** supplied set of pre-computed replicate weight columns and scale
** factor(s), reusing Milestone 26's quaidsFit() weight argument for both
** the full-sample point estimate and every replicate refit.
**
** Checks, in order: (1) the base (full-sample) point estimate matches a
** direct quaidsFit() call exactly, both unweighted and weighted; (2)
** shape/finiteness/non-negativity/symmetry of se/v, and the reshape/
** cell-position regression guard (the exact bug class Milestone 18 found
** twice already); (3) scaleFactor accepted as either a scalar or an Rx1
** vector, and echoed correctly either way; (4) the method label is
** echoed verbatim; (5) an EXACT zero-variance identity check -- when
** every replicate weight column is IDENTICAL to the base weight, every
** replicate refit must reproduce the full-sample bestB exactly, so v/se
** must be exactly zero, not just small -- the strongest possible
** correctness check of the variance formula itself, contrasted with (6)
** a genuinely different (real JK1-style) replicate design producing
** genuinely nonzero se; (7) deterministic repeatability; (8) rOut.b/
** rOut.v are directly usable as (b, v) inputs to quaidsSharesFit()/
** quaidsElasFit()/quaidsWelfareFit() with no expansion step, unlike
** quaidsRobustFit()'s own reduced basis; (9) partial-failure handling --
** a pathological low-effective-sample-size replicate column is dropped
** (nFailed>0, nCompleted<nRequested) via the effNMin pre-check rather
** than crashing the whole job (see src/quaidsreplicate.src's own header
** for the real, empirically-reproduced G0058 crash this pre-check
** avoids -- a non-trappable failure mode, confirmed directly, not
** assumed, the same class already documented for eighv()/glm()).
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_replicate_test.e
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
#include ../src/quaidsreplicate.src
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

/* --- JK1-style delete-one-cluster replicate weights from 20 pseudo-clusters. --- */
nClusters = 20;
idx = seqa(1, 1, tobs);
clusterId = ceil(idx/(tobs/nClusters));
replicateWeights = zeros(tobs, nClusters);
g = 1;
do while g <= nClusters;
    replicateWeights[., g] = (clusterId ./= g) * nClusters/(nClusters-1);
    g = g + 1;
endo;
scaleFactorJK1 = (nClusters-1)/nClusters;

struct quaidsReplicateOut rOut;
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactorJK1, method="JK1");

call check(rOut.weighted == 0, "unweighted base fit reports weighted=0");
call check(maxc(maxc(abs(rOut.b - qOutBase.bestB))) == 0,
    "base point estimate matches a direct unweighted quaidsFit() call exactly");
call check(rOut.nRequested == nClusters and rOut.nCompleted == nClusters and rOut.nFailed == 0,
    "all JK1 replicates converge on the clean fixture");
call check(rOut.method $== "JK1", "method label is echoed verbatim");

call check(rows(rOut.se) == rows(rOut.b) and cols(rOut.se) == cols(rOut.b), "se shape matches b shape");
call check(sumc(sumc(rOut.se .== rOut.se)) == rows(rOut.se)*cols(rOut.se), "se contains no NaN/missing values");
call check(minc(minc(rOut.se)) >= 0, "se are all non-negative");
call check(maxc(maxc(abs(rOut.v - rOut.v'))) < 1e-8, "v is symmetric");

/* Reshape/cell-position regression guard -- written from day one. */
seCheck = reshape(sqrt(diag(rOut.v)), cols(rOut.b), rows(rOut.b))';
call check(maxc(maxc(abs(rOut.se - seCheck))) == 0,
    "se's cells are correctly positioned relative to v (reshape/transpose regression guard)");

call check(rows(rOut.bReplicate) == rOut.nCompleted and cols(rOut.bReplicate) == rows(rOut.b)*cols(rOut.b),
    "bReplicate shape matches nCompleted x vec(b) dimensions");

call check(rows(rOut.scaleFactor) == nClusters and maxc(abs(rOut.scaleFactor - scaleFactorJK1)) == 0,
    "scalar scaleFactor is broadcast to an Rx1 vector correctly");

/* --- Rx1 (BRR-style) scale factor vector. --- */
scaleFactorVec = (1/nClusters)*ones(nClusters, 1);
struct quaidsReplicateOut rOutVecScale;
rOutVecScale = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactorVec, method="BRR");
call check(maxc(abs(rOutVecScale.scaleFactor - scaleFactorVec)) == 0,
    "Rx1 scaleFactor vector is echoed unchanged");
call check(rOutVecScale.method $== "BRR", "second method label is echoed verbatim");

/* --- Weighted base fit. --- */
rndseed 55;
wgtBase = 1 + rndu(tobs, 1)*3;
struct quaidsOut qOutW;
qOutW = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, wgtBase);
struct quaidsReplicateOut rOutWeighted;
rOutWeighted = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactorJK1, weight=wgtBase, method="JK1");
call check(rOutWeighted.weighted == 1, "weighted base fit reports weighted=1");
call check(maxc(maxc(abs(rOutWeighted.b - qOutW.bestB))) == 0,
    "base point estimate matches a direct weighted quaidsFit() call exactly");

/* --- Exact zero-variance identity: every replicate weight column
   IDENTICAL to the base weight must reproduce the full-sample bestB
   exactly on every replicate, so v/se must be EXACTLY zero, not just
   small -- the strongest possible correctness check of the variance
   formula itself. --- */
identicalReplicateWeights = ones(tobs, 5);
struct quaidsReplicateOut rOutIdentical;
rOutIdentical = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, identicalReplicateWeights, 1, method="identical");
call check(maxc(maxc(abs(rOutIdentical.v))) == 0,
    "replicate weights identical to the base weight give an EXACTLY zero covariance");
call check(maxc(maxc(rOutIdentical.se)) == 0,
    "replicate weights identical to the base weight give EXACTLY zero standard errors");
call check(rOutIdentical.nCompleted == 5 and rOutIdentical.nFailed == 0,
    "identical-weight replicates all converge (same design as the base fit)");

/* Contrast: a genuinely different (real JK1) replicate design gives
   genuinely nonzero se -- confirms the zero-variance case above is not
   simply a vacuous always-zero formula. */
call check(maxc(maxc(rOut.se)) > 1e-6,
    "a genuinely different replicate design gives non-vacuous (nonzero) standard errors");

/* --- Deterministic repeatability. --- */
struct quaidsReplicateOut rOutRepeat;
rOutRepeat = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeights, scaleFactorJK1, method="JK1");
call check(maxc(maxc(abs(rOut.se - rOutRepeat.se))) == 0 and maxc(maxc(abs(rOut.v - rOutRepeat.v))) == 0,
    "repeated evaluation with identical inputs is deterministic");

/* --- rOut.b/rOut.v are directly usable as (b, v) inputs to
   quaidsSharesFit()/quaidsElasFit()/quaidsWelfareFit() with no expansion
   step, unlike quaidsRobustFit()'s own reduced basis. --- */
n = qOutBase.n;
nint = qOutBase.nint;
m_ = meanc(qOutBase.intcptFull~prices~totexp);
intcptPt = m_[1:1+nint];
pricesPt = m_[1+nint+1:1+nint+n];
totexpPt = m_[1+nint+n+1];

struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(rOut.b, rOut.v, intcptPt, pricesPt, totexpPt, aCtl);
call check(sumc(sharesOut.w .== sharesOut.w) == rows(sharesOut.w) and abs(sumc(sharesOut.w) - 1) < 1e-8,
    "rOut.b/rOut.v feed quaidsSharesFit() directly and produce adding-up shares");
call check(sumc(sharesOut.se .== sharesOut.se) == rows(sharesOut.se) and minc(sharesOut.se) >= 0,
    "quaidsSharesFit() propagated SE from rOut.v are finite and non-negative");

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(rOut.b, rOut.v, intcptPt, pricesPt, totexpPt, aCtl);
call check(sumc(elasOut.er .== elasOut.er) == rows(elasOut.er),
    "rOut.b/rOut.v feed quaidsElasFit() directly and produce finite income elasticities");

pricesPt1 = pricesPt;
pricesPt1[1] = pricesPt1[1] + ln(1.05);
struct quaidsWelfareOut welfareOut;
welfareOut = quaidsWelfareFit(rOut.b, rOut.v, intcptPt, pricesPt, pricesPt1, totexpPt, aCtl);
call check(welfareOut.seCV .== welfareOut.seCV and welfareOut.seCV >= 0,
    "rOut.b/rOut.v feed quaidsWelfareFit() directly and produce a finite, non-negative welfare SE");

/* --- Partial-failure handling: a pathological low-effN replicate is
   dropped via the effNMin pre-check rather than crashing the job. --- */
replicateWeightsMixed = ones(tobs, 3);
replicateWeightsMixed[1:500, 1] = 0.5*ones(500, 1);
replicateWeightsMixed[501:1000, 1] = 1.5*ones(500, 1);
replicateWeightsMixed[1:500, 2] = 1.5*ones(500, 1);
replicateWeightsMixed[501:1000, 2] = 0.5*ones(500, 1);
replicateWeightsMixed[., 3] = zeros(tobs, 1);
replicateWeightsMixed[1:5, 3] = 200*ones(5, 1);

struct quaidsReplicateOut rOutPartial;
rOutPartial = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, replicateWeightsMixed, 0.5, method="mixed");
call check(rOutPartial.nRequested == 3 and rOutPartial.nCompleted == 2 and rOutPartial.nFailed == 1,
    "a pathological low-effective-sample-size replicate is dropped, not crashed");
call check(sumc(sumc(rOutPartial.se .== rOutPartial.se)) == rows(rOutPartial.se)*cols(rOutPartial.se),
    "se remains finite after a partial-failure run");

call printQuaidsReplicateWeight(rOut);
call check(1, "printQuaidsReplicateWeight() runs without error");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "REPLICATE-WEIGHT TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "REPLICATE-WEIGHT TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
