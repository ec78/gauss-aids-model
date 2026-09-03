/*
** 08_robust_standard_errors.e
**
** Every other covariance in this library rests on a pooled,
** homoskedastic sandwich. quaidsRobustFit() generalizes that to a
** heteroskedasticity-robust or cluster-robust sandwich for an
** already-fitted quaidsFit() result, given the raw sample. This example
** uses data with a genuine within-cluster-correlated shock (e.g.
** households sharing a regional shock) so cluster-robust standard
** errors are measurably, honestly larger than naive ones -- not just a
** formula that runs. See docs/USAGE_GUIDE.md's "Robust and
** Cluster-Robust Standard Errors" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 08_robust_standard_errors.e
*/

new;
library quaids;
#include example_data.src

nClusters = 40;
{ w, intcpt, prices, totexp, instr, clusterId } = quaidsExampleClusterData(1000, 204, nClusters);

aCtl = quaidsControlCreate();
aCtl.linear = 1;         // AIDS -- keeps this example's focus on the SE
                         // formula, not iteration/convergence details
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

if not qOut.converged;
    print "quaidsFit() did not converge; quaidsRobustFit() requires it.";
    end;
endif;

/* ---------------------------------------------------------------------
** 1. Heteroskedasticity-robust (clusterId omitted -- every observation
** is its own "cluster")
** --------------------------------------------------------------------- */

print "=== Heteroskedasticity-robust ===";
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl);
call printQuaidsRobust(rOut);

/* ---------------------------------------------------------------------
** 2. Cluster-robust, using the true cluster labels
** --------------------------------------------------------------------- */

print "";
print "=== Cluster-robust (" nClusters "clusters) ===";
rOutCluster = quaidsRobustFit(qOut, w, prices, totexp, aCtl, clusterId=clusterId);
call printQuaidsRobust(rOutCluster);

print "";
print "Cluster-robust se is larger than heteroskedasticity-robust se here";
print "because the data has a genuine shared within-cluster shock --";
print "ignoring it (treating every row as independent) understates the";
print "true sampling variability:";
print "mean se, robust:   " meanc(vec(rOut.se));
print "mean se, cluster:  " meanc(vec(rOutCluster.se));

/* ---------------------------------------------------------------------
** A real, documented caveat: quaidsRobustFit()'s closed-form sandwich
** uses a SIMPLIFIED bread (not quaidsFit()'s own nonlinear-feedback-
** corrected Jacobian), which makes it dramatically more conservative
** than qOut's own classical SE. Confirmed as an expected consequence of
** comparing a simple equation-by-equation sandwich against the full,
** cross-equation-efficient FGLS system, not a formula bug.
** --------------------------------------------------------------------- */

print "";
print "For comparison, qOut's own classical (non-robust) coefficient SE:";
print "mean se, classical:" meanc(vec(qOut.homogSE));
print "(often 10-100x smaller than the robust se above -- a documented,";
print "expected property of the simplified-bread design, not a bug --";
print "see docs/USAGE_GUIDE.md's own Robust SE section.)";

/* ---------------------------------------------------------------------
** 3. quaidsRobustCovariance(): expand into qOut.bestB's full basis
**
** rOut.se above only covers the n1 independently-estimated equations,
** in a reduced regressor-aligned basis -- to feed shares/elasticities/
** welfare, expand it into qOut.bestB's full basis first.
** --------------------------------------------------------------------- */

{ bFull, vFull } = quaidsRobustCovariance(qOut, rOutCluster, aCtl);

pricesPt = meanc(prices);
totexpPt = meanc(totexp);
intcptPt = meanc(qOut.intcptFull);

sharesRobust = quaidsSharesFit(bFull, vFull, intcptPt, pricesPt, totexpPt, aCtl);
print "";
print "Predicted shares' cluster-robust SE at the sample mean:" sharesRobust.se';

/* ---------------------------------------------------------------------
** 4. A cluster-aware bootstrap alternative
**
** Resamples whole clusters and refits quaidsFit() on each resample --
** typically closer to qOut's own SE than the closed-form sandwich above,
** since it resamples the ACTUAL, cross-equation-efficient estimator.
** A small B here keeps this example fast; see
** docs/USAGE_GUIDE.md#limitations for realistic B choices.
** --------------------------------------------------------------------- */

print "";
print "=== Cluster-aware bootstrap (B=30, for comparison) ===";
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 30, clusterId=clusterId, seed=42);
call printQuaidsRobustBootstrap(rbOut);
