/*
** 09_replicate_weights.e
**
** Many household-expenditure surveys ship pre-computed REPLICATE weight
** columns (jackknife, BRR, ...) instead of expecting you to implement
** your own resampling scheme. quaidsReplicateWeightFit() refits
** quaidsFit() once per caller-supplied replicate weight column and
** combines the results via the standard linearized replication-variance
** formula. This example builds a small, explicit JK1 (delete-one-
** cluster jackknife) design by hand, so every replicate weight's
** construction is visible rather than hidden inside a survey package.
** See docs/USAGE_GUIDE.md's "Replicate-Weight (Jackknife/BRR) Standard
** Errors" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 09_replicate_weights.e
*/

new;
library quaids;
#include example_data.src

/* A small number of clusters (here, "regions") keeps a delete-one-
   cluster design easy to follow -- 10 replicates, one per region. */
nRegions = 10;
{ w, intcpt, prices, totexp, instr, regionId } = quaidsExampleClusterData(1000, 204, nRegions);

baseWeight = ones(rows(w), 1);

/* ---------------------------------------------------------------------
** Build the JK1 design by hand: replicate r drops region r entirely
** (weight 0) and reweights every remaining region's rows by
** nRegions/(nRegions-1), so the replicate's total weight is unchanged.
** scaleFactorJK1 = (R-1)/R is the standard JK1 prescription for this
** design (see docs/USAGE_GUIDE.md).
** --------------------------------------------------------------------- */

replicateWeights = zeros(rows(w), nRegions);
r = 1;
do while r <= nRegions;
    keepThisRegion = 1 - (regionId .== r);
    replicateWeights[., r] = baseWeight .* keepThisRegion * (nRegions/(nRegions-1));
    r = r + 1;
endo;

scaleFactorJK1 = (nRegions-1)/nRegions;

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl,
    replicateWeights, scaleFactorJK1, weight=baseWeight, method="JK1");

call printQuaidsReplicateWeight(rOut);

print "";
print "Replicates completed:" rOut.nCompleted "  failed:" rOut.nFailed;
print "rOut.b/rOut.se are already in quaidsFit()'s own full coefficient";
print "basis -- unlike quaidsRobustFit(), no separate expansion step is";
print "needed before feeding them to quaidsSharesFit()/quaidsElasFit()/";
print "quaidsWelfareFit().";

/* quaidsReplicateOut has no intcptFull field (unlike quaidsOut), so the
   constant-prepended intercept point is built by hand: a leading 1 for
   the constant term, then the mean of the actual demographic shifter. */
pricesPt = meanc(prices);
totexpPt = meanc(totexp);
intcptPt = 1 | meanc(intcpt);

sharesOut = quaidsSharesFit(rOut.b, rOut.v, intcptPt, pricesPt, totexpPt, aCtl);
print "";
print "Predicted shares' replicate-weight SE at the sample mean:" sharesOut.se';
