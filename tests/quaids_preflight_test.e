/*
** quaids_preflight_test.e
**
** Milestone 23: validates quaidsPreflight() as a silent, estimator-free
** diagnostics layer. These tests check data/design flags and metrics, not
** econometric estimates.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_preflight_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidsdiagnostics.src

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

tobs = 80;
N = 3;
idx = seqa(1, 1, tobs);
s1 = 123;
s2 = 124;
s3 = 125;
s4 = 126;
s5 = 127;
s6 = 128;

w1 = .23 + .01*rndns(tobs, 1, s1);
w2 = .32 + .01*rndns(tobs, 1, s2);
w3 = 1 - w1 - w2;
w = w1~w2~w3;
intcpt = 2 + rndns(tobs, 1, s3);
prices = rndns(tobs, N, s4);
instr = 5 + .2*idx + rndns(tobs, 1, s5);
totexp = .85*instr + .01*rndns(tobs, 1, s6);
clusterId = ceil(idx/10);

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 50;
aCtl.homogenous = 1;

pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, 0);

call check(pOut.ok == 1 and pOut.nErrors == 0,
    "clean fixture passes preflight with no hard errors");
call check(pOut.dimensionsOk == 1 and pOut.finiteOk == 1 and pOut.designInvOk == 1,
    "clean fixture dimensions, finiteness, and design checks pass");
call check(pOut.shareAddOk == 1 and pOut.maxShareSumDev < 1e-12,
    "clean fixture share adding-up is recognized");
call check(pOut.zeroShareCount == 0 and pOut.negativeShareCount == 0,
    "clean fixture has no zero or negative shares");
call check(pOut.ivDiagnosticsValid == 1 and pOut.weakIV == 0 and pOut.ivFstat > 10,
    "clean fixture first-stage diagnostics are valid and strong");
call check(pOut.clusterValid == 1 and pOut.nClusters == tobs and pOut.clusterWarning == 0,
    "clusterId=0 is treated as one observation per robust cluster");

call printQuaidsPreflight(pOut);

pCluster = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, clusterId, 0);
call check(pCluster.clusterValid == 1 and pCluster.nClusters == 8 and pCluster.minClusterSize == 10,
    "explicit cluster labels produce correct cluster counts");
call check(pCluster.clusterWarning == 1 and pCluster.ok == 1,
    "few-cluster condition is a warning, not a hard preflight error");

wZero = w;
wZero[1, 3] = wZero[1, 3] + wZero[1, 1];
wZero[1, 1] = 0;
pZero = quaidsPreflight(wZero, intcpt, prices, totexp, instr, aCtl, 0, 0);
call check(pZero.ok == 1 and pZero.zeroShareCount == 1 and pZero.nWarnings > 0,
    "zero shares are flagged as a warning while preserving ok status");

wNeg = w;
wNeg[1, 1] = -.10;
wNeg[1, 3] = 1 - wNeg[1, 1] - wNeg[1, 2];
pNeg = quaidsPreflight(wNeg, intcpt, prices, totexp, instr, aCtl, 0, 0);
call check(pNeg.ok == 0 and pNeg.negativeShareCount == 1 and pNeg.convergenceRisk == 2,
    "negative shares are a hard preflight error");

wAdd = w;
wAdd[1, 1] = wAdd[1, 1] + .01;
pAdd = quaidsPreflight(wAdd, intcpt, prices, totexp, instr, aCtl, 0, 0);
call check(pAdd.ok == 0 and pAdd.shareAddOk == 0 and pAdd.maxShareSumDev > .009,
    "share adding-up violations are hard preflight errors");

pricesLow = prices;
pricesLow[., 2] = ones(tobs, 1);
pLow = quaidsPreflight(w, intcpt, pricesLow, totexp, instr, aCtl, 0, 0);
call check(pLow.lowPriceVariation == 1 and pLow.nWarnings > 0,
    "low price variation is flagged as a warning");

pDim = quaidsPreflight(w, intcpt, prices[1:tobs-1, .], totexp, instr, aCtl, 0, 0);
call check(pDim.ok == 0 and pDim.dimensionsOk == 0 and pDim.nErrors == 1,
    "dimension mismatch returns a diagnostic struct instead of fitting");

/* Milestone 26: weight validation -- a genuinely unequal, valid weight is
   accepted (weightValid=1) and its Kish's effective sample size is
   correctly reduced below the raw weight sum; a weight with a non-finite
   entry is a hard preflight error, the same tier as clusterId. */
rndseed 88;
wgtUnequal = 1 + rndu(tobs, 1)*4;
pWeighted = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, wgtUnequal);
call check(pWeighted.weightValid == 1 and pWeighted.weightSum == sumc(wgtUnequal),
    "a genuinely unequal valid weight vector is accepted and its sum is reported");
call check(pWeighted.effN < pWeighted.weightSum,
    "unequal weighting reduces Kish's effective sample size below the raw weight sum");

wgtNonFinite = ones(tobs, 1);
wgtNonFinite[1] = miss(1, 1);
pWeightBad = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0, wgtNonFinite);
call check(pWeightBad.weightValid == 0 and pWeightBad.ok == 0,
    "a non-finite weight entry is a hard preflight error");

print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "PREFLIGHT TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "PREFLIGHT TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
