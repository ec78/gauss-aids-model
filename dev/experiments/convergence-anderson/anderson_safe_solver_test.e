new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Isolated validation for the safeguarded Anderson additions. */

/* Rank-deficient least squares: col 3 = col 1 + col 2. */
diffMat = (1~2~3)
    |(0~1~1)
    |(1~3~4)
    |(0~0~0);
targetResid = 1|2|3|4;
{ mixCoef, rankUse } = _andersonSvdMix(targetResid, diffMat, 1e-10);
pinvCoef = pinv(diffMat)*targetResid;
normalEqErr = maxc(abs(diffMat'*(targetResid-diffMat*mixCoef)));
fitDiff = maxc(abs(diffMat*mixCoef-diffMat*pinvCoef));

print "rank-deficient SVD solve: rank" rankUse
    "normal-equation residual" normalEqErr "fit diff vs pinv" fitDiff;

if rankUse == 2 and normalEqErr < 1e-10 and fitDiff < 1e-10;
    print "PASS: rank-truncated SVD solve is consistent.";
else;
    print "FAIL: rank-truncated SVD solve is not trustworthy.";
endif;

/* Real QUAIDS smoke case through the exact safeguarded prototype path. */
seed = 7;
tobs = 3000;
aCtl = quaidsControlCreate;
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.err = .0001;
aCtl.homogenous = 1;

{ w, intcpt, prices, totexp, instr, trueParams } =
    _quaidsSyntheticDGP(tobs, seed, 1, 1);
qOutSafe = _quaidsFitAnderson(w, intcpt, prices, totexp, instr,
    aCtl, 0, 8, 1, 1, 1e-10, 10, 10);
recErrSafe = maxc(maxc(abs(qOutSafe.bS-trueParams)));

print "safe QUAIDS seed 7: converged" qOutSafe.converged
    "iterations" qOutSafe.iterations "fixed residual" qOutSafe.finalErr
    "recErr" recErrSafe;

if qOutSafe.converged and qOutSafe.finalErr <= aCtl.err and recErrSafe < 1;
    print "PASS: safeguarded real-data path converged to the known good basin.";
else;
    print "FAIL: safeguarded real-data path did not pass the smoke check.";
endif;
