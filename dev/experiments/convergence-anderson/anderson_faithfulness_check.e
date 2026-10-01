new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

seed = 7;
tobs = 3000;
q = 1;  // QUAIDS

aCtl = quaidsControlCreate;
aCtl.linear = 1 - q;
aCtl.maxiter = 100;
aCtl.err = .0001;
aCtl.homogenous = 1;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, q, 1);

qOutReal = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
qOutProto = _quaidsFitAnderson(w, intcpt, prices, totexp, instr, aCtl, 0, 0);   // mDepth=0

print "real quaidsFit():   converged" qOutReal.converged "iters" qOutReal.iterations "finalErr" qOutReal.finalErr;
print "prototype mDepth=0: converged" qOutProto.converged "iters" qOutProto.iterations "finalErr" qOutProto.finalErr;
print "max abs diff in bS:" maxc(maxc(abs(qOutReal.bS - qOutProto.bS)));
print "max abs diff in vS:" maxc(maxc(abs(qOutReal.vS - qOutProto.vS)));

if maxc(maxc(abs(qOutReal.bS - qOutProto.bS))) < 1e-10;
    print "FAITHFUL: prototype with mDepth=0 exactly reproduces real quaidsFit().";
else;
    print "MISMATCH: prototype does NOT reproduce real quaidsFit() -- do not trust the sweep comparison yet.";
endif;

print;
print "=== Anderson depth sweep on this same seed ===";
recErrBase = maxc(maxc(abs(qOutReal.bS - trueParams)));
print "baseline (relax=1, no Anderson): iters" qOutReal.iterations "converged" qOutReal.converged "recErr" recErrBase;

d = 2;
do while d <= 12;
    qOutA = _quaidsFitAnderson(w, intcpt, prices, totexp, instr, aCtl, 0, d);
    recErrA = maxc(maxc(abs(qOutA.bS - trueParams)));
    print "mDepth" d ": converged" qOutA.converged "iters" qOutA.iterations "finalErr" qOutA.finalErr "recErr" recErrA;
    d = d + 2;
endo;
