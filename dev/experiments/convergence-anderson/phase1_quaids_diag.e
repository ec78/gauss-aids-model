new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include direct_objective.src;

/* Scratch diagnostic: is the QUAIDS gradient sane (does a small step in
** -gradient direction actually decrease the objective), or is optmt's
** line search itself the problem? Not part of the Phase 1 record. */

tobs = 500;
seed = 4242;
modelFlag = 1;

{ w, intcpt, prices, totexp, instr, trueParams } =
    _quaidsNoiselessDGP(tobs, seed, modelFlag);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 1-modelFlag;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsControl aCtlStone;
aCtlStone = aCtl;
aCtlStone.maxiter = 1;
struct quaidsOut qOutStone;
qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);

nGoods = cols(prices);
nReduced = nGoods - 1;
nIntFull = qOutStone.nint + 1;
ng = qOutStone.ng;
u = qOutStone.u;
intcptFull = qOutStone.intcptFull;

bRawTrue = _directTrueRaw(trueParams, nIntFull, nReduced, nGoods);
scaleVec = _directObjScale(qOutStone.homogB, nReduced);

rndseed 900001;
bRawStart = bRawTrue[., 1:nReduced].*(1+.01*rndn(ng, nReduced));
thetaStart = _directObjPack(bRawStart, nReduced)./scaleVec;

indFn = 1|0|0;
struct PV p0;
p0 = pvCreate;
p0 = pvPacki(p0, thetaStart, "theta", 1);
struct modelResults mm0;
mm0 = _directObjFn(p0, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, indFn);
f0 = mm0.function;
print "f0 (at start)" f0;

grad = _directObjFdGrad(thetaStart, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
print "||grad||_inf" maxc(abs(grad)) "||grad||_2" sqrt(grad'grad);

alpha = 1e-8;
do while alpha <= 1;
    thetaTry = thetaStart - alpha*grad;
    struct PV pTry;
    pTry = pvCreate;
    pTry = pvPacki(pTry, thetaTry, "theta", 1);
    struct modelResults mmTry;
    mmTry = _directObjFn(pTry, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, indFn);
    print "alpha" alpha "f(theta - alpha*grad)" mmTry.function;
    alpha = alpha*10;
endo;

print "diag complete.";
