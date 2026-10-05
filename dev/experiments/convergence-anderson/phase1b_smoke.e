new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include direct_objective.src;

/* Scratch smoke test for the Phase 1b additions, before wiring the full
** audit script. Not part of the audit record itself. */

tobs = 500;
seed = 4242;
modelFlag = 1;  /* QUAIDS -- the harder case */

print "=== stationary fixture ===";
{ w, intcpt, prices, totexp, instr, trueParams } =
    _quaidsStationaryFixtureDGP(tobs, seed, modelFlag);

struct quaidsControl aCtlStone;
aCtlStone = quaidsControlCreate;
aCtlStone.linear = 1-modelFlag;
aCtlStone.homogenous = 1;
aCtlStone.maxiter = 1;
struct quaidsOut qOutStone;
qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);

nGoods = cols(prices);
nReduced = nGoods-1;
nIntFull = qOutStone.nint+1;
ng = qOutStone.ng;
u = qOutStone.u;
intcptFull = qOutStone.intcptFull;

bRawTrue = _directTrueRaw(trueParams, nIntFull, nReduced, nGoods);
scaleVec = _directObjScale(bRawTrue, nReduced);
thetaTrue = _directObjPack(bRawTrue, nReduced)./scaleVec;

struct quaidsControl aCtl;
aCtl = aCtlStone;
aCtl.maxiter = 0;

{ fTrue, feasTrue, residTrue } = _directObjCore(thetaTrue, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
print "feasible at truth:" feasTrue "objective (minimized) at truth:" fTrue;

gradFd = _directObjFdGrad(thetaTrue, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
print "||FD grad||_inf at truth:" maxc(abs(gradFd));

print;
print "=== semi-analytic gradient/Hessian cross-check (at truth) ===";
J = _directResidJacobianFD(thetaTrue, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
residCov = residTrue[.,1:nReduced]'residTrue[.,1:nReduced]/tobs;
Sinv = invpd(residCov);
gGN = _directGNGradient(residTrue, J, Sinv, tobs, nReduced);
print "||GN grad||_inf at truth:" maxc(abs(gGN));
print "max abs diff (FD grad vs GN grad) at truth:" maxc(abs(gradFd-gGN));

B = _directGNHessian(J, Sinv, tobs, nReduced);
print "GN Hessian: rows" rows(B) "cols" cols(B) "min eig" minc(eigh(B)) "max eig" maxc(eigh(B));

print;
print "=== cross-check at a PERTURBED (non-stationary) point ===";
rndseed 777;
thetaPert = thetaTrue + .2*rndn(rows(thetaTrue), 1);
gradFdP = _directObjFdGrad(thetaPert, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
{ fP, feasP, residP } = _directObjCore(thetaPert, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
print "feasible at perturbed point:" feasP "objective:" fP;
if feasP;
    JP = _directResidJacobianFD(thetaPert, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
    residCovP = residP[.,1:nReduced]'residP[.,1:nReduced]/tobs;
    SinvP = invpd(residCovP);
    gGNP = _directGNGradient(residP, JP, SinvP, tobs, nReduced);
    print "||FD grad||_inf:" maxc(abs(gradFdP)) "||GN grad||_inf:" maxc(abs(gGNP));
    print "max abs diff (FD vs GN) at perturbed point:" maxc(abs(gradFdP-gGNP));
endif;

print;
print "=== _directLM smoke run from a 5% perturbation ===";
rndseed 888;
bRawStart = bRawTrue[.,1:nReduced].*(1+.05*rndn(ng,nReduced));
thetaStart = _directObjPack(bRawStart, nReduced)./scaleVec;
{ thetaHat, fHat, retcode, gnorm, iterUsed, rejectedSteps } = _directLM(thetaStart, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 200, 1e-3, 1e-6, 0);
print "retcode" retcode "gnorm" gnorm "fHat" fHat "fTrue" fTrue;
bRawHat = _directObjUnpack(thetaHat.*scaleVec, ng, nReduced);
dummy = _directBlockErrors(bRawHat, bRawTrue, nIntFull, nReduced, aCtl.linear);

print;
print "=== SVD/rank/weak-direction diagnostics (Jacobian at truth) ===";
dummy = _directJacobianSvdReport(J, 1e-10);
{ leftVecs, singMat, rightVecs } = svd2(J);
singVals = diag(singMat);
dummy = _directWeakDirectionReport(rightVecs, singVals, nIntFull, nReduced, aCtl.linear, ng, 3);

print;
print "=== line profile: start -> truth ===";
dummy = _directLineProfile(thetaStart, thetaTrue, 7, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);

print;
print "phase1b_smoke.e: run complete.";
