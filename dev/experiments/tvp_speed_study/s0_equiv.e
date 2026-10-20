/*
** s0_equiv.e -- harness-equivalence gate. On the Stage 6 fixture (as
** ../tvp_conv_study/conv_equiv.e: n = 3, tobs = 120, seed 7, draft truth
** Q, H = 2.5e-5*eye(2), shipped q0 rule, caps 500 iterations x 2
** minutes) the numerical arm of s0_lib.src's _s0Chain() (useGrad = 0) must
** reproduce public quaidsTVPFit()'s Qfit, bFinal and retcode EXACTLY
** (difference == 0). The analytic arm (useGrad = 1) is run on the same
** fixture and only reported here.
*/
new;
library cmlmt, tsmt, sslib;
#include s0_includes.inc

n1 = 2;
tobs = 120;
state0 = _convState0(n1);
trueQ = _convTrueQ(n1, 1);
H = _convTrueH(n1);
q0 = 0.1*meanc(diag(H))*ones(rows(trueQ), 1);
capIters = 500;
capTime = 2;
{ w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, 7, state0, trueQ, H);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = capIters;
tvpCtl.maxTime = capTime;
struct quaidsTVPOut tvOut;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);

Zarr = _s0Zarr(w, prices, totexp);
{ qN, rcN, itN, fctN, gN, bN, sN, cN, uN, dN } = _s0Chain(w, prices, totexp, H, q0, capIters, capTime, 0);
llN = _convLL(Zarr, qN, H, w[., 1:n1], n1);
{ Qfit, rc, itA, fctA, gA, bFinal, sA, cA, uA, dA } = _s0Chain(w, prices, totexp, H, q0, capIters, capTime, 1);
llA = _convLL(Zarr, Qfit, H, w[., 1:n1], n1);

format /re 24,16;
print "S0_EQ_QFIT_MAXDIFF" maxc(abs(tvOut.Qfit - qN));
print "S0_EQ_BFINAL_MAXDIFF" maxc(maxc(abs(tvOut.bFinal - bN)));
print "S0_EQ_RC" tvOut.mleRetcode rcN;
eqOk = maxc(abs(tvOut.Qfit - qN)) == 0 and maxc(maxc(abs(tvOut.bFinal - bN))) == 0 and tvOut.mleRetcode == rcN;
print "S0_EQ_PASS" eqOk;
print "S0_EQ_N_ITERS" itN;
print "S0_EQ_N_LL" llN;
print "S0_EQ_A_RC" rc;
print "S0_EQ_A_ITERS" itA;
print "S0_EQ_N_COUNTS" uN dN;
print "S0_EQ_A_COUNTS" uA dA;
print "S0_EQ_A_LL" llA;
print "S0_EQ_A_DLL" llA - llN;
print "S0_EQ_A_QFIT_MAXDIFF" maxc(abs(Qfit - qN));
print "S0_EQ_A_BFINAL_MAXDIFF" maxc(maxc(abs(bFinal - bN)));
print "S0_EQ_DONE";
