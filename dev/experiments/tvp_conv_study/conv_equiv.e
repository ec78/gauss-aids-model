/*
** conv_equiv.e -- harness-equivalence gate and Fct scaling calibration.
** Runs the Stage 6 fixture (tests/quaidstvpfit_test.e: n = 3, tobs =
** 120, seed 7, the draft truth Q, H = 2.5e-5*eye(2), q0 = 0.1*mean(diag(H)),
** caps 500 iterations x 2 minutes) through public quaidsTVPFit() and
** through conv_chain.inc, and requires Qfit, bFinal and the retcode to
** be EXACTLY equal (difference == 0). Also reports how CMLMT's stored Fct
** relates to the total filter log-likelihood at Qfit.
*/
new;
library cmlmt, tsmt, sslib;
#include conv_includes.inc

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

#include conv_chain.inc

llKF = sumc(sOut.kfResults.loglikelihood);
llFit = _convLL(Zarr, Qfit, H, w[., 1:n1], n1);
fct = sOut.mleResults.Fct;
format /re 24,16;
print "CONV_EQ_QFIT_MAXDIFF" maxc(abs(tvOut.Qfit - Qfit));
print "CONV_EQ_BFINAL_MAXDIFF" maxc(maxc(abs(tvOut.bFinal - bFinal)));
print "CONV_EQ_RC" tvOut.mleRetcode rc;
print "CONV_EQ_LLKF" llKF;
print "CONV_EQ_LLFIT" llFit;
print "CONV_EQ_FCT" fct;
print "CONV_EQ_NOBS" nobs;
print "CONV_EQ_FCT_OVER_SUM" fct/llKF;
print "CONV_EQ_FCT_TIMES_NOBS_MINUS_SUM" fct*nobs - llKF;
print "CONV_EQ_FCT_MINUS_SUM" fct - llKF;
eqOk = maxc(abs(tvOut.Qfit - Qfit)) == 0 and maxc(maxc(abs(tvOut.bFinal - bFinal))) == 0 and tvOut.mleRetcode == rc and llKF == llFit;
print "CONV_EQ_PASS" eqOk;
print "CONV_EQ_DONE";
