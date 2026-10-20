/*
** conv_fit.template.e -- one study fit per job. run_conv.ps1 substitutes
** 1 (1 = known-truth DGP cell, 2 = P6 anchor cell), 2,
** 200, 1001, 1 (1 = oracle true H, 2 = static H),
** 10 (start multiplier on the shipped q0 rule), 500,
** 600 and 0.01, and runs the generated copy from this
** directory. See PROTOCOL.md.
*/
new;
library cmlmt, tsmt, sslib;
#include conv_includes.inc

mode = 1;
n1Cell = 2;
tobs = 200;
seed = 1001;
harm = 1;
qmult = 10;
capIters = 500;
capTime = 600;
qscale = 0.01;

kCell = 2*n1Cell + n1Cell*(n1Cell+1)/2;
if mode == 1;
    trueQ = _convTrueQ(n1Cell, qscale);
    trueH = _convTrueH(n1Cell);
    { w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, seed, _convState0(n1Cell), trueQ, trueH);
    if harm == 1;
        H = trueH;
    else;
        H = _convStaticH(w, prices, totexp, n1Cell);
    endif;
else;
    /* P6 anchor: examples/14_tvp_aids_estimation.e's data and H rule
       with the first n1+1 goods (dev/experiments/tvp_p6_scale). */
    { wFull, intcpt, pricesFull, totexp, instr } = quaidsExampleData(tobs, seed);
    wRaw = wFull[., 1:n1Cell+1];
    w = wRaw ./ sumc(wRaw');
    prices = pricesFull[., 1:n1Cell+1];
    aCtl = quaidsControlCreate();
    aCtl.linear = 1;
    aCtl.homogenous = 1;
    screen off;
    qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
    screen on;
    H = diagrv(eye(n1Cell), qOut.homogSse[1:n1Cell]/qOut.nobs);
    trueQ = miss(0, 0)*ones(kCell, 1);
    print "CONV_STATIC_CONVERGED" qOut.converged;
endif;
qShipped = 0.1*meanc(diag(H))*ones(kCell, 1);
q0 = qmult*qShipped;

format /re 24,16;
print "CONV_CONFIG mode=" mode " n1=" n1Cell " tobs=" tobs " seed=" seed " harm=" harm " qmult=" qmult " maxIters=" capIters " maxTime=" capTime " qscale=" qscale;
print "CONV_MINW" minc(minc(w));
print "CONV_MAXW" maxc(maxc(w));
print "CONV_H";
print diag(H)';

t0 = date;
#include conv_chain.inc
el = ethsec(t0, date)/100;

llKF = sumc(sOut.kfResults.loglikelihood);
llFit = _convLL(Zarr, Qfit, H, w[., 1:n1], n1);
if mode == 1;
    llTrue = _convLL(Zarr, trueQ, H, w[., 1:n1], n1);
    qRef = trueQ;
else;
    llTrue = miss(0, 0);
    qRef = qShipped;
endif;
llStart = _convLL(Zarr, q0, H, w[., 1:n1], n1);
tol = 1e-6*nobs;

print "CONV_ELAPSED_SEC" el;
print "CONV_RC" rc;
print "CONV_ITERS" sOut.mleResults.NumIterations;
print "CONV_CMLMT_MIN" sOut.mleResults.ElapsedTime;
print "CONV_FCT" sOut.mleResults.Fct;
print "CONV_LLKF" llKF;
print "CONV_LLFIT" llFit;
print "CONV_LLTRUE" llTrue;
print "CONV_LLSTART" llStart;
print "CONV_TOL" tol;
print "CONV_MAXGRAD" maxc(abs(sOut.mleResults.Gradient));
print "CONV_NBOUNDARY" sumc(Qfit .< 1e-3*qRef);
print "CONV_QFIT_SUM" sumc(Qfit);
print "CONV_Q0";
print q0';
print "CONV_QTRUE";
print trueQ';
print "CONV_QFIT";
print Qfit';
print "CONV_BFINAL";
print bFinal;
if mode == 1 and rc == 0 and llFit < llTrue - tol;
    print "CONV_FAILC3 1";
else;
    print "CONV_FAILC3 0";
endif;
print "CONV_DONE";
