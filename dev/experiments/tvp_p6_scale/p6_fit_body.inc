/*
** p6_fit_body.inc -- shared body of the TVP-AIDS P6 scale re-measurement.
** #include'd by run_p6.ps1-generated copies of p6_scale_fit.template.e, which set nSub, tobs, seed,
** capIters and capTime first. Mirrors examples/14_tvp_aids_estimation.e
** exactly (same data generator, H from a static homogeneous LA-AIDS fit,
** q0 = 0.1*mean(diag(H)), smooth = 1), changing only the number of goods.
** One quaidsTVPFit() call per job; the outer hard timeout is applied by
** run_p6.ps1.
*/

{ wFull, intcpt, pricesFull, totexp, instr } = quaidsExampleData(tobs, seed);

wRaw = wFull[., 1:nSub];
w = wRaw ./ sumc(wRaw');
prices = pricesFull[., 1:nSub];
goodNames = 0$+quaidsExampleGoodNames();
goodNames = goodNames[1:nSub];

n = cols(prices);
n1 = n - 1;
k_states = 2*n1 + n1*(n1+1)/2;

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

H = diagrv(eye(n1), qOut.homogSse[1:n1]/qOut.nobs);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.othnam = goodNames;
tvpCtl.maxIters = capIters;
tvpCtl.maxTime = capTime;

print "P6_CONFIG nSub=" nSub " n1=" n1 " k_states=" k_states " tobs=" tobs " seed=" seed " maxIters=" capIters " maxTime_min=" capTime;
print "P6_STATIC_CONVERGED" qOut.converged;

t0 = date;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
el = ethsec(t0, date)/100;

format /re 24,16;
print "P6_ELAPSED_SEC" el;
print "P6_MLE_RETCODE" tvOut.mleRetcode;
print "P6_QFIT";
print tvOut.Qfit';
print "P6_QFIT_SUM" sumc(tvOut.Qfit);
print "P6_BFINAL";
print tvOut.bFinal;
format /rd 10,4;
print "";
call printQuaidsTVP(tvOut);
print "P6_DONE";
