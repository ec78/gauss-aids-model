new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc

n1 = 3;
tobs = 500;
seed = 204;
{ wFull, intcpt, pricesFull, totexp, instr } =
    quaidsExampleData(tobs, seed);
wRaw = wFull[., 1:n1+1];
w = wRaw ./ sumc(wRaw');
prices = pricesFull[., 1:n1+1];

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
screen off;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
screen on;
H = diagrv(eye(n1), qOut.homogSse[1:n1]/qOut.nobs);
q0 = 0.1*meanc(diag(H))*ones(12, 1);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 2000;
tvpCtl.maxTime = 600;
struct quaidsTVPMultiStartOut dOut;
screen off;
dOut = quaidsTVPMultiStartDiagnostics(w, prices, totexp, H, q0, tvpCtl);
screen on;

priorLL = -6841.8654840494210|-6841.8655741471166|-6841.8654911245685;
priorQSum = 0.0013269734775710979|0.0013406784859562100|0.0013243386530964532;
qSums = sumc(dOut.qfits);
priorOk = maxc(abs(dOut.loglikelihoods - priorLL)) <= 1e-9 and
    maxc(abs(qSums - priorQSum)) <= 1e-12;
accept = qOut.converged == 1 and sumc(dOut.retcodes .== 0) == 3 and
    sumc(dOut.fitOk .== 1) == 3 and sumc(dOut.qValid .== 1) == 3 and
    sumc(dOut.llValid .== 1) == 3 and dOut.status == 5 and
    dOut.statusLabel $== "SCALE-CAUTION" and dOut.startStable == 1 and
    dOut.anyScaleCaution == 1 and dOut.bestStartIndex == 1 and
    abs(dOut.bestMultiplier - 0.1) <= 1e-15 and
    sumc(dOut.nScaleCollapse .== 12) == 3 and
    dOut.maxLLGap <= dOut.tolerance and abs(dOut.tolerance - 5e-4) <= 1e-15 and
    priorOk;

format /re 24,16;
print "TVPDIAG_D5_CONFIG" n1 tobs seed tvpCtl.maxIters tvpCtl.maxTime;
print "TVPDIAG_D5_STATIC_CONVERGED" qOut.converged;
print "TVPDIAG_D5_H" diag(H)';
print "TVPDIAG_D5_MULT" dOut.multipliers';
print "TVPDIAG_D5_RC" dOut.retcodes';
print "TVPDIAG_D5_FIT_Q_LL_VALID" dOut.fitOk' dOut.qValid' dOut.llValid';
print "TVPDIAG_D5_LL" dOut.loglikelihoods';
print "TVPDIAG_D5_QSUM" qSums';
print "TVPDIAG_D5_COLLAPSE" dOut.nScaleCollapse';
print "TVPDIAG_D5_GAP_TOL" dOut.maxLLGap dOut.tolerance;
print "TVPDIAG_D5_STATUS" dOut.status dOut.statusLabel;
print "TVPDIAG_D5_BEST" dOut.bestStartIndex dOut.bestMultiplier;
print "TVPDIAG_D5_STABILITY_CAUTION" dOut.startStable dOut.anyScaleCaution;
print "TVPDIAG_D5_PRIOR_MATCH" priorOk;
print "TVPDIAG_D5_ACCEPT" accept;
print "TVPDIAG_D5_DONE";
