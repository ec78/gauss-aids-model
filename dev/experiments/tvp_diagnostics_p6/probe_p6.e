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
struct _quaidsTVPDiagMultiOut dOut;
screen off;
dOut = _quaidsTVPDiagMultiStart(w, prices, totexp, H, q0, tvpCtl);
screen on;

priorLL = -6841.8654840494210|-6841.8655741471166|-6841.8654911245685;
priorQSum = 0.0013269734775710979|0.0013406784859562100|0.0013243386530964532;
qSums = sumc(dOut.qfits);
priorOk = maxc(abs(dOut.loglikelihoods - priorLL)) <= 1e-9 and
    maxc(abs(qSums - priorQSum)) <= 1e-12;
accept = qOut.converged == 1 and sumc(dOut.retcodes .== 0) == 3 and
    dOut.status == 0 and dOut.startStable == 1 and dOut.bestStart == 1 and
    sumc(dOut.nScaleCollapse .== 12) == 3 and
    dOut.maxLLGap <= dOut.tolerance and abs(dOut.tolerance - 5e-4) <= 1e-15 and
    priorOk;

format /re 24,16;
print "TVPDIAG_P6_CONFIG" n1 tobs seed tvpCtl.maxIters tvpCtl.maxTime;
print "TVPDIAG_P6_STATIC_CONVERGED" qOut.converged;
print "TVPDIAG_P6_H" diag(H)';
print "TVPDIAG_P6_MULT" dOut.multipliers';
print "TVPDIAG_P6_RC" dOut.retcodes';
print "TVPDIAG_P6_LL" dOut.loglikelihoods';
print "TVPDIAG_P6_QSUM" qSums';
print "TVPDIAG_P6_COLLAPSE" dOut.nScaleCollapse';
print "TVPDIAG_P6_GAP_TOL" dOut.maxLLGap dOut.tolerance;
print "TVPDIAG_P6_STATUS" dOut.status dOut.bestStart dOut.startStable;
i = 1;
do while i <= rows(dOut.qfits);
    j = 1;
    do while j <= cols(dOut.qfits);
        print "TVPDIAG_P6_Q" i j dOut.qfits[i, j];
        j = j + 1;
    endo;
    i = i + 1;
endo;
print "TVPDIAG_P6_PRIOR_MATCH" priorOk;
print "TVPDIAG_P6_ACCEPT" accept;
print "TVPDIAG_P6_DONE";
