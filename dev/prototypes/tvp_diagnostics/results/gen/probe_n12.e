new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc
#include ../../experiments/tvp_conv_study/conv_lib.src

n1 = 2;
tobs = 200;
trueQ = _convTrueQ(n1, 0.01);
H = _convTrueH(n1);
q0 = 0.1*meanc(diag(H))*ones(rows(trueQ), 1);
{ w, prices, totexp, trueStateArr, lxTrue } =
    _quaidsTVPRawSyntheticDGP(tobs, 1001, _convState0(n1), trueQ, H);
tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 600;
struct _quaidsTVPDiagMultiOut dOut;
screen off;
dOut = _quaidsTVPDiagMultiStart(w, prices, totexp, H, q0, tvpCtl);
screen on;
format /re 24,16;
print "TVPDIAG_PROBE_MULT" dOut.multipliers';
print "TVPDIAG_PROBE_RC" dOut.retcodes';
print "TVPDIAG_PROBE_LL" dOut.loglikelihoods';
print "TVPDIAG_PROBE_COLLAPSE" dOut.nScaleCollapse';
print "TVPDIAG_PROBE_GAP_TOL" dOut.maxLLGap dOut.tolerance;
print "TVPDIAG_PROBE_STATUS" dOut.status dOut.startStable;
print "TVPDIAG_PROBE_DONE";
