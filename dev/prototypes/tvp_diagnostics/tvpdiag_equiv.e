new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc
#include ../../experiments/tvp_conv_study/conv_lib.src

n1 = 2;
tobs = 120;
trueQ = _convTrueQ(n1, 1);
H = _convTrueH(n1);
q0 = 0.1*meanc(diag(H))*ones(rows(trueQ), 1);
{ w, prices, totexp, trueStateArr, lxTrue } =
    _quaidsTVPRawSyntheticDGP(tobs, 7, _convState0(n1), trueQ, H);
tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 2;
struct quaidsTVPOut tvOut;
struct _quaidsTVPDiagMultiOut dOut;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
dOut = _quaidsTVPDiagMultiStart(w, prices, totexp, H, q0, tvpCtl);

pricesRel = prices[., 1:n1] - prices[., n1+1];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n1+1);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
ll1 = _convLL(Zarr, tvOut.Qfit, H, w[., 1:n1], n1);
qdiff = maxc(abs(tvOut.Qfit - dOut.qfits[., 2]));
lldiff = abs(ll1 - dOut.loglikelihoods[2]);
pass = qdiff == 0 and tvOut.mleRetcode == dOut.retcodes[2]
    and lldiff == 0 and dOut.startStable == 1;
format /re 24,16;
print "TVPDIAG_EQ_QDIFF" qdiff;
print "TVPDIAG_EQ_LLDIFF" lldiff;
print "TVPDIAG_EQ_RC" tvOut.mleRetcode dOut.retcodes[2];
print "TVPDIAG_EQ_GAP_TOL" dOut.maxLLGap dOut.tolerance;
print "TVPDIAG_EQ_PASS" pass;
if not pass;
    end;
endif;
