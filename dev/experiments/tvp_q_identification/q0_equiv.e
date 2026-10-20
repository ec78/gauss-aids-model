/*
** q0_equiv.e -- profile-harness gate on the Stage 6 fixture (n1 = 2,
** tobs 120, seed 7, draft Q, H = 2.5e-5*eye(2), shipped start, caps 500
** iterations x 2 minutes). Requires: (1) _q0MLEFitReduced with no fixed
** element reproduces quaidsTVPFit()'s Qfit and retcode EXACTLY; (2) with
** element 1 fixed at a probe value, that element is returned exactly and
** CMLMT's Fct equals the total filter LL at the rebuilt full Q exactly;
** (3) the fixed-element fit's LL is not below its start's LL.
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

n1 = 2;
tobs = 120;
trueQ = _convTrueQ(n1, 1);
H = _convTrueH(n1);
k_states = rows(trueQ);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);
{ w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, 7, _convState0(n1), trueQ, H);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 2;
struct quaidsTVPOut tvOut;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);

n = cols(prices);
pricesRel = prices[., 1:n1] - prices[., n];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
y = w[., 1:n1];

struct ssOut rOut;
rOut = _q0MLEFitReduced(Zarr, H, y, q0, n1, zeros(k_states, 1), 500, 2);
d1 = maxc(abs(rOut.final_params - tvOut.Qfit));
rc1 = rOut.mleResults.retcode;

fm = zeros(k_states, 1);
fm[1] = 1;
qs = tvOut.Qfit;
qs[1] = trueQ[1];
qs = maxc((qs~(1e-12*trueQ))');
llStart = _convLL(Zarr, qs, H, y, n1);
rOut = _q0MLEFitReduced(Zarr, H, y, qs, n1, fm, 500, 2);
qp = _q0FullQ(rOut.final_params, qs, fm);
llp = _convLL(Zarr, qp, H, y, n1);
format /re 24,16;
print "Q0_EQ_QFIT_MAXDIFF" d1;
print "Q0_EQ_RC" tvOut.mleRetcode rc1;
print "Q0_EQ_FIXDIFF" abs(qp[1] - trueQ[1]);
print "Q0_EQ_FCT_MINUS_LL" rOut.mleResults.Fct - llp;
print "Q0_EQ_PROF_MINUS_START" llp - llStart;
eqOk = d1 == 0 and tvOut.mleRetcode == rc1 and qp[1] == trueQ[1] and rOut.mleResults.Fct == llp and llp >= llStart;
print "Q0_EQ_PASS" eqOk;
print "Q0_EQ_DONE";
