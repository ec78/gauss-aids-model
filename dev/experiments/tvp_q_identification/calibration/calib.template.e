/*
** calib.template.e -- PRE-FREEZE calibration only (reserve / non-study
** data). Substitute __MODE__, __N1__, __TOBS__, __SEED__, __HARM__,
** __MAXITERS__ and run from the study directory. Checks: reduced-fit
** equivalence with no fixed element, save/load round trip, floor ladder,
** Hessian step ladder, one grouped profile point (two starts) and
** per-evaluation cost.
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

mode = __MODE__;
n1Cell = __N1__;
tobs = __TOBS__;
seed = __SEED__;
harm = __HARM__;
qscale = 0.01;
capIters = __MAXITERS__;
capTime = 600;
#include q0_data.inc
format /re 24,16;
print "CAL_CONFIG" mode n1Cell tobs seed harm;

q0 = qShipped;
t0 = date;
#include ../tvp_conv_study/conv_chain.inc
print "CAL_FIT_SEC" ethsec(t0, date)/100;
print "CAL_RC" rc;
print "CAL_ITERS" sOut.mleResults.NumIterations;
qhat = Qfit;
llhat = _convLL(Zarr, qhat, H, w[., 1:n1], n1);
print "CAL_LLHAT" llhat;
print "CAL_QHAT_OVER_QREF";
print (qhat ./ qRef)';

/* 1. reduced fit with no fixed element == shipped _quaidsTVPMLEFit (skipped when skipRed = 1) */
skipRed = __SKIPRED__;
y = w[., 1:n1];
struct ssOut rOut;
if skipRed == 0;
    t0 = date;
    rOut = _q0MLEFitReduced(Zarr, H, y, q0, n1, zeros(k_states, 1), capIters, capTime);
    print "CAL_RED_SEC" ethsec(t0, date)/100;
    print "CAL_RED_QDIFF" maxc(abs(rOut.final_params - Qfit));
    print "CAL_RED_RCITERS" rOut.mleResults.retcode rOut.mleResults.NumIterations;
endif;

/* 2. save/load round trip */
qfile = "results/_cal_qhat.fmt";
save ^qfile = qhat;
load qback = ^qfile;
print "CAL_SAVELOAD_DIFF" maxc(abs(qback - qhat));

/* 3. floor ladder: all elements at c*qRef; and boundary/each element */
ex = seqa(-2, -2, 7);
i = 1;
do while i <= rows(ex);
    print "CAL_FLOOR_ALL" ex[i] _convLL(Zarr, (10^ex[i])*qRef, H, y, n1);
    print "CAL_FLOOR_QHATMAX" ex[i] _convLL(Zarr, maxc((qhat~(10^ex[i])*qRef)'), H, y, n1) - llhat;
    i = i + 1;
endo;

/* 4. Hessian step ladder in log-Q coordinates */
hs = { 0.1, 0.03, 0.01, 0.003, 0.001 };
t0 = date;
i = 1;
do while i <= rows(hs);
    { f0, g, Hm } = _q0HessLog(Zarr, qhat, H, y, n1, hs[i]);
    va = eigh(-Hm);
    print "CAL_HESS_H" hs[i];
    print "CAL_HESS_G";
    print g';
    print "CAL_HESS_DIAG";
    print diag(-Hm)';
    print "CAL_HESS_EIG";
    print va';
    i = i + 1;
endo;
print "CAL_HESS_SEC_TOTAL" ethsec(t0, date)/100;

/* 5. one grouped profile point: gamma block at qRef level, two starts */
gi = _q0GroupIdx(n1, 2);
fm = zeros(k_states, 1);
fm[gi] = ones(rows(gi), 1);
qs1 = qhat;
qs1[gi] = qRef[gi];
llslice = _convLL(Zarr, qs1, H, y, n1);
qs2 = qShipped;
qs2[gi] = qRef[gi];
t0 = date;
rOut = _q0MLEFitReduced(Zarr, H, y, qs1, n1, fm, capIters, capTime);
qp1 = _q0FullQ(rOut.final_params, qs1, fm);
print "CAL_PROF1_SEC" ethsec(t0, date)/100;
print "CAL_PROF1_RCITERS" rOut.mleResults.retcode rOut.mleResults.NumIterations;
print "CAL_PROF1_FIXDIFF" maxc(abs(qp1[gi] - qRef[gi]));
print "CAL_PROF1_LL_MINUS_SLICE" _convLL(Zarr, qp1, H, y, n1) - llslice;
print "CAL_PROF1_FCT_MINUS_LL" rOut.mleResults.Fct - _convLL(Zarr, qp1, H, y, n1);
t0 = date;
rOut = _q0MLEFitReduced(Zarr, H, y, qs2, n1, fm, capIters, capTime);
qp2 = _q0FullQ(rOut.final_params, qs2, fm);
print "CAL_PROF2_SEC" ethsec(t0, date)/100;
print "CAL_PROF2_RCITERS" rOut.mleResults.retcode rOut.mleResults.NumIterations;
print "CAL_PROF2_LL_MINUS_PROF1" _convLL(Zarr, qp2, H, y, n1) - _convLL(Zarr, qp1, H, y, n1);
print "CAL_SLICE_LOSS" llhat - llslice;
print "CAL_PROF_LOSS" llhat - maxc(_convLL(Zarr, qp1, H, y, n1)|_convLL(Zarr, qp2, H, y, n1));

/* 6. per-evaluation cost */
t0 = date;
i = 1;
do while i <= 50;
    call _convLL(Zarr, qhat, H, y, n1);
    i = i + 1;
endo;
print "CAL_SEC_PER_EVAL" ethsec(t0, date)/100/50;
print "CAL_DONE";
