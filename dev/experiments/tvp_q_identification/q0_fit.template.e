/*
** q0_fit.template.e -- one start fit per job, exactly as Phase C's
** conv_fit.template.e (same data, H rule, start rule, caps and private
** chain ../tvp_conv_study/conv_chain.inc), plus a binary save of Qfit for
** the geometry and profile jobs. run_q0.ps1 substitutes __MODE__,
** __N1__, __TOBS__, __SEED__, __HARM__, __QMULT__, __MAXITERS__,
** __MAXTIME__, __QSCALE__ and __QFILE__. See PROTOCOL.md.
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

mode = __MODE__;
n1Cell = __N1__;
tobs = __TOBS__;
seed = __SEED__;
harm = __HARM__;
qmult = __QMULT__;
capIters = __MAXITERS__;
capTime = __MAXTIME__;
qscale = __QSCALE__;
qfile = "__QFILE__";
#include q0_data.inc
q0 = qmult*qShipped;

format /re 24,16;
print "Q0F_CONFIG" mode n1Cell tobs seed harm qmult capIters capTime qscale;

t0 = date;
#include ../tvp_conv_study/conv_chain.inc
el = ethsec(t0, date)/100;

llFit = _convLL(Zarr, Qfit, H, w[., 1:n1], n1);
save ^qfile = Qfit;
load qBack = ^qfile;

print "Q0F_ELAPSED_SEC" el;
print "Q0F_RC" rc;
print "Q0F_ITERS" sOut.mleResults.NumIterations;
print "Q0F_CMLMT_MIN" sOut.mleResults.ElapsedTime;
print "Q0F_FCT" sOut.mleResults.Fct;
print "Q0F_LLFIT" llFit;
print "Q0F_MAXGRAD" maxc(abs(sOut.mleResults.Gradient));
print "Q0F_NBOUNDARY" sumc(Qfit .< 1e-3*qRef);
print "Q0F_QFIT_SUM" sumc(Qfit);
print "Q0F_SAVE_DIFF" maxc(abs(qBack - Qfit));
print "Q0F_QFIT";
print Qfit';
print "Q0F_DONE";
