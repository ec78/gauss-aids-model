/*
** q0_fit.template.e -- one start fit per job, exactly as Phase C's
** conv_fit.template.e (same data, H rule, start rule, caps and private
** chain ../tvp_conv_study/conv_chain.inc), plus a binary save of Qfit for
** the geometry and profile jobs. run_q0.ps1 substitutes 2,
** 3, 500, 204, 2, 0.1, 2000,
** 600, 0.01 and C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q0p1.fmt. See PROTOCOL.md.
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

mode = 2;
n1Cell = 3;
tobs = 500;
seed = 204;
harm = 2;
qmult = 0.1;
capIters = 2000;
capTime = 600;
qscale = 0.01;
qfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q0p1.fmt";
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
