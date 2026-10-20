/*
** q0_prof.template.e -- TRUE (reoptimized) grouped profiles for one
** variance group at one endpoint. Group 2 (1 alpha, 2 gamma,
** 3 beta) is held fixed at the common level 10^c * qRef on the levels
** -12|-2|-1|0|1; every other Q element is reoptimized by _q0MLEFitReduced
** (same objective, transform and CMLMT settings as the shipped fit), from
** start a = qhat with the group replaced (the grouped-slice point) and,
** for levels listed in -12|-2|-1|0|1, also from start b = the shipped start
** with the group replaced. Caps 500 x 600 minutes.
** Dataset placeholders as q0_geom.template.e; C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/n12-t200-s1001-Horacle-qhat.fmt is the
** geometry job's saved qhat. See PROTOCOL.md.
**
** Tagged numeric rows:
**   Q0P_ROW  g log10c start rc iters llslice llprof llprof-llslice fixdiff nfloored sec
**   Q0P_Q    g log10c start q_1..q_k   (full reoptimized Q)
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

mode = 1;
n1Cell = 2;
tobs = 200;
seed = 1001;
harm = 1;
qscale = 0.01;
phi = -12;
grpId = 2;
levels = -12|-2|-1|0|1;
levelsB = -12|-2|-1|0|1;
capIters = 500;
capTime = 600;
#include q0_data.inc
#include q0_design.inc

format /re 24,16;
qhatfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/n12-t200-s1001-Horacle-qhat.fmt";
load qhat = ^qhatfile;
llhat = _convLL(Zarr, qhat, H, y, n1);
print "Q0P_CONFIG" mode n1Cell tobs seed harm qscale phi grpId capIters capTime;
print "Q0P_LLHAT" llhat;

qfl = (10^phi)*qRef;
gi = _q0GroupIdx(n1, grpId);
fm = zeros(k_states, 1);
fm[gi] = ones(rows(gi), 1);
struct ssOut rOut;
i = 1;
do while i <= rows(levels);
    c = levels[i];
    sb = 1;
    if sumc(abs(levelsB - c) .< 1e-12) > 0;
        sb = 2;
    endif;
    st = 1;
    do while st <= sb;
        if st == 1;
            qs = qhat;
        else;
            qs = qShipped;
        endif;
        qs[gi] = (10^c)*qRef[gi];
        llslice = _convLL(Zarr, qs, H, y, n1);
        nfl = sumc(qs .<= 0);
        qs = maxc((qs~qfl)');
        t0 = date;
        rOut = _q0MLEFitReduced(Zarr, H, y, qs, n1, fm, capIters, capTime);
        el = ethsec(t0, date)/100;
        qp = _q0FullQ(rOut.final_params, qs, fm);
        llp = _convLL(Zarr, qp, H, y, n1);
        print "Q0P_ROW" grpId c st rOut.mleResults.retcode rOut.mleResults.NumIterations llslice llp llp-llslice maxc(abs(qp[gi] - (10^c)*qRef[gi])) nfl el;
        print "Q0P_Q" grpId c st qp';
        st = st + 1;
    endo;
    i = i + 1;
endo;
print "Q0P_DONE";
