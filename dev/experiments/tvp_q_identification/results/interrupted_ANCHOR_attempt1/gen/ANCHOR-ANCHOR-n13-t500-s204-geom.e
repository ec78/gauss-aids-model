/*
** q0_geom.template.e -- fixed-point likelihood geometry at one endpoint
** (no optimization). run_q0.ps1 substitutes the dataset placeholders
** (2, 3, 500, 204, 2, 0.01), the three
** start-fit files C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q1.fmt/C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q10.fmt/C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q0p1.fmt (starts 1x, 10x, 0.1x
** the shipped rule), C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-qhat.fmt (output), -12 (log10 floor),
** 0.01/0.03 (log-Q Hessian steps). See PROTOCOL.md.
**
** Every LL is the total exact-diffuse filter log-likelihood
** (_convLL, Phase C). "loss" = LL(qhat) - LL(point). Slices hold every
** other element at qhat: they are FIXED-PATH SLICES, not profiles.
**
** Tagged numeric rows (parsed by run_q0.ps1):
**   Q0G_START   s qmult ll ll-llbest nboundary maxdlog10_nonbnd
**   Q0G_ELEM    j group qhat qref qtrue qhat/qref bnd belowfloor
**   Q0G_POINT   code ll loss            (codes: PROTOCOL.md)
**   Q0G_SLICE   kind idx log10c ll loss (kinds: PROTOCOL.md)
**   Q0G_HEIG    step i eigval
**   Q0G_HVEC    step i v_1..v_k        (eigenvector i, ascending eigval)
**   Q0G_HDIAG   step j -H_jj g_j
**   Q0G_HINT    step nint mineig maxeig cond
**   Q0G_CORR    step i j c_ij          (interior block, i < j)
*/
new;
library cmlmt, tsmt, sslib;
#include q0_includes.inc

mode = 2;
n1Cell = 3;
tobs = 500;
seed = 204;
harm = 2;
qscale = 0.01;
phi = -12;
hsteps = 0.01|0.03;
#include q0_data.inc
#include q0_design.inc

format /re 24,16;
print "Q0G_CONFIG" mode n1Cell tobs seed harm qscale phi hsteps';

/* ---- start endpoints and qhat (best LL; ties -> earliest start) ---- */
qmults = 1|10|0.1;
qfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q1.fmt";
load qa = ^qfile;
qfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q10.fmt";
load qb = ^qfile;
qfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-q0p1.fmt";
load qc = ^qfile;
qs = qa~qb~qc;
lls = zeros(3, 1);
s = 1;
do while s <= 3;
    lls[s] = _convLL(Zarr, qs[., s], H, y, n1);
    s = s + 1;
endo;
ib = 1;
s = 2;
do while s <= 3;
    if lls[s] > lls[ib];
        ib = s;
    endif;
    s = s + 1;
endo;
qhat = qs[., ib];
llhat = lls[ib];
qhatfile = "C:/Users/eclow/Documents/GitHub/gauss-aids-model/dev/experiments/tvp_q_identification/results/work/ANCHOR-n13-t500-s204-qhat.fmt";
save ^qhatfile = qhat;
bnd = qhat .< 1e-3*qRef;
qfl = (10^phi)*qRef;
print "Q0G_BEST" ib llhat;
s = 1;
do while s <= 3;
    nb = sumc(qs[., s] .< 1e-3*qRef);
    nbAll = bnd .or (qs[., s] .< 1e-3*qRef);
    if sumc(nbAll .== 0) > 0;
        dl = maxc(abs(log(selif(qs[., s], nbAll .== 0)) - log(selif(qhat, nbAll .== 0))))/ln(10);
    else;
        dl = 0;
    endif;
    print "Q0G_START" s qmults[s] lls[s] lls[s]-llhat nb dl;
    s = s + 1;
endo;

/* ---- element table ---- */
grp = zeros(k_states, 1);
g = 1;
do while g <= 3;
    grp[_q0GroupIdx(n1, g)] = g*ones(rows(_q0GroupIdx(n1, g)), 1);
    g = g + 1;
endo;
j = 1;
do while j <= k_states;
    print "Q0G_ELEM" j grp[j] qhat[j] qRef[j] trueQ[j] qhat[j]/qRef[j] bnd[j] (qhat[j] < qfl[j]);
    j = j + 1;
endo;

/* ---- special points ---- */
print "Q0G_POINT" 1 llhat 0;
if mode == 1;
    v = _convLL(Zarr, trueQ, H, y, n1);
    print "Q0G_POINT" 2 v llhat-v;
endif;
v = _convLL(Zarr, qShipped, H, y, n1);
print "Q0G_POINT" 3 v llhat-v;
v = _convLL(Zarr, qfl, H, y, n1);
print "Q0G_POINT" 4 v llhat-v;
v = _convLL(Zarr, maxc((qhat~qfl)'), H, y, n1);
print "Q0G_POINT" 5 v llhat-v;
if sumc(bnd) > 0;
    qq = qhat;
    qq[selif(seqa(1, 1, k_states), bnd)] = qRef[selif(seqa(1, 1, k_states), bnd)];
    v = _convLL(Zarr, qq, H, y, n1);
    print "Q0G_POINT" 6 v llhat-v;
endif;

/* ---- fixed-path slices ---- */
cj = { -3, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 3 };
i = 1;
do while i <= rows(cj);
    v = _convLL(Zarr, (10^cj[i])*qhat, H, y, n1);
    print "Q0G_SLICE" 1 0 cj[i] v llhat-v;
    i = i + 1;
endo;
cr = phi|{ -4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2 };
i = 1;
do while i <= rows(cr);
    v = _convLL(Zarr, (10^cr[i])*qRef, H, y, n1);
    print "Q0G_SLICE" 2 0 cr[i] v llhat-v;
    i = i + 1;
endo;
j = 1;
do while j <= k_states;
    i = 1;
    do while i <= rows(cr);
        qq = qhat;
        qq[j] = (10^cr[i])*qRef[j];
        v = _convLL(Zarr, qq, H, y, n1);
        print "Q0G_SLICE" 3 j cr[i] v llhat-v;
        i = i + 1;
    endo;
    j = j + 1;
endo;
g = 1;
do while g <= 3;
    gi = _q0GroupIdx(n1, g);
    i = 1;
    do while i <= rows(cr);
        qq = qhat;
        qq[gi] = (10^cr[i])*qRef[gi];
        v = _convLL(Zarr, qq, H, y, n1);
        print "Q0G_SLICE" 4 g cr[i] v llhat-v;
        i = i + 1;
    endo;
    g = g + 1;
endo;

/* ---- local curvature in log-Q coordinates, two steps ---- */
t0 = date;
si = 1;
do while si <= 2;
    { f0, gr, Hm } = _q0HessLog(Zarr, qhat, H, y, n1, hsteps[si]);
    { va, ve } = eighv(-Hm);
    i = 1;
    do while i <= k_states;
        print "Q0G_HEIG" si i va[i];
        print "Q0G_HVEC" si i ve[., i]';
        print "Q0G_HDIAG" si i -Hm[i, i] gr[i];
        i = i + 1;
    endo;
    ii = selif(seqa(1, 1, k_states), bnd .== 0);
    if rows(ii) >= 1 and not scalmiss(ii);
        hi = -Hm[ii, ii];
        vi = eigh(hi);
        print "Q0G_HINT" si rows(ii) minc(vi) maxc(vi) maxc(vi)/minc(vi);
        if minc(diag(hi)) > 0;
            dd = 1 ./ sqrt(diag(hi));
            cc = dd .* hi .* dd';
            i = 1;
            do while i <= rows(ii);
                j = i + 1;
                do while j <= rows(ii);
                    print "Q0G_CORR" si ii[i] ii[j] cc[i, j];
                    j = j + 1;
                endo;
                i = i + 1;
            endo;
        endif;
    endif;
    si = si + 1;
endo;
print "Q0G_HESS_SEC" ethsec(t0, date)/100;
print "Q0G_DONE";
