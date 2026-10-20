/*
** s0_fit.template.e -- one Phase S0 pilot cell per job: the same Phase C
** known-truth dataset and start fitted by the numerical-gradient arm (N,
** shipped path) and the analytic-gradient arm (A), three timed
** repetitions each in the fixed order N A | A N | N A, after one untimed
** warm-up evaluation per arm. run_s0.ps1 substitutes 3, 200,
** 1001, 1 (1 oracle true H, 2 static H), 10,
** 500, 600 and 0.01. Data, H rule, start and caps
** are ../tvp_conv_study/conv_fit.template.e's (mode 1). See PROTOCOL.md.
*/
new;
library cmlmt, tsmt, sslib;
#include s0_includes.inc

n1Cell = 3;
tobs = 200;
seed = 1001;
harm = 1;
qmult = 10;
capIters = 500;
capTime = 600;
qscale = 0.01;

kCell = 2*n1Cell + n1Cell*(n1Cell+1)/2;
trueQ = _convTrueQ(n1Cell, qscale);
trueH = _convTrueH(n1Cell);
{ w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, seed, _convState0(n1Cell), trueQ, trueH);
if harm == 1;
    H = trueH;
else;
    H = _convStaticH(w, prices, totexp, n1Cell);
endif;
qShipped = 0.1*meanc(diag(H))*ones(kCell, 1);
q0 = qmult*qShipped;
n1 = n1Cell;
nobs = rows(w);
Zarr = _s0Zarr(w, prices, totexp);
tol = 1e-6*nobs;

format /re 24,16;
print "S0F_CONFIG n1=" n1Cell " tobs=" tobs " seed=" seed " harm=" harm " qmult=" qmult " maxIters=" capIters " maxTime=" capTime " qscale=" qscale;
print "S0F_MINW" minc(minc(w));
print "S0F_MAXW" maxc(maxc(w));

/* Warm-up: one untimed objective evaluation per arm at the start. */
struct ssControl ssW;
ssW = _s0Ctl(Zarr, H, q0, n1, capIters, capTime, 1);
{ llW, dumW } = _s0Eval(sqrt(q0), w[., 1:n1], ssW, 0);
{ llW, dumW } = _s0Eval(sqrt(q0), w[., 1:n1], ssW, 1);

/* Timed repetitions; arms[i] = 0 numerical (N), 1 analytic (A). */
arms = 0|1|1|0|0|1;
reps = 1|1|2|2|3|3;
nrun = rows(arms);
runTab = zeros(nrun, 12);
qN1 = 0; bN1 = 0; qA1 = 0; bA1 = 0;
detN = 1; detA = 1;
ir = 1;
do while ir <= nrun;
    { Qfit, rc, iters, fct, mxg, bFinal, mleSec, chainSec, nUpd, nGrad } = _s0Chain(w, prices, totexp, H, q0, capIters, capTime, arms[ir]);
    llFit = _convLL(Zarr, Qfit, H, w[., 1:n1], n1);
    runTab[ir, .] = arms[ir]~reps[ir]~rc~iters~nUpd~nGrad~mleSec~chainSec~fct~llFit~mxg~sumc(Qfit);
    print "S0F_RUN" runTab[ir, .];
    if reps[ir] == 1;
        if arms[ir] == 0;
            qN1 = Qfit; bN1 = bFinal; rcN = rc; itN = iters;
        else;
            qA1 = Qfit; bA1 = bFinal; rcA = rc; itA = iters;
        endif;
    else;
        if arms[ir] == 0;
            detN = detN and maxc(abs(Qfit - qN1)) == 0 and maxc(maxc(abs(bFinal - bN1))) == 0 and rc == rcN and iters == itN;
        else;
            detA = detA and maxc(abs(Qfit - qA1)) == 0 and maxc(maxc(abs(bFinal - bA1))) == 0 and rc == rcA and iters == itA;
        endif;
    endif;
    ir = ir + 1;
endo;

llN = _convLL(Zarr, qN1, H, w[., 1:n1], n1);
llA = _convLL(Zarr, qA1, H, w[., 1:n1], n1);
llTrue = _convLL(Zarr, trueQ, H, w[., 1:n1], n1);
dqAbs = maxc(abs(qA1 - qN1));
print "S0F_TOL" tol;
print "S0F_LLTRUE" llTrue;
print "S0F_LLN" llN;
print "S0F_LLA" llA;
print "S0F_DLL" llA - llN;
print "S0F_DQ_ABS" dqAbs;
print "S0F_DQ_REL" dqAbs/maxc(abs(qN1));
print "S0F_DB_ABS" maxc(maxc(abs(bA1 - bN1)));
print "S0F_DET_N" detN;
print "S0F_DET_A" detA;
print "S0F_QFIT_N";
print qN1';
print "S0F_QFIT_A";
print qA1';
print "S0F_BFINAL_N";
print bN1;
print "S0F_BFINAL_A";
print bA1;
print "S0F_DONE";
