/*
** s0_deriv.template.e -- derivative-correctness gate and per-evaluation
** profile for one dataset. run_s0.ps1 substitutes 2, 200,
** 1001, 1 (1 oracle true H, 2 static H) and 0.01.
** 0 = 1 instead uses the Stage 6 fixture (n1 = 2, tobs 120,
** seed 7, draft Q, true H) for pre-freeze calibration only.
**
** At each prespecified raw point p (q = p^2): the analytic total gradient
** of the objective CMLMT sees (sslib _tvpLFNImpl, gradient filter,
** positive_vars chain rule) against central finite differences of the
** plain-filter objective with step h_j = 1e-4*p_j (gate) and 1e-3*p_j
** (reported), and the two filters' total LL against each other. See
** PROTOCOL.md for the frozen gate.
*/
new;
library cmlmt, tsmt, sslib;
#include s0_includes.inc

fixture = 0;
n1Cell = 2;
tobs = 200;
seed = 1001;
harm = 1;
qscale = 0.01;

if fixture;
    n1Cell = 2;
    tobs = 120;
    seed = 7;
    harm = 1;
    qscale = 1;
endif;
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
alt = 1 + 99*(seqa(1, 1, kCell) % 2);

n1 = n1Cell;
n = n1 + 1;
nobs = rows(w);
pricesRel = prices[., 1:n1] - prices[., n];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
y = w[., 1:n1];

struct ssControl ssCtl;
ssCtl = _s0Ctl(Zarr, H, qShipped, n1, 500, 600, 1);

/* Prespecified points (natural units). */
qPts = qShipped~(10*qShipped)~(0.1*qShipped)~trueQ~(trueQ .* alt);
ptNames = "SHIP1" $| "SHIP10" $| "SHIP01" $| "TRUEQ" $| "ALTQ";

format /re 24,16;
print "S0D_CONFIG fixture=" fixture " n1=" n1Cell " tobs=" tobs " seed=" seed " harm=" harm " qscale=" qscale;
print "S0D_MINW" minc(minc(w));
allOk = 1;
ip = 1;
do while ip <= cols(qPts);
    p = sqrt(qPts[., ip]);
    { llP, dum } = _s0Eval(p, y, ssCtl, 0);
    { llG, gA } = _s0Eval(p, y, ssCtl, 1);
    gF4 = _s0FDGrad(p, y, ssCtl, 1e-4);
    gF3 = _s0FDGrad(p, y, ssCtl, 1e-3);
    gnorm = maxc(abs(gF4'));
    err4 = maxc(abs(gA - gF4)')/maxc(gnorm|1);
    err3 = maxc(abs(gA - gF3)')/maxc(gnorm|1);
    dll = abs(llG - llP)/maxc(abs(llP)|1);
    ok = err4 <= 1e-5 and dll <= 1e-12;
    allOk = allOk and ok;
    print "S0D_POINT " $+ ptNames[ip];
    print "S0D_LLPLAIN" llP;
    print "S0D_LLGRAD" llG;
    print "S0D_DLL_REL" dll;
    print "S0D_GNORM" gnorm;
    print "S0D_ERR4" err4;
    print "S0D_ERR3" err3;
    print "S0D_OK" ok;
    print "S0D_GA";
    print gA;
    print "S0D_GFD4";
    print gF4;
    ip = ip + 1;
endo;
print "S0D_ALLOK" allOk;

/* Per-evaluation profile at the shipped start (descriptive). */
p = sqrt(qShipped);
nr = 5;
t0 = hsec;
ir = 1;
do while ir <= nr;
    { llP, dum } = _s0Eval(p, y, ssCtl, 0);
    ir = ir + 1;
endo;
secPlain = (hsec - t0)/100/nr;
t0 = hsec;
ir = 1;
do while ir <= nr;
    { llG, gA } = _s0Eval(p, y, ssCtl, 1);
    ir = ir + 1;
endo;
secGrad = (hsec - t0)/100/nr;
struct tvpModel tvpProf;
tvpProf = ssCtl.tvp;
nu = 20;
t0 = hsec;
ir = 1;
do while ir <= nu;
    _s0CountQUpdate(&tvpProf, qShipped);
    ir = ir + 1;
endo;
secUpd = (hsec - t0)/100/nu;
print "S0D_SEC_PLAIN_EVAL" secPlain;
print "S0D_SEC_GRAD_EVAL" secGrad;
print "S0D_SEC_QUPDATE" secUpd;
print "S0D_GRAD_OVER_PLAIN" secGrad/secPlain;
print "S0D_DONE";
