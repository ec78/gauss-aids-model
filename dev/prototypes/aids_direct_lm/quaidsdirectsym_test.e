/*
** quaidsdirectsym_test.e
**
** Implementation checks I1-I5, I7 and I8 of SYMMETRY_PROTOCOL.md v1
** (this directory) for the PRIVATE in-objective symmetry prototype
** quaidsdirectsym.src and its harness sym_harness.src. Fixture data only
** (protocol section 5: seeds 601 and 602); no development or
** confirmation seed is touched. I6 is test_sym_analyze.py; I9 and the
** byte-identity half of I5 are run_sym_checks.ps1.
**
** Run with the repository's tests/ directory as the working directory:
**   tgauss -b -x ../dev/prototypes/aids_direct_lm/quaidsdirectsym_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirectsym.src;
#include quaidsfixtures.src;
#include ../dev/prototypes/aids_direct_lm/sym_harness.src;

nfail = 0;
ncheck = 0;

proc (0) = check(cond, label);
    local i;
    i = ncheck + 1;
    ncheck = i;
    if cond;
        print "PASS  " $+ label;
    else;
        print "FAIL  " $+ label;
        nfail = nfail + 1;
    endif;
endp;

/* exact, missing-aware equality of two matrices */
proc (1) = sameMat(a, b);
    if rows(a) /= rows(b) or cols(a) /= cols(b);
        retp(0);
    endif;
    retp(missrv(a, -1e300) == missrv(b, -1e300));
endp;

/* reduced raw truth (ng x n1) from a trueParams (bS-layout) matrix */
proc (1) = trueRawOf(tp, nIntFull, n1);
    local n, r;
    n = n1 + 1;
    r = tp[1:nIntFull, .] | tp[nIntFull+1:nIntFull+n1, .] | tp[nIntFull+n+1:rows(tp), .];
    retp(r[., 1:n1]);
endp;

/* Verbatim copy of src/quaids.src:561-573 (the symmetry X), with
   nint = nIntFull - 1 and nendog + nu = 2 (AIDS rows beta and u). */
proc (1) = quaidsFitSymX(nIntFull, n1);
    local nint, nendog, nu, ng, ic, ig, i, R1, k, X;
    nint = nIntFull - 1;
    nendog = 1;
    nu = 1;
    ng = 1+nint+n1+nendog+nu;
        ic = seqa(1, 1, 1+nint)|seqa(1+nint+n1+1, 1, nendog+nu);
        ig = seqa(1+nint+1, 1, n1);
        ic = vec(ic+seqa(0, ng, n1)');
        ig = vec(ig+seqa(0, ng, n1)');
        i = ic|ig;
        i = sortind(i);

        R1 = design(vec(xpnd(seqa(1, 1, n1*(n1+1)/2))));
        /* vec(G) = R1*vech(G) (if G=G') */
        k = (ng-n1)*n1;
        X = ( eye(k)          ~zeros(k, cols(R1)) )|
            ( zeros(rows(R1), k)~     R1           );
        X = X[i, .];
    retp(X);
endp;

struct quaidsControl aCtl, ctlStone;
struct quaidsOut stone;
struct _quaidsDirectOut dS, dS2, dH, dH2, dNone, dTie;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
ctlStone = aCtl;
ctlStone.maxiter = 1;

/* ==========================================================================
   I1. Index map equals quaidsFit()'s symmetry X
   ========================================================================== */

cases = { 1 2, 1 3, 2 4 };
ci = 1;
do while ci <= rows(cases);
    nifC = cases[ci, 1];
    n1C = cases[ci, 2];
    ngC = nifC + n1C + 2;
    idxC = _quaidsDirectSymIndex(ngC, n1C);
    tC = design(idxC);
    xC = quaidsFitSymX(nifC, n1C);
    call check(sameMat(tC, xC) and cols(tC) == ngC*n1C - n1C*(n1C-1)/2
        and rows(tC) == ngC*n1C,
        "I1: design(idx) == quaidsFit X and kS = ng*n1 - n1(n1-1)/2 for (nIntFull, n1) = ("
        $+ ftocv(nifC, 1, 0) $+ ", " $+ ftocv(n1C, 1, 0) $+ ")");
    ci = ci + 1;
endo;

/* ==========================================================================
   I2. Round trips and bitwise symmetry
   ========================================================================== */

sd7001 = 7001;
sd7002 = 7002;
sd7003 = 7003;
sd7004 = 7004;
ng8 = 8;
n14 = 4;
idx8 = _quaidsDirectSymIndex(ng8, n14);
kS8 = maxc(idx8);
bR = rndns(ng8, n14, sd7001);
gS = bR[3:6, .];
bR[3:6, .] = (gS + gS')/2;
bR[3:6, .] = xpnd(vech(bR[3:6, .]));
vSym = vec(bR);
thR = _quaidsDirectSymProj(vSym, idx8);
call check(thR[idx8] == vSym, "I2a: for symmetric G, proj(v)[idx] == v exactly");
thRand = rndns(kS8, 1, sd7002);
rawRand = _quaidsDirectSymRaw(thRand, idx8, ng8, n14);
call check(_quaidsDirectSymAsym(rawRand, 2) == 0 and rows(rawRand) == ng8 and cols(rawRand) == n14,
    "I2b: a random theta reconstructs to a bitwise-symmetric G");
call check(_quaidsDirectSymProj(vec(rawRand), idx8) == thRand,
    "I2c: proj(theta[idx]) == theta exactly");
gA = rndns(4, 4, sd7003);
bA = zeros(ng8, n14);
bA[3:6, .] = gA;
thA = _quaidsDirectSymProj(vec(bA), idx8);
rawA = _quaidsDirectSymRaw(thA, idx8, ng8, n14);
call check(rawA[3:6, .] == (gA + gA')/2,
    "I2d: the projection of an asymmetric G is exactly (G+G')/2");

/* ==========================================================================
   I3. Exact recovery at the T3 stationary truth (fixture 601)
   ========================================================================== */

{ w0, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, 601, 0, 1);
stone = quaidsFit(w0, intcpt, prices, totexp, instr, ctlStone);
intcptFull = stone.intcptFull;
u = stone.u;
ng = stone.ng;
nGoods = cols(prices);
n1 = nGoods - 1;
nIntFull = cols(intcptFull);
tobs = rows(w0);

/* T3 construction, verbatim from quaidsdirect_test.e section 1 */
bStar = trueRawOf(tp, nIntFull, n1) ~ zeros(ng, 1);
scaleStar = _quaidsDirectScale(bStar, n1);
thetaStar = vec(bStar[., 1:n1])./scaleStar;
fitted = -_quaidsDirectResid(zeros(tobs, nGoods), bStar, intcptFull, prices, totexp, u, aCtl);
Jstar = _quaidsDirectJacobian(thetaStar, 1e-6, zeros(tobs, nGoods), intcptFull, prices, totexp, u, aCtl, ng, n1, scaleStar);
Zall = {};
jj = 1;
do while jj <= n1;
    Zall = Zall ~ Jstar[(jj-1)*tobs+1:jj*tobs, .];
    jj = jj + 1;
endo;
{ Uz, Sz, Vz } = svd2(Zall);
sz = diag(Sz);
Qbas = Uz[., 1:sumc(sz .> 1e-10*sz[1])];
seedE = 6011;
E = .05*rndns(tobs, n1, seedE);
E = E - Qbas*(Qbas'E);
wStar = fitted[., 1:n1] + E;
wStar = wStar ~ (1 - sumc(wStar'));

idx = _quaidsDirectSymIndex(ng, n1);
scaleSymStar = _quaidsDirectSymScale(scaleStar, idx);
call check(_quaidsDirectSymAsym(bStar[., 1:n1], nIntFull) == 0,
    "I3a: the fixture truth is symmetric (F2)");
thStarS = _quaidsDirectSymProj(vec(bStar[., 1:n1]), idx)./scaleSymStar;
{ thE, fE, rcE, gnE, itE, rjE, dgE } = _quaidsDirectSymLM(thStarS, wStar,
    intcptFull, prices, totexp, u, aCtl, ng, n1, scaleSymStar, idx, 500, 1e-3, 1e-6);
print "  stationary fixture: S gnorm at truth" gnE;
call check(rcE == 0 and itE == 1 and gnE < 1e-8,
    "I3b: S's LM started at the truth stops at iteration 1 (retcode 0, gnorm < 1e-8)");
call check(thE == thStarS, "I3c: S's LM started at the truth returns it exactly");
dS = _quaidsDirectSymFit(wStar, intcpt, prices, totexp, instr, aCtl);
errStar = maxc(maxc(abs(dS.bRaw[1:ng-1, .] - bStar[1:ng-1, 1:n1])));
print "  stationary fixture: S selected start" dS.selStart "max |bRaw - truth| (excl. u)" errStar;
call check(dS.ok == 1 and errStar < 1e-6,
    "I3d: _quaidsDirectSymFit's production starts recover the stationary truth to < 1e-6");

/* ==========================================================================
   I4. SymCore / Core equivalence and the reduced Jacobian (fixture 602)
   ========================================================================== */

{ w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, 602, 0, 1);
stone = quaidsFit(w, intcpt, prices, totexp, instr, ctlStone);
bStone = stone.homogB;
intcptFull = stone.intcptFull;
u = stone.u;
ng = stone.ng;
n1 = cols(prices) - 1;
nIntFull = cols(intcptFull);
idx = _quaidsDirectSymIndex(ng, n1);
tm = design(idx);
scaleVec = _quaidsDirectScale(bStone, n1);
scaleSym = _quaidsDirectSymScale(scaleVec, idx);
call check(rows(idx) == 32 and maxc(idx) == 26 and rows(scaleSym) == 26,
    "I4a: E2 layout: ng*n1 = 32 full, kS = 26 free coordinates");
thStone = _quaidsDirectSymProj(vec(bStone[., 1:n1]), idx)./scaleSym;
pts = thStone ~ (thStone + .1*rndns(rows(thStone), 1, sd7004));
ptI = 1;
do while ptI <= 2;
    thS = pts[., ptI];
    { fA, feasA, resA } = _quaidsDirectSymCore(thS, w, intcptFull, prices, totexp, u, aCtl, ng, n1, scaleSym, idx);
    thU = thS.*scaleSym;
    vFull = thU[idx];
    { fB, feasB, resB } = _quaidsDirectCore(vFull, w, intcptFull, prices, totexp, u, aCtl, ng, n1, ones(ng*n1, 1));
    call check(resA == resB and fA == fB,
        "I4b: point " $+ ftocv(ptI, 1, 0) $+ ": SymCore equals Core at the equivalent full point (unit scale), exactly");
    jSym = _quaidsDirectSymJacobian(thS, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, n1, scaleSym, idx);
    jV = _quaidsDirectJacobian(vFull, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, n1, ones(ng*n1, 1));
    jRef = (jV*tm).*scaleSym';
    relJ = maxc(maxc(abs(jSym - jRef)))/maxc(maxc(abs(jRef)));
    print "  point" ptI "J_sym relative error" relJ;
    call check(relJ <= 1e-7 and rows(jSym) == rows(w)*n1 and cols(jSym) == 26,
        "I4c: point " $+ ftocv(ptI, 1, 0) $+ ": J_sym == J_v*T*diag(scaleSym) to <= 1e-7 relative");
    ptI = ptI + 1;
endo;

/* paired starts (protocol 3.3) */
stS = _quaidsDirectSymStarts(bStone, nIntFull, n1, ng, idx, scaleSym, 0);
stF = _quaidsDirectAllStarts(bStone, nIntFull, n1);
okSt = rows(stS) == 16 and cols(stS) == 26;
s = 1;
do while s <= 16;
    okSt = okSt and stS[s, .]' == _quaidsDirectSymProj(vec(stF[(s-1)*ng+1:s*ng, 1:n1]), idx)./scaleSym;
    s = s + 1;
endo;
call check(okSt and stS[1, .]' == thStone,
    "I4d: S's 16 starts are the projected H starts (same draws); start 1 is the symmetrized Stone start");

/* ==========================================================================
   I5. H through the harness; S's Stone inputs equal H's
   ========================================================================== */

dH = _symFitH(w, intcpt, prices, totexp, instr, aCtl);
dH2 = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
call check(sameMat(dH.bRaw, dH2.bRaw) and sameMat(dH.startTable, dH2.startTable)
    and sameMat(dH.startRaw, dH2.startRaw) and dH.selStart == dH2.selStart
    and sameMat(dH.homogCrit, dH2.homogCrit),
    "I5a: H through the harness equals a direct _quaidsDirectFit call");
rndseed 424242;
r1 = rndn(5, 1);
rndseed 424242;
dS = _quaidsDirectSymFit(w, intcpt, prices, totexp, instr, aCtl);
r2 = rndn(5, 1);
call check(r1 == r2, "I5b: the caller's rndseed/rndn stream is unaffected by an S fit");
call check(sameMat(dS.bStone, dH.bStone) and sameMat(dS.u, dH.u)
    and sameMat(dS.intcptFull, dH.intcptFull)
    and sameMat(dS.scaleVec, _quaidsDirectSymScale(dH.scaleVec, idx)),
    "I5c: S's Stone fit, u, intercepts and scale derive bitwise from H's");

/* S output shape, exact symmetry, internal consistency, descriptive gap */
okSym = 1;
s = 1;
do while s <= dS.nStarts;
    okSym = okSym and (_quaidsDirectSymAsym(dS.startRaw[(s-1)*ng+1:s*ng, .], nIntFull) == 0
        or ismiss(dS.startRaw[(s-1)*ng+1:s*ng, .]));
    s = s + 1;
endo;
print "  seed 602: S selected start" dS.selStart "eligible" dS.nEligible "critS" dS.homogCrit
    "critH" dH.homogCrit "N1" (dS.homogCrit - dH.homogCrit);
call check(dS.ok == 1 and rows(dS.bRaw) == ng and cols(dS.bRaw) == n1
    and rows(dS.bFull) == ng+1 and cols(dS.bFull) == n1+1
    and rows(dS.startTable) == 16 and cols(dS.startTable) == 13
    and rows(dS.startRaw) == 16*ng and rows(dS.scaleVec) == 26,
    "I5d: S returns the unchanged _quaidsDirectOut layout (scaleVec holds scaleSym)");
call check(okSym and _quaidsDirectSymAsym(dS.bRaw, nIntFull) == 0,
    "I5e: every S endpoint has a bitwise-symmetric G");
cIndep = _symCritFromFull(dS.bFull, w, dS.intcptFull, prices, totexp, dS.u, aCtl.alpha0);
call check(abs(cIndep - dS.homogCrit) < 1e-9,
    "I5f: S's criterion matches the independent absolute-price recomputation (<= 1e-9)");
{ asymS, addS, homS } = _symStructErrs(dS.bRaw, nIntFull);
call check(addS <= 1e-12 and homS <= 1e-12,
    "I5g: S's bFull satisfies adding-up and homogeneity (<= 1e-12)");
elig = dS.startTable[., 3];
call check(dS.nEligible == sumc(elig) and elig[dS.selStart] == 1
    and dS.homogCrit == maxc(selif(dS.startTable[., 4], elig)),
    "I5h: S's selected start is eligible with the largest criterion");
call check(dS.homogCrit <= dH.homogCrit + .01,
    "I5i: nesting holds (critS <= critH + .01)");
dS2 = _quaidsDirectSymFit(w, intcpt, prices, totexp, instr, aCtl);
call check(sameMat(dS.bRaw, dS2.bRaw) and sameMat(dS.startTable, dS2.startTable)
    and sameMat(dS.startRaw, dS2.startRaw),
    "I5j: two S fits on the same data are identical");

/* ==========================================================================
   I7. P comparator
   ========================================================================== */

vH = vec(dH.bRaw);
B = _quaidsDirectSymPWeight(vH, w, dH.intcptFull, prices, totexp, dH.u, aCtl, ng, n1, dH.scaleVec);
call check(rows(B) == ng*n1 and cols(B) == ng*n1 and maxc(maxc(abs(B - B'))) <= 1e-12*maxc(maxc(abs(B))),
    "I7a: the unscaled weight B is (ng*n1) x (ng*n1) and symmetric");
thH = _quaidsDirectSymProj(vH, idx);
vHsym = thH[idx];
{ okI, vI } = _quaidsDirectSymPSolve(vHsym, B, idx);
errI = maxc(abs(vI - vHsym));
print "  P idempotence max |vP - v|" errI;
call check(okI == 1 and errI <= 1e-10, "I7b: P of an already-symmetric v returns v (<= 1e-10)");
{ okZ, vZ } = _quaidsDirectSymPSolve(vH, zeros(ng*n1, ng*n1), idx);
call check(okZ == 0 and scalmiss(vZ), "I7c: a non-positive-definite T'BT gives no result");
{ okM, vM } = _quaidsDirectSymPSolve(vH, miss(0, 0), idx);
call check(okM == 0, "I7d: a failed weight gives no result");
{ okP, pRaw, critP } = _quaidsDirectSymP(dH, w, prices, totexp, aCtl);
vP = vec(pRaw);
normEq = maxc(abs(tm'B*(vH - vP)))/maxc(abs(tm'B*vH));
print "  P: critP" critP "N2 = critP - critS" (critP - dS.homogCrit) "normal-equation residual" normEq;
call check(okP == 1 and rows(pRaw) == ng and cols(pRaw) == n1
    and _quaidsDirectSymAsym(pRaw, nIntFull) == 0,
    "I7e: P returns an ng x n1 endpoint with a bitwise-symmetric G");
call check(normEq <= 1e-8, "I7f: P satisfies the minimum-distance normal equations T'B(vH - vP) = 0");
{ fP, feasP, resP } = _quaidsDirectCore(vP./dH.scaleVec, w, dH.intcptFull, prices, totexp, dH.u, aCtl, ng, n1, dH.scaleVec);
call check(abs(-fP - critP) < 1e-10, "I7g: critP (unit scale) agrees with H-scaled evaluation");
dNone = dH;
dNone.ok = 0;
{ okN, pN, cN } = _quaidsDirectSymP(dNone, w, prices, totexp, aCtl);
call check(okN == 0 and scalmiss(pN) and scalmiss(cN), "I7h: P is no result when H has no result");

/* ==========================================================================
   I8. Selection on S through the injected-start seam
   ========================================================================== */

bad = miss(zeros(rows(bStone), cols(bStone)), 0);
dNone = _quaidsDirectSymFitStarts(w, intcpt, prices, totexp, instr, aCtl, bad|bad);
call check(dNone.ok == 0 and dNone.status == 1 and dNone.nEligible == 0
    and scalmiss(dNone.bRaw) and scalmiss(dNone.bFull) and scalmiss(dNone.homogCrit)
    and scalmiss(dNone.selStart) and dNone.nStarts == 2
    and sumc(dNone.startTable[., 3]) == 0 and minc(dNone.startTable[., 2]) > 0,
    "I8a: S with no eligible start returns ok = 0, missing estimates, no fallback");
dTie = _quaidsDirectSymFitStarts(w, intcpt, prices, totexp, instr, aCtl, bad|bStone|bStone);
call check(dTie.ok == 1 and dTie.selStart == 2 and dTie.exactTies == 1
    and dTie.startTable[1, 3] == 0,
    "I8b: S injected [failed, Stone, Stone]: ineligible never selected; tie to the lower start, counted");
call check(dTie.homogCrit == dS.startTable[1, 4] and sameMat(dTie.bRaw, dS.startRaw[1:ng, .]),
    "I8c: an injected Stone start reproduces S's internal start 1 exactly");

print;
if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
