/*
** quaidsdirect_test.e
**
** Focused tests for the PRIVATE AIDS direct-LM prototype
** (dev/prototypes/aids_direct_lm/quaidsdirect.src; design: dev/AIDS_DIRECT_LM_DESIGN.md). Needs no
** optional package and nothing under dev/. The cross-implementation
** checks against the frozen research code and the read-only E2
** reproduction are in quaidsdirect_evidence_test.e.
**
** Checks, in order:
**   1. Exact recovery on a constructed fixture whose truth is an exactly
**      stationary point of the objective (equation noise projected off
**      every Jacobian column at the truth), plus near-recovery from the
**      production starts.
**   2. Noisy synthetic recovery (_quaidsSyntheticDGP, seeds outside the
**      research seeds 1-200) at the protocol's pfErr <= 1 cutoff, and an
**      independent absolute-price recomputation of the returned
**      criterion.
**   3. Selection rules (_quaidsDirectSelect) and the end-to-end
**      no-result / tie / ineligible paths via injected starts.
**   4. Determinism and isolation of the caller's RNG stream.
**   5. Default-estimator compatibility: quaidsControlCreate() defaults
**      and quaidsFit() output are unchanged by including and calling
**      the prototype.
** Input guards are separate scripts in dev/prototypes/aids_direct_lm/guard_error_cases/.
**
** Run with the repository's tests/ directory as the working directory
** (its includes and quaidsfixtures.src resolve from there):
**   tgauss -b -x ../dev/prototypes/aids_direct_lm/quaidsdirect_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
#include quaidsfixtures.src;

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

/* independent criterion: absolute-price form, from the full recovered
   (ng+1) x n matrix. Uses p'Gp with G's zero row/column sums instead of
   the solver's relative-price index. */
proc (1) = critFromFull(bFull, w, intcptFull, prices, totexp, u, alpha0);
    local nif, n, n1, alphaT, G, a_p, lx, fitted, e, S;
    nif = cols(intcptFull);
    n = cols(prices);
    n1 = n - 1;
    alphaT = intcptFull*bFull[1:nif, .];
    G = bFull[nif+1:nif+n, .];
    a_p = alpha0 + sumc((prices.*alphaT)') + .5*sumc(((prices*G).*prices)');
    lx = totexp - a_p;
    fitted = alphaT + prices*G + lx*bFull[nif+n+1, .] + u*bFull[nif+n+2, .];
    e = w[., 1:n1] - fitted[., 1:n1];
    S = e'e/rows(w);
    retp(-ln(det(S)));
endp;

struct quaidsControl aCtl, ctlStone, ctlFP;
struct _quaidsDirectOut dStar, dOut, dNone, dTie, dA, dB, dC;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;

/* ==========================================================================
   1. Exact recovery at a constructed, exactly stationary truth
   ========================================================================== */

{ w0, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, 601, 0, 1);
ctlStone = aCtl;
ctlStone.maxiter = 1;
stone = quaidsFit(w0, intcpt, prices, totexp, instr, ctlStone);
intcptFull = stone.intcptFull;
u = stone.u;
ng = stone.ng;
nGoods = cols(prices);
n1 = nGoods - 1;
nIntFull = cols(intcptFull);
tobs = rows(w0);

bStar = trueRawOf(tp, nIntFull, n1) ~ zeros(ng, 1);
scaleStar = _quaidsDirectScale(bStar, n1);
thetaStar = vec(bStar[., 1:n1])./scaleStar;

/* fitted shares at the truth, and the residual Jacobian there (both
   independent of w) */
fitted = -_quaidsDirectResid(zeros(tobs, nGoods), bStar, intcptFull, prices, totexp, u, aCtl);
Jstar = _quaidsDirectJacobian(thetaStar, 1e-6, zeros(tobs, nGoods), intcptFull, prices, totexp, u, aCtl, ng, n1, scaleStar);

/* equation noise orthogonal to every equation block of the Jacobian, so
   every term of the Sinv-weighted gradient vanishes at the truth */
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

{ thE, fE, rcE, gnE, itE, rjE, dgE } = _quaidsDirectLM(thetaStar, wStar,
    intcptFull, prices, totexp, u, aCtl, ng, n1, scaleStar, 500, 1e-3, 1e-6);
print "  stationary fixture: gnorm at truth" gnE;
call check(rcE == 0 and itE == 1 and gnE < 1e-8,
    "1a: LM started at the constructed truth converges immediately (retcode 0, gnorm < 1e-8)");
call check(thE == thetaStar,
    "1b: LM started at the truth returns the truth exactly");

dStar = _quaidsDirectFit(wStar, intcpt, prices, totexp, instr, aCtl);
errStar = maxc(maxc(abs(dStar.bRaw[1:ng-1, .] - bStar[1:ng-1, 1:n1])));
print "  stationary fixture: selected start" dStar.selStart "max |bRaw - truth| (excl. u)" errStar;
call check(dStar.ok == 1 and errStar < 1e-6,
    "1c: production starts recover the stationary truth to < 1e-6");
critStar = -fE;
call check(dStar.homogCrit >= critStar - 1e-8,
    "1d: selected criterion is not below the criterion at the truth");

/* ==========================================================================
   2. Noisy synthetic recovery, and an independent criterion check
   ========================================================================== */

seeds = 501|502;
si = 1;
do while si <= rows(seeds);
    { w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, seeds[si], 0, 1);
    dOut = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
    pf = maxc(maxc(abs(dOut.bFull[1:rows(tp)-1, .] - tp[1:rows(tp)-1, .])));
    print "  seed" seeds[si] "selected start" dOut.selStart "eligible" dOut.nEligible
        "pfErr" pf "nearTie" dOut.nearTie;
    call check(dOut.ok == 1 and dOut.status == 0 and pf <= 1,
        "2a: seed " $+ ftocv(seeds[si], 1, 0) $+ " returns an estimate with pfErr <= 1");
    call check(rows(dOut.bFull) == dOut.ng+1 and cols(dOut.bFull) == dOut.n
        and rows(dOut.bRaw) == dOut.ng and cols(dOut.bRaw) == dOut.n1
        and rows(dOut.startTable) == 16 and cols(dOut.startTable) == 13,
        "2b: seed " $+ ftocv(seeds[si], 1, 0) $+ " output shapes");
    cIndep = critFromFull(dOut.bFull, w, dOut.intcptFull, prices, totexp, dOut.u, aCtl.alpha0);
    call check(abs(cIndep - dOut.homogCrit) < 1e-9,
        "2c: seed " $+ ftocv(seeds[si], 1, 0) $+ " criterion matches an independent absolute-price recomputation");
    call check(maxc(abs(sumc(dOut.bFull') - (1|zeros(rows(dOut.bFull)-1, 1)))) < 1e-12
        and maxc(abs(sumc(dOut.bFull[dOut.nIntFull+1:dOut.nIntFull+dOut.n, .]))) < 1e-12,
        "2d: seed " $+ ftocv(seeds[si], 1, 0) $+ " recovery satisfies adding-up and homogeneity");
    elig = dOut.startTable[., 3];
    call check(dOut.nEligible == sumc(elig) and elig[dOut.selStart] == 1
        and dOut.homogCrit == maxc(selif(dOut.startTable[., 4], elig)),
        "2e: seed " $+ ftocv(seeds[si], 1, 0) $+ " selected start is eligible with the largest criterion");
    /* descriptive only: the shipped fixed-point estimator on the same data */
    ctlFP = aCtl;
    ctlFP.maxiter = 100;
    qFP = quaidsFit(w, intcpt, prices, totexp, instr, ctlFP);
    print "  (descriptive) quaidsFit converged" qFP.converged "homogCrit" qFP.homogCrit
        "direct homogCrit" dOut.homogCrit;
    si = si + 1;
endo;

/* ==========================================================================
   3. Selection rules and the no-result path
   ========================================================================== */

rawsA = zeros(9, 2);
{ b, t, m, nt } = _quaidsDirectSelect(0|1|1, 10|5|6, rawsA, 3);
call check(b == 3, "3a: an ineligible start is never selected, even with the largest criterion");
{ b, t, m, nt } = _quaidsDirectSelect(1|1|1, 5|7|7, rawsA, 3);
call check(b == 2 and t == 1, "3b: exact ties go to the lower start index and are counted");
{ b, t, m, nt } = _quaidsDirectSelect(0|0|0, 5|7|7, rawsA, 3);
call check(b == 0, "3c: no eligible start gives no selection");
rawsB = rawsA;
rawsB[3, 1] = 1;           /* start 1's u row (row 3) differs: same basin */
rawsB[4, 1] = .06;         /* start 2 differs by .06 in a non-u row: other basin */
{ b, t, m, nt } = _quaidsDirectSelect(1|1|0, 5|4.995|0, rawsB, 3);
call check(b == 1 and abs(m - .005) < 1e-12 and nt == 1,
    "3d: basin margin uses non-u rows with tolerance .05; near-tie when margin < .01");
{ b, t, m, nt } = _quaidsDirectSelect(1|0|1, 5|0|4, rawsB, 3);
call check(b == 1 and scalmiss(m) and nt == 0,
    "3e: a difference only in the u row is the same basin (no margin)");

{ w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, 503, 0, 1);
stone = quaidsFit(w, intcpt, prices, totexp, instr, ctlStone);
bS1 = stone.homogB;
bad = miss(zeros(rows(bS1), cols(bS1)), 0);
dNone = _quaidsDirectFitStarts(w, intcpt, prices, totexp, instr, aCtl, bad|bad);
call check(dNone.ok == 0 and dNone.status == 1 and dNone.nEligible == 0
    and scalmiss(dNone.bRaw) and scalmiss(dNone.bFull) and scalmiss(dNone.homogCrit)
    and scalmiss(dNone.selStart) and dNone.nStarts == 2,
    "3f: no eligible start returns ok = 0 with missing estimates (no crash, no fallback)");
call check(sumc(dNone.startTable[., 3]) == 0 and minc(dNone.startTable[., 2]) > 0,
    "3g: no-result per-start table records the failures");
dTie = _quaidsDirectFitStarts(w, intcpt, prices, totexp, instr, aCtl, bad|bS1|bS1);
call check(dTie.ok == 1 and dTie.selStart == 2 and dTie.exactTies == 1
    and dTie.startTable[1, 3] == 0,
    "3h: injected [failed, Stone, Stone]: start 2 selected, exact tie counted, failed start ineligible");

/* ==========================================================================
   4. Determinism and RNG isolation
   ========================================================================== */

rndseed 424242;
r1 = rndn(5, 1);
rndseed 424242;
dA = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
r2 = rndn(5, 1);
call check(r1 == r2, "4a: the caller's rndseed/rndn stream is unaffected by a fit");
dB = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
call check(sameMat(dA.bRaw, dB.bRaw) and sameMat(dA.startTable, dB.startTable)
    and dA.selStart == dB.selStart,
    "4b: two fits on the same data are identical");
st1 = _quaidsDirectAllStarts(bS1, nIntFull, n1);
st2 = _quaidsDirectAllStarts(bS1, nIntFull, n1);
call check(st1 == st2 and rows(st1) == 16*rows(bS1) and st1[1:rows(bS1), .] == bS1,
    "4c: internal starts are reproducible, 16 of them, start 1 = Stone");

/* ==========================================================================
   5. Default-estimator compatibility
   ========================================================================== */

defCtl = quaidsControlCreate();
call check(defCtl.maxiter == 50 and defCtl.alpha0 == 0 and defCtl.linear == 0
    and defCtl.homogenous == 1 and defCtl.err == .0001 and defCtl.b0 == 0 and defCtl.relax == 1,
    "5a: quaidsControlCreate() defaults unchanged");
ctlFP = aCtl;
ctlFP.maxiter = 100;
q1 = quaidsFit(w, intcpt, prices, totexp, instr, ctlFP);
dC = _quaidsDirectFit(w, intcpt, prices, totexp, instr, ctlFP);
q2 = quaidsFit(w, intcpt, prices, totexp, instr, ctlFP);
call check(sameMat(q1.homogB, q2.homogB) and sameMat(q1.homogV, q2.homogV)
    and sameMat(q1.bS, q2.bS) and sameMat(q1.vS, q2.vS)
    and sameMat(q1.homogCrit, q2.homogCrit)
    and q1.iterations == q2.iterations and q1.converged == q2.converged,
    "5b: quaidsFit() output is identical before and after a direct fit");
call check(ctlFP.maxiter == 100 and ctlFP.linear == 1 and ctlFP.b0 == 0,
    "5c: the caller's control struct is not modified");

print;
if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
