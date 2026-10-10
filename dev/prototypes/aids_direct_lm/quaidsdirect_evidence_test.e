/*
** quaidsdirect_evidence_test.e
**
** Cross-implementation checks for the PRIVATE AIDS direct-LM prototype
** (dev/prototypes/aids_direct_lm/quaidsdirect.src) against the FROZEN research code and the
** committed E2 evidence in dev/experiments/convergence-anderson/.
**
** READ-ONLY with respect to that directory: it #includes the frozen
** .src files and reads e2_out/cells_m0.csv and
** confirm_e2_out/selected.csv. It writes nothing anywhere.
** Requires the optmt package (direct_objective.src's comparison-only
** optmt wrapper uses its PV/modelResults types). Gated in
** run_source_tests.ps1 by -SkipDirectEvidence.
**
** Every comparison is EXACT (zero difference): the prototype is a copy of
** the frozen arithmetic, and the E2 CSVs store full double precision
** (%.17g from GAUSS, shortest round-trip repr from Python).
**
**   1. Arithmetic equivalence on E2 held-out seed 31 data: scale, unpack,
**      residual, objective, FD Jacobian, Gauss-Newton gradient and
**      Hessian at several points, full LM solves from several E2 starts,
**      and the deterministic start recipes.
**   2. E2 reproduction for AIDS seeds 31, 115, 200 (the E2 determinism
**      seeds) and 47, 90, 104, 123 (the four E2 near-tie seeds): feeding
**      E2's own 16 starts through the private start seam reproduces every
**      stored cold-LM cell (retcode, eligibility, critExact, iterations,
**      rejected steps, gnorm, diagnostic counters, raw endpoint) and the
**      L16 selection (start, criterion, ties, margin, near-tie flag,
**      eligible count, pfErr).
**
** Run with the repository's tests/ directory as the working directory
** (its includes and quaidsfixtures.src resolve from there):
**   tgauss -b -x ../dev/prototypes/aids_direct_lm/quaidsdirect_evidence_test.e
*/

new;
library optmt;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/experiments/convergence-anderson/quaidsfit_anderson_prototype.src;
#include quaidsfixtures.src;
#include ../dev/experiments/convergence-anderson/direct_objective.src;
#include ../dev/experiments/convergence-anderson/pilot_lib.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;

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

/* a CSV text field -> number; empty field -> missing */
proc (1) = fieldNum(s);
    if strlen(s) == 0;
        retp(miss(0, 0));
    endif;
    retp(stof(s));
endp;

evDir = "../dev/experiments/convergence-anderson/";

struct quaidsControl ctl, ctlStone;
struct quaidsOut stone;
struct _quaidsDirectOut dOut;

/* the E2 harness's own control settings */
ctl = quaidsControlCreate;
ctl.linear = 1;
ctl.maxiter = 100;
ctl.err = .0001;
ctl.homogenous = 1;
ctlStone = ctl;
ctlStone.maxiter = 1;

/* ==========================================================================
   1. Arithmetic equivalence with the frozen procs (seed 31)
   ========================================================================== */

{ w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, 31, 0, 1);
stone = quaidsFit(w, intcpt, prices, totexp, instr, ctlStone);
bStone = stone.homogB;
intcptFull = stone.intcptFull;
u = stone.u;
ng = stone.ng;
n1 = cols(prices) - 1;
nIntFull = cols(intcptFull);
nobs = rows(w);
sv = _directObjScale(bStone, n1);
call check(sameMat(sv, _quaidsDirectScale(bStone, n1)), "1a: scale vector identical");
call check(sameMat(_directObjUnpack(vec(bStone[., 1:n1]), ng, n1),
    _quaidsDirectUnpack(vec(bStone[., 1:n1]), ng, n1)), "1b: unpack identical");

startsE2 = _pilotAllStarts(bStone, nIntFull, n1, 0, 31, 16);
pts = 1|2|11|13;
pk = 1;
do while pk <= rows(pts);
    bP = startsE2[(pts[pk]-1)*ng+1:pts[pk]*ng, .];
    th = vec(bP[., 1:n1])./sv;
    lbl = "start " $+ ftocv(pts[pk], 1, 0);
    call check(sameMat(_directObjResidHomotopy(w, _directObjUnpack(th.*sv, ng, n1), intcptFull, prices, totexp, u, ctl, 1),
        _quaidsDirectResid(w, _quaidsDirectUnpack(th.*sv, ng, n1), intcptFull, prices, totexp, u, ctl)),
        "1c: residual identical at E2 " $+ lbl);
    { f1, fe1, r1 } = _directObjCoreHomotopy(th, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv, 1);
    { f2, fe2, r2 } = _quaidsDirectCore(th, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv);
    call check(sameMat(f1, f2) and fe1 == fe2 and sameMat(r1, r2),
        "1d: objective/feasibility/residual identical at E2 " $+ lbl);
    J1 = _directResidJacobianFDHomotopy(th, 1e-6, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv, 1);
    J2 = _quaidsDirectJacobian(th, 1e-6, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv);
    call check(sameMat(J1, J2), "1e: FD Jacobian identical at E2 " $+ lbl);
    Sinv = invpd(r1[., 1:n1]'r1[., 1:n1]/nobs);
    call check(sameMat(_directGNGradient(r1, J1, Sinv, nobs, n1), _quaidsDirectGNGradient(r1, J1, Sinv, nobs, n1))
        and sameMat(_directGNHessian(J1, Sinv, nobs, n1), _quaidsDirectGNHessian(J1, Sinv, nobs, n1)),
        "1f: GN gradient and Hessian identical at E2 " $+ lbl);
    { a1, a2, a3, a4, a5, a6, a7 } = _directLMHomotopyDiag(th, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv, 500, 1e-3, 1e-6, 0, 1, 0, 0, 0);
    { c1, c2, c3, c4, c5, c6, c7 } = _quaidsDirectLM(th, w, intcptFull, prices, totexp, u, ctl, ng, n1, sv, 500, 1e-3, 1e-6);
    call check(sameMat(a1, c1) and sameMat(a2, c2) and a3 == c3 and a4 == c4
        and a5 == c5 and a6 == c6 and sameMat(a7, c7),
        "1g: full LM solve identical (all 7 outputs) from E2 " $+ lbl);
    pk = pk + 1;
endo;

/* deterministic recipes are identical; noise recipes perturb the same rows */
sameDet = 1;
samePattern = 1;
rndseed 1;
ix = 1;
do while ix <= 16;
    bF = _pilotBenchmarkStart(bStone, ix, nIntFull, n1, 0);
    { bQ, stDummy } = _quaidsDirectStart(bStone, ix, nIntFull, n1, 1);
    if ix == 1 or ix == 11 or ix == 13 or ix == 16;
        sameDet = sameDet and sameMat(bF, bQ);
    else;
        samePattern = samePattern and
            sameMat(sumc(((bF ./= bStone))'), sumc(((bQ ./= bStone))'));
    endif;
    ix = ix + 1;
endo;
call check(sameDet, "1h: deterministic start recipes 1, 11, 13, 16 identical to the frozen recipes");
call check(samePattern, "1i: noise start recipes 2-10, 12, 14, 15 perturb the same rows as the frozen recipes");

/* ==========================================================================
   2. Read-only E2 reproduction (AIDS, cold LM, configuration L16)
   ========================================================================== */

cells = csvReadM(evDir $+ "e2_out/cells_m0.csv", 1);
sel = csvReadSA(evDir $+ "confirm_e2_out/selected.csv", 2);

seeds = 31|115|200|47|90|104|123;
si = 1;
do while si <= rows(seeds);
    sd = seeds[si];
    lbl = "seed " $+ ftocv(sd, 1, 0);
    { w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(3000, sd, 0, 1);
    stone = quaidsFit(w, intcpt, prices, totexp, instr, ctlStone);
    bStone = stone.homogB;
    ng = stone.ng;
    nIntFull = cols(stone.intcptFull);
    startsE2 = _pilotAllStarts(bStone, nIntFull, n1, 0, sd, 16);
    dOut = _quaidsDirectFitStarts(w, intcpt, prices, totexp, instr, ctl, startsE2);

    cs = selif(cells, (cells[., 1] .== 0) .and (cells[., 2] .== sd) .and (cells[., 3] .== 2));
    cs = sortc(cs, 4);
    tbl = dOut.startTable;
    call check(rows(cs) == 16 and cs[., 4] == seqa(1, 1, 16), "2a: " $+ lbl $+ " has 16 stored cold-LM cells");
    call check(sameMat(cs[., 6], tbl[., 2]) and sameMat(cs[., 5], tbl[., 3]),
        "2b: " $+ lbl $+ " retcode and eligibility reproduce every cell");
    call check(sameMat(cs[., 8], tbl[., 4]), "2c: " $+ lbl $+ " critExact reproduces every cell exactly");
    call check(sameMat(cs[., 21], tbl[., 6]) and sameMat(cs[., 22], tbl[., 7])
        and sameMat(cs[., 24], tbl[., 5]) and sameMat(cs[., 27:32], tbl[., 8:13]),
        "2d: " $+ lbl $+ " iterations, rejected steps, gnorm and LM counters reproduce every cell");
    rawOk = 1;
    s = 1;
    do while s <= 16;
        rawOk = rawOk and sameMat(cs[s, 59:58+ng*n1]', vec(dOut.startRaw[(s-1)*ng+1:s*ng, .]));
        s = s + 1;
    endo;
    call check(rawOk, "2e: " $+ lbl $+ " every raw endpoint reproduces exactly");

    r = 1;
    srow = 0;
    do while r <= rows(sel);
        if sel[r, 1] $== "0" and sel[r, 2] $== ftocv(sd, 1, 0) and sel[r, 3] $== "L16";
            srow = r;
        endif;
        r = r + 1;
    endo;
    if srow == 0;
        call check(0, "2f: " $+ lbl $+ " L16 row found in selected.csv");
    else;
        pf = maxc(maxc(abs(dOut.bFull[1:rows(tp)-1, .] - tp[1:rows(tp)-1, .])));
        print "  " $+ lbl $+ ": selected start" dOut.selStart "nearTie" dOut.nearTie "pfErr" pf;
        call check(dOut.ok == 1 and dOut.selStart == fieldNum(sel[srow, 7])
            and dOut.nEligible == fieldNum(sel[srow, 5])
            and dOut.homogCrit == fieldNum(sel[srow, 8]),
            "2f: " $+ lbl $+ " selected start, eligible count and criterion match confirm_e2_out/selected.csv");
        call check(dOut.exactTies == fieldNum(sel[srow, 15])
            and sameMat(dOut.marginOtherBasin, fieldNum(sel[srow, 16]))
            and dOut.nearTie == (sel[srow, 17] $== "True"),
            "2g: " $+ lbl $+ " exact ties, basin margin and near-tie flag match");
        call check(pf == fieldNum(sel[srow, 9]) and pf <= 1,
            "2h: " $+ lbl $+ " pfErr matches the stored value and is <= 1");
    endif;
    si = si + 1;
endo;

print;
if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
