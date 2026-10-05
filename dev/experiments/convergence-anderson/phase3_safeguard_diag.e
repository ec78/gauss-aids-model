new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;
#include direct_objective.src;

/*
** phase3_safeguard_diag.e
**
** Phase 3 hardening item 5: targeted validation that an invpd() failure
** inside _directLMHomotopyDiag() is trapped, counted, never reuses an
** old inverse, is never accepted as a successful LM step, and returns
** or retries with the documented status. Failures are INJECTED through
** the solver's own fault hook (faultMode 1/2), which hands invpd() a
** genuinely non-PD matrix (-I) so the real trap path runs; no safeguard
** is relaxed. Small, cheap fixture: iterated AIDS, _quaidsSyntheticDGP
** tobs=500 seed=1, Stone start.
*/

nfail = 0;
ncheck = 0;

proc (1) = _sgCheck(cond, label);
    if cond;
        print "  PASS: " $+ label;
    else;
        print "  FAIL: " $+ label;
    endif;
    retp(cond);
endp;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(500, 1, 0, 1);
struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 1;
aCtl.homogenous = 1;
aCtl.maxiter = 1;
struct quaidsOut qOutStone;
qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
nReduced = cols(prices)-1;
ng = qOutStone.ng;
u = qOutStone.u;
intcptFull = qOutStone.intcptFull;
bStone = qOutStone.homogB;
scaleVec = _directObjScale(bStone, nReduced);
theta0 = _directObjPack(bStone, nReduced)./scaleVec;
{ f0, feas0, resid0 } = _directObjCore(theta0, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);

/* ---- T0: baselines (no injection) ---- */
print "T0: baselines (faultMode 0)";
{ thB, fBase, rcB, gnB, itB, rjB, dgB } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3,
    1e-6, 0, 1, 0, 0, 0);
print "  full: retcode" rcB "iter" itB "rejected" rjB "gnorm" gnB;
print "  diag [cov step infeas rho accepted lambda]:" dgB';
{ thB1, fB1, rcB1, gnB1, itB1, rjB1, dgB1 } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 1, 1e-3,
    1e-6, 0, 1, 0, 0, 0);
print "  maxIter=1: retcode" rcB1 "iter" itB1 "accepted" dgB1[5];
{ thW, fW, rcW, gnW, itW, rjW } = _directLMHomotopy(theta0, w, intcptFull,
    prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3, 1e-6, 0, 1);
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(rcB == 0 and itB >= 3,
    "baseline converges and needs >= 3 iterations (so iteration-2 injection is meaningful)");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(maxc(abs(thW-thB)) == 0 and fW == fBase and rcW == rcB and itW == itB and rjW == rjB,
    "wrapper _directLMHomotopy == Diag(faultMode 0), bit-identical");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(sumc(dgB[1:2]) == 0 and rjB == sumc(dgB[2:4]),
    "no inversion failures without injection; rejectedSteps == step+infeas+rho counters");
print;

/* ---- T1: covariance inversion fails at iteration 2 (a valid Sinv
   from iteration 1 exists and must NOT be reused) ---- */
print "T1: faultMode 1, covariance invpd() fails at outer iteration 2";
{ th1, f1, rc1, gn1, it1, rj1, dg1 } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3,
    1e-6, 0, 1, 1, 2, 0);
print "  retcode" rc1 "iter" it1 "gnorm" gn1 "diag:" dg1';
{ fChk1, feasChk1, rChk1 } = _directObjCore(th1, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(dg1[1] == 1, "failure trapped (no crash) and covariance-failure counter == 1");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(rc1 == 5 and it1 == 2, "returns documented retcode 5 at iteration 2 (no further progress)");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(maxc(abs(th1-thB1)) == 0 and f1 == fB1,
    "returned theta/fVal == the iteration-1 accepted point exactly (no step taken with a stale inverse)");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(gn1 == -1, "gnorm reported -1 (not a stale gradient norm)");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(fChk1 == f1 and feasChk1 == 1, "returned fVal is the returned point's own objective");
print;

/* ---- T2: covariance inversion fails at iteration 1 ---- */
print "T2: faultMode 1, covariance invpd() fails at outer iteration 1";
{ th2, f2, rc2, gn2, it2, rj2, dg2 } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3,
    1e-6, 0, 1, 1, 1, 0);
print "  retcode" rc2 "iter" it2 "gnorm" gn2 "diag:" dg2';
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(rc2 == 5 and dg2[1] == 1 and dg2[5] == 0 and maxc(abs(th2-theta0)) == 0 and f2 == f0,
    "retcode 5, counter 1, zero accepted steps, theta/fVal unchanged from start");
print;

/* ---- T3: every damped-step inversion fails ---- */
print "T3: faultMode 2, EVERY damped-step invpd() fails";
{ th3, f3, rc3, gn3, it3, rj3, dg3 } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3,
    1e-6, 0, 1, 2, 0, 1e9);
print "  retcode" rc3 "iter" it3 "rejected" rj3 "diag:" dg3';
nExpect = ceil(ln(1e12/1e-3)/ln(4));
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(dg3[2] == nExpect and rj3 == nExpect,
    "every failure trapped and counted (expected " $+ ftocv(nExpect,1,0) $+ " = damping 1e-3 -> 1e12 by x4)");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(dg3[5] == 0 and maxc(abs(th3-theta0)) == 0 and f3 == f0,
    "no failed inversion was accepted as a step: theta/fVal unchanged, 0 accepted");
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(rc3 == 4 and dg3[3] == 0 and dg3[4] == 0,
    "returns documented retcode 4 (stalled), no objective trial was even evaluated");
print;

/* ---- T4: first 3 damped-step inversions fail, then the solver
   retries with larger damping and recovers ---- */
print "T4: faultMode 2, first 3 damped-step invpd() calls fail";
{ th4, f4, rc4, gn4, it4, rj4, dg4 } = _directLMHomotopyDiag(theta0, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-3,
    1e-6, 0, 1, 2, 0, 3);
print "  retcode" rc4 "iter" it4 "rejected" rj4 "gnorm" gn4 "diag:" dg4';
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(dg4[2] == 3 and rj4 == sumc(dg4[2:4]),
    "exactly 3 trapped step-inversion failures counted as rejected steps");
{ fChk4, feasChk4, rChk4 } = _directObjCore(th4, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
ncheck = ncheck+1;
nfail = nfail + 1-_sgCheck(rc4 == 0 and gn4 < 1e-3 and feasChk4 == 1 and fChk4 == f4,
    "solver retries with increased damping and reaches genuine convergence (retcode 0, gnorm < 1e-3, feasible, own objective)");
/* Informational, NOT gated: where the endpoint differs from the
   un-injected baseline. An earlier version of this script gated on
   |f4 - fBase| < 1e-6, which failed (|df| ~1.2e-3): both runs stop at
   the gradTol=1e-3 criterion, but the injected failures change the
   damping path, so they stop at different points. The breakdown below
   shows whether that difference is in the weakly identified u row
   (Phase 1b audit) or in the structural rows. */
bRaw4 = _directObjUnpack(th4.*scaleVec, ng, nReduced);
bRawBase = _directObjUnpack(thB.*scaleVec, ng, nReduced);
print "  info: f4 - fBase (signed):" f4-fBase " max|theta4 - thetaBase| (scaled):" maxc(abs(th4-thB));
print "  info: max|raw diff| structural rows (excl. u):"
    maxc(maxc(abs(bRaw4[1:ng-1,1:nReduced]-bRawBase[1:ng-1,1:nReduced])))
    " u row:" maxc(abs(bRaw4[ng,1:nReduced]-bRawBase[ng,1:nReduced])');
{ th4b, f4b, rc4b, gn4b, it4b, rj4b, dg4b } = _directLMHomotopyDiag(th4, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-6,
    1e-6, 0, 1, 0, 0, 0);
{ thBb, fBb, rcBb, gnBb, itBb, rjBb, dgBb } = _directLMHomotopyDiag(thB, w,
    intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec, 500, 1e-6,
    1e-6, 0, 1, 0, 0, 0);
print "  info: both endpoints re-polished at gradTol 1e-6 (no injection): retcodes" rc4b rcBb
    " |f diff|:" abs(f4b-fBb) " max|scaled theta diff|:" maxc(abs(th4b-thBb));
print;

if nfail == 0;
    print "ALL " $+ ftocv(ncheck,1,0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail,1,0) $+ " CHECKS FAILED";
endif;
