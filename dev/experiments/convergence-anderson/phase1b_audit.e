new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include direct_objective.src;

/*
** phase1b_audit.e
**
** Phase 1b audit, requested after Phase 1's own QUAIDS recovery-gate
** failure, BEFORE proceeding to Phase 3 continuation/homotopy work. Does
** NOT conclude direct optimization cannot work or that QUAIDS is weakly
** identified -- it tests exactly that question, separating two
** confounds Phase 1 left entangled: (1) a fixture whose first-stage
** residual regressor was itself made weakly identified by the same
** noise-shrinking that kept the GLS criterion well-conditioned, and (2)
** an optimizer (generic optmt, plus a flat 1e8 infeasibility penalty)
** with no real safeguards. See SUMMARY.md's Phase 1b section for the
** full writeup; this file is the record of how each number there was
** produced.
**
** Item numbers below match the audit's own numbered work items.
*/

tobs = 500;
seed = 4242;

modelFlag = 0;
do while modelFlag <= 1;
    if modelFlag == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;
    print "================================================================";
    print "" $+ modelName $+ ": Phase 1b audit";
    print "================================================================";

    /* ---- Item 2: corrected fixture, built + validated ---- */
    { w, intcpt, prices, totexp, instr, trueParams } =
        _quaidsStationaryFixtureDGP(tobs, seed, modelFlag);

    struct quaidsControl aCtlStone;
    aCtlStone = quaidsControlCreate;
    aCtlStone.linear = 1-modelFlag;
    aCtlStone.homogenous = 1;
    aCtlStone.maxiter = 1;
    struct quaidsOut qOutStone;
    qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);

    nGoods = cols(prices);
    nReduced = nGoods-1;
    nIntFull = qOutStone.nint+1;
    ng = qOutStone.ng;
    u = qOutStone.u;
    intcptFull = qOutStone.intcptFull;

    bRawTrue = _directTrueRaw(trueParams, nIntFull, nReduced, nGoods);
    scaleVec = _directObjScale(bRawTrue, nReduced);
    thetaTrue = _directObjPack(bRawTrue, nReduced)./scaleVec;

    struct quaidsControl aCtl;
    aCtl = aCtlStone;
    aCtl.maxiter = 0;

    { fTrue, feasTrue, residTrue } = _directObjCore(thetaTrue, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
    gradTrue = _directObjFdGrad(thetaTrue, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
    print "Item 2: fixture validation";
    print "  feasible at truth:" feasTrue "  objective at truth:" fTrue;
    print "  ||FD gradient||_inf at truth (fixture is NEAR-, not exactly,";
    print "  stationary -- this is the actual observed gap, not ~0):"
        maxc(abs(gradTrue));
    print "  structural recovery at truth (should be exactly 0, sanity check):";
    dummy = _directBlockErrors(bRawTrue, bRawTrue, nIntFull, nReduced, aCtl.linear);
    print;

    /* ---- Item 3: local identification diagnostics at truth ---- */
    print "Item 3: local identification diagnostics at truth";
    Jtrue = _directResidJacobianFD(thetaTrue, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);

    residCovT = residTrue[.,1:nReduced]'residTrue[.,1:nReduced]/tobs;
    SinvT = invpd(residCovT);

    /* Phase 1b cleanup item 4: the semi-analytic (Sinv-weighted
       Gauss-Newton) gradient vs. pure finite-difference cross-check,
       promoted from the scratch phase1b_smoke.e into this canonical,
       saved-output script. */
    gGNTrue = _directGNGradient(residTrue, Jtrue, SinvT, tobs, nReduced);
    print "  semi-analytic (GN) vs. finite-difference gradient cross-check,";
    print "  at truth: max abs diff" maxc(abs(gradTrue-gGNTrue));

    rndseed 700000+modelFlag;
    thetaPertChk = thetaTrue + .2*rndn(rows(thetaTrue), 1);
    { fPertChk, feasPertChk, residPertChk } = _directObjCore(thetaPertChk, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
    if feasPertChk;
        gradPertChk = _directObjFdGrad(thetaPertChk, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
        JPertChk = _directResidJacobianFD(thetaPertChk, 1e-6, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
        residCovPertChk = residPertChk[.,1:nReduced]'residPertChk[.,1:nReduced]/tobs;
        SinvPertChk = invpd(residCovPertChk);
        gGNPertChk = _directGNGradient(residPertChk, JPertChk, SinvPertChk, tobs, nReduced);
        print "  same cross-check at a 0.2-scaled-unit perturbed (non-stationary)";
        print "  point: max abs diff" maxc(abs(gradPertChk-gGNPertChk));
    else;
        print "  (perturbed-point cross-check skipped: infeasible draw)";
    endif;
    print;

    dummy = _directJacobianSvdReport(Jtrue, 1e-10);
    { leftVecsT, singMatT, rightVecsT } = svd2(Jtrue);
    singValsT = diag(singMatT);
    dummy = _directWeakDirectionReport(rightVecsT, singValsT, nIntFull, nReduced, aCtl.linear, ng, 3);

    BT = _directGNHessian(Jtrue, SinvT, tobs, nReduced);
    eigBT = sortc(eigh(BT), 1);
    print "  GN-Hessian eigenvalues at truth: min" eigBT[1] "max"
        eigBT[rows(eigBT)] "condition number" eigBT[rows(eigBT)]/eigBT[1];
    print;

    /* ---- Items 4-5: fixed LM optimizer, local recovery gate ---- */
    print "Items 4-5: LM recovery gate, several 1%/5% perturbations";
    perturbIdx = 1;
    do while perturbIdx <= 2;
        if perturbIdx == 1;
            pScale = .01;
        else;
            pScale = .05;
        endif;

        trial = 1;
        do while trial <= 2;
            rndseed 500000+1000*modelFlag+100*perturbIdx+trial;
            bRawStart = bRawTrue[.,1:nReduced].*(1+pScale*rndn(ng, nReduced));
            thetaStart = _directObjPack(bRawStart, nReduced)./scaleVec;

            { fStart, feasStart, residStart } = _directObjCore(thetaStart, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
            print "  --- pScale" pScale "trial" trial "---";
            print "  start: feasible" feasStart "objective" fStart;
            print "  start blockwise recovery error:";
            dummy = _directBlockErrors(bRawStart, bRawTrue, nIntFull, nReduced, aCtl.linear);

            { thetaHat, fHat, retcode, gnorm, iterUsed, rejectedSteps } = _directLM(thetaStart, w,
                intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec,
                300, 1e-3, 1e-6, 0);
            bRawHat = _directObjUnpack(thetaHat.*scaleVec, ng, nReduced);

            print "  LM result: retcode" retcode "(0=converged,2=maxiter,";
            print "  3=infeasible start,4=stalled) gnorm" gnorm "fHat" fHat
                "fTrue" fTrue "objective gap" abs(fHat-fTrue);
            print "  LM endpoint blockwise recovery error:";
            dummy = _directBlockErrors(bRawHat, bRawTrue, nIntFull, nReduced, aCtl.linear);

            structErr = 0;
            betaRowChk = nIntFull+nReduced+1;
            if modelFlag == 0;
                uRowStartChk = betaRowChk+1;
            else;
                uRowStartChk = betaRowChk+2;
            endif;
            structErr = maxc(maxc(abs(bRawHat[1:uRowStartChk-1,1:nReduced]
                - bRawTrue[1:uRowStartChk-1,1:nReduced])));

            pass = (retcode==0) and feasStart and (abs(fHat-fTrue)<1) and (structErr<.05);
            if pass;
                print "  VERDICT: PASS (feasible+converged, objective near";
                print "  truth, structural blocks recovered; u-row reported";
                print "  separately above, not gated).";
            else;
                print "  VERDICT: FAIL -- retcode" retcode "structErr" structErr
                    "objective gap" abs(fHat-fTrue);
            endif;
            print;

            trial = trial+1;
        endo;
        perturbIdx = perturbIdx+1;
    endo;

    /* ---- Item 3 (continued): line profiling ---- */
    print "Item 3 (continued): objective line profiles";
    print "  start(5%,trial1) -> truth:";
    rndseed 500000+1000*modelFlag+100*2+1;
    bRawStart = bRawTrue[.,1:nReduced].*(1+.05*rndn(ng, nReduced));
    thetaStart = _directObjPack(bRawStart, nReduced)./scaleVec;
    dummy = _directLineProfile(thetaStart, thetaTrue, 7, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);

    print "  weakest Jacobian direction from truth (+/- 2 scaled units):";
    weakDir = rightVecsT[., cols(rightVecsT)];
    dummy = _directLineProfile(thetaTrue-2*weakDir, thetaTrue+2*weakDir, 9, w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec);
    print;

    modelFlag = modelFlag+1;
endo;

/* ---- Item 1: retrospective blockwise diagnosis of Phase 1's own saved
   QUAIDS result (the ORIGINAL _quaidsNoiselessDGP fixture + the warm
   start used in phase1_objective_validation.e), to identify exactly
   which block drove the 0.225 recErr Phase 1 reported. ---- */
print "================================================================";
print "Item 1: retrospective blockwise diagnosis of Phase 1's QUAIDS result";
print "================================================================";
{ w, intcpt, prices, totexp, instr, trueParams } =
    _quaidsNoiselessDGP(tobs, seed, 1);

struct quaidsControl aCtlOld;
aCtlOld = quaidsControlCreate;
aCtlOld.linear = 0;
aCtlOld.homogenous = 1;
aCtlOld.err = .0001;

struct quaidsControl aCtlStoneOld;
aCtlStoneOld = aCtlOld;
aCtlStoneOld.maxiter = 1;
struct quaidsOut qOutStoneOld;
qOutStoneOld = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStoneOld);

nGoods = cols(prices);
nReduced = nGoods-1;
nIntFull = qOutStoneOld.nint+1;
ng = qOutStoneOld.ng;
u = qOutStoneOld.u;
intcptFull = qOutStoneOld.intcptFull;

bRawTrue = _directTrueRaw(trueParams, nIntFull, nReduced, nGoods);
scaleVec = _directObjScale(qOutStoneOld.homogB, nReduced);

rndseed 900001;
bRawStart = bRawTrue[., 1:nReduced].*(1+.05*rndn(ng, nReduced));
thetaStart = _directObjPack(bRawStart, nReduced)./scaleVec;

{ thetaWarm, fWarm } = _directObjGradDescent(thetaStart, w, intcptFull,
    prices, totexp, u, aCtlOld, ng, nReduced, scaleVec, 400, 1e-6);
bRawWarm = _directObjUnpack(thetaWarm.*scaleVec, ng, nReduced);

print "Reproducing phase1_objective_validation.e's QUAIDS warm-start";
print "endpoint (same seeds/fixture/method) for blockwise diagnosis:";
warmRecoveryErr = maxc(maxc(abs(bRawWarm[.,1:nReduced] - bRawTrue[.,1:nReduced])));
print "  warm-start-only recErr (should match the ~0.1-0.3 range Phase 1";
print "  reported across runs):" warmRecoveryErr;
print "  blockwise breakdown:";
dummy = _directBlockErrors(bRawWarm, bRawTrue, nIntFull, nReduced, aCtlOld.linear);

print;
print "phase1b_audit.e: run complete.";
