new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;
#include direct_objective.src;

/*
** phase3_smoke_gate.e
**
** Phase 3 (continuation/homotopy) smoke gate. Per CODEX_HANDOFF.md's own
** stage-gate discipline and the user's explicit Phase 3 instructions:
** a TIGHTLY BOUNDED test on one known-difficult seed per model, BEFORE
** any 30-seed pilot. Uses the REAL noisy _quaidsSyntheticDGP() fixture
** (tests/quaidsfixtures.src), NOT Phase 1b's near-stationary fixtures.
**
** Known-difficult seeds: chosen from
** anderson_full_multistart_benchmark_output.txt's own "no safeguarded
** start converged" lines -- seeds where ALL 16 structured starts x
** safeguarded Anderson(8) failed to converge in the completed 200-seed
** benchmark. AIDS seed=6, QUAIDS seed=9, tobs=3000.
**
** Compares, on identical data and starts: direct LM without continuation
** (_directLM, homotopy=1 throughout), continuation+LM (this script's own
** driver around _directLMHomotopyDiag), and safe Anderson
** (_quaidsFitAnderson, safeMode=1, depth 8).
**
** Phase 3 hardening (2026-10-05):
**   - the driver now SOLVES the homotopy=0 problem first (stage zero)
**     and only advances after a feasible, genuinely converged stage-zero
**     endpoint exists; the pre-hardening schedule (first stage at
**     homotopy=step0) is still run as "legacy" so the two endpoints can
**     be compared directly;
**   - per-failure-mode bookkeeping (anyStageInfeasible renamed
**     anyStageFailedOnce, plus separate counters);
**   - structural error now uses the CORRECT row range (all rows except
**     the single u row). The pre-hardening slice
**     1:nIntFull+nReduced+2-modelFlag INCLUDED the u row for AIDS and
**     EXCLUDED the lambda row for QUAIDS; it is still printed as
**     "pre-fix slice" for traceability against the earlier output;
**   - explicit endpoint-vector differences across starts;
**   - absolute correctness threshold structErr <= 1 in the verdict, with
**     threshold sensitivity at .5/1/2/5;
**   - every computed gate predicate enters the verdict.
*/

/* ---- Step 0: faithfulness check -- homotopy=1 IS the original
   objective (by construction, and verified numerically here), and
   homotopy=0 really is a different objective (so the stage-zero solve
   is not a no-op). ---- */
print "================================================================";
print "Step 0: continuation faithfulness check (homotopy=1 == original)";
print "================================================================";
{ wChk, intcptChk, pricesChk, totexpChk, instrChk, trueChk } =
    _quaidsSyntheticDGP(200, 1, 1, 1);
struct quaidsControl aCtlChk;
aCtlChk = quaidsControlCreate;
aCtlChk.linear = 0;
aCtlChk.homogenous = 1;
nGoodsChk = cols(pricesChk);
nRedChk = nGoodsChk-1;
intcptFullChk = ones(200,1)~intcptChk;
ngChk = cols(intcptFullChk)+nRedChk+2+1;
rndseed 123456;
thetaChk = rndn(ngChk*nRedChk, 1);
scaleChk = ones(ngChk*nRedChk, 1);
uChk = .05*rndn(200, 1);
{ fA, feasA, residA } = _directObjCore(thetaChk, wChk, intcptFullChk, pricesChk, totexpChk, uChk, aCtlChk, ngChk, nRedChk, scaleChk);
{ fB, feasB, residB } = _directObjCoreHomotopy(thetaChk, wChk, intcptFullChk, pricesChk, totexpChk, uChk, aCtlChk, ngChk, nRedChk, scaleChk, 1);
{ fZ, feasZ, residZ } = _directObjCoreHomotopy(thetaChk, wChk, intcptFullChk, pricesChk, totexpChk, uChk, aCtlChk, ngChk, nRedChk, scaleChk, 0);
print "  max abs diff, objective:" abs(fA-fB) "  feasible match:" (feasA==feasB);
print "  max abs diff, residual:" maxc(maxc(abs(residA-residB)));
print "  homotopy=0 vs 1 residual max abs diff (must be > 0):" maxc(maxc(abs(residA-residZ)));
if abs(fA-fB)==0 and maxc(maxc(abs(residA-residB)))==0 and maxc(maxc(abs(residA-residZ))) > 0;
    print "  PASS: homotopy=1 is bit-identical to the original objective;";
    print "        homotopy=0 is a genuinely different objective.";
else;
    print "  FAIL: homotopy endpoint check failed -- STOP.";
endif;
print;

/* ---- Error decomposition with the CORRECT row layout:
   rows 1:nIntFull alpha | gamma n1 rows | beta | [lambda] | u (1 row).
   Returns structural (all rows except u), lambda (-1 when linear), u,
   overall, and the pre-fix slice for traceability. ---- */
proc (5) = _phase3Errs(bHat, bTrue, nIntFull, n1, linear);
    local nStruct, bH, bT, structErr, lamErr, uErr, overallErr, preFixErr,
        preFixRows;
    nStruct = nIntFull+n1+1+(1-linear);
    bH = bHat[., 1:n1];
    bT = bTrue[., 1:n1];
    if rows(bH) /= nStruct+1;
        print "  LAYOUT ERROR: rows(bHat)=" rows(bH) " expected" nStruct+1;
        end;
    endif;
    structErr = maxc(maxc(abs(bH[1:nStruct, .] - bT[1:nStruct, .])));
    if linear;
        lamErr = -1;
    else;
        lamErr = maxc(maxc(abs(bH[nStruct, .] - bT[nStruct, .])));
    endif;
    uErr = maxc(maxc(abs(bH[nStruct+1, .] - bT[nStruct+1, .])));
    overallErr = maxc(maxc(abs(bH - bT)));
    preFixRows = nIntFull+n1+2-(1-linear);
    preFixErr = maxc(maxc(abs(bH[1:preFixRows, .] - bT[1:preFixRows, .])));
    retp(structErr, lamErr, uErr, overallErr, preFixErr);
endp;

/* ---- Continuation driver.
   doStageZero=1 (hardened): solve homotopy=0 from theta0 first; stop
   (no advancing) if that stage is not genuinely converged. Then advance
   by `step`, warm-starting every stage from the last ACCEPTED endpoint;
   on any nonzero stage retcode, do not advance, halve the step and retry
   from the last accepted endpoint; stop honestly (minStepExhausted=1) if
   the step falls below minStep. The last stage is snapped to exactly
   homotopy=1.0, i.e. exactly the original objective.
   doStageZero=0 (legacy, pre-hardening schedule): first stage at
   homotopy=step0, otherwise identical.

   A stage is ACCEPTED on retcode==0; its feasibility is re-evaluated
   independently and recorded (allAcceptedFeasible), so that predicate is
   a real check, not a tautology.

   info (20x1):
     1 reachedOne           2 stage0Ok (1/0; -1 = not attempted, legacy)
     3 stagesAttempted      4 stagesAccepted
     5 totalIter            6 totalRejectedSteps
     7 anyStageFailedOnce   8 stages retcode 3 (infeasible start)
     9 stages retcode 5 (covariance-inversion failure)
    10 stages retcode 2 (maxIter)
    11 stages retcode 4 (stalled: damping exhausted)
    12 adaptive step reductions
    13 minStepExhausted     14 allAcceptedConverged (retcode 0)
    15 allAcceptedFeasible (independent re-evaluation)
    16 sum covariance-inversion failures (inner events)
    17 sum damped-step inversion failures (inner events)
    18 sum infeasible trial points
    19 sum rho-rejected feasible trials
    20 final step size ---- */
proc (5) = _phase3Continuation(theta0, w, intcptFull, prices, totexp, u, struct quaidsControl aCtl, ng, n1, scaleVec, maxIterPerStage, gradTol, h, step0, minStep, doStageZero);
    local homotopy, step, theta, info, fFinal, tStart, wallTime, stopNow,
        firstStage, homotopyTry, thetaTry, fTry, retcodeTry, gnormTry,
        iterTry, rejTry, diagTry, fChk, feasChk, residChk, feasOk;

    tStart = hsec;
    homotopy = 0;
    step = step0;
    theta = theta0;
    fFinal = miss(0, 0);
    info = zeros(20, 1);
    info[2] = -1;
    info[14] = 1;
    info[15] = 1;
    firstStage = doStageZero;
    stopNow = 0;

    do while not stopNow;
        if firstStage;
            homotopyTry = 0;
        else;
            if homotopy >= 1;
                break;
            endif;
            if step < minStep;
                info[13] = 1;
                print "    continuation STOPPED: step" step "< minStep" minStep
                    "at homotopy" homotopy "-- reported as NOT reaching 1.";
                break;
            endif;
            homotopyTry = homotopy+step;
            if homotopyTry >= 1-1e-12;
                homotopyTry = 1;
            endif;
        endif;

        { thetaTry, fTry, retcodeTry, gnormTry, iterTry, rejTry, diagTry } =
            _directLMHomotopyDiag(theta, w, intcptFull, prices, totexp, u,
            aCtl, ng, n1, scaleVec, maxIterPerStage, gradTol, h, 0,
            homotopyTry, 0, 0, 0);
        { fChk, feasChk, residChk } = _directObjCoreHomotopy(thetaTry, w,
            intcptFull, prices, totexp, u, aCtl, ng, n1, scaleVec,
            homotopyTry);
        feasOk = feasChk and not scalmiss(fChk);

        info[3] = info[3]+1;
        info[5] = info[5]+iterTry;
        info[6] = info[6]+rejTry;
        info[16:19] = info[16:19]+diagTry[1:4];

        print "    stage" info[3] "homotopy->" homotopyTry "retcode" retcodeTry
            "feasible" feasOk "f" fTry "gnorm" gnormTry "iter" iterTry
            "rejected" rejTry;
        print "      [covInvFail stepInvFail infeasTrial rhoReject]"
            diagTry[1:4]';

        if retcodeTry == 0;
            theta = thetaTry;
            fFinal = fTry;
            homotopy = homotopyTry;
            info[4] = info[4]+1;
            info[15] = info[15] and feasOk;
            if firstStage;
                info[2] = 1;
            endif;
            if homotopy >= 1;
                info[1] = 1;
            endif;
        else;
            info[7] = 1;
            if retcodeTry == 3;
                info[8] = info[8]+1;
            elseif retcodeTry == 5;
                info[9] = info[9]+1;
            elseif retcodeTry == 2;
                info[10] = info[10]+1;
            elseif retcodeTry == 4;
                info[11] = info[11]+1;
            endif;
            if firstStage;
                info[2] = 0;
                print "    stage ZERO failed -- not advancing; continuation stops.";
                stopNow = 1;
            else;
                print "    stage REJECTED at homotopy" homotopyTry
                    "-- halving step from" step "to" step/2
                    "and retrying from the last ACCEPTED endpoint.";
                step = step/2;
                info[12] = info[12]+1;
            endif;
        endif;
        firstStage = 0;
    endo;

    if info[4] == 0;
        info[14] = 0;
        info[15] = 0;
    endif;
    info[20] = step;
    wallTime = (hsec-tStart)/100;
    retp(theta, fFinal, homotopy, wallTime, info);
endp;

proc (1) = _phase3PrintInfo(info, homotopyFinal, wallTime);
    print "      reachedHomotopy=1" info[1] "finalHomotopy" homotopyFinal
        "stage0Ok" info[2] "stagesAttempted" info[3] "stagesAccepted"
        info[4] "wallTime(s)" wallTime;
    print "      totalIter" info[5] "totalRejectedSteps" info[6]
        "anyStageFailedOnce" info[7] "stepReductions" info[12]
        "minStepExhausted" info[13];
    print "      failed stages by retcode [3 infeasStart, 5 covInv, 2 maxIter, 4 stalled]:"
        info[8:11]';
    print "      inner events [covInvFail stepInvFail infeasTrial rhoReject]:"
        info[16:19]';
    print "      allAcceptedConverged" info[14] "allAcceptedFeasible" info[15];
    retp(1);
endp;

/* ---- Structured ("far") perturbed start (unchanged from the
   pre-hardening gate, so the same start is reproduced). ---- */
proc (1) = _phase3FarStart(bStone, nIntFull, nReduced);
    local bStart, priceRows;
    priceRows = seqa(nIntFull+1, 1, nReduced);
    bStart = bStone + .8*abs(bStone).*rndn(rows(bStone), cols(bStone));
    bStart[priceRows, .] = -bStone[priceRows, .];
    retp(bStart);
endp;

tobs = 3000;
anDepth = 8;
rankTol = 1e-10;
residualGrowth = 10;
stepFactor = 10;
primaryThr = 1;
thrGrid = .5|1|2|5;
regTol = .05;
eqTol = 1e-6;
nGatePass = 0;
nGateCases = 0;

modelFlag = 0;
do while modelFlag <= 1;
    if modelFlag == 0;
        modelName = "Iterated AIDS (linear)";
        hardSeed = 6;
    else;
        modelName = "QUAIDS";
        hardSeed = 9;
    endif;

    print "================================================================";
    print "" $+ modelName $+ ": Phase 3 smoke gate, known-difficult seed"
        $+ ftocv(hardSeed, 1, 0);
    print "================================================================";

    { w, intcpt, prices, totexp, instr, trueParams } =
        _quaidsSyntheticDGP(tobs, hardSeed, modelFlag, 1);

    struct quaidsControl aCtl;
    aCtl = quaidsControlCreate;
    aCtl.linear = 1-modelFlag;
    aCtl.homogenous = 1;
    aCtl.maxiter = 100;
    aCtl.err = .0001;

    struct quaidsControl aCtlStone;
    aCtlStone = aCtl;
    aCtlStone.maxiter = 1;
    struct quaidsOut qOutStone;
    qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);

    nGoods = cols(prices);
    nReduced = nGoods-1;
    nIntFull = qOutStone.nint+1;
    ng = qOutStone.ng;
    u = qOutStone.u;
    intcptFull = qOutStone.intcptFull;
    bStone = qOutStone.homogB;
    linear = aCtl.linear;
    nStruct = nIntFull+nReduced+1+(1-linear);

    print "  layout: ng" ng "nIntFull" nIntFull "n1" nReduced
        "structural rows (excl. u)" nStruct;
    if ng /= nStruct+1;
        print "  LAYOUT MISMATCH -- STOP.";
        end;
    endif;

    trueRaw = _directTrueRaw(trueParams, nIntFull, nReduced, nGoods);
    scaleVec = _directObjScale(bStone, nReduced);

    rndseed 300000+modelFlag;
    bFar = _phase3FarStart(bStone, nIntFull, nReduced);

    startIdx = 1;
    do while startIdx <= 2;
        if startIdx == 1;
            startName = "Stone start";
            bStart = bStone;
        else;
            startName = "far/structured start (0.8x noise + gamma sign flip)";
            bStart = bFar;
        endif;

        print "----------------------------------------------------------";
        print "" $+ modelName $+ "/" $+ startName;
        print "----------------------------------------------------------";

        thetaStart = _directObjPack(bStart, nReduced)./scaleVec;
        startRecErr = maxc(maxc(abs(bStart[.,1:nReduced] - trueRaw[.,1:nReduced])));
        print "  start overall recErr vs true DGP params:" startRecErr;

        /* ---- (a) direct LM, no continuation ---- */
        tA = hsec;
        { thetaA, fA2, retcodeA, gnormA, iterA, rejA } = _directLM(thetaStart,
            w, intcptFull, prices, totexp, u, aCtl, ng, nReduced, scaleVec,
            500, 1e-3, 1e-6, 0);
        wallA = (hsec-tA)/100;
        bRawA = _directObjUnpack(thetaA.*scaleVec, ng, nReduced);
        { structErrA, lamErrA, uErrA, overallErrA, preFixA } =
            _phase3Errs(bRawA, trueRaw, nIntFull, nReduced, linear);
        homogCritA = -fA2;
        print "  (a) direct LM (no continuation): retcode" retcodeA
            "gnorm" gnormA "iter" iterA "rejected" rejA "wallTime(s)" wallA;
        print "      homogCrit" homogCritA "structErr(excl.u)" structErrA
            "lambdaErr" lamErrA "uErr" uErrA "overallErr" overallErrA
            "pre-fix slice" preFixA;

        /* ---- (b) continuation + LM, hardened (stage zero) ---- */
        print "  (b) continuation+LM (hardened: stage zero at homotopy=0):";
        { thetaB, fB2, homotopyFinalB, wallB, infoB } =
            _phase3Continuation(thetaStart, w, intcptFull, prices, totexp,
            u, aCtl, ng, nReduced, scaleVec, 500, 1e-3, 1e-6, .25, .01, 1);
        call _phase3PrintInfo(infoB, homotopyFinalB, wallB);
        bRawB = _directObjUnpack(thetaB.*scaleVec, ng, nReduced);
        { structErrB, lamErrB, uErrB, overallErrB, preFixB } =
            _phase3Errs(bRawB, trueRaw, nIntFull, nReduced, linear);
        homogCritB = -fB2;
        print "      homogCrit" homogCritB "structErr(excl.u)" structErrB
            "lambdaErr" lamErrB "uErr" uErrB "overallErr" overallErrB
            "pre-fix slice" preFixB;
        print "      approx total evals (iter*(2k+1)+rejected, k="
            $+ ftocv(ng*nReduced,1,0) $+"):"
            infoB[5]*(2*ng*nReduced+1)+infoB[6];

        /* ---- (b') continuation+LM, legacy schedule (no stage zero) ---- */
        print "  (b') continuation+LM (legacy: first stage at homotopy=.25):";
        { thetaL, fL2, homotopyFinalL, wallL, infoL } =
            _phase3Continuation(thetaStart, w, intcptFull, prices, totexp,
            u, aCtl, ng, nReduced, scaleVec, 500, 1e-3, 1e-6, .25, .01, 0);
        call _phase3PrintInfo(infoL, homotopyFinalL, wallL);
        bRawL = _directObjUnpack(thetaL.*scaleVec, ng, nReduced);
        { structErrL, lamErrL, uErrL, overallErrL, preFixL } =
            _phase3Errs(bRawL, trueRaw, nIntFull, nReduced, linear);
        homogCritL = -fL2;
        print "      homogCrit" homogCritL "structErr(excl.u)" structErrL
            "lambdaErr" lamErrL "uErr" uErrL "overallErr" overallErrL
            "pre-fix slice" preFixL;
        dThetaBL = maxc(abs(thetaB-thetaL));
        dStructBL = maxc(maxc(abs(bRawB[1:nStruct,1:nReduced]-bRawL[1:nStruct,1:nReduced])));
        print "  stage-zero vs legacy endpoint: max|dTheta(scaled)|" dThetaBL
            "max|dRaw structural|" dStructBL "|dObjective|" abs(fB2-fL2);

        /* ---- (c) safe Anderson ---- */
        struct quaidsControl aCtlAnd;
        aCtlAnd = aCtl;
        /* bStart already has nGoods (n) columns (homogB is ng x n). */
        aCtlAnd.b0 = bStart;
        tC = hsec;
        struct quaidsOut qOutAnd;
        qOutAnd = _quaidsFitAnderson(w, intcpt, prices, totexp, instr,
            aCtlAnd, 0, anDepth, 1, 1, rankTol, residualGrowth, stepFactor);
        wallC = (hsec-tC)/100;
        { structErrC, lamErrC, uErrC, overallErrC, preFixC } =
            _phase3Errs(qOutAnd.homogB, trueRaw, nIntFull, nReduced, linear);
        print "  (c) safe Anderson: converged" qOutAnd.converged
            "iterations" qOutAnd.iterations "finalErr" qOutAnd.finalErr
            "wallTime(s)" wallC;
        print "      homogCrit" qOutAnd.homogCrit "structErr(excl.u)"
            structErrC "lambdaErr" lamErrC "uErr" uErrC "overallErr"
            overallErrC "pre-fix slice" preFixC;

        /* ---- threshold sensitivity (structErr excl. u) ---- */
        print "  structErr <= threshold [.5 1 2 5]:";
        print "    cold LM            " (structErrA .<= thrGrid)';
        print "    continuation (b)   " (structErrB .<= thrGrid)';
        print "    continuation (b')  " (structErrL .<= thrGrid)';
        print "    safe Anderson      " (structErrC .<= thrGrid)';

        /* ---- strengthened verdict: every predicate below is used ---- */
        pReached = (infoB[1] == 1) and (homotopyFinalB == 1);
        pStageZero = (infoB[2] == 1);
        pAcceptedConverged = (infoB[14] == 1);
        pAcceptedFeasible = (infoB[15] == 1);
        pCorrect = structErrB <= primaryThr;
        pNoRegCold = structErrB <= structErrA+regTol;
        pNoRegAnd = structErrB <= structErrC+regTol;
        pClean = pReached and pStageZero and pAcceptedConverged
            and pAcceptedFeasible and pCorrect and pNoRegCold and pNoRegAnd;
        print "  VERDICT predicates:";
        print "    reached homotopy=1 exactly           " pReached;
        print "    stage zero genuinely converged       " pStageZero;
        print "    all accepted stages retcode 0        " pAcceptedConverged;
        print "    all accepted stages feasible         " pAcceptedFeasible;
        print "    correct: structErr <= 1              " pCorrect;
        print "    no regression vs cold LM (+.05)      " pNoRegCold;
        print "    no regression vs safe Anderson (+.05)" pNoRegAnd;
        nGateCases = nGateCases+1;
        if pClean;
            nGatePass = nGatePass+1;
            print "  VERDICT: CLEAN PASS";
        else;
            print "  VERDICT: NOT A CLEAN PASS";
        endif;
        print;

        if startIdx == 1;
            thetaBS = thetaB; fBS = fB2; bRawBS = bRawB;
            thetaLS = thetaL; fLS = fL2; bRawLS = bRawL;
            thetaAS = thetaA; fAS = fA2; bRawAS = bRawA;
        else;
            thetaBF = thetaB; fBF = fB2; bRawBF = bRawB;
            thetaLF = thetaL; fLF = fL2; bRawLF = bRawL;
            thetaAF = thetaA; fAF = fA2; bRawAF = bRawA;
        endif;

        startIdx = startIdx+1;
    endo;

    /* ---- endpoint equality across the two starts ---- */
    print "----------------------------------------------------------";
    print "" $+ modelName $+ ": endpoint differences, Stone vs far start";
    print "----------------------------------------------------------";
    dT = maxc(abs(thetaBS-thetaBF));
    dR = maxc(maxc(abs(bRawBS[1:nStruct,1:nReduced]-bRawBF[1:nStruct,1:nReduced])));
    dF = abs(fBS-fBF);
    print "  continuation (b, stage zero): max|dTheta(scaled)|" dT
        "max|dRaw structural|" dR "|dObjective|" dF;
    if dT == 0;
        print "    -> bit-identical endpoints";
    elseif dT <= eqTol and dR <= eqTol;
        print "    -> numerically equivalent (both diffs <= 1e-6), NOT identical";
    else;
        print "    -> DIFFERENT endpoints at tolerance 1e-6";
    endif;
    dT = maxc(abs(thetaLS-thetaLF));
    dR = maxc(maxc(abs(bRawLS[1:nStruct,1:nReduced]-bRawLF[1:nStruct,1:nReduced])));
    dF = abs(fLS-fLF);
    print "  continuation (b', legacy):    max|dTheta(scaled)|" dT
        "max|dRaw structural|" dR "|dObjective|" dF;
    if dT == 0;
        print "    -> bit-identical endpoints";
    elseif dT <= eqTol and dR <= eqTol;
        print "    -> numerically equivalent (both diffs <= 1e-6), NOT identical";
    else;
        print "    -> DIFFERENT endpoints at tolerance 1e-6";
    endif;
    dT = maxc(abs(thetaAS-thetaAF));
    dR = maxc(maxc(abs(bRawAS[1:nStruct,1:nReduced]-bRawAF[1:nStruct,1:nReduced])));
    dF = abs(fAS-fAF);
    print "  cold LM (context):            max|dTheta(scaled)|" dT
        "max|dRaw structural|" dR "|dObjective|" dF;
    print;

    modelFlag = modelFlag+1;
endo;

print "================================================================";
print "Smoke gate summary: clean passes" nGatePass "of" nGateCases;
if nGatePass == nGateCases;
    print "SMOKE GATE: PASS (all scenarios clean)";
else;
    print "SMOKE GATE: NOT PASSED";
endif;
print "phase3_smoke_gate.e: run complete.";
