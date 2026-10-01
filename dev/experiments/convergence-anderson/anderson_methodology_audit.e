new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Methodology audit for the handoff results.  This intentionally checks:
**   1. Anderson convergence by the true fixed-point residual, not merely
**      the accelerated step size;
**   2. both directions of homogCrit selection, because quaidsFit defines
**      homogCrit as -ln(det(Sigma)), so LARGER is the better fit;
**   3. a directly recomputed residual-covariance criterion on final bS;
**   4. basin diversity with and without Anderson under identical starts.
**
** This was smoke-tested first at nSeeds=3 before the 30-seed pilot.
*/

proc (1) = _structuralFitCrit(w, bFinal, intcptFull, prices, totexp, u, struct quaidsControl aCtl);
    local nIntFull, nGoods, nEndog, alpha, priceCoef, expendCoef,
        quadCoef, controlCoef, priceIndex, realExp, expScale, realExp2,
        fitted, resid, residCov, critDet;

    nIntFull = cols(intcptFull);
    nGoods = cols(prices);
    nEndog = 2 - aCtl.linear;

    alpha = intcptFull*bFinal[1:nIntFull, .];
    priceCoef = bFinal[nIntFull+1:nIntFull+nGoods, .];
    expendCoef = bFinal[nIntFull+nGoods+1, .];

    priceIndex = aCtl.alpha0 + sumc((prices.*alpha)')
        + .5*sumc(((prices*priceCoef).*prices)');
    realExp = totexp - priceIndex;
    fitted = alpha + prices*priceCoef + realExp*expendCoef;

    if not aCtl.linear;
        quadCoef = bFinal[nIntFull+nGoods+2, .];
        expScale = exp(prices*expendCoef');
        realExp2 = (realExp^2)./expScale;
        fitted = fitted + realExp2*quadCoef;
    endif;

    controlCoef = bFinal[nIntFull+nGoods+nEndog+1:rows(bFinal), .];
    fitted = fitted + u*controlCoef;
    resid = w - fitted;
    residCov = resid'resid/rows(w);
    critDet = det(residCov[1:nGoods-1, 1:nGoods-1]);
    if critDet <= 0;
        retp(miss(0, 0));
    endif;
    retp(-ln(critDet));
endp;

proc (1) = _rawHomogFitCrit(w, bRaw, intcptFull, prices, totexp, u, struct quaidsControl aCtl);
    local nIntFull, nGoods, nReduced, alpha, priceCoef, expendCoef,
        relativePrices, priceIndex, realExp, expScale, realExp2,
        designMat, fitted, resid, residCov, critDet;

    nIntFull = cols(intcptFull);
    nGoods = cols(prices);
    nReduced = nGoods - 1;
    relativePrices = (prices[., 1:nReduced] - prices[., nGoods])~prices[., nGoods];

    alpha = intcptFull*bRaw[1:nIntFull, 1:nReduced];
    priceCoef = bRaw[nIntFull+1:nIntFull+nReduced, 1:nReduced]
        |zeros(1, nReduced);
    priceIndex = aCtl.alpha0 + prices[., nGoods]
        + sumc((relativePrices[., 1:nReduced].*alpha)')
        + .5*sumc(((relativePrices*priceCoef)
        .*relativePrices[., 1:nReduced])');
    realExp = totexp - priceIndex;

    if aCtl.linear;
        designMat = intcptFull~relativePrices[., 1:nReduced]~realExp~u;
    else;
        expendCoef = bRaw[nIntFull+nReduced+1, 1:nReduced];
        expScale = exp(relativePrices[., 1:nReduced]*expendCoef');
        realExp2 = (realExp^2)./expScale;
        designMat = intcptFull~relativePrices[., 1:nReduced]~realExp~realExp2~u;
    endif;

    fitted = designMat*bRaw;
    resid = w - fitted;
    residCov = resid'resid/rows(w);
    critDet = det(residCov[1:nReduced, 1:nReduced]);
    if critDet <= 0;
        retp(miss(0, 0));
    endif;
    retp(-ln(critDet));
endp;

proc (2) = _appendDistinct(representatives, candidate, distinctCount, distinctTol);
    local distances;
    if distinctCount == 0;
        representatives = candidate;
        distinctCount = 1;
    else;
        distances = maxc(abs(representatives - candidate));
        if minc(distances) > distinctTol;
            representatives = representatives~candidate;
            distinctCount = distinctCount + 1;
        endif;
    endif;
    retp(representatives, distinctCount);
endp;

nSeeds = 30;
tobs = 3000;
wrongCutoff = 1.0;
anDepth = 8;
numStarts = 8;
pertScale = 0.3;
distinctTol = 1e-4;

struct quaidsControl aCtlStone;
struct quaidsControl aCtlStart;

modelFlag = 0;
do while modelFlag <= 1;
    if modelFlag == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    maxHomogCorrect = 0;
    minHomogCorrect = 0;
    maxStructCorrect = 0;
    oracleCorrect = 0;
    wrongBeatsTruthSeeds = 0;
    wrongRawBeatsTruthSeeds = 0;
    legacyOnlyConverged = 0;
    residualOnlyConverged = 0;
    sameConvergenceFlag = 0;
    sumDistinctAnderson = 0;
    sumDistinctPlain = 0;
    sumConvergedAnderson = 0;
    sumConvergedPlain = 0;
    maxCritReplicaDiff = 0;
    maxCritDiffCorrect = 0;
    maxCritDiffWrong = 0;
    maxRawCritReplicaDiff = 0;

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1 - modelFlag;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } =
            _quaidsSyntheticDGP(tobs, seed, modelFlag, 1);

        aCtlStone = aCtl;
        aCtlStone.maxiter = 1;
        qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);
        bStone = qOutStone.homogB;
        trueCrit = _structuralFitCrit(w, trueParams, qOutStone.intcptFull,
            prices, totexp, qOutStone.u, aCtl);
        nIntFull = cols(qOutStone.intcptFull);
        nGoods = cols(prices);
        trueRaw = trueParams[1:nIntFull, .]
            |trueParams[nIntFull+1:nIntFull+nGoods-1, .]
            |trueParams[nIntFull+nGoods+1:rows(trueParams), .];
        trueRawCrit = _rawHomogFitCrit(w, trueRaw, qOutStone.intcptFull,
            prices, totexp, qOutStone.u, aCtl);
        zeroFrac = meanc(vec(abs(bStone) .< 1e-8));

        rndseed 1000*modelFlag + seed;
        foundAnderson = 0;
        foundPlain = 0;
        distinctAnderson = 0;
        distinctPlain = 0;
        repsAnderson = 0;
        repsPlain = 0;
        foundWrongBetterTruth = 0;
        foundWrongRawBetterTruth = 0;

        startIdx = 1;
        do while startIdx <= numStarts;
            if startIdx == 1;
                bStart = bStone;
            else;
                bStart = bStone + pertScale*abs(bStone)
                    .*rndn(rows(bStone), cols(bStone));
            endif;

            aCtlStart = aCtl;
            aCtlStart.b0 = bStart;

            qOutAnderson = _quaidsFitAnderson(w, intcpt, prices, totexp,
                instr, aCtlStart, 0, anDepth, 1);
            qOutPlain = _quaidsFitAnderson(w, intcpt, prices, totexp,
                instr, aCtlStart, 0, 0, 0);

            if startIdx == 1;
                qOutLegacy = _quaidsFitAnderson(w, intcpt, prices, totexp,
                    instr, aCtlStart, 0, anDepth, 0);
                if qOutLegacy.converged and not qOutAnderson.converged;
                    legacyOnlyConverged = legacyOnlyConverged + 1;
                elseif not qOutLegacy.converged and qOutAnderson.converged;
                    residualOnlyConverged = residualOnlyConverged + 1;
                else;
                    sameConvergenceFlag = sameConvergenceFlag + 1;
                endif;
            endif;

            if qOutAnderson.converged;
                sumConvergedAnderson = sumConvergedAnderson + 1;
                recErrAnderson = maxc(maxc(abs(qOutAnderson.bS - trueParams)));
                structCritAnderson = _structuralFitCrit(w, qOutAnderson.bS,
                    qOutAnderson.intcptFull, prices, totexp, qOutAnderson.u, aCtl);
                rawCritAnderson = _rawHomogFitCrit(w, qOutAnderson.homogB,
                    qOutAnderson.intcptFull, prices, totexp, qOutAnderson.u, aCtl);
                maxRawCritReplicaDiff = maxc(maxRawCritReplicaDiff
                    |abs(rawCritAnderson - qOutAnderson.homogCrit));
                maxCritReplicaDiff = maxc(maxCritReplicaDiff
                    |abs(structCritAnderson - qOutAnderson.symcCrit));
                if recErrAnderson <= wrongCutoff;
                    maxCritDiffCorrect = maxc(maxCritDiffCorrect
                        |abs(structCritAnderson - qOutAnderson.symcCrit));
                else;
                    maxCritDiffWrong = maxc(maxCritDiffWrong
                        |abs(structCritAnderson - qOutAnderson.symcCrit));
                endif;
                { repsAnderson, distinctAnderson } = _appendDistinct(
                    repsAnderson, vec(qOutAnderson.bS), distinctAnderson, distinctTol);

                if not foundAnderson;
                    maxHomog = qOutAnderson.homogCrit;
                    minHomog = qOutAnderson.homogCrit;
                    maxStruct = structCritAnderson;
                    maxHomogRecErr = recErrAnderson;
                    minHomogRecErr = recErrAnderson;
                    maxStructRecErr = recErrAnderson;
                    oracleRecErr = recErrAnderson;
                    foundAnderson = 1;
                else;
                    if qOutAnderson.homogCrit > maxHomog;
                        maxHomog = qOutAnderson.homogCrit;
                        maxHomogRecErr = recErrAnderson;
                    endif;
                    if qOutAnderson.homogCrit < minHomog;
                        minHomog = qOutAnderson.homogCrit;
                        minHomogRecErr = recErrAnderson;
                    endif;
                    if structCritAnderson > maxStruct;
                        maxStruct = structCritAnderson;
                        maxStructRecErr = recErrAnderson;
                    endif;
                    if recErrAnderson < oracleRecErr;
                        oracleRecErr = recErrAnderson;
                    endif;
                endif;

                if recErrAnderson > wrongCutoff and structCritAnderson > trueCrit;
                    foundWrongBetterTruth = 1;
                endif;
                if recErrAnderson > wrongCutoff and rawCritAnderson > trueRawCrit;
                    foundWrongRawBetterTruth = 1;
                endif;
            endif;

            if qOutPlain.converged;
                sumConvergedPlain = sumConvergedPlain + 1;
                { repsPlain, distinctPlain } = _appendDistinct(
                    repsPlain, vec(qOutPlain.bS), distinctPlain, distinctTol);
                foundPlain = 1;
            endif;

            startIdx = startIdx + 1;
        endo;

        sumDistinctAnderson = sumDistinctAnderson + distinctAnderson;
        sumDistinctPlain = sumDistinctPlain + distinctPlain;
        wrongBeatsTruthSeeds = wrongBeatsTruthSeeds + foundWrongBetterTruth;
        wrongRawBeatsTruthSeeds = wrongRawBeatsTruthSeeds + foundWrongRawBetterTruth;

        if foundAnderson;
            maxHomogCorrect = maxHomogCorrect + (maxHomogRecErr <= wrongCutoff);
            minHomogCorrect = minHomogCorrect + (minHomogRecErr <= wrongCutoff);
            maxStructCorrect = maxStructCorrect + (maxStructRecErr <= wrongCutoff);
            oracleCorrect = oracleCorrect + (oracleRecErr <= wrongCutoff);
            print "seed" seed "model" modelName "zeroFrac" zeroFrac
                "convA" foundAnderson "distinctA" distinctAnderson
                "convP" foundPlain "distinctP" distinctPlain
                "recErr maxHomog/minHomog/maxStruct/oracle"
                maxHomogRecErr minHomogRecErr maxStructRecErr oracleRecErr
                "trueCrit" trueCrit "bestStructCrit" maxStruct;
        else;
            print "seed" seed "model" modelName
                "no residual-validated Anderson start converged"
                "convP" foundPlain "distinctP" distinctPlain;
        endif;

        seed = seed + 1;
    endo;

    print;
    print "==============================================================";
    print "METHODOLOGY AUDIT:" modelName "seeds" nSeeds "starts" numStarts;
    print "==============================================================";
    print "correct selections -- MAX homogCrit:" maxHomogCorrect
        "MIN homogCrit:" minHomogCorrect "MAX structural crit:" maxStructCorrect
        "oracle:" oracleCorrect;
    print "seeds where a wrong solution's structural criterion beats truth:"
        wrongBeatsTruthSeeds;
    print "seeds where a wrong solution's raw homogeneity criterion beats truth:"
        wrongRawBeatsTruthSeeds;
    print "unperturbed legacy-step-only converged:" legacyOnlyConverged
        "residual-only converged:" residualOnlyConverged
        "same convergence flag:" sameConvergenceFlag;
    print "total converged starts -- Anderson:" sumConvergedAnderson
        "plain:" sumConvergedPlain;
    print "sum of within-seed distinct converged endpoints -- Anderson:"
        sumDistinctAnderson "plain:" sumDistinctPlain;
    print "max abs diff: recomputed structural criterion vs qOut.symcCrit:"
        maxCritReplicaDiff;
    print "  max diff among correct endpoints:" maxCritDiffCorrect
        "among wrong endpoints:" maxCritDiffWrong;
    print "max abs diff: recomputed raw criterion vs qOut.homogCrit:"
        maxRawCritReplicaDiff;
    print;

    modelFlag = modelFlag + 1;
endo;

print "anderson_methodology_audit.e: run complete.";
