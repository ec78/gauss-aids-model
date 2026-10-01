new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Full recommended benchmark:
**   - residual-validated, rank-truncated SVD Anderson(8);
**   - 16 deterministic start recipes spanning relative, additive,
**     block-targeted, sign-flip, and joint perturbations;
**   - selection by MAXIMUM homogCrit (-ln(det(Sigma)));
**   - oracle and known-truth-criterion checks for synthetic diagnosis.
**
** This file is experiment-only.  It does not change the shipped package.
** It was first smoke-tested with nSeeds=3 before the full 200-seed run.
*/

proc (1) = _rawBenchmarkCrit(w, bRaw, intcptFull, prices, totexp, u, struct quaidsControl aCtl);
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

proc (1) = _addScaledNoise(block, noiseScale);
    local blockRms, scaleFloor;
    blockRms = sqrt(meanc(vec(block^2)));
    scaleFloor = maxc(blockRms|0.05);
    retp(block + noiseScale*(abs(block)+scaleFloor)
        .*rndn(rows(block), cols(block)));
endp;

proc (1) = _benchmarkStart(bStone, startIdx, nIntFull, nReduced, quadratic);
    local bStart, alphaRows, priceRows, betaRow, quadRow, coreLast;

    bStart = bStone;
    alphaRows = seqa(1, 1, nIntFull);
    priceRows = seqa(nIntFull+1, 1, nReduced);
    betaRow = nIntFull+nReduced+1;
    quadRow = betaRow+1;
    coreLast = betaRow+quadratic;

    if startIdx == 1;
        bStart = bStone;
    elseif startIdx == 2;
        bStart = bStone + .3*abs(bStone).*rndn(rows(bStone), cols(bStone));
    elseif startIdx == 3;
        bStart = bStone + .8*abs(bStone).*rndn(rows(bStone), cols(bStone));
    elseif startIdx == 4;
        bStart = _addScaledNoise(bStone, .1);
    elseif startIdx == 5;
        bStart = _addScaledNoise(bStone, .3);
    elseif startIdx == 6;
        bStart = _addScaledNoise(bStone, .6);
    elseif startIdx == 7;
        bStart[alphaRows, .] = _addScaledNoise(bStone[alphaRows, .], .3);
    elseif startIdx == 8;
        bStart[alphaRows, .] = _addScaledNoise(bStone[alphaRows, .], .8);
    elseif startIdx == 9;
        bStart[priceRows, .] = _addScaledNoise(bStone[priceRows, .], .3);
    elseif startIdx == 10;
        bStart[priceRows, .] = _addScaledNoise(bStone[priceRows, .], .8);
    elseif startIdx == 11;
        bStart[priceRows, .] = -bStone[priceRows, .];
    elseif startIdx == 12;
        bStart[betaRow, .] = _addScaledNoise(bStone[betaRow, .], .5);
    elseif startIdx == 13;
        bStart[betaRow, .] = -bStone[betaRow, .];
    elseif startIdx == 14;
        bStart[1:coreLast, .] = _addScaledNoise(bStone[1:coreLast, .], .3);
    elseif startIdx == 15;
        bStart[1:coreLast, .] = _addScaledNoise(bStone[1:coreLast, .], .8);
    else;
        if quadratic;
            bStart[quadRow, .] = _addScaledNoise(bStone[quadRow, .], .8);
        else;
            bStart[priceRows, .] = -bStone[priceRows, .];
            bStart[betaRow, .] = -bStone[betaRow, .];
        endif;
    endif;

    retp(bStart);
endp;

proc (2) = _appendBenchmarkDistinct(representatives, candidate, distinctCount, distinctTol);
    local distances;
    if distinctCount == 0;
        representatives = candidate;
        distinctCount = 1;
    else;
        distances = maxc(abs(representatives-candidate));
        if minc(distances) > distinctTol;
            representatives = representatives~candidate;
            distinctCount = distinctCount+1;
        endif;
    endif;
    retp(representatives, distinctCount);
endp;

nSeeds = 200;
tobs = 3000;
wrongCutoff = 1.0;
anDepth = 8;
numStarts = 16;
distinctTol = 1e-4;
rankTol = 1e-10;
residualGrowth = 10;
stepFactor = 10;

struct quaidsControl aCtlStone;
struct quaidsControl aCtlStart;
struct quaidsOut qOutStart;

print "START RECIPES:";
print "1 Stone; 2-3 relative all (.3,.8); 4-6 standardized all (.1,.3,.6);";
print "7-8 alpha (.3,.8); 9-10 gamma (.3,.8); 11 gamma sign flip;";
print "12 beta (.5); 13 beta sign flip; 14-15 core joint (.3,.8);";
print "16 lambda (.8) for QUAIDS, gamma+beta sign flip for AIDS.";
print;

modelFlag = 0;
do while modelFlag <= 1;
    if modelFlag == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    singleNever = 0;
    singleWrong = 0;
    singleCorrect = 0;
    multiNever = 0;
    multiWrong = 0;
    multiCorrect = 0;
    oracleWrong = 0;
    oracleCorrect = 0;
    selectorMiss = 0;
    wrongBeatsTruthSeeds = 0;
    multiEndpointSeeds = 0;
    totalConvergedStarts = 0;
    maxCritReplicaDiff = 0;
    selectedWithin05 = 0;
    selectedWithin1 = 0;
    selectedWithin2 = 0;
    selectedWithin5 = 0;
    oracleWithin05 = 0;
    oracleWithin1 = 0;
    oracleWithin2 = 0;
    oracleWithin5 = 0;
    winnerCounts = zeros(numStarts, 1);
    convergedByStart = zeros(numStarts, 1);

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1-modelFlag;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } =
            _quaidsSyntheticDGP(tobs, seed, modelFlag, 1);

        aCtlStone = aCtl;
        aCtlStone.maxiter = 1;
        qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);
        bStone = qOutStone.homogB;
        nIntFull = cols(qOutStone.intcptFull);
        nGoods = cols(prices);
        nReduced = nGoods-1;
        trueRaw = trueParams[1:nIntFull, .]
            |trueParams[nIntFull+1:nIntFull+nReduced, .]
            |trueParams[nIntFull+nGoods+1:rows(trueParams), .];
        trueRawCrit = _rawBenchmarkCrit(w, trueRaw, qOutStone.intcptFull,
            prices, totexp, qOutStone.u, aCtl);

        rndseed 100000*modelFlag + 100*seed;
        bestFound = 0;
        oracleFound = 0;
        selectedCrit = 0;
        oracleErr = 0;
        foundWrongBetterTruth = 0;
        distinctCount = 0;
        representatives = 0;
        convergedThisSeed = 0;

        startIdx = 1;
        do while startIdx <= numStarts;
            bStart = _benchmarkStart(bStone, startIdx, nIntFull,
                nReduced, modelFlag);
            aCtlStart = aCtl;
            aCtlStart.b0 = bStart;

            qOutStart = _quaidsFitAnderson(w, intcpt, prices, totexp,
                instr, aCtlStart, 0, anDepth, 1, 1, rankTol,
                residualGrowth, stepFactor);

            if qOutStart.converged;
                convergedThisSeed = convergedThisSeed+1;
                totalConvergedStarts = totalConvergedStarts+1;
                convergedByStart[startIdx] = convergedByStart[startIdx]+1;
                recErrStart = maxc(maxc(abs(qOutStart.bS-trueParams)));
                rawCritStart = _rawBenchmarkCrit(w, qOutStart.homogB,
                    qOutStart.intcptFull, prices, totexp, qOutStart.u, aCtl);
                maxCritReplicaDiff = maxc(maxCritReplicaDiff
                    |abs(rawCritStart-qOutStart.homogCrit));
                { representatives, distinctCount } = _appendBenchmarkDistinct(
                    representatives, vec(qOutStart.bS), distinctCount, distinctTol);

                if startIdx == 1;
                    if recErrStart <= wrongCutoff;
                        singleCorrect = singleCorrect+1;
                    else;
                        singleWrong = singleWrong+1;
                    endif;
                endif;

                if not bestFound or qOutStart.homogCrit > selectedCrit;
                    bestFound = 1;
                    selectedCrit = qOutStart.homogCrit;
                    bestRecErr = recErrStart;
                    bestStartIdx = startIdx;
                endif;
                if not oracleFound or recErrStart < oracleErr;
                    oracleFound = 1;
                    oracleErr = recErrStart;
                endif;
                if recErrStart > wrongCutoff and rawCritStart > trueRawCrit;
                    foundWrongBetterTruth = 1;
                endif;
            elseif startIdx == 1;
                singleNever = singleNever+1;
            endif;

            startIdx = startIdx+1;
        endo;

        if distinctCount > 1;
            multiEndpointSeeds = multiEndpointSeeds+1;
        endif;
        wrongBeatsTruthSeeds = wrongBeatsTruthSeeds+foundWrongBetterTruth;

        if not bestFound;
            multiNever = multiNever+1;
            print "seed" seed "model" modelName
                "no safeguarded start converged";
        else;
            winnerCounts[bestStartIdx] = winnerCounts[bestStartIdx]+1;
            if bestRecErr <= wrongCutoff;
                multiCorrect = multiCorrect+1;
            else;
                multiWrong = multiWrong+1;
            endif;
            if oracleErr <= wrongCutoff;
                oracleCorrect = oracleCorrect+1;
            else;
                oracleWrong = oracleWrong+1;
            endif;
            if bestRecErr > wrongCutoff and oracleErr <= wrongCutoff;
                selectorMiss = selectorMiss+1;
            endif;

            selectedWithin05 = selectedWithin05+(bestRecErr <= .5);
            selectedWithin1 = selectedWithin1+(bestRecErr <= 1);
            selectedWithin2 = selectedWithin2+(bestRecErr <= 2);
            selectedWithin5 = selectedWithin5+(bestRecErr <= 5);
            oracleWithin05 = oracleWithin05+(oracleErr <= .5);
            oracleWithin1 = oracleWithin1+(oracleErr <= 1);
            oracleWithin2 = oracleWithin2+(oracleErr <= 2);
            oracleWithin5 = oracleWithin5+(oracleErr <= 5);

            print "seed" seed "model" modelName "convergedStarts"
                convergedThisSeed "distinct" distinctCount "selectedStart"
                bestStartIdx "selectedRecErr" bestRecErr "oracleRecErr"
                oracleErr "selectedCrit" selectedCrit "trueCrit" trueRawCrit;
        endif;

        seed = seed+1;
    endo;

    print;
    print "================================================================";
    print "FULL SAFEGUARDED MULTI-START:" modelName
        "seeds" nSeeds "starts" numStarts;
    print "================================================================";
    print "SAFE ANDERSON STONE-ONLY: never/wrong/correct"
        singleNever singleWrong singleCorrect;
    print "MAX-homogCrit MULTI-START: never/wrong/correct"
        multiNever multiWrong multiCorrect;
    print "ORACLE among converged starts: wrong/correct"
        oracleWrong oracleCorrect;
    print "selector misses where oracle had a correct endpoint:" selectorMiss;
    print "seeds with multiple distinct converged endpoints:" multiEndpointSeeds;
    print "total converged starts:" totalConvergedStarts
        "of" nSeeds*numStarts;
    print "seeds where any wrong endpoint beat truth criterion:"
        wrongBeatsTruthSeeds;
    print "selected recErr <= .5/1/2/5:"
        selectedWithin05 selectedWithin1 selectedWithin2 selectedWithin5;
    print "oracle recErr <= .5/1/2/5:"
        oracleWithin05 oracleWithin1 oracleWithin2 oracleWithin5;
    print "max abs raw-criterion replica difference:" maxCritReplicaDiff;
    print "converged count by start recipe:";
    print convergedByStart';
    print "selected winner count by start recipe:";
    print winnerCounts';
    print;

    modelFlag = modelFlag+1;
endo;

print "anderson_full_multistart_benchmark.e: run complete.";
