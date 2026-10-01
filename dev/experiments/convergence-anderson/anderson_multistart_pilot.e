new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Pilot: multi-start on top of Anderson(8), numStarts perturbed starts around
** the deterministic Stone/LA-AIDS starting point, selected by
** qOut.homogCrit (a real, already-shipped, no-truth-needed criterion --
** NOT recErr, which is only available here because the DGP is
** synthetic). Reports both the realistic (homogCrit-selected) outcome
** AND an oracle upper bound (best-of-numStarts by recErr) to see whether
** homogCrit is actually a trustworthy selector or just noise.
**
** Piloted on a SMALLER seed subset first (cheaper) before committing to
** the full 200-seed sweep.
*/

nSeeds = 30;
tobs = 3000;
structTol = 0.10;
wrongMult = 10;
anDepth = 8;
numStarts = 8;
pertScale = 0.3;   // relative perturbation scale around bStone

// Plain struct-to-struct copy assignments (aCtlStone = aCtl; aCtlK = aCtl;
// below) need these pre-declared -- CLAUDE.md's documented gotcha:
// inference does not propagate through a plain struct-to-struct copy.
struct quaidsControl aCtlStone;
struct quaidsControl aCtlK;

q = 0;
do while q <= 1;
    if q == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    msNeverConv = 0;
    msConvWrong = 0;
    msConvCorrect = 0;
    oracleConvWrong = 0;
    oracleConvCorrect = 0;
    anyConvergedCount = 0;

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1 - q;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, q, 1);

        // Deterministic Stone/LA-AIDS starting point, reused as the
        // perturbation center -- cheap (maxiter=1, no iteration).
        aCtlStone = aCtl;
        aCtlStone.maxiter = 1;
        qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);
        bStone = qOutStone.homogB;

        rndseed 1000*q + seed;

        bestCrit = 0;
        bestFound = 0;
        oracleRecErr = 0;
        oracleFound = 0;
        anyConverged = 0;

        k = 1;
        do while k <= numStarts;
            if k == 1;
                bStart = bStone;   // always include the unperturbed start
            else;
                bStart = bStone + pertScale*abs(bStone).*rndn(rows(bStone), cols(bStone));
            endif;

            aCtlK = aCtl;
            aCtlK.b0 = bStart;
            qOutK = _quaidsFitAnderson(w, intcpt, prices, totexp, instr, aCtlK, 0, anDepth);

            if qOutK.converged;
                anyConverged = 1;
                recErrK = maxc(maxc(abs(qOutK.bS - trueParams)));

                // Realistic (no-truth) selection: lowest homogCrit among converged.
                if not bestFound or qOutK.homogCrit < bestCrit;
                    bestCrit = qOutK.homogCrit;
                    bestRecErr = recErrK;
                    bestFound = 1;
                endif;

                // Oracle upper bound: lowest recErr among converged (uses truth).
                if not oracleFound or recErrK < oracleRecErr;
                    oracleRecErr = recErrK;
                    oracleFound = 1;
                endif;
            endif;

            k = k + 1;
        endo;

        if anyConverged;
            anyConvergedCount = anyConvergedCount + 1;
        endif;

        if not bestFound;
            msNeverConv = msNeverConv + 1;
            bucketMs = "never-converged";
        else;
            if bestRecErr > wrongMult*structTol;
                msConvWrong = msConvWrong + 1;
                bucketMs = "converged-but-wrong";
            else;
                msConvCorrect = msConvCorrect + 1;
                bucketMs = "converged-correctly";
            endif;
        endif;

        if oracleFound;
            if oracleRecErr > wrongMult*structTol;
                oracleConvWrong = oracleConvWrong + 1;
            else;
                oracleConvCorrect = oracleConvCorrect + 1;
            endif;
        endif;

        print "seed" seed "model" modelName "homogCrit-picked bucket" bucketMs
            "  oracle recErr" oracleRecErr "  picked recErr" bestRecErr;

        seed = seed + 1;
    endo;

    print;
    print "===================================================================";
    print "MULTI-START PILOT SUMMARY:" modelName "  (" nSeeds "seeds, numStarts=" numStarts ", anDepth=" anDepth ")";
    print "===================================================================";
    print "homogCrit-selected (realistic, no-truth-needed):";
    print "  never-converged:      " msNeverConv "  (" 100*msNeverConv/nSeeds "%)";
    print "  converged-but-wrong:   " msConvWrong "  (" 100*msConvWrong/nSeeds "%)";
    print "  converged-correctly:   " msConvCorrect "  (" 100*msConvCorrect/nSeeds "%)";
    print "ORACLE best-of-numStarts by recErr (upper bound, uses truth):";
    print "  converged-but-wrong:   " oracleConvWrong "  (" 100*oracleConvWrong/nSeeds "%)";
    print "  converged-correctly:   " oracleConvCorrect "  (" 100*oracleConvCorrect/nSeeds "%)";
    print "at least one of numStarts starts converged: " anyConvergedCount "  (" 100*anyConvergedCount/nSeeds "%)";
    print;

    q = q + 1;
endo;

print "anderson_multistart_pilot.e: pilot run complete.";
