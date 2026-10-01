new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Head-to-head: same 200 seeds x 2 models as
** tests/quaids_convergence_sweep.e, same classification rule, run
** through BOTH the real quaidsFit() (baseline, relax=1) and
** _quaidsFitAnderson() with mDepth=anDepth, in one pass so the
** comparison is on identical data/seeds/aCtl settings. */

nSeeds = 200;
tobs = 3000;
structTol = 0.10;
wrongMult = 10;
anDepth = 8;

q = 0;
do while q <= 1;
    if q == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    baseNeverConv = 0;
    baseConvWrong = 0;
    baseConvCorrect = 0;
    anNeverConv = 0;
    anConvWrong = 0;
    anConvCorrect = 0;
    firstConv = 1;
    firstConvAn = 1;

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1 - q;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, q, 1);

        qOutBase = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
        qOutAn = _quaidsFitAnderson(w, intcpt, prices, totexp, instr, aCtl, 0, anDepth);

        recErrBase = maxc(maxc(abs(qOutBase.bS - trueParams)));
        recErrAn = maxc(maxc(abs(qOutAn.bS - trueParams)));

        if qOutBase.converged == 0;
            baseNeverConv = baseNeverConv + 1;
            bucketBase = "never-converged";
        else;
            if firstConv;
                iterListBase = qOutBase.iterations;
                firstConv = 0;
            else;
                iterListBase = iterListBase | qOutBase.iterations;
            endif;
            if recErrBase > wrongMult*structTol;
                baseConvWrong = baseConvWrong + 1;
                bucketBase = "converged-but-wrong";
            else;
                baseConvCorrect = baseConvCorrect + 1;
                bucketBase = "converged-correctly";
            endif;
        endif;

        if qOutAn.converged == 0;
            anNeverConv = anNeverConv + 1;
            bucketAn = "never-converged";
        else;
            if firstConvAn;
                iterListAn = qOutAn.iterations;
                firstConvAn = 0;
            else;
                iterListAn = iterListAn | qOutAn.iterations;
            endif;
            if recErrAn > wrongMult*structTol;
                anConvWrong = anConvWrong + 1;
                bucketAn = "converged-but-wrong";
            else;
                anConvCorrect = anConvCorrect + 1;
                bucketAn = "converged-correctly";
            endif;
        endif;

        print "seed" seed "model" modelName
            "  BASE conv" qOutBase.converged "iters" qOutBase.iterations "bucket" bucketBase
            "  ANDERSON conv" qOutAn.converged "iters" qOutAn.iterations "bucket" bucketAn;

        seed = seed + 1;
    endo;

    print;
    print "===================================================================";
    print "SUMMARY:" modelName "  (" nSeeds "seeds, tobs=" tobs ", anDepth=" anDepth ")";
    print "===================================================================";
    print "BASELINE (relax=1, no acceleration):";
    print "  never-converged:      " baseNeverConv "  (" 100*baseNeverConv/nSeeds "%)";
    print "  converged-but-wrong:   " baseConvWrong "  (" 100*baseConvWrong/nSeeds "%)";
    print "  converged-correctly:   " baseConvCorrect "  (" 100*baseConvCorrect/nSeeds "%)";
    print "ANDERSON(" $+ ftocv(anDepth,1,0) $+ "):";
    print "  never-converged:      " anNeverConv "  (" 100*anNeverConv/nSeeds "%)";
    print "  converged-but-wrong:   " anConvWrong "  (" 100*anConvWrong/nSeeds "%)";
    print "  converged-correctly:   " anConvCorrect "  (" 100*anConvCorrect/nSeeds "%)";
    print;

    q = q + 1;
endo;

print "anderson_sweep_compare.e: diagnostic run complete.";
