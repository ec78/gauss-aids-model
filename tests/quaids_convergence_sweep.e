/*
** quaids_convergence_sweep.e
**
** Milestone 12: a real, committed seed-sweep diagnostic for the iterated
** estimator's convergence reliability -- replaces the informal 8-seed
** probe referenced throughout CLAUDE.md/GOLD_STANDARD_TODO.md since
** Milestone 3, which never survived as a committed artifact (confirmed by
** searching this repo's full git history -- it doesn't exist anywhere).
**
** This is a DIAGNOSTIC REPORT GENERATOR, not a pass/fail test: it prints
** one row per (seed, model) combination plus a summary block, but does
** NOT print an "ALL N CHECKS PASSED" line and is deliberately NOT wired
** into tests/run_source_tests.ps1 -- there is no convergence guarantee to
** gate on, so treating this as a hard pass/fail check would just be a
** flaky test around a known, documented property of the estimator.
** tests/quaids_reliability_regression_test.e is the real automated test;
** this file exists to make the underlying failure rate re-measurable and
** diffable across code changes (that's the whole point -- see the
** "Two symptoms conflated" note below).
**
** Two symptoms conflated under "roughly half of seeds" in the prior
** docs are separated here into three buckets:
**   - never-converged:     qOut.converged == 0 (exhausted aCtl.maxiter)
**   - converged-but-wrong:  qOut.converged == 1, but the recovered
**                           coefficients are far (>10x the normal
**                           structural tolerance) from the true DGP
**                           parameters -- i.e. the fixed-point iteration
**                           settled into a self-consistent but WRONG
**                           answer, a materially different failure mode
**                           from simply running out of iterations.
**   - converged-correctly:  converged and within the normal tolerance
**                           band used by quaids_synthetic_validation_test.e.
**
** Near-zero-true-coefficient hypothesis: quaidsFit() does not expose its
** intermediate per-iteration coefficient vector (b0 inside the loop is
** private to src/quaids.src), so this sweep cannot directly observe
** whether b0 hits an exact zero mid-iteration. As a practical proxy, it
** instead checks whether the TRUE DGP parameters themselves (returned by
** _quaidsSyntheticDGP(), which rounds every coefficient to the nearest
** 0.1) contain a near-zero entry -- if the iteration is converging toward
** the truth at all, an exactly-zero true coefficient implies b0 will pass
** through near-zero values on the way there, which is exactly the
** scenario the guarded convergence-check formula (src/quaids.src) is
** meant to protect against.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_convergence_sweep.e
** or via tests/run_convergence_sweep.ps1, which captures the output to
** tests/convergence_sweep_report.txt (gitignored -- regenerate on demand).
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include quaidsfixtures.src;

/* Tunable sweep parameters -- tobs/aCtl.err/aCtl.maxiter match the
   settings the original informal 8-seed probe used (CLAUDE.md's
   "Seed sensitivity" note), so this sweep is a direct, comparable
   superset of it, not a different experiment. */
nSeeds = 200;
tobs = 3000;
structTol = 0.10;
wrongMult = 10;

struct quaidsControl aCtl;
struct quaidsOut qOut;

q = 0;
do while q <= 1;
    if q == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    neverConv = 0;
    convWrong = 0;
    convCorrect = 0;
    nearZeroTrueCount = 0;
    firstConv = 1;

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1 - q;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, q, 1);
        qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

        recErr = maxc(maxc(abs(qOut.bS - trueParams)));
        nearZeroTrue = sumc(sumc(abs(trueParams) .< 1e-8)) > 0;
        if nearZeroTrue;
            nearZeroTrueCount = nearZeroTrueCount + 1;
        endif;

        if qOut.converged == 0;
            neverConv = neverConv + 1;
            bucket = "never-converged";
        else;
            if firstConv;
                iterList = qOut.iterations;
                firstConv = 0;
            else;
                iterList = iterList | qOut.iterations;
            endif;
            if recErr > wrongMult*structTol;
                convWrong = convWrong + 1;
                bucket = "converged-but-wrong";
            else;
                convCorrect = convCorrect + 1;
                bucket = "converged-correctly";
            endif;
        endif;

        print "seed" seed "model" modelName "converged" qOut.converged "iters" qOut.iterations
            "finalErr" qOut.finalErr "recErr" recErr "bucket" bucket "nearZeroTrue" nearZeroTrue;

        seed = seed + 1;
    endo;

    print;
    print "===================================================================";
    print "SUMMARY:" modelName "  (" nSeeds "seeds, tobs=" tobs ")";
    print "===================================================================";
    print "never-converged:      " neverConv "  (" 100*neverConv/nSeeds "%)";
    print "converged-but-wrong:   " convWrong "  (" 100*convWrong/nSeeds "%)";
    print "converged-correctly:   " convCorrect "  (" 100*convCorrect/nSeeds "%)";
    print "seeds w/ near-zero true coefficient: " nearZeroTrueCount "  (" 100*nearZeroTrueCount/nSeeds "%)";
    if not firstConv;
        /* Manual median (not quantile()) -- quantile() errors on very
           short vectors (observed: type mismatch in GAUSS's own
           quantile.src for a 1-row input), which is a real possibility
           here since iterList can be as short as a single converged
           seed. */
        sortedIter = sortc(iterList, 1);
        nIter = rows(sortedIter);
        halfIter = trunc(nIter/2);
        if nIter - 2*halfIter == 0;
            medIter = (sortedIter[halfIter] + sortedIter[halfIter+1])/2;
        else;
            medIter = sortedIter[halfIter+1];
        endif;
        print "iterations among converged runs -- min:" minc(iterList) " median:" medIter " max:" maxc(iterList);
    else;
        print "iterations among converged runs -- N/A (no seed converged)";
    endif;
    print;

    q = q + 1;
endo;

print "quaids_convergence_sweep.e: diagnostic run complete (no pass/fail gate).";
