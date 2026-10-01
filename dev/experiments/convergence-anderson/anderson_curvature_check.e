new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Does curvature (Slutzky negative semidefiniteness) discriminate
** "converged-correctly" from "converged-but-wrong" fixed points, where
** the model's own homogCrit fit statistic did NOT?
**
** Same 30-seed x 2-model x K=8-perturbed-Anderson-start setup as
** anderson_multistart_pilot.e, but for EVERY converged start (not just
** the best), records recErr (needs truth, synthetic-only) AND the max
** Slutzky eigenvalue (a real, theory-grounded, no-truth-needed
** diagnostic -- >0 means curvature is violated somewhere in-sample).
*/

proc (1) = _slutzkyMaxEig(b, intcpt, prices, totexp, struct quaidsControl aCtl);
    local nint, n, alpha, gama, _beta, lambda, w, wepc, mu, a_p,
        lx, lx2, b_p, nobs, i, va, maxEig;

    nint = cols(intcpt);
    n = cols(prices);

    alpha = intcpt*b[1:nint, .];
    gama = b[nint+1:nint+n, .];
    _beta = b[nint+n+1, .];

    a_p = aCtl.alpha0 + sumc((prices.*alpha)') + .5*sumc(((prices*gama).*prices)');
    lx = totexp - a_p;

    w = alpha + prices*gama + lx*_beta;
    if not aCtl.linear;
        b_p = exp(prices*_beta');
        lambda = b[nint+n+2, .];
        lx2 = (lx^2)./b_p;
        w = w + lx2*lambda;
    endif;

    nobs = rows(w);
    maxEig = -1e300;
    i = 1;
    do while i <= nobs;
        if not aCtl.linear;
            mu = lambda*lx[i]/b_p[i];
        else;
            mu = 0;
        endif;

        wepc = -diagrv(eye(n), w[i, .]') + w[i, .]'w[i, .] + gama
            + (_beta'_beta + _beta'mu + mu'_beta + 2*mu'mu)*lx[i];

        va = eigh(wepc);
        if maxc(va) > maxEig;
            maxEig = maxc(va);
        endif;
        i = i + 1;
    endo;

    retp(maxEig);
endp;

nSeeds = 30;
tobs = 3000;
structTol = 0.10;
wrongMult = 10;
anDepth = 8;
numStarts = 8;
pertScale = 0.3;

struct quaidsControl aCtlStone;
struct quaidsControl aCtlK;

// Cross-tab counters: [curvatureHolds x correctness] over EVERY converged
// (seed, model, start) triple -- the broadest, most statistically useful
// version of the check, not just the 5 disagreement cases.
nHoldsCorrect = 0;
nHoldsWrong = 0;
nViolatesCorrect = 0;
nViolatesWrong = 0;

q = 0;
do while q <= 1;
    if q == 0;
        modelName = "Iterated AIDS (linear)";
    else;
        modelName = "QUAIDS";
    endif;

    seed = 1;
    do while seed <= nSeeds;
        aCtl = quaidsControlCreate;
        aCtl.linear = 1 - q;
        aCtl.maxiter = 100;
        aCtl.err = .0001;
        aCtl.homogenous = 1;

        { w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, q, 1);

        aCtlStone = aCtl;
        aCtlStone.maxiter = 1;
        qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);
        bStone = qOutStone.homogB;

        rndseed 1000*q + seed;

        k = 1;
        do while k <= numStarts;
            if k == 1;
                bStart = bStone;
            else;
                bStart = bStone + pertScale*abs(bStone).*rndn(rows(bStone), cols(bStone));
            endif;

            aCtlK = aCtl;
            aCtlK.b0 = bStart;
            qOutK = _quaidsFitAnderson(w, intcpt, prices, totexp, instr, aCtlK, 0, anDepth);

            if qOutK.converged;
                recErrK = maxc(maxc(abs(qOutK.bS - trueParams)));
                isCorrect = recErrK <= wrongMult*structTol;

                maxEigK = _slutzkyMaxEig(qOutK.bS, qOutK.intcptFull, prices, totexp, aCtl);
                curvHolds = maxEigK <= 1e-6;   // small tolerance for numerical noise at exactly 0

                if curvHolds and isCorrect;
                    nHoldsCorrect = nHoldsCorrect + 1;
                elseif curvHolds and not isCorrect;
                    nHoldsWrong = nHoldsWrong + 1;
                elseif not curvHolds and isCorrect;
                    nViolatesCorrect = nViolatesCorrect + 1;
                else;
                    nViolatesWrong = nViolatesWrong + 1;
                endif;

                print "seed" seed "model" modelName "k" k "recErr" recErrK "isCorrect" isCorrect
                    "maxSlutzkyEig" maxEigK "curvHolds" curvHolds;
            endif;

            k = k + 1;
        endo;

        seed = seed + 1;
    endo;

    q = q + 1;
endo;

print;
print "===================================================================";
print "CURVATURE vs. CORRECTNESS CROSS-TAB (all converged starts, both models pooled)";
print "===================================================================";
print "curvature HOLDS   & correct:   " nHoldsCorrect;
print "curvature HOLDS   & wrong:     " nHoldsWrong;
print "curvature VIOLATED & correct:  " nViolatesCorrect;
print "curvature VIOLATED & wrong:    " nViolatesWrong;
print;
totalCorrect = nHoldsCorrect + nViolatesCorrect;
totalWrong = nHoldsWrong + nViolatesWrong;
print "Among CORRECT fits, % where curvature holds: " 100*nHoldsCorrect/totalCorrect;
print "Among WRONG fits,   % where curvature holds: " 100*nHoldsWrong/totalWrong;
print;
print "anderson_curvature_check.e: run complete.";
