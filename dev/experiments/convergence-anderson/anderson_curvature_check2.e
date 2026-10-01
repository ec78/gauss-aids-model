new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Corrected curvature-vs-correctness check, using
** _quaidsCurvatureSyntheticDGP() -- whose TRUE gamma genuinely IS
** curvature-consistent at its own sample mean by construction (unlike
** _quaidsSyntheticDGP, which is not curvature-consistent at all, the
** confound found in the first attempt). AIDS/linear only (this fixture
** has no QUAIDS version); one fixed internal seed, so diversity comes
** from many perturbed Anderson-accelerated starts on this ONE dataset
** rather than from many seeds.
**
** Known reference point from this fixture's own header: a normal,
** single-start unconstrained fit on this data recovers structure
** reasonably (max abs diff ~0.16 from truth) but still shows a SMALL
** curvature violation (~+0.17 max eigenvalue) from ordinary sampling
** noise -- so "curvature holds" should not be judged by a strict <=0
** cutoff; the real question is whether GROSS violations (orders of
** magnitude larger) concentrate in the wrong-fixed-point bucket.
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

tobs = 3000;
structTol = 0.10;
wrongMult = 10;
anDepth = 8;
numStarts = 40;
pertScale = 0.3;

struct quaidsControl aCtlStone;
struct quaidsControl aCtlK;

aCtl = quaidsControlCreate;
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.err = .0001;
aCtl.homogenous = 1;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsCurvatureSyntheticDGP(tobs);

// Reference point check: a normal single-start (Stone-seeded, no
// perturbation, no acceleration) fit's own recErr/curvature, matching
// this fixture's own documented ~0.16 / ~+0.17 reference numbers.
qOutRef = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
refRecErr = maxc(maxc(abs(qOutRef.bS - trueParams)));
refMaxEig = _slutzkyMaxEig(qOutRef.bS, qOutRef.intcptFull, prices, totexp, aCtl);
print "REFERENCE (normal single-start fit, no acceleration): converged" qOutRef.converged
    "recErr" refRecErr "maxSlutzkyEig" refMaxEig;
print "(fixture's own documented reference: recErr ~0.16, maxEig ~+0.17)";
print;

aCtlStone = aCtl;
aCtlStone.maxiter = 1;
qOutStone = quaidsFit(w, intcpt, prices, totexp, instr, aCtlStone);
bStone = qOutStone.homogB;

rndseed 4242;

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

        print "k" k "converged 1  recErr" recErrK "isCorrect" isCorrect "maxSlutzkyEig" maxEigK;
    else;
        print "k" k "converged 0 (never-converged)";
    endif;

    k = k + 1;
endo;

print;
print "anderson_curvature_check2.e: run complete.";
