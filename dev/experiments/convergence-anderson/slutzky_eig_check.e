new;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;

/* Numeric (returns-a-value) replica of quaidsSlutzky()'s own eigenvalue
** computation -- validated directly against quaidsSlutzky()'s own
** printed min/max before being trusted in the bulk wrong-convergence
** check. */
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

/* --- Validate against a real quaidsFit() + quaidsSlutzky() printed run --- */
seed = 7;
tobs = 500;
aCtl = quaidsControlCreate;
aCtl.linear = 1;
aCtl.maxiter = 1;
aCtl.homogenous = 1;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, seed, 0, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

print "=== quaidsSlutzky() printed output (check its own Maximum column) ===";
call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);

myMax = _slutzkyMaxEig(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);
print "=== _slutzkyMaxEig() replica result ===";
print "my computed max eigenvalue (should match quaidsSlutzky's own printed Maximum, last row):" myMax;
