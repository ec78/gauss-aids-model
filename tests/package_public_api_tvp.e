new;

/* Installed-artifact smoke for the optional TVP surface. Source files are
   resolved from the installed package search path, not ../src. */
library cmlmt, tsmt, sslib, quaids;
#include quaids.sdf
#include quaidstvp.src
#include quaidstvpkalman.src
#include quaidstvpmle.src
#include quaidstvpsmooth.src
#include quaidstvpelas.src
#include quaidstvpfit.src

proc (0) = assert_true(ok, msg);
    if not ok;
        errorlog "package_public_api_tvp.e failed: " $+ msg;
        end;
    endif;
endp;

/* Self-contained raw-input n1=2 fixture matching the frozen D1/Stage-6
   integration case. */
proc (3) = tvpInstalledFixture(tobs, seed, state0, trueQdiag, trueH);
    local n1, ngamma, k_states, gidx, cnt, i, j, t, pricesRel, lx, pref,
        stateArr, wsub, w, state_t, gmat, Hchol, Qsd, eps, prices, totexp;

    rndseed seed;
    n1 = rows(trueH);
    ngamma = n1*(n1+1)/2;
    k_states = n1 + ngamma + n1;
    gidx = zeros(n1, n1);
    cnt = 1;
    i = 1;
    do while i <= n1;
        j = i;
        do while j <= n1;
            gidx[i, j] = cnt;
            gidx[j, i] = cnt;
            cnt = cnt + 1;
            j = j + 1;
        endo;
        i = i + 1;
    endo;

    pricesRel = .3*rndn(tobs, n1);
    lx = .3*rndn(tobs, 1);
    pref = .1*rndn(tobs, 1);
    stateArr = zeros(k_states, tobs);
    wsub = zeros(tobs, n1);
    state_t = state0;
    Qsd = sqrt(trueQdiag);
    Hchol = chol(trueH)';

    t = 1;
    do while t <= tobs;
        if t > 1;
            state_t = state_t + Qsd.*rndn(k_states, 1);
        endif;
        stateArr[., t] = state_t;
        gmat = zeros(n1, n1);
        i = 1;
        do while i <= n1;
            j = 1;
            do while j <= n1;
                gmat[i, j] = state_t[n1 + gidx[i, j]];
                j = j + 1;
            endo;
            i = i + 1;
        endo;
        eps = Hchol*rndn(n1, 1);
        i = 1;
        do while i <= n1;
            wsub[t, i] = state_t[i] + pricesRel[t, .]*gmat[i, .]' +
                state_t[n1+ngamma+i]*lx[t] + eps[i];
            i = i + 1;
        endo;
        t = t + 1;
    endo;

    w = wsub~(1-sumc(wsub'));
    prices = (pricesRel+pref)~pref;
    totexp = lx + prices*meanc(w);
    retp(w, prices, totexp);
endp;

n1 = 2;
tobs = 120;
k_states = 7;
state0 = .40|.30|.05|-.02|.04|-.04|.03;
trueQ = 1e-5*ones(n1, 1)|1e-4*ones(k_states-n1, 1);
H = 2.5e-5*eye(n1);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);
{ w, prices, totexp } = tvpInstalledFixture(tobs, 7, state0, trueQ, H);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.smooth = 0;
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 2;

struct quaidsTVPOut fit;
struct quaidsTVPDiagOut dOut;
struct quaidsTVPMultiStartOut msOut;
screen off;
fit = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
msOut = quaidsTVPMultiStartDiagnostics(w, prices, totexp, H, q0, tvpCtl);
screen on;
dOut = quaidsTVPDiagnostics(fit);

call assert_true(fit.mleRetcode == 0 and dOut.fitOk == 1 and dOut.qValid == 1,
    "cheap diagnostic did not accept the converged installed-package fit");
call assert_true(maxc(abs(fit.Qfit-msOut.qfits[., 2])) == 0 and
    fit.mleRetcode == msOut.retcodes[2],
    "fixed-grid 1x arm did not reproduce quaidsTVPFit");
call assert_true(msOut.startStable == 1 and (msOut.status == 0 or msOut.status == 5),
    "frozen n1=2 case did not reproduce observed start stability");

print "package_public_api_tvp.e: PASS";
