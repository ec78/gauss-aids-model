new;

/*
** package_public_api.e
**
** Milestone 7: release-gate smoke test for the installed GAUSS package
** catalog. Unlike every other test in tests/ (which #include the source
** tree directly), this loads `library quaids;` and verifies that the
** procedures registered in package.json's src array are callable from the
** installed package -- catching load-order bugs, missing src entries, or
** a stale .lcg catalog that source-tree testing cannot catch. Adapted from
** gauss-qardl's tests/package_public_api.e.
**
** Builds its own small synthetic dataset inline (mirroring
** examples/quaids_example.e's DGP) rather than reusing
** tests/quaidsfixtures.src's private _quaidsSyntheticDGP() -- that helper
** is tests/-only source, not part of the installed package, and the point
** of this test is to exercise exactly what an installed-package consumer
** actually has available.
**
** pubtable_quaids.src is deliberately NOT exercised here: it is not in
** package.json's src array (see CLAUDE.md's "Milestone 6: reporting via
** pubtable" section), so `library quaids;` does not load it -- a caller
** who wants it still #includes src/pubtable_quaids.src directly, same as
** in the source tree. tests/quaids_pubtable_test.e already covers it.
**
** quaidsCurvatureFit()/printQuaidsCurvature() (Milestone 10) ARE in
** package.json's src array (required public API), but have a hard
** compile-time dependency on the optmt package's struct types (struct
** PV/optmtControl/optmtResults) -- confirmed empirically that `library
** quaids;` alone does NOT resolve these (errors with "Undefined
** structure 'modelResults'"), so this file now loads `library optmt,
** quaids;`, matching what every curvature-related doc page already
** documents as the required usage pattern.
**
** The main inline DGP below (seed=11) is known (CLAUDE.md's Milestone 3
** notes) to be one of the seeds for which the iterated estimator does not
** converge cleanly -- fine for exercising quaidsFit/quaids/quaidsFull/
** quaidsElasFit/quaidsSlutzky/quaidsWelfareFit/the hypothesis tests, none
** of which require convergence, but not suitable for quaidsCurvatureFit
** (which needs an already-converged, homogeneity+symmetry-constrained
** AIDS starting fit).
** A second, separate inline dataset below (seed=500, with a self-
** consistent fixed-point construction so its true gamma is curvature-
** consistent at its own sample mean) mirrors
** tests/quaidsfixtures.src's _quaidsCurvatureSyntheticDGP() -- duplicated
** here rather than reused, for the same "no dependency on tests/-only
** fixture code" reason as the main dataset above.
**
** Milestone 12 added aCtl.relax (an opt-in damping field on the existing
** quaidsControl struct, not a new proc) -- exercised here by asserting
** its default value and by running one fit with it set to a non-default
** value, confirming the installed package's .sdf/.lcg picked up the new
** field, not just that it compiles source-tree-side.
**
** Milestone 13 extended quaidsCurvatureFit() to also accept QUAIDS
** (aCtl.linear=0) fits, not just AIDS. A third inline dataset below
** (seed=204, quadratic=1/endogenous=1, mirroring
** tests/quaidsfixtures.src's _quaidsSyntheticDGP()) exercises the QUAIDS
** path -- requires aCtl.relax=.25, since QUAIDS's curvature outer loop
** is measurably less stable than AIDS's own (undamped runs on this
** fixture diverge, confirmed directly -- see src/quaidscurvature.src's
** header and tests/quaids_curvature_test.e for the full story).
**
** Milestone 21 adds quaidsWorkflowFit()/quaidsWorkflowScenarioFit(), both
** exercised against the converged seed=500 AIDS fixture used by the
** curvature block below. Milestone 22 adds
** quaidsRobustCovariance()/quaidsRobustBootstrapCovariance(), exercising
** robust/cluster covariance propagation from the installed package.
** Milestone 23 adds quaidsPreflight()/printQuaidsPreflight(); Milestone 24
** echoes a compact preflight summary from quaidsWorkflowFit(). Milestone 25
** adds quaidsSurveyWorkflowFit() for sampling-weighted workflow evaluation
** points.
**
** Run this after building/installing the package (see
** scripts/run_release_verification.ps1 -InstallArtifact).
*/

library optmt, quaids;

proc (0) = assert_true(ok, msg);
    if not ok;
        errorlog "package_public_api.e failed: " $+ msg;
        end;
    endif;
endp;

/* Same 5-good synthetic DGP shape as examples/quaids_example.e (homogeneity
   and adding-up true by construction), inlined so this test has no
   dependency on tests/-only fixture code. */
seed = 11;
tobs = 1000;
N = 5;

al = round(rndns(1, N-1, seed)*10)/10;
al = al~(1-sumc(al'));
al1 = .5*round(rndns(1, N-1, seed)*10)/10;
al1 = al1~(-sumc(al1'));
ga = round(rndns(N-1, N-1, seed)*10)/10;
ga = xpnd(vech(ga));
ga = ga|(-sumc(ga)');
ga = ga~(-sumc(ga'));
be = .5*round(rndns(1, N-1, seed)*10)/10;
be = be~(-sumc(be'));
la = .01*round(rndns(1, N-1, seed)*10)/10;
la = la~(-sumc(la'));
ro = round(rndns(1, N-1, seed)*10)/10;

prices = 1+rndns(tobs, N, seed);
instr = 5+5*rndns(tobs, 1, seed);
intcpt = 2+2*rndns(tobs, 1, seed);
u = .1*rndns(tobs, 1, seed);
totexp = .85*instr + u;
e = 2*rndns(tobs, N-1, seed) + u*ro;
e = e~(-sumc(e'));

a_p = sumc((prices.*(al+intcpt*al1))') + .5*sumc(((prices*ga).*prices)');
lx = totexp - a_p;
b_p = prices*be';
lx2 = (lx^2)./exp(b_p);

w = al + prices*ga + lx*be + e + intcpt*al1 + lx2*la;


/* --- quaidsControlCreate() / getDefaultQuaidsControl() --- */

struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
call assert_true(aCtl.maxiter > 1 and aCtl.homogenous == 1, "quaidsControlCreate defaults look wrong");
call assert_true(aCtl.relax == 1, "quaidsControlCreate: aCtl.relax (Milestone 12) default should be 1 (no damping)");

struct quaidsControl aCtlAlias;
aCtlAlias = getDefaultQuaidsControl();
call assert_true(aCtlAlias.maxiter == aCtl.maxiter and aCtlAlias.homogenous == aCtl.homogenous,
    "getDefaultQuaidsControl does not match quaidsControlCreate");

aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .001;


/* --- quaidsFit() / quaids() --- */

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(qOut.model $== "QUAIDS" and qOut.n == N and qOut.nobs == tobs,
    "quaidsFit metadata invalid");
call assert_true(rows(qOut.bestB) > 0 and rows(qOut.bestV) == rows(qOut.bestB)*N and cols(qOut.bestV) == rows(qOut.bestB)*N,
    "quaidsFit bestB/bestV shape invalid");

/* aCtl.relax (Milestone 12): confirm the installed package's .sdf/.lcg
   actually picked up the new quaidsControl field, not just that it
   compiles source-tree-side. */
struct quaidsOut qOutRelax;
struct quaidsControl aCtlRelax;
aCtlRelax = aCtl;
aCtlRelax.relax = .75;
qOutRelax = quaidsFit(w, intcpt, prices, totexp, instr, aCtlRelax);
call assert_true(rows(qOutRelax.bestB) == rows(qOut.bestB), "quaidsFit with aCtl.relax=.75 (installed package) produced a valid fit");

struct quaidsOut qOutLA;
struct quaidsControl aCtlLA;
aCtlLA = quaidsControlCreate;
aCtlLA.linear = 1;
aCtlLA.maxiter = 1;
qOutLA = quaidsFit(w, intcpt, prices, totexp, instr, aCtlLA);
call assert_true(qOutLA.model $== "LA-AIDS", "quaidsFit LA-AIDS model tag invalid");

{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(maxc(maxc(abs(b2 - qOut.bS))) == 0 and maxc(maxc(abs(v2 - qOut.vS))) == 0,
    "legacy quaids() wrapper output does not match quaidsFit()");

call printQuaids(qOut);


/* --- quaidsFull(): dataframe entry point --- */

xnames = "W1"$|"W2"$|"W3"$|"W4"$|"W5";
pnames = "P1"$|"P2"$|"P3"$|"P4"$|"P5";

data = asDF(w[., 1], "W1");
data = dfaddcol(data, "W2", w[., 2]);
data = dfaddcol(data, "W3", w[., 3]);
data = dfaddcol(data, "W4", w[., 4]);
data = dfaddcol(data, "W5", w[., 5]);
data = dfaddcol(data, "P1", prices[., 1]);
data = dfaddcol(data, "P2", prices[., 2]);
data = dfaddcol(data, "P3", prices[., 3]);
data = dfaddcol(data, "P4", prices[., 4]);
data = dfaddcol(data, "P5", prices[., 5]);
data = dfaddcol(data, "TOTEXP", totexp);
data = dfaddcol(data, "Z1", instr);
data = dfaddcol(data, "X1", intcpt);

struct quaidsOut qOutFull;
qOutFull = quaidsFull(data, xnames, pnames, "TOTEXP", "Z1", "X1", aCtl);
call assert_true(maxc(maxc(abs(qOutFull.bestB - qOut.bestB))) == 0,
    "quaidsFull does not match quaidsFit on the same data");


/* --- quaidsElasFit() / quaidsElas() / printQuaidsElas() --- */

nint = qOut.nint;
n = qOut.n;
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call assert_true(rows(elasOut.er) == N and rows(elasOut.ep) == N and cols(elasOut.ep) == N,
    "quaidsElasFit output shape invalid");

/* Engel aggregation (sum_i w_i*er_i == 1) is a property of the MODEL-
   IMPLIED share at the evaluation point, not the noisy empirical mean
   share meanc(w). Milestone 16: uses the real, installed
   quaidsSharesFit() rather than a fourth hand-inlined copy of the share
   formula (this file used to recompute it inline here, duplicating
   quaidsElas_() and tests/quaids_elasticities_test.e's own former
   modelShareAt() helper). */
struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call assert_true(rows(sharesOut.w) == N and rows(sharesOut.se) == N and rows(sharesOut.v) == N and cols(sharesOut.v) == N,
    "quaidsSharesFit output shape invalid");
call assert_true(abs(sumc(sharesOut.w) - 1) < 1e-6,
    "quaidsSharesFit predicted shares fail adding-up (sum to 1)");
call assert_true(abs(sumc(sharesOut.w.*elasOut.er) - 1) < 1e-6,
    "quaidsElasFit income elasticities fail Engel aggregation");

call printQuaidsShares(sharesOut);
call printQuaidsElas(elasOut);

/* quaidsElas() is a thin fit-then-print wrapper (see
   src/quaidselas.src) -- calling it here, not just its split
   quaidsElasFit()/printQuaidsElas() halves, closes a real gap this test
   previously had: quaidsElas() and printQuaids() are both required,
   package.json-listed public procs that were never actually exercised
   through `library quaids;` until Milestone 9's final integration pass
   found it. */
call quaidsElas(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);


/* --- quaidsSlutzky() (prints a report, no return value) --- */

call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);


/* --- quaidsWelfareFit() / printQuaidsWelfare() (Milestone 11) ---
   Needs no convergence, just a fitted bestB/bestV, so the seed=11
   dataset above (known non-converging for the iterated estimator) is
   fine here -- unlike quaidsCurvatureFit below. */

pricesPt1 = pricesMean;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

struct quaidsWelfareOut wOut;
wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, pricesPt1, totexpMean, aCtl);
call assert_true(wOut.seCV >= 0 and wOut.seEV >= 0, "quaidsWelfareFit standard errors invalid");

struct quaidsWelfareOut wOutZero;
wOutZero = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, pricesMean, totexpMean, aCtl);
call assert_true(wOutZero.cv == 0 and wOutZero.ev == 0,
    "quaidsWelfareFit zero-price-change identity failed");

call printQuaidsWelfare(wOut);


/* --- quaidsHomogeneityTest() / quaidsJointTest() (need an unconstrained fit) --- */

struct quaidsControl aCtlU;
aCtlU = quaidsControlCreate;
aCtlU.linear = 0;
aCtlU.maxiter = 100;
aCtlU.homogenous = 0;
aCtlU.err = .001;

struct quaidsOut qOutU;
qOutU = quaidsFit(w, intcpt, prices, totexp, instr, aCtlU);

{ statH, pvalH, dfH } = quaidsHomogeneityTest(qOutU);
call assert_true(pvalH >= 0 and pvalH <= 1 and dfH == N-1, "quaidsHomogeneityTest output invalid");

{ statJ, pvalJ, dfJ } = quaidsJointTest(qOutU);
call assert_true(pvalJ >= 0 and pvalJ <= 1 and dfJ == (N-1) + (N-1)*(N-2)/2, "quaidsJointTest output invalid");

/* quaidsQuadraticTest() (Milestone 17) -- reuses qOutU, already an
   unconstrained QUAIDS (aCtlU.linear=0) fit, exactly what this test
   requires. */
{ statQ2, pvalQ2, dfQ2 } = quaidsQuadraticTest(qOutU);
call assert_true(pvalQ2 >= 0 and pvalQ2 <= 1 and dfQ2 == N-1, "quaidsQuadraticTest output invalid");


/* --- quaidsCurvatureFit() / printQuaidsCurvature() (Milestone 10) ---
   Needs its own dataset: the seed=11 fixture above is a known non-
   converging seed for the iterated estimator (fine for everything above,
   since none of it requires convergence, but quaidsCurvatureFit needs an
   already-converged AIDS starting fit). Mirrors
   tests/quaidsfixtures.src's _quaidsCurvatureSyntheticDGP() -- see that
   file's header comment for why the self-consistent fixed-point
   construction and seed=500 are necessary, not arbitrary. */

Nc = 5;
n1c = Nc - 1;
seedC = 500;
aScaleC = 0.15;
nroundsC = 8;
tobsC = 3000;

AtrueC = aScaleC * lowmat(rndns(n1c, n1c, seedC));

alC = round(rndns(1, Nc-1, seedC)*10)/10;
alC = alC~(1-sumc(alC'));
beC = .5*round(rndns(1, Nc-1, seedC)*10)/10;
beC = beC~(-sumc(beC'));
roC = round(rndns(1, Nc-1, seedC)*10)/10;

pricesC = 1+rndns(tobsC, Nc, seedC);
instrC = 5+5*rndns(tobsC, 1, seedC);
uC = .1*rndns(tobsC, 1, seedC);
totexpC = .85*instrC + uC;
eC = 2*rndns(tobsC, Nc-1, seedC) + uC*roC;
eC = eC~(-sumc(eC'));

wbarC = ones(Nc, 1)/Nc;
lxbarC = 3.25;
rC = 1;
do while rC <= nroundsC;
    K0C = -diagrv(eye(Nc), wbarC) + wbarC*wbarC' + beC'beC*lxbarC;
    CC = -AtrueC*AtrueC' - K0C[1:n1c, 1:n1c];
    gaC = CC | (-sumc(CC)');
    gaC = gaC ~ (-sumc(gaC'));

    a_pC = sumc((pricesC.*alC)') + .5*sumc(((pricesC*gaC).*pricesC)');
    lxC = totexpC - a_pC;
    wC = alC + pricesC*gaC + lxC*beC + eC;

    wbarC = meanc(wC);
    lxbarC = meanc(lxC);
    rC = rC + 1;
endo;

struct quaidsControl aCtlC;
aCtlC = quaidsControlCreate();
aCtlC.linear = 1;
aCtlC.maxiter = 100;
aCtlC.homogenous = 1;
aCtlC.err = .0001;

struct quaidsOut qOutC;
qOutC = quaidsFit(wC, 0, pricesC, totexpC, instrC, aCtlC);
call assert_true(qOutC.converged == 1, "quaidsCurvatureFit prerequisite AIDS fit did not converge");

/* --- Milestone 26: quaidsFit()'s own optional estimator-level sampling
   weight argument, against the real installed package. An explicit
   uniform weight is used (not a fresh non-uniform draw) so this stays a
   release gate exercising the new argument end-to-end, not a second
   convergence/recovery validation -- that already lives in
   tests/quaids_survey_test.e. */
weightC = ones(tobsC, 1);
struct quaidsOut qOutCWeighted;
qOutCWeighted = quaidsFit(wC, 0, pricesC, totexpC, instrC, aCtlC, weightC);
call assert_true(qOutCWeighted.weighted == 1 and qOutCWeighted.weightSum == tobsC and qOutCWeighted.effN == tobsC,
    "quaidsFit: weight diagnostics invalid for an explicit uniform weight");
call assert_true(maxc(maxc(abs(qOutCWeighted.bestB - qOutC.bestB))) == 0,
    "quaidsFit: an explicit uniform weight did not reproduce the unweighted fit exactly");

/* --- quaidsPreflight() / printQuaidsPreflight() (Milestone 23) --- */
struct quaidsPreflightOut pOutC;
pOutC = quaidsPreflight(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, 0);
call assert_true(pOutC.dimensionsOk == 1 and pOutC.designInvOk == 1,
    "quaidsPreflight: converged seed=500 fixture has invalid dimensions or design");
call assert_true(pOutC.ivDiagnosticsValid == 1 and rows(pOutC.ivFstat) == 1,
    "quaidsPreflight: IV diagnostics invalid");
call assert_true(pOutC.clusterValid == 1 and pOutC.nClusters == tobsC,
    "quaidsPreflight: clusterId=0 robust cluster summary invalid");
call assert_true(pOutC.weightValid == 1 and pOutC.weightSum == tobsC and pOutC.effN == tobsC,
    "quaidsPreflight: weight=0 (unweighted) diagnostics invalid");

struct quaidsPreflightOut pOutCWeighted;
pOutCWeighted = quaidsPreflight(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, weightC);
call assert_true(pOutCWeighted.weightValid == 1 and pOutCWeighted.weightSum == tobsC and pOutCWeighted.effN == tobsC,
    "quaidsPreflight: weight diagnostics invalid for an explicit uniform weight");

call printQuaidsPreflight(pOutC);

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOutC, wC, pricesC, totexpC, aCtlC);
call assert_true(cOut.converged == 1, "quaidsCurvatureFit did not converge");
call assert_true(maxc(cOut.eigenvalues) < 1e-3,
    "quaidsCurvatureFit: Slutzky matrix at the reference point is not negative semidefinite");

call printQuaidsCurvature(cOut);

/* --- quaidsCurvatureBootstrapFit() / printQuaidsCurvatureBootstrap()
   (Milestone 15) --- B kept tiny (this is a release gate, not a
   statistical validation -- see tests/quaids_curvature_bootstrap_test.e
   for the real check) since it refits the whole pipeline B times. */
struct quaidsCurvBootOut bootOut;
bootOut = quaidsCurvatureBootstrapFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 2, 42);
call assert_true(bootOut.nCompleted >= 1, "quaidsCurvatureBootstrapFit: no replications completed");
call assert_true(rows(bootOut.seBoot) == rows(cOut.b) and cols(bootOut.seBoot) == cols(cOut.b),
    "quaidsCurvatureBootstrapFit: seBoot shape does not match cOut.b");

call printQuaidsCurvatureBootstrap(bootOut);

/* quaidsCurvatureBootstrapCI() (Milestone 18) -- percentile CIs from
   bootOut's already-computed raw draws, no new resampling. */
{ ciLowerBoot, ciUpperBoot } = quaidsCurvatureBootstrapCI(bootOut, 0.05);
call assert_true(rows(ciLowerBoot) == rows(cOut.b) and cols(ciLowerBoot) == cols(cOut.b),
    "quaidsCurvatureBootstrapCI: ciLower shape does not match cOut.b");
call assert_true(minc(minc(ciUpperBoot - ciLowerBoot)) >= 0,
    "quaidsCurvatureBootstrapCI: ciUpper must be >= ciLower elementwise");


/* --- quaidsCurvatureFit() on QUAIDS (Milestone 13) ---
   Needs aCtl.relax (Milestone 12) -- QUAIDS's curvature outer loop is
   measurably less stable than AIDS's own (see src/quaidscurvature.src's
   header); undamped runs on this fixture diverge, confirmed directly.
   Mirrors tests/quaidsfixtures.src's _quaidsSyntheticDGP(seed=204,
   quadratic=1, endogenous=1) -- duplicated inline rather than reused,
   per this file's own "no dependency on tests/-only fixture code"
   principle, same as the AIDS block above. */

Nq = 5;
seedQ = 204;
tobsQ = 3000;

alQ = round(rndns(1, Nq-1, seedQ)*10)/10;
alQ = alQ~(1-sumc(alQ'));
al1Q = .5*round(rndns(1, Nq-1, seedQ)*10)/10;
al1Q = al1Q~(-sumc(al1Q'));
gaQ0 = round(rndns(Nq-1, Nq-1, seedQ)*10)/10;
gaQ0 = xpnd(vech(gaQ0));
gaQ0 = gaQ0|(-sumc(gaQ0)');
gaQ0 = gaQ0~(-sumc(gaQ0'));
beQ = .5*round(rndns(1, Nq-1, seedQ)*10)/10;
beQ = beQ~(-sumc(beQ'));
laQ = .01*round(rndns(1, Nq-1, seedQ)*10)/10;
laQ = laQ~(-sumc(laQ'));
roQ = round(rndns(1, Nq-1, seedQ)*10)/10;

pricesQ = 1+rndns(tobsQ, Nq, seedQ);
instrQ = 5+5*rndns(tobsQ, 1, seedQ);
intcptQ = 2+2*rndns(tobsQ, 1, seedQ);
uQ = .1*rndns(tobsQ, 1, seedQ);
totexpQ = .85*instrQ + uQ;
eNoiseQ = 2*rndns(tobsQ, Nq-1, seedQ) + uQ*roQ;
eNoiseQ = eNoiseQ~(-sumc(eNoiseQ'));

a_pQ = sumc((pricesQ.*(alQ+intcptQ*al1Q))') + .5*sumc(((pricesQ*gaQ0).*pricesQ)');
lxQ = totexpQ - a_pQ;
b_pQ = pricesQ*beQ';
lx2Q = (lxQ^2)./exp(b_pQ);

wQ = alQ + pricesQ*gaQ0 + lxQ*beQ + eNoiseQ + intcptQ*al1Q + lx2Q*laQ;

struct quaidsControl aCtlQ;
aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 1;
aCtlQ.err = .0001;

struct quaidsOut qOutQ;
qOutQ = quaidsFit(wQ, intcptQ, pricesQ, totexpQ, instrQ, aCtlQ);
call assert_true(qOutQ.converged == 1, "quaidsCurvatureFit (QUAIDS) prerequisite fit did not converge");

struct quaidsControl aCtlCurvQ;
aCtlCurvQ = aCtlQ;
aCtlCurvQ.relax = .25;

struct quaidsCurvOut cOutQ;
cOutQ = quaidsCurvatureFit(qOutQ, wQ, pricesQ, totexpQ, aCtlCurvQ);
call assert_true(cOutQ.converged == 1, "quaidsCurvatureFit (QUAIDS) did not converge");
call assert_true(maxc(cOutQ.eigenvalues) < 1e-3,
    "quaidsCurvatureFit (QUAIDS): Slutzky matrix at the reference point is not negative semidefinite");
call assert_true(rows(cOutQ.b) == 1+qOutQ.nint+Nq+2, "quaidsCurvatureFit (QUAIDS): cOut.b missing the lambda row");

call printQuaidsCurvature(cOutQ);


/* --- quaidsZeroFit() / printQuaidsZero() (Milestone 19) ---
   Needs a dataset with a genuine zero-budget-share censoring pattern.
   Mirrors tests/quaidsfixtures.src's _quaidsZeroSyntheticDGP(3000, 1) --
   duplicated inline rather than reused, per this file's own "no
   dependency on tests/-only fixture code" principle, same as the other
   dedicated blocks above. See that fixture's own header for why the
   structural coefficients are scaled down (0.25x) and why seed=1 was
   chosen (direct screening, not arbitrary). */

Nz = 5;
seedZ = 1;
tobsZ = 3000;

alZ = round(rndns(1, Nz-1, seedZ)*10)/10;
alZ = alZ~(1-sumc(alZ'));
al1Z = .25*.5*round(rndns(1, Nz-1, seedZ)*10)/10;
al1Z = al1Z~(-sumc(al1Z'));
gaZ = .25*round(rndns(Nz-1, Nz-1, seedZ)*10)/10;
gaZ = xpnd(vech(gaZ));
gaZ = gaZ|(-sumc(gaZ)');
gaZ = gaZ~(-sumc(gaZ'));
beZ = .25*.5*round(rndns(1, Nz-1, seedZ)*10)/10;
beZ = beZ~(-sumc(beZ'));
laZ = .25*.01*round(rndns(1, Nz-1, seedZ)*10)/10;
laZ = laZ~(-sumc(laZ'));
roZ = round(rndns(1, Nz-1, seedZ)*10)/10;

pricesZ = 1+rndns(tobsZ, Nz, seedZ);
instrZ = 5+5*rndns(tobsZ, 1, seedZ);
intcptZ = 2+2*rndns(tobsZ, 1, seedZ);
uZ = .1*rndns(tobsZ, 1, seedZ);
totexpZ = .85*instrZ + uZ;
eZ = 0.04*rndns(tobsZ, Nz-1, seedZ) + uZ*roZ;
eZ = eZ~(-sumc(eZ'));

a_pZ = sumc((pricesZ.*(alZ+intcptZ*al1Z))') + .5*sumc(((pricesZ*gaZ).*pricesZ)');
lxZ = totexpZ - a_pZ;
b_pZ = pricesZ*beZ';
lx2Z = (lxZ^2)./exp(b_pZ);

wLatentZ = alZ + pricesZ*gaZ + lxZ*beZ + eZ + intcptZ*al1Z + lx2Z*laZ;

wZ = wLatentZ .* (wLatentZ .> 0);
wZ = wZ ./ sumc(wZ');

struct quaidsControl aCtlZ;
aCtlZ = quaidsControlCreate();
aCtlZ.linear = 0;
aCtlZ.maxiter = 100;
aCtlZ.homogenous = 0;
aCtlZ.err = .0001;

struct quaidsZeroOut zOut;
zOut = quaidsZeroFit(wZ, intcptZ, pricesZ, totexpZ, instrZ, aCtlZ);
call assert_true(zOut.converged == 1, "quaidsZeroFit did not converge");
call assert_true(rows(zOut.probitB) > 0 and cols(zOut.probitB) == Nz, "quaidsZeroFit probitB shape invalid");
call assert_true(sumc(zOut.probitConverged) == Nz, "quaidsZeroFit: not all first-stage probits converged");

call printQuaidsZero(zOut);


/* --- quaidsRobustFit() / printQuaidsRobust() / quaidsRobustBootstrapFit()
   / printQuaidsRobustBootstrap() (Milestone 20) plus full-basis robust
   covariance expansion helpers (Milestone 22) --- both robust paths need
   a converged base fit: quaidsRobustFit() requires the caller's qOut to
   have converged, and quaidsRobustBootstrapFit() requires its internal
   base quaidsFit() call to converge. This block therefore reuses the
   already-converging seed=500 AIDS fixture (qOutC/wC/pricesC/totexpC/
   instrC/aCtlC) from the quaidsCurvatureFit block above. clusterId=0
   (heteroskedasticity-robust) is exercised here -- cluster-robust
   behavior itself is already thoroughly validated in
   tests/quaids_robust_test.e/quaids_robust_bootstrap_test.e, this is a
   release gate, not a re-validation. B kept tiny, same reasoning as the
   quaidsCurvatureBootstrapFit block above. */

struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOutC, wC, pricesC, totexpC, aCtlC, 0);
call assert_true(rOut.nClusters == qOutC.nobs, "quaidsRobustFit: nClusters == nobs when clusterId=0");
call assert_true(rows(rOut.se) == rows(rOut.b) and cols(rOut.se) == cols(rOut.b), "quaidsRobustFit: se shape does not match b");

{ bRFull, vRFull } = quaidsRobustCovariance(qOutC, rOut, aCtlC);
call assert_true(rows(bRFull) == rows(qOutC.bestB) and cols(bRFull) == cols(qOutC.bestB),
    "quaidsRobustCovariance: expanded b shape does not match qOutC.bestB");
call assert_true(rows(vRFull) == rows(qOutC.bestV) and cols(vRFull) == cols(qOutC.bestV),
    "quaidsRobustCovariance: expanded v shape does not match qOutC.bestV");

call printQuaidsRobust(rOut);

/* Milestone 26: quaidsRobustFit()'s own optional weight argument -- an
   explicit uniform weight must reproduce the unweighted sandwich exactly. */
struct quaidsRobustOut rOutWeighted;
rOutWeighted = quaidsRobustFit(qOutC, wC, pricesC, totexpC, aCtlC, 0, weightC);
call assert_true(maxc(maxc(abs(rOutWeighted.se - rOut.se))) == 0,
    "quaidsRobustFit: an explicit uniform weight did not reproduce the unweighted sandwich exactly");

struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, 2, 42);
call assert_true(rbOut.nCompleted >= 2, "quaidsRobustBootstrapFit: fewer than two replications completed");
call assert_true(rows(rbOut.seBoot) == rows(rbOut.b) and cols(rbOut.seBoot) == cols(rbOut.b),
    "quaidsRobustBootstrapFit: seBoot shape does not match rbOut.b");

/* Milestone 26: same optional weight argument on the bootstrap variant. */
struct quaidsRobustBootOut rbOutWeighted;
rbOutWeighted = quaidsRobustBootstrapFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, 2, 42, weightC);
call assert_true(maxc(maxc(abs(rbOutWeighted.b - rbOut.b))) == 0 and maxc(maxc(abs(rbOutWeighted.seRobust - rbOut.seRobust))) == 0,
    "quaidsRobustBootstrapFit: an explicit uniform weight did not reproduce the unweighted base fit/sandwich exactly");

{ bRBFull, vRBFull } = quaidsRobustBootstrapCovariance(qOutC, rbOut, aCtlC);
call assert_true(rows(bRBFull) == rows(qOutC.bestB) and cols(bRBFull) == cols(qOutC.bestB),
    "quaidsRobustBootstrapCovariance: expanded b shape does not match qOutC.bestB");
call assert_true(rows(vRBFull) == rows(qOutC.bestV) and cols(vRBFull) == cols(qOutC.bestV),
    "quaidsRobustBootstrapCovariance: expanded v shape does not match qOutC.bestV");

call printQuaidsRobustBootstrap(rbOut);


/* --- quaidsWorkflowFit() / quaidsWorkflowScenarioFit() (Milestone 21) ---
   Uses the same converged seed=500 AIDS fixture as the curvature and
   robust-bootstrap blocks, so post-estimation and robust outputs should
   be valid. */

struct quaidsWorkflowOut wfOut;
wfOut = quaidsWorkflowFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0);
call assert_true(wfOut.postValid == 1 and wfOut.robustValid == 1,
    "quaidsWorkflowFit: post/robust outputs were not computed");
call assert_true(wfOut.postRobustValid == 1 and rows(wfOut.sharesRobustSE) == Nc,
    "quaidsWorkflowFit: robust post-estimation outputs were not computed");
call assert_true(rows(wfOut.shares) == Nc and rows(wfOut.incomeElas) == Nc,
    "quaidsWorkflowFit: post-estimation output shapes invalid");
call assert_true(wfOut.preflightOk == pOutC.ok and wfOut.preflightDesignInvOk == pOutC.designInvOk,
    "quaidsWorkflowFit: preflight status does not match quaidsPreflight");
call assert_true(wfOut.preflightIVFstat == pOutC.ivFstat and wfOut.preflightNClusters == pOutC.nClusters,
    "quaidsWorkflowFit: preflight IV/cluster summary does not match quaidsPreflight");

/* Milestone 26: quaidsWorkflowFit()'s own optional estimator-level weight
   argument -- an explicit uniform weight must reproduce the unweighted
   workflow fit exactly. */
struct quaidsWorkflowOut wfOutWeighted;
wfOutWeighted = quaidsWorkflowFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, weightC);
call assert_true(wfOutWeighted.weighted == 1 and wfOutWeighted.weightSum == tobsC,
    "quaidsWorkflowFit: weight diagnostics invalid for an explicit uniform weight");
call assert_true(maxc(maxc(abs(wfOutWeighted.bestB - wfOut.bestB))) == 0,
    "quaidsWorkflowFit: an explicit uniform weight did not reproduce the unweighted workflow fit exactly");

mWF = meanc(qOutC.intcptFull~pricesC~totexpC);
intcptWF = mWF[1:1+qOutC.nint];
pricesWF0 = mWF[1+qOutC.nint+1:1+qOutC.nint+qOutC.n];
totexpWF0 = mWF[1+qOutC.nint+qOutC.n+1];
pricesWF1 = pricesWF0;
pricesWF1[1] = pricesWF1[1] + ln(1.02);

weightWF = seqa(1, 1, tobsC);
weightWF = weightWF/sumc(weightWF)*tobsC;
struct quaidsWorkflowOut wfSurvey;
wfSurvey = quaidsSurveyWorkflowFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, weightWF);
call assert_true(wfSurvey.surveyWeighted == 1 and wfSurvey.surveyWeightValid == 1,
    "quaidsSurveyWorkflowFit: survey weight diagnostics invalid");
call assert_true(wfSurvey.postValid == 1 and rows(wfSurvey.shares) == Nc and rows(wfSurvey.incomeElas) == Nc,
    "quaidsSurveyWorkflowFit: weighted post-estimation outputs invalid");
call assert_true(abs(wfSurvey.surveyWeightSum - tobsC) < 1e-8 and wfSurvey.surveyWeightNPositive == tobsC,
    "quaidsSurveyWorkflowFit: survey weight summary does not match input");

struct quaidsWorkflowOut wfScenario;
wfScenario = quaidsWorkflowScenarioFit(wC, 0, pricesC, totexpC, instrC, aCtlC, 0, intcptWF, pricesWF0, pricesWF1, totexpWF0);
call assert_true(wfScenario.welfareValid == 1 and wfScenario.seCV >= 0 and wfScenario.seEV >= 0,
    "quaidsWorkflowScenarioFit: welfare outputs invalid");
call assert_true(wfScenario.welfareRobustValid == 1 and wfScenario.seCVRobust >= 0 and wfScenario.seEVRobust >= 0,
    "quaidsWorkflowScenarioFit: robust welfare outputs invalid");


/* --- quaidsReplicateWeightFit() / printQuaidsReplicateWeight()
   (Milestone 27) --- reuses the same converged seed=500 AIDS fixture. A
   small (3-column) JK1-style delete-one-cluster replicate design is
   built inline; an explicit uniform weight (weightC, already defined
   above) is used as the base weight for a real, non-vacuous release-gate
   exercise rather than the scalar-0 unweighted path already covered by
   tests/quaids_replicate_test.e. */
nRepC = 3;
clusterIdC = ceil(seqa(1, 1, tobsC)/(tobsC/nRepC));
replicateWeightsC = zeros(tobsC, nRepC);
gC = 1;
do while gC <= nRepC;
    replicateWeightsC[., gC] = (clusterIdC ./= gC) * nRepC/(nRepC-1);
    gC = gC + 1;
endo;

struct quaidsReplicateOut rOutC;
rOutC = quaidsReplicateWeightFit(wC, 0, pricesC, totexpC, instrC, aCtlC, weightC, replicateWeightsC, (nRepC-1)/nRepC, "JK1");
call assert_true(rOutC.weighted == 1 and rOutC.nRequested == nRepC,
    "quaidsReplicateWeightFit: weighted/nRequested diagnostics invalid");
call assert_true(rows(rOutC.se) == rows(rOutC.b) and cols(rOutC.se) == cols(rOutC.b),
    "quaidsReplicateWeightFit: se shape does not match b");

call printQuaidsReplicateWeight(rOutC);


print "package_public_api.e: PASS";
