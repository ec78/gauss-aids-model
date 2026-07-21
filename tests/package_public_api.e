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
** Run this after building/installing the package (see
** scripts/run_release_verification.ps1 -InstallArtifact).
*/

library quaids;

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
call assert_true(abs(sumc(meanc(w).*elasOut.er) - 1) < 1e-6,
    "quaidsElasFit income elasticities fail Engel aggregation");

call printQuaidsElas(elasOut);


/* --- quaidsSlutzky() (prints a report, no return value) --- */

call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);


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


print "package_public_api.e: PASS";
