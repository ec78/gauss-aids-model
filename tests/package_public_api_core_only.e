new;

/*
** package_public_api_core_only.e
**
** Public release roadmap PR-101 acceptance evidence: "A clean machine
** without QUAIDS or optmt can install and run the core quick start using
** the documented procedure." tests/package_public_api.e (the original
** installed-package release gate) loads `library optmt, quaids;` for its
** own curvature-testing block, so it cannot demonstrate core
** independence from optmt on its own -- this file exists specifically to
** prove that independence, with `optmt` never loaded or referenced
** anywhere in this script.
**
** A genuinely separate file, not a refactor splitting
** package_public_api.e in two: an earlier attempt to issue `library
** quaids;` alone, call quaidsFit() for real, and only THEN add a second
** `library optmt;` statement later in the SAME script (immediately before
** an isolated curvature block) was tested directly and found to break
** compilation of quaids.src's own cross-file references (`error G0025:
** Undefined symbol '_quaidsIVFirstStage'`/`'quaidsElas'`/`'quaidsSlutzky'`)
** -- a real, reproducible GAUSS quirk where a second `library` statement
** appearing anywhere in a script disrupts symbol resolution for an
** earlier `library`-loaded package's own cross-file calls, not something
** to work around by careful within-file ordering.
**
** Deliberately covers a representative slice of core functionality, not
** an exhaustive re-test of package_public_api.e's own ~600 assertions --
** the point is proving `library quaids;` alone is sufficient, not
** re-validating every proc's behavior a second time.
**
** Run this after building/installing the package (see
** scripts/run_release_verification.ps1 -InstallArtifact).
*/

library quaids;

proc (0) = assert_true(ok, msg);
    if not ok;
        errorlog "package_public_api_core_only.e failed: " $+ msg;
        end;
    endif;
endp;

/* Same 5-good synthetic DGP shape as examples/01_basic_estimation.e. */
seed = 204;
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


/* --- Control struct, including the compatibility setter/getter --- */

aCtl = quaidsControlCreate();
aCtlAlias = getDefaultQuaidsControl();
call assert_true(aCtlAlias.maxiter == aCtl.maxiter, "getDefaultQuaidsControl does not match quaidsControlCreate");
call assert_true(quaidsGetHomogeneity(aCtl) == aCtl.homogenous, "quaidsGetHomogeneity does not match aCtl.homogenous");
aCtl = quaidsSetHomogeneity(aCtl, 1);
call assert_true(aCtl.homogenous == 1, "quaidsSetHomogeneity(aCtl, 1) did not take effect");

aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.err = .001;


/* --- Preflight, estimation, the legacy wrapper --- */

pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(pOut.n == N and pOut.nobs == tobs, "quaidsPreflight metadata invalid");

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(qOut.model $== "QUAIDS" and qOut.converged == 1, "quaidsFit did not converge");
call assert_true(qOut.homogeneous == qOut.homogenous, "qOut.homogeneous does not match the deprecated alias");
call printQuaids(qOut);

{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(maxc(maxc(abs(b2 - qOut.bS))) == 0, "legacy quaids() wrapper does not match quaidsFit()");

shareVars = "W1"$|"W2"$|"W3"$|"W4"$|"W5";
priceVars = "P1"$|"P2"$|"P3"$|"P4"$|"P5";
data = asDF(w[.,1], "W1");
data = dfaddcol(data, "W2", w[.,2]);
data = dfaddcol(data, "W3", w[.,3]);
data = dfaddcol(data, "W4", w[.,4]);
data = dfaddcol(data, "W5", w[.,5]);
data = dfaddcol(data, "P1", prices[.,1]);
data = dfaddcol(data, "P2", prices[.,2]);
data = dfaddcol(data, "P3", prices[.,3]);
data = dfaddcol(data, "P4", prices[.,4]);
data = dfaddcol(data, "P5", prices[.,5]);
data = dfaddcol(data, "TOTEXP", totexp);
data = dfaddcol(data, "Z1", instr);
data = dfaddcol(data, "X1", intcpt);

qOutFull = quaidsFull(data, shareVars, priceVars, "TOTEXP", "Z1", "X1", aCtl);
call assert_true(maxc(maxc(abs(qOutFull.bestB - qOut.bestB))) == 0,
    "quaidsFull does not match quaidsFit under library quaids; alone");


/* --- Elasticities, shares, Slutzky, welfare --- */

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptPt = m_[1:1+nint];
pricesPt = m_[1+nint+1:1+nint+n];
totexpPt = m_[1+nint+n+1];

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call assert_true(rows(elasOut.er) == N, "quaidsElasFit shape invalid");
call printQuaidsElas(elasOut);

{ erDeprecated, epDeprecated, epcDeprecated } = quaidsElas_(qOut.bestB, intcptPt, pricesPt, totexpPt, aCtl);
call assert_true(maxc(maxc(abs(erDeprecated - elasOut.er))) == 0, "quaidsElas_ does not match quaidsElasFit");

sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call assert_true(abs(sumc(sharesOut.w) - 1) < 1e-8, "quaidsSharesFit adding-up failed");
call printQuaidsShares(sharesOut);

call quaidsSlutzky(qOut.bestB, qOut.intcptFull, prices, totexp, aCtl);

wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, pricesPt, totexpPt, aCtl);
call assert_true(wOut.cv == 0 and wOut.ev == 0, "quaidsWelfareFit zero-price-change identity failed");
call printQuaidsWelfare(wOut);


/* --- Hypothesis tests (need an unconstrained fit) --- */

aCtlUnc = quaidsControlCreate();
aCtlUnc.linear = 0;
aCtlUnc.maxiter = 100;
aCtlUnc = quaidsSetHomogeneity(aCtlUnc, 0);
aCtlUnc.err = .001;
qOutUnc = quaidsFit(w, intcpt, prices, totexp, instr, aCtlUnc);

{ statH, pvalH, dfH } = quaidsHomogeneityTest(qOutUnc);
call assert_true(dfH == N-1, "quaidsHomogeneityTest df invalid");
{ statJ, pvalJ, dfJ } = quaidsJointTest(qOutUnc);
call assert_true(dfJ > dfH, "quaidsJointTest df invalid");
{ statQ, pvalQ, dfQ } = quaidsQuadraticTest(qOutUnc);
call assert_true(dfQ == N-1, "quaidsQuadraticTest df invalid");


/* --- Zero-share correction, robust SE, replicate weights, workflow --- */

zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(zOut.converged == 1, "quaidsZeroFit did not converge");
call printQuaidsZero(zOut);

rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl);
call assert_true(rows(rOut.se) > 0, "quaidsRobustFit shape invalid");
call printQuaidsRobust(rOut);

nRep = 5;
repWeights = ones(tobs, nRep);
i = 1;
do while i <= nRep;
    repWeights[((i-1)*(tobs/nRep)+1):(i*(tobs/nRep)), i] = 0;
    i = i + 1;
endo;
repOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, repWeights, (nRep-1)/nRep);
call assert_true(repOut.nCompleted >= 1, "quaidsReplicateWeightFit: no replicates completed");
call printQuaidsReplicateWeight(repOut);

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);
call assert_true(wfOut.converged == 1, "quaidsWorkflowFit did not converge");

pricesPt1 = wfOut.evalPrices;
pricesPt1[1] = pricesPt1[1] + ln(1.05);
wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);
call assert_true(wfScenario.welfareValid == 1, "quaidsWorkflowScenarioFit welfare block invalid");

sampwt = ones(tobs, 1);
wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, sampwt);
call assert_true(wfSurvey.converged == 1, "quaidsSurveyWorkflowFit did not converge");

print "";
print "package_public_api_core_only.e: ALL CHECKS PASSED (optmt never loaded)";
