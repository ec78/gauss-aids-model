/* Public TVP diagnostic D3: pure contract checks plus Stage-6 parity. */
new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src;
#include ../src/quaidsiv.src;
#include ../src/quaidselas.src;
#include ../src/quaidsslutzky.src;
#include ../src/quaids.src;
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;
#include ../src/quaidstvpsmooth.src;
#include ../src/quaidstvpelas.src;
#include ../src/quaidstvpfit.src;
#include quaidsfixtures.src;

ncheck = 0;
nfail = 0;
proc (0) = check(cond, label);
    ncheck = ncheck + 1;
    if cond;
        print "PASS  " $+ label;
    else;
        nfail = nfail + 1;
        print "FAIL  " $+ label;
    endif;
endp;

/* Cheap tier: strict threshold, missing semantics and precedence. */
struct quaidsTVPOut tv;
struct quaidsTVPDiagOut d;
tv.mleRetcode = 0;
tv.q0 = ones(3, 1);
tv.Qfit = 1|1e-4|1e-3;
d = quaidsTVPDiagnostics(tv);
call check(d.status == 5 and d.statusLabel $== "SCALE-CAUTION" and
    d.scaleCaution == 1 and d.nScaleCollapse == 1,
    "cheap: strict threshold counts collapse but not equality");
call check(d.minQOverQ0 == 1e-4 and d.maxQOverQ0 == 1,
    "cheap: ratio extrema are reported");

tv.mleRetcode = 2;
d = quaidsTVPDiagnostics(tv);
call check(d.status == 1 and d.statusLabel $== "FIT-FAILED",
    "cheap: fit failure has precedence over scale caution");
tv.mleRetcode = 0;
tv.Qfit = 1|miss(0, 0)|1;
d = quaidsTVPDiagnostics(tv);
call check(d.status == 2 and d.qValid == 0 and
    scalmiss(d.scaleCaution) and scalmiss(d.nScaleCollapse) and
    scalmiss(d.minQOverQ0) and scalmiss(d.maxQOverQ0),
    "cheap: invalid Q suppresses all ratio summaries");
tv.Qfit = 1|2|3;
d = quaidsTVPDiagnostics(tv);
call check(d.status == 0 and d.statusLabel $== "OK" and
    d.scaleCaution == 0 and d.nScaleCollapse == 0,
    "cheap: valid noncollapsed Q is OK");

/* Finiteness and Q/LL predicates never summarize invalid values. */
call check(_quaidsTVPDiagQValid({}, 2) == 0 and
    _quaidsTVPDiagQValid(ones(3, 1), 2) == 0 and
    _quaidsTVPDiagQValid(ones(2, 2), 2) == 0 and
    _quaidsTVPDiagQValid(1|0, 2) == 0 and
    _quaidsTVPDiagQValid(1|-1, 2) == 0,
    "Q predicate rejects empty/wrong shapes, zero and negative values");
call check(_quaidsTVPDiagQValid(1|miss(0, 0), 2) == 0 and
    _quaidsTVPDiagQValid(1|(__INFp-__INFp), 2) == 0 and
    _quaidsTVPDiagQValid(1|__INFp, 2) == 0 and
    _quaidsTVPDiagQValid(1|-__INFp, 2) == 0,
    "Q predicate rejects missing, NaN and infinity");
call check(_quaidsTVPDiagLLValid(miss(0, 0)) == 0 and
    _quaidsTVPDiagLLValid(__INFp-__INFp) == 0 and
    _quaidsTVPDiagLLValid(__INFp) == 0 and
    _quaidsTVPDiagLLValid(1|2) == 0 and
    _quaidsTVPDiagLLValid(-12.5) == 1,
    "likelihood predicate accepts only a finite scalar");
{ qNormValid, qNorm } = _quaidsTVPDiagNormalizeQ(1|2|3, 2);
call check(qNormValid == 0 and rows(qNorm) == 2 and cols(qNorm) == 1 and
    ismiss(qNorm),
    "Q normalization makes wrong-shaped output safely conformable and missing");
{ llNormValid, llNorm } = _quaidsTVPDiagNormalizeLL(__INFp-__INFp);
call check(llNormValid == 0 and scalmiss(llNorm),
    "likelihood normalization replaces invalid output with missing");
{ armCaution1, armCollapse1 } =
    _quaidsTVPDiagScaleCaution(1e-4|5e-5, 0.1|0.1, 1e-3);
{ armCaution2, armCollapse2 } =
    _quaidsTVPDiagScaleCaution(1e-3|5e-4, 1|1, 1e-3);
{ armCaution3, armCollapse3 } =
    _quaidsTVPDiagScaleCaution(1e-2|5e-3, 10|10, 1e-3);
call check(armCaution1 == 1 and armCollapse1 == 1 and
    armCaution2 == 1 and armCollapse2 == 1 and
    armCaution3 == 1 and armCollapse3 == 1,
    "multi: all three arms use a strict threshold relative to their own scaled start");
call check(_quaidsTVPDiagScaleRatioValid(1e-3) == 1 and
    _quaidsTVPDiagScaleRatioValid(0) == 0 and
    _quaidsTVPDiagScaleRatioValid(1) == 0 and
    _quaidsTVPDiagScaleRatioValid(-1) == 0 and
    _quaidsTVPDiagScaleRatioValid(miss(0, 0)) == 0 and
    _quaidsTVPDiagScaleRatioValid(__INFp-__INFp) == 0 and
    _quaidsTVPDiagScaleRatioValid(__INFp) == 0 and
    _quaidsTVPDiagScaleRatioValid(0.1|0.2) == 0 and
    _quaidsTVPDiagScaleRatioValid("bad") == 0,
    "scale ratio accepts only a finite numeric scalar strictly inside (0,1)");

/* Pure multi-start precedence and comparison contract. */
tol = 1e-6*200;
call check(_quaidsTVPDiagMultiStatus(0|2|0, 0|1|1, 0|1|1,
    miss(0, 0), miss(0, 0), tol) == 1,
    "multi: FIT-FAILED has highest precedence without reading missing summaries");
call check(_quaidsTVPDiagMultiStatus(0|0|0, 1|0|1, 1|0|1,
    miss(0, 0), miss(0, 0), tol) == 2,
    "multi: INVALID-Q precedes invalid LL without reading missing summaries");
call check(_quaidsTVPDiagMultiStatus(0|0|0, 1|1|1, 1|0|1,
    miss(0, 0), miss(0, 0), tol) == 3,
    "multi: INVALID-LOGLIKELIHOOD precedes warnings without reading missing summaries");
call check(_quaidsTVPDiagMultiStatus(0|0|0, 1|1|1, 1|1|1, 1,
    tol+1, tol) == 4, "multi: START-SENSITIVE precedes scale caution");
call check(_quaidsTVPDiagMultiStatus(0|0|0, 1|1|1, 1|1|1, 1,
    tol, tol) == 5, "multi: stable likelihood plus collapse is SCALE-CAUTION");
call check(_quaidsTVPDiagMultiStatus(0|0|0, 1|1|1, 1|1|1, 0,
    tol, tol) == 0, "multi: clean stable result is tier success");
{ stableCode5, bestIdxCode5, bestLLCode5, maxGapCode5, tolCode5 } =
    _quaidsTVPDiagCompareLL(0|0|0, 10|10|10, 200);
statusCode5 = _quaidsTVPDiagMultiStatus(0|0|0, 1|1|1, 1|1|1, 1,
    maxGapCode5, tolCode5);
call check(statusCode5 == 5 and stableCode5 == 1,
    "multi: SCALE-CAUTION code 5 coexists with startStable=1");

{ stable, bestIdx, bestLL, maxGap, tolOut } =
    _quaidsTVPDiagCompareLL(0|0|0, tol|0|tol/2, 200);
call check(stable == 1 and bestIdx == 1 and bestLL == tol and
    maxGap == tolOut and tolOut == tol,
    "multi: tolerance is inclusive and exact ties select first maximum");
{ stable, bestIdx, bestLL, maxGap, tolOut } =
    _quaidsTVPDiagCompareLL(0|0|0, tol+1e-9|0|tol/2, 200);
call check(stable == 0, "multi: gap above tolerance is start-sensitive");
{ stable, bestIdx, bestLL, maxGap, tolOut } =
    _quaidsTVPDiagCompareLL(0|2|0, 10|10|10, 200);
call check(stable == 0 and bestIdx == 1,
    "multi: nonzero retcode prevents stability and tie remains first-arm");

/* Derived-start pure guard: extreme scaling must be detected pre-fit. */
call check(_quaidsTVPDiagQValid((1e308*10)*ones(2, 1), 2) == 0 and
    _quaidsTVPDiagQValid((1e-323*0.1)*ones(2, 1), 2) == 0,
    "derived starts detect overflow and underflow before optimization");

/* Stage-6 fixture: the fixed-grid 1x arm is exactly quaidsTVPFit(). */
n = 3;
n1 = 2;
tobs = 120;
k_states = 7;
state0 = .40|.30|.05|-.02|.04|-.04|.03;
trueQ = 1e-5*ones(n1, 1)|1e-4*ones(k_states-n1, 1);
H = 2.5e-5*eye(n1);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);
{ w, prices, totexp, trueStateArr, lxTrue } =
    _quaidsTVPRawSyntheticDGP(tobs, 7, state0, trueQ, H);
tvpCtl = quaidsTVPControlCreate();
tvpCtl.smooth = 0;
tvpCtl.maxIters = 500;
tvpCtl.maxTime = 2;

struct quaidsTVPOut fit;
struct quaidsTVPMultiStartOut md;
screen off;
fit = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
md = quaidsTVPMultiStartDiagnostics(w, prices, totexp, H, q0, tvpCtl);
screen on;
pricesRel = prices[., 1:n1] - prices[., n];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
ll1 = _quaidsTVPDiagTotalLL(Zarr, fit.Qfit, H, w[., 1:n1], n1);
call check(maxc(abs(fit.Qfit - md.qfits[., 2])) == 0 and
    fit.mleRetcode == md.retcodes[2] and
    abs(ll1 - md.loglikelihoods[2]) == 0,
    "integration: 1x arm exactly matches quaidsTVPFit Q, retcode and total LL");
call check(rows(md.multipliers) == 3 and rows(md.qfits) == k_states and
    cols(md.qfits) == 3 and rows(md.scaleCaution) == 3 and
    md.startStable == 1,
    "integration: public multi-start result has frozen shapes and observed stability");
call check(rows(md.status)*cols(md.status) == 1 and
    rows(md.scaleRatio)*cols(md.scaleRatio) == 1 and
    rows(md.anyScaleCaution)*cols(md.anyScaleCaution) == 1 and
    rows(md.bestStartIndex)*cols(md.bestStartIndex) == 1 and
    rows(md.bestMultiplier)*cols(md.bestMultiplier) == 1 and
    rows(md.bestLL)*cols(md.bestLL) == 1 and
    rows(md.maxLLGap)*cols(md.maxLLGap) == 1 and
    rows(md.tolerance)*cols(md.tolerance) == 1 and
    rows(md.startStable)*cols(md.startStable) == 1 and
    rows(md.multipliers) == 3 and cols(md.multipliers) == 1 and
    rows(md.retcodes) == 3 and cols(md.retcodes) == 1 and
    rows(md.fitOk) == 3 and cols(md.fitOk) == 1 and
    rows(md.qValid) == 3 and cols(md.qValid) == 1 and
    rows(md.llValid) == 3 and cols(md.llValid) == 1 and
    rows(md.scaleCaution) == 3 and cols(md.scaleCaution) == 1 and
    rows(md.loglikelihoods) == 3 and cols(md.loglikelihoods) == 1 and
    rows(md.nScaleCollapse) == 3 and cols(md.nScaleCollapse) == 1 and
    rows(md.qfits) == k_states and cols(md.qfits) == 3,
    "integration: every numeric multi-start field has its frozen shape");
call check(md.status == 0 or md.status == 5,
    "integration: stable fixture is success or explicit scale caution");

print "quaidstvpdiag_test.e: " $+ ftocv(ncheck, 1, 0) $+ " checks, "
    $+ ftocv(nfail, 1, 0) $+ " failed";
if nfail;
    end;
endif;
print "quaidstvpdiag_test.e: ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
