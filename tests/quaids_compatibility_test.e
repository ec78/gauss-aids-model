/*
** quaids_compatibility_test.e
**
** Public release roadmap PR-003: tests for the compatibility surface
** introduced when preparing this library for public release --
** quaidsSetHomogeneity()/quaidsGetHomogeneity() (the correctly spelled
** control accessors), the homogeneous/homogenous dual fields on
** quaidsOut/quaidsZeroOut/quaidsWorkflowOut, and the quaidsElas_()
** deprecated compatibility wrapper. See docs/public-api.json and
** README.md's "Compatibility Policy" section.
**
** quaidsSetHomogeneity()'s own invalid-input guard is tested separately
** in guard_error_cases/quaids_set_homogeneity_invalid.e, since a real
** guard error aborts the whole batch job (GAUSS's `errorlog`+`end` has
** no "confirm it errors, then keep going" idiom -- see CLAUDE.md's
** Milestone 17 finding).
**
** Run from the tests/ directory so the relative #includes resolve:
**   cd tests
**   tgauss -b -x quaids_compatibility_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidszerocorrect.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstests.src
#include ../src/quaidswelfare.src
#include ../src/quaidsrobust.src
#include ../src/quaidsdiagnostics.src
#include ../src/quaidsworkflow.src

nfail = 0;
ncheck = 0;

proc (0) = check(cond, label);
    local i;
    i = ncheck + 1;
    ncheck = i;
    if cond;
        print "PASS  " $+ label;
    else;
        print "FAIL  " $+ label;
        nfail = nfail + 1;
    endif;
endp;

proc (0) = checkEqual(a, b, label);
    call check(maxc(maxc(abs(a - b))) == 0, label);
endp;


/* Deterministic synthetic 5-good dataset, same DGP as examples/01_basic_estimation.e */

seed = 204;
tobs = 1000;
N = 5;

al = round(rndns(1,N-1,seed)*10)/10;
al = al~(1-sumc(al'));
al1 = .5*round(rndns(1,N-1,seed)*10)/10;
al1 = al1~(-sumc(al1'));
ga = round(rndns(N-1,N-1,seed)*10)/10;
ga = xpnd(vech(ga));
ga = ga|(-sumc(ga)');
ga = ga~(-sumc(ga'));
be = .5*round(rndns(1,N-1,seed)*10)/10;
be = be~(-sumc(be'));
la = .01*round(rndns(1,N-1,seed)*10)/10;
la = la~(-sumc(la'));
ro = round(rndns(1,N-1,seed)*10)/10;

prices = 1+rndns(tobs,N,seed);
instr = 5+5*rndns(tobs,1,seed);
intcpt = 2+2*rndns(tobs,1,seed);
u =  .1*rndns(tobs,1,seed);
totexp = .85*instr + u;
e = 2*rndns(tobs,N-1,seed) + u*ro;
e = e~(-sumc(e'));

a_p = sumc( (prices.*(al+intcpt*al1) )') + .5*sumc(((prices*ga).*prices)');
lx = totexp -a_p;
b_p = prices*be';
lx2 = (lx^2)./exp(b_p);

w = al +  prices*ga + lx*be + e +intcpt*al1 + lx2*la ;


/* --- quaidsSetHomogeneity()/quaidsGetHomogeneity() -- valid inputs --- */

aCtl = quaidsControlCreate();
call check(quaidsGetHomogeneity(aCtl) == 1, "quaidsGetHomogeneity reads quaidsControlCreate()'s default (1)");

aCtl = quaidsSetHomogeneity(aCtl, 0);
call check(aCtl.homogenous == 0, "quaidsSetHomogeneity(aCtl, 0) sets the underlying field");
call check(quaidsGetHomogeneity(aCtl) == 0, "quaidsGetHomogeneity reads back 0 after quaidsSetHomogeneity(aCtl, 0)");

aCtl = quaidsSetHomogeneity(aCtl, 1);
call check(aCtl.homogenous == 1, "quaidsSetHomogeneity(aCtl, 1) sets the underlying field");
call check(quaidsGetHomogeneity(aCtl) == 1, "quaidsGetHomogeneity reads back 1 after quaidsSetHomogeneity(aCtl, 1)");


/* --- homogeneous/homogenous dual fields on returned structs --- */

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl = quaidsSetHomogeneity(aCtl, 1);
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "quaidsFit() prerequisite fit converged");
call check(qOut.homogeneous == qOut.homogenous, "qOut.homogeneous matches the deprecated qOut.homogenous alias exactly");
call check(qOut.homogeneous == 1, "qOut.homogeneous reflects the aCtl setting used to fit");

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);
call check(wfOut.homogeneous == wfOut.homogenous, "wfOut.homogeneous matches the deprecated wfOut.homogenous alias exactly");
call check(wfOut.homogeneous == qOut.homogeneous, "wfOut.homogeneous matches the underlying qOut.homogeneous");

zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call check(zOut.homogeneous == zOut.homogenous, "zOut.homogeneous matches the deprecated zOut.homogenous alias exactly");
call check(zOut.homogeneous == 1, "zOut.homogeneous reflects the aCtl setting used to fit");


/* --- quaidsElas_() deprecated compatibility wrapper --- */

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptPt = m_[1:1+nint];
pricesPt = m_[1+nint+1:1+nint+n];
totexpPt = m_[1+nint+n+1];

{ erDeprecated, epDeprecated, epcDeprecated } = quaidsElas_(qOut.bestB, intcptPt, pricesPt, totexpPt, aCtl);

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

call checkEqual(erDeprecated, elasOut.er, "quaidsElas_() income elasticities match quaidsElasFit()'s er exactly");
call checkEqual(epDeprecated, elasOut.ep, "quaidsElas_() uncompensated price elasticities match quaidsElasFit()'s ep exactly");
call checkEqual(epcDeprecated, elasOut.epc, "quaidsElas_() compensated price elasticities match quaidsElasFit()'s epc exactly");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "COMPATIBILITY TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "COMPATIBILITY TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
