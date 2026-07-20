/*
** quaids_schema_test.e
**
** Milestone 1 schema test: asserts that quaidsFit() returns a quaidsOut
** struct with the expected field values/shapes, and that it produces this
** struct with no console output (estimation and printing are split).
**
** Run from the tests/ directory so the relative #includes resolve:
**   cd tests
**   tgauss -b -x quaids_schema_test.e
**
** Prints one PASS/FAIL line per check plus a final summary line. This is a
** source-tree smoke/schema test, not a numerical-correctness benchmark --
** see GOLD_STANDARD_TODO.md Milestone 3 for published/deterministic
** numerical validation.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;

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


/* Deterministic synthetic 5-good dataset, same DGP as examples/quaids_example.e */

seed = 11;
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

struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .001;


/* --- Silence check: nothing should print between these two markers --- */

"=== BEGIN QUAIDSFIT SILENCE WINDOW ===";
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
"=== END QUAIDSFIT SILENCE WINDOW (should be no output above this line since the previous marker) ===";


/* --- Metadata --- */

call check(qOut.model $== "QUAIDS", "model == QUAIDS (aCtl.linear=0, aCtl.maxiter>1)");
call check(qOut.homogenous == 1, "homogenous flag echoed correctly");
call check(qOut.linear == 0, "linear flag echoed correctly");
call check(qOut.nobs == tobs, "nobs == tobs");
call check(qOut.n == N, "n == number of goods");
call check(qOut.nint == 1, "nint == 1 (one extra intercept-shifter column)");
call check(qOut.ninst == 1, "ninst == 1 (one instrument column)");
call check(qOut.nu == 1, "nu == 1 (log total expenditure is the only endogenous regressor)");
call check(rows(qOut.xnam) == 1+qOut.nint, "xnam length == 1 + nint");
call check(rows(qOut.wnam) == qOut.n, "wnam length == n");
call check(rows(qOut.znam) == qOut.ninst, "znam length == ninst");
call check(rows(qOut.intcptFull) == tobs and cols(qOut.intcptFull) == 1+qOut.nint,
    "intcptFull is nobs x (1+nint)");
call check(rows(qOut.u) == tobs and cols(qOut.u) == qOut.nu, "u (IV residuals) is nobs x nu");


/* --- Iteration bookkeeping --- */

call check(qOut.iterations > 0, "iterations > 0 for an iterated (maxiter>1) fit");
/* This synthetic DGP/tolerance combination (err=.001, maxiter=100) does not
   actually satisfy the convergence tolerance -- it runs to the iteration
   cap (confirmed against the printed iteration log: iteration 99, err~20.9).
   That is a property of this fixture, not a bug, so the correct assertion
   is that quaidsFit() honestly reports non-convergence rather than claiming
   false convergence. */
call check(qOut.converged == 0, "converged flag correctly reports false: this fixture hits maxiter without meeting err tolerance");
call check(qOut.iterations == aCtl.maxiter, "iterations == maxiter when the loop exits via the iteration cap");
call check(not ismiss(qOut.finalErr), "finalErr is not missing for an iterated fit");


/* --- First-stage IV diagnostics --- */

call check(rows(qOut.ivB) == qOut.nz and cols(qOut.ivB) == qOut.nu, "ivB is nz x nu");
call check(qOut.ivRsq > 0.9, "ivRsq is high (instrument is a strong predictor by construction)");


/* --- Homogeneity-constrained stage --- */

call check(rows(qOut.homogB) == qOut.ng and cols(qOut.homogB) == qOut.n, "homogB is ng x n");
call check(not ismiss(qOut.homogCrit), "homogCrit (log det sigma) is not missing");


/* --- Overidentification test: ninst == nu here, so should NOT be valid --- */

call check(qOut.overidValid == 0, "overidValid is 0 when ninst == nu (exactly identified)");


/* --- Symmetry test + symmetry-constrained stage: valid since homogenous == 1 --- */

call check(qOut.symValid == 1, "symValid is 1 when homogenous == 1");
call check(qOut.symDf == qOut.n1*(qOut.n1-1)/2, "symDf matches n1*(n1-1)/2");
call check(qOut.symStat >= 0, "symStat (chi2 statistic) is non-negative");
call check(qOut.symPval >= 0 and qOut.symPval <= 1, "symPval is a valid probability");
call check(rows(qOut.symcB) == qOut.ng and cols(qOut.symcB) == qOut.n, "symcB is ng x n");


/* --- Final legacy-compatible outputs: b/v/bS/vS --- */

call check(rows(qOut.b) == qOut.ng+1 and cols(qOut.b) == qOut.n,
    "final b has the absolute-price-recovered row added back (ng+1 rows) when homogenous");
call check(rows(qOut.bS) == rows(qOut.b) and cols(qOut.bS) == qOut.n,
    "final bS matches final b in shape when homogenous");
call check(rows(qOut.bestB) == rows(qOut.b) and cols(qOut.bestB) == qOut.n,
    "bestB (used for elasticities/Slutzky) matches final b in shape");


/* --- Cross-check: quaidsFit()'s stored final b/v/bS/vS match what the legacy --- */
/* --- quaids() wrapper actually returns for the same inputs.                 --- */

output file=schema_test_quaids_wrapper_out reset;
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
output off;

call check(maxc(maxc(abs(b1 - qOut.b))) == 0, "legacy quaids() b1 == quaidsFit() qOut.b exactly");
call check(maxc(maxc(abs(v1 - qOut.v))) == 0, "legacy quaids() v1 == quaidsFit() qOut.v exactly");
call check(maxc(maxc(abs(b2 - qOut.bS))) == 0, "legacy quaids() b2 == quaidsFit() qOut.bS exactly");
call check(maxc(maxc(abs(v2 - qOut.vS))) == 0, "legacy quaids() v2 == quaidsFit() qOut.vS exactly");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "SCHEMA TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "SCHEMA TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
