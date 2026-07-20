/*
** quaids_formula_parity_test.e
**
** Milestone 2 formula-vs-matrix parity test: builds the same synthetic
** 5-good dataset both as plain matrices and as a named-column dataframe,
** estimates with quaidsFit(w, intcpt, prices, totexp, instr, aCtl) and with
** quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars,
** aCtl), and asserts the two produce numerically identical output.
**
** Run from the tests/ directory so the relative #includes resolve:
**   cd tests
**   tgauss -b -x quaids_formula_parity_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsformula.src;

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


/* --- Matrix-API estimate --- */

struct quaidsOut qOutMatrix;
qOutMatrix = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);


/* --- Assemble the identical data as a named-column dataframe --- */

shareVars = "W1"$|"W2"$|"W3"$|"W4"$|"W5";
priceVars = "P1"$|"P2"$|"P3"$|"P4"$|"P5";
totexpVar = "TOTEXP";
instrVars = "Z1";
extraVars = "X1";

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

print "dataframe columns:" getcolnames(data)';
print "dataframe rows:" rows(data);
print;


/* --- Dataframe-API estimate --- */

struct quaidsOut qOutFormula;
qOutFormula = quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars, aCtl);


/* --- Parity checks: matrix API vs. dataframe API on identical data --- */

call check(qOutFormula.model $== qOutMatrix.model, "model matches");
call check(qOutFormula.nobs == qOutMatrix.nobs, "nobs matches");
call check(qOutFormula.n == qOutMatrix.n, "n matches");
call check(qOutFormula.nint == qOutMatrix.nint, "nint matches");
call check(qOutFormula.ninst == qOutMatrix.ninst, "ninst matches");
call checkEqual(qOutFormula.u, qOutMatrix.u, "IV first-stage residuals (u) match exactly");
call checkEqual(qOutFormula.ivB, qOutMatrix.ivB, "IV first-stage coefficients (ivB) match exactly");
call checkEqual(qOutFormula.homogB, qOutMatrix.homogB, "homogeneity-constrained coefficients (homogB) match exactly");
call checkEqual(qOutFormula.homogV, qOutMatrix.homogV, "homogeneity-constrained covariance (homogV) match exactly");
call checkEqual(qOutFormula.symcB, qOutMatrix.symcB, "symmetry-constrained coefficients (symcB) match exactly");
call check(qOutFormula.symStat == qOutMatrix.symStat, "symmetry test statistic matches exactly");
call checkEqual(qOutFormula.b, qOutMatrix.b, "final b (legacy b1) matches exactly");
call checkEqual(qOutFormula.v, qOutMatrix.v, "final v (legacy v1) matches exactly");
call checkEqual(qOutFormula.bS, qOutMatrix.bS, "final bS (legacy b2) matches exactly");
call checkEqual(qOutFormula.vS, qOutMatrix.vS, "final vS (legacy v2) matches exactly");


/* --- extraVars == 0 (no extra intercept shifters) path --- */

struct quaidsOut qOutNoExtra;
qOutNoExtra = quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, 0, aCtl);
call check(qOutNoExtra.nint == 0, "quaidsFull with extraVars=0 yields nint == 0");

struct quaidsOut qOutNoExtraMatrix;
qOutNoExtraMatrix = quaidsFit(w, 0, prices, totexp, instr, aCtl);
call checkEqual(qOutNoExtra.b, qOutNoExtraMatrix.b,
    "quaidsFull(extraVars=0) matches quaidsFit(intcpt=0) exactly");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "FORMULA PARITY TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "FORMULA PARITY TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
