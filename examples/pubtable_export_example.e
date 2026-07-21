/*
** pubtable_export_example.e
**
** Milestone 6: demonstrates exporting quaidsFit()/quaidsElasFit() results
** to publication-quality LaTeX/Markdown/CSV tables via the pubtable
** adapter (src/pubtable_quaids.src). Manual, eyeball-comparison example
** (no assertions) -- mirrors quaids_example.e's style, but for the
** reporting layer instead of the estimator itself. See
** tests/quaids_pubtable_test.e for the automated version.
**
** Requires the pubtable package installed (this machine has it at
** c:\gauss26\pkgs\pubtable -- see CLAUDE.md's "What GAUSS Already
** Provides" section). Run from the examples/ directory:
**   tgauss -b -x pubtable_export_example.e
*/

new;
library pubtable;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/pubtable_quaids.src;

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

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

/* Coefficient table: one comparison column per good, symmetry-constrained
   estimates (qOut.bestB/qOut.bestV), with significance stars and SEs. */
struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);

call ptExport(coefTbl, "quaids_coefficients.tex");
call ptExport(coefTbl, "quaids_coefficients.md");
call ptExport(coefTbl, "quaids_coefficients.csv");

/* Elasticities at the sample mean. */
n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);

/* Income elasticities, uncompensated price elasticities, compensated
   price elasticities -- three tables, one export call each per format. */
struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);

call ptExport(elasTbls[1], "quaids_income_elasticities.md");
call ptExport(elasTbls[2], "quaids_uncompensated_elasticities.tex");
call ptExport(elasTbls[3], "quaids_compensated_elasticities.csv");

print "Exported quaids_coefficients.{tex,md,csv},";
print "quaids_income_elasticities.md, quaids_uncompensated_elasticities.tex,";
print "quaids_compensated_elasticities.csv to the current working directory.";
