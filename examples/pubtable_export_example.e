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
**
** src/pubtable_quaids.src is not part of the installed quaids package
** (see CLAUDE.md's "Milestone 6: reporting via pubtable" section), so it
** is still #included from the source tree. The bare "#include quaids.sdf"
** below (rather than a full "library quaids;"-only setup) is required
** so that pubtable_quaids.src's #ifDef QUAIDS_SDF_INCLUDED guard is
** active -- `library quaids;` alone lazily loads procs on demand and does
** not run quaids.sdf's #define, confirmed empirically; explicitly
** including the .sdf (resolved via the installed package's search path)
** is what activates the guard, matching how pubtable's own bundled
** pubtable_qardl.src documents the identical requirement for qardl.sdf.
*/

new;
library pubtable, quaids;
#include quaids.sdf
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
