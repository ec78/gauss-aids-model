/*
** workflow_example.e
**
** Minimal applied workflow example using the installed quaids package.
** Run from the repository root or examples/ after installing the package:
**   tgauss -b -x examples/workflow_example.e
*/

new;
library quaids;

seed = 204;
tobs = 1000;
n = 5;

al = round(rndns(1, n-1, seed)*10)/10;
al = al~(1-sumc(al'));
al1 = .5*round(rndns(1, n-1, seed)*10)/10;
al1 = al1~(-sumc(al1'));
ga = round(rndns(n-1, n-1, seed)*10)/10;
ga = xpnd(vech(ga));
ga = ga|(-sumc(ga)');
ga = ga~(-sumc(ga'));
be = .5*round(rndns(1, n-1, seed)*10)/10;
be = be~(-sumc(be'));
ro = round(rndns(1, n-1, seed)*10)/10;

prices = 1+rndns(tobs, n, seed);
instr = 5+5*rndns(tobs, 1, seed);
intcpt = 2+2*rndns(tobs, 1, seed);
u = .1*rndns(tobs, 1, seed);
totexp = .85*instr + u;
e = 2*rndns(tobs, n-1, seed) + u*ro;
e = e~(-sumc(e'));

a_p = sumc((prices.*(al+intcpt*al1))') + .5*sumc(((prices*ga).*prices)');
lx = totexp - a_p;
w = al + prices*ga + lx*be + e + intcpt*al1;

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);

print "Model:" wfOut.model;
print "Converged:" wfOut.converged;
print "Predicted shares at sample mean:";
print wfOut.shares;
print "Income elasticities at sample mean:";
print wfOut.incomeElas;

pricesPt1 = wfOut.evalPrices;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);

if wfScenario.welfareValid;
    print "CV / seCV:" wfScenario.cv wfScenario.seCV;
    print "EV / seEV:" wfScenario.ev wfScenario.seEV;
endif;
