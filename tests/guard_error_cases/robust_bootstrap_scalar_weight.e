/*
** Expected-failure guard test: quaidsRobustBootstrapFit() must reject
** scalar nonzero weights before any resampling loop starts.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
#include quaidsfixtures.src;


tobs = 1000;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, 204, 1, 1);
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 2, weight=1);
