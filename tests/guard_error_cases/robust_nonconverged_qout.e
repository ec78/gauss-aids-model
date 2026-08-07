/*
** Expected-failure guard test: quaidsRobustFit() must reject a non-converged
** quaidsFit() result with a clear diagnostic.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
#include quaidsfixtures.src;


aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 43, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);
