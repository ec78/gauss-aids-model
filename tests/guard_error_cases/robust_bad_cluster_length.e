/*
** Expected-failure guard test: quaidsRobustFit() must reject cluster labels
** that do not match the fitted sample length.
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

struct quaidsControl aCtl;
struct quaidsOut qOut;
struct quaidsRobustOut rOut;

tobs = 1000;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, 204, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, seqa(1, 1, tobs-1));
