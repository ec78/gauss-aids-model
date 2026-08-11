/*
** Expected-failure guard test: quaidsWorkflowFit() must reject scalar
** nonzero weights instead of silently routing them to the unweighted path.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaidstests.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
#include ../src/quaidsdiagnostics.src
#include ../src/quaidswelfare.src
#include ../src/quaidsworkflow.src
#include quaidsfixtures.src;


tobs = 1000;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, 204, 1, 1);
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, 0, 1);
