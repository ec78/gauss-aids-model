/*
** Expected-failure guard test: quaidsReplicateWeightFit() must reject a
** non-positive scaleFactor.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsreplicate.src
#include quaidsfixtures.src;

struct quaidsControl aCtl;
struct quaidsReplicateOut rOut;

tobs = 1000;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, 204, 1, 1);
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl, ones(tobs, 5), -0.5, method="JK1");
