/*
** Expected-failure guard test: quaidsFit() must reject scalar nonzero
** weights. Scalar 0 is the only unweighted sentinel; real sampling weights
** must be Tx1 vectors.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include quaidsfixtures.src;


tobs = 1000;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(tobs, 204, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl, 1);
