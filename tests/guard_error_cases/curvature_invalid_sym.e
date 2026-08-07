/*
** Expected-failure guard test: quaidsCurvatureFit() must reject a qOut whose
** symmetry-constrained estimate has been marked invalid.
*/

new;
library optmt;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidscurvature.src;
#include quaidsfixtures.src;


tobs = 3000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsCurvatureSyntheticDGP(tobs);

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
qOut.symValid = 0;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
