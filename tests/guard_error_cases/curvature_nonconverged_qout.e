/*
** Expected-failure guard test: quaidsCurvatureFit() must reject a
** non-converged quaidsFit() result before imposing curvature.
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

struct quaidsControl aCtl;
struct quaidsOut qOut;
struct quaidsCurvOut cOut;

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 43, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
