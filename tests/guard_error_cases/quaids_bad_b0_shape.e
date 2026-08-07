/*
** Expected-failure guard test: quaidsFit() must reject a supplied aCtl.b0
** matrix whose shape cannot be used as the reduced raw coefficient block.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include quaidsfixtures.src;


aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 1;
aCtl.homogenous = 1;
aCtl.b0 = zeros(2, 2);

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(1000, 204, 1, 1);
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
