/*
** Expected-failure guard test: quaidsTrendFit() must reject
** aCtl.homogenous /= 1 -- the exact adding-up/homogeneity guarantee this
** diagnostic relies on (see quaidstrend.src's Remarks) only holds when
** the reference good's price is dropped from the shared design, exactly
** as quaidsFit()'s own homogeneous branch does.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidstrend.src
#include quaidsfixtures.src;

{ w, intcpt, prices, totexp, instr } = _quaidsTrendSyntheticDGP(400, 12345, 0);

aCtl = quaidsControlCreate();
aCtl.homogenous = 0;

tOut = quaidsTrendFit(w, intcpt, prices, totexp, instr, aCtl);
