/*
** Expected-failure guard test: _quaidsTVPElasFit() must reject a state
** vector whose length does not match k_states implied by n1 -- TVP-AIDS
** initiative, Stage 5.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidselas.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpelas.src;
#include quaidsfixtures.src;

n1 = 2;
tobs = 20;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 1, n1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;

prices = zeros(n1+1, 1);
badState = trueState|0;  /* one element too many */
{ er, ep, epc } = _quaidsTVPElasFit(badState, n1, prices, 0, aCtl);
