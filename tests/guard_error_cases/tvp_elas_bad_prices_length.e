/*
** Expected-failure guard test: _quaidsTVPElasFit() must reject a prices
** vector whose length is not n = n1+1 (e.g. a caller mistakenly passing
** the n1-length RELATIVE prices the state space itself was built from,
** instead of the full n-good ABSOLUTE prices _quaidsElas() needs) --
** TVP-AIDS initiative, Stage 5.
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

badPrices = zeros(n1, 1);  /* should be n1+1 */
{ er, ep, epc } = _quaidsTVPElasFit(trueState, n1, badPrices, 0, aCtl);
