/*
** Expected-failure guard test: quaidsTVPElasFit() must reject
** aCtl.linear == 0 -- TVP-AIDS initiative, Stage 5. This initiative's own
** Stage 1 scope (Stone index, no quadratic term) means the recovered
** state has no lambda row for _quaidsElas() to read if aCtl.linear is
** left at its non-TVP default (0).
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
/* aCtl.linear left at its default 0 -- deliberately wrong for this call. */

prices = zeros(n1+1, 1);
{ er, ep, epc } = quaidsTVPElasFit(trueState, n1, prices, 0, aCtl);
