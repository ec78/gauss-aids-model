/*
** Expected-failure guard test: quaidsTVPFit() must reject a non-positive tvpCtl.maxTime (-1 minutes)
** before any estimation work starts (_quaidsTVPResolveCaps,
** src/quaidstvpfit.src) -- TVP-AIDS promotion P3, finite optimizer caps.
** Every other input is valid (a well-posed DGP from
** _quaidsTVPRawSyntheticDGP), so only the cap can trigger the error.
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;
#include ../src/quaidstvpsmooth.src;
#include ../src/quaidstvpelas.src;
#include ../src/quaidstvpfit.src;
#include quaidsfixtures.src;

n1 = 2;
k_states = 7;
state0 = .40|.30|.05|-.02|.04|-.04|.03;
trueQ = 1e-5*ones(n1, 1)|1e-4*ones(k_states-n1, 1);
H = 2.5e-5*eye(n1);
q0 = 0.1*meanc(diag(H))*ones(k_states, 1);
{ w, prices, totexp, trueStateArr, lx } = _quaidsTVPRawSyntheticDGP(40, 7, state0, trueQ, H);

tvpCtl = quaidsTVPControlCreate();
tvpCtl.maxTime = -1;
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
