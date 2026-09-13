/*
** Expected-failure guard test: _quaidsTVPMLEFit() must reject a q0 with
** a non-positive element -- TVP-AIDS initiative, Stage 3 (every q0
** element's own sqrt is the raw CMLMT starting value, so a non-positive
** entry has no valid unconstrained starting point).
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;
#include quaidsfixtures.src;

n1 = 3;
tobs = 20;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 1, n1);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);

ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

H = eye(n1);
q0 = ones(k_states, 1);
q0[1] = 0;          /* not strictly positive */

sOut = _quaidsTVPMLEFit(Zarr, H, w, q0, n1);
