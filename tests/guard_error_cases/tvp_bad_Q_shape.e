/*
** Expected-failure guard test: _quaidsTVPBuildModel() must reject a Q
** that is not k_states x k_states (k_states implied by Zarr's own page
** width) -- TVP-AIDS initiative, Stage 2.
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include quaidsfixtures.src;

n1 = 3;
tobs = 20;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 1, n1);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);

Q = eye(n1);            /* wrong size: should be k_states x k_states */
H = eye(n1);

tvpm = _quaidsTVPBuildModel(Zarr, Q, H, n1);
