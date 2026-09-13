/*
** Expected-failure guard test: _quaidsTVPSmoothFit() must reject a
** kalmanResult whose filtered_state column count (nobs) does not match
** the tvpModel it is paired with -- TVP-AIDS initiative, Stage 4.
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpsmooth.src;
#include quaidsfixtures.src;

n1 = 2;
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

tobs = 20;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 1, n1);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
tvpm = _quaidsTVPBuildModel(Zarr, 0.01*eye(k_states), 0.01*eye(n1), n1);

/* A second model with the SAME n1 but a DIFFERENT tobs, so its
   rslt.filtered_state has the wrong number of columns for tvpm above. */
tobsb = 30;
{ wb, pricesRelb, lxb, trueStateb, trueGammaFullb } = _quaidsTVPStaticSyntheticDGP(tobsb, 2, n1);
Zarrb = _quaidsTVPBuildZ(pricesRelb, lxb, n1);
tvpmb = _quaidsTVPBuildModel(Zarrb, 0.01*eye(k_states), 0.01*eye(n1), n1);
rsltb = _quaidsTVPKalmanFit(tvpmb, wb', 1);

{ aTS, pTS } = _quaidsTVPSmoothFit(tvpm, rsltb);
