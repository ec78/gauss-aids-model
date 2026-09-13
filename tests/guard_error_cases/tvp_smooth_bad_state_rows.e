/*
** Expected-failure guard test: _quaidsTVPSmoothFit() must reject a
** kalmanResult whose filtered_state row count (k_states) does not match
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
tobs = 20;
{ w, pricesRel, lx, trueState, trueGammaFull } = _quaidsTVPStaticSyntheticDGP(tobs, 1, n1);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

tvpm = _quaidsTVPBuildModel(Zarr, 0.01*eye(k_states), 0.01*eye(n1), n1);

/* A second model with a DIFFERENT n1 (hence different k_states), so its
   rslt.filtered_state has the wrong number of rows for tvpm above. */
n1b = 3;
{ wb, pricesRelb, lxb, trueStateb, trueGammaFullb } = _quaidsTVPStaticSyntheticDGP(tobs, 2, n1b);
Zarrb = _quaidsTVPBuildZ(pricesRelb, lxb, n1b);
ngammab = n1b*(n1b+1)/2;
k_statesb = n1b + ngammab + n1b;
tvpmb = _quaidsTVPBuildModel(Zarrb, 0.01*eye(k_statesb), 0.01*eye(n1b), n1b);
rsltb = _quaidsTVPKalmanFit(tvpmb, wb', 1);

{ aTS, pTS } = _quaidsTVPSmoothFit(tvpm, rsltb);
