/*
** Expected-failure guard test: _quaidsTVPMLEFit() must reject a y that
** is not nobs x n1 -- TVP-AIDS initiative, Stage 3. Note this is this
** library's OWN usual Txn shares convention (NOT the transposed
** n1 x nobs form _quaidsTVPKalmanFit() takes -- see quaidstvpmle.src's
** own header for why ssFitTVP()'s convention differs).
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

sOut = _quaidsTVPMLEFit(Zarr, H, w', q0, n1);          /* wrong shape: transposed */
