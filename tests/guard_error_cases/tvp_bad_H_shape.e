/*
** Expected-failure guard test: _quaidsTVPBuildModel() must reject an H
** that is not n1 x n1 -- TVP-AIDS initiative, Stage 2.
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

ngamma = n1*(n1+1)/2;
k_states = n1 + ngamma + n1;

Q = eye(k_states);
H = eye(n1+1);          /* wrong size: should be n1 x n1 */

tvpm = _quaidsTVPBuildModel(Zarr, Q, H, n1);
