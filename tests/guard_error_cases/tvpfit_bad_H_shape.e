/*
** Expected-failure guard test: quaidsTVPFit() must reject an H that is
** not n1 x n1 -- reaches the same guard as tvp_mle_bad_H_shape.e
** (_quaidsTVPMLEFit's own shape check), but via the public quaidsTVPFit()
** wrapper -- TVP-AIDS initiative, Stage 6.
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

rndseed 1;

n = 3;
n1 = n - 1;
tobs = 20;
k_states = 2*n1 + n1*(n1+1)/2;

prices = 1 + 0.1*rndn(tobs, n);
totexp = 5 + 0.2*rndn(tobs, 1);
raw = 0.1 + abs(rndn(tobs, n));
w = raw ./ sumc(raw');

H = eye(n1+1);          /* wrong size: should be n1 x n1 */
q0 = ones(k_states, 1);

tvpCtl = quaidsTVPControlCreate();
tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
