/*
** conv_chain.inc -- the private component chain quaidsTVPFit() runs,
** written out line for line (src/quaidstvpfit.src) so the study can read
** CMLMT's own result fields (iterations, Fct, gradient, time), which
** quaidsTVPOut does not expose. The equivalence job (conv_equiv.e) checks
** that this chain reproduces quaidsTVPFit()'s Qfit and bFinal exactly.
**
** Inputs set by the including job: w, prices, totexp, H, q0, capIters,
** capTime. Outputs: n, n1, nobs, k_states, lx, pricesRel, Zarr, sOut,
** Qfit, rc, aTS, pTS, bFinal. Smoothing is always on (quaidsTVPFit's
** default).
*/
n = cols(prices);
n1 = n - 1;
nobs = rows(w);
k_states = 2*n1 + n1*(n1+1)/2;
pricesRel = prices[., 1:n1] - prices[., n];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
struct ssOut sOut;
sOut = _quaidsTVPMLEFit(Zarr, H, w[., 1:n1], q0, n1, capIters, capTime);
Qfit = sOut.final_params;
rc = sOut.mleResults.retcode;
aTS = 0;
pTS = 0;
struct tvpModel tvpFinalChain;
tvpFinalChain = sOut.tvpFinal;
{ aTS, pTS } = _quaidsTVPSmoothFit(tvpFinalChain, sOut.kfResults);
bFinal = quaidsTVPStateToFullB(aTS[., nobs], n1);
