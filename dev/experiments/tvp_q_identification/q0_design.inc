/*
** q0_design.inc -- the first lines of Phase C's conv_chain.inc (no fit):
** the measurement design Zarr for the dataset q0_data.inc built.
*/
n = cols(prices);
n1 = n - 1;
nobs = rows(w);
k_states = 2*n1 + n1*(n1+1)/2;
pricesRel = prices[., 1:n1] - prices[., n];
lx = _quaidsTVPStoneIndex(w, prices, totexp, n);
Zarr = _quaidsTVPBuildZ(pricesRel, lx, n1);
y = w[., 1:n1];
