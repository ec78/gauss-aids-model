/*
** export_tvp_loglik_fixture.e
**
** TVP-AIDS promotion P5: writes the committed deterministic fixture
** tvp_loglik_fixture.txt used by the independent fixed-Q log-likelihood
** cross-check (statsmodels_reference.py, quaidstvp_loglik_xcheck_test.e).
**
** The data are the P4 Stage 6 DGP exactly as tests/quaidstvpfit_test.e
** builds it: _quaidsTVPRawSyntheticDGP(), n = 3 goods (n1 = 2,
** k_states = 7), tobs = 120, seed 7, the same state0/trueQ/H. No sslib
** dependency (pure #include of the fixture generator).
**
** Columns (single-space separated, %.17g, so every double round-trips
** exactly): t, w1, w2, w3, prel1, prel2, lx, where
** prel_j = prices[., j] - prices[., 3] (the relative log prices
** _quaidsTVPBuildZ() takes) and lx is the Stone-deflated regressor used to
** simulate (lxTrue).
**
** The committed .txt is the authoritative input. This script only
** documents (and can regenerate) where it came from; the GAUSS check
** reads the .txt and never re-simulates.
**
** Run from this directory (tests/fixtures/tvp_loglik/):
**   tgauss -b -x export_tvp_loglik_fixture.e
*/

new;
#include ../../quaidsfixtures.src;

n = 3;
n1 = n - 1;
tobs = 120;
k_states = 2*n1 + n1*(n1+1)/2;

state0 = .40|.30|.05|-.02|.04|-.04|.03;
trueQ = 1e-5*ones(n1, 1)|1e-4*ones(k_states-n1, 1);
H = 2.5e-5*eye(n1);

{ w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, 7, state0, trueQ, H);

prel = prices[., 1:n1] - prices[., n];
dat = seqa(1, 1, tobs)~w~prel~lxTrue;

/* One line per period, "%.17g" per value (round-trips a double
   exactly), single-space separated, no trailing whitespace. */
fmt = "%.17g %.17g %.17g %.17g %.17g %.17g %.17g";
lines = sprintf(fmt, dat[., 1], dat[., 2], dat[., 3], dat[., 4], dat[., 5], dat[., 6], dat[., 7]);
fh = fopen("tvp_loglik_fixture.txt", "w");
call fputs(fh, "t w1 w2 w3 prel1 prel2 lx" $+ chrs(10));
i = 1;
do while i <= rows(lines);
    call fputs(fh, lines[i] $+ chrs(10));
    i = i + 1;
endo;
call close(fh);
print "wrote tvp_loglik_fixture.txt: " $+ ftocv(rows(dat), 1, 0) $+ " rows";
