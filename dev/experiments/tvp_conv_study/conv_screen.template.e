/*
** conv_screen.template.e -- share-admissibility screen. run_conv.ps1
** substitutes __N1__, __TOBS__, __SEEDS__ and __QSCALE__. For each seed, simulates
** the study DGP (no fitting) and reports whether every share of all n
** goods lies strictly inside (0,1).
*/
new;
library cmlmt, tsmt, sslib;
#include conv_includes.inc

n1 = __N1__;
tobs = __TOBS__;
seeds = __SEEDS__;
qscale = __QSCALE__;

state0 = _convState0(n1);
trueQ = _convTrueQ(n1, qscale);
trueH = _convTrueH(n1);
format /rd 12,6;
i = 1;
do while i <= rows(seeds);
    { w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, seeds[i], state0, trueQ, trueH);
    ok = minc(minc(w)) > 0 and maxc(maxc(w)) < 1;
    print "CONV_SCREEN qscale=" qscale " n1=" n1 " tobs=" tobs " seed=" seeds[i] " minw=" minc(minc(w)) " maxw=" maxc(maxc(w)) " admissible=" ok;
    i = i + 1;
endo;
print "CONV_SCREEN_DONE";
