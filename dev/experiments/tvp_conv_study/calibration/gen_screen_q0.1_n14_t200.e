/*
** conv_screen.template.e -- share-admissibility screen. run_conv.ps1
** substitutes 4, 200, seqa(1101,1,10) and 0.1. For each seed, simulates
** the study DGP (no fitting) and reports whether every share of all n
** goods lies strictly inside (0,1).
*/
new;
library cmlmt, tsmt, sslib;
#include conv_includes.inc

n1 = 4;
tobs = 200;
seeds = seqa(1101,1,10);
qscale = 0.1;

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
