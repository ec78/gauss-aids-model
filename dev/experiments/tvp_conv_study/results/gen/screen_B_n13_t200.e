/*
** conv_screen.template.e -- share-admissibility screen. run_conv.ps1
** substitutes 3, 200, {1001,1002,1003,1004,1005,1006,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110} and 0.01. For each seed, simulates
** the study DGP (no fitting) and reports whether every share of all n
** goods lies strictly inside (0,1).
*/
new;
library cmlmt, tsmt, sslib;
#include conv_includes.inc

n1 = 3;
tobs = 200;
seeds = {1001,1002,1003,1004,1005,1006,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110};
qscale = 0.01;

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
