/*
** q0_data.inc -- builds one endpoint's dataset exactly as Phase C's
** conv_fit.template.e does. Inputs set by the including job: mode
** (1 = known-truth DGP, 2 = P6 anchor), n1Cell, tobs, seed, harm
** (1 = oracle H, 2 = static H), qscale. Outputs: w, prices, totexp, H,
** trueQ (missing for the anchor), qShipped, qRef (true Q for DGP cells,
** the shipped start for the anchor), kCell.
*/
kCell = 2*n1Cell + n1Cell*(n1Cell+1)/2;
if mode == 1;
    trueQ = _convTrueQ(n1Cell, qscale);
    trueH = _convTrueH(n1Cell);
    { w, prices, totexp, trueStateArr, lxTrue } = _quaidsTVPRawSyntheticDGP(tobs, seed, _convState0(n1Cell), trueQ, trueH);
    if harm == 1;
        H = trueH;
    else;
        H = _convStaticH(w, prices, totexp, n1Cell);
    endif;
else;
    { wFull, intcpt, pricesFull, totexp, instr } = quaidsExampleData(tobs, seed);
    wRaw = wFull[., 1:n1Cell+1];
    w = wRaw ./ sumc(wRaw');
    prices = pricesFull[., 1:n1Cell+1];
    aCtl = quaidsControlCreate();
    aCtl.linear = 1;
    aCtl.homogenous = 1;
    screen off;
    qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
    screen on;
    H = diagrv(eye(n1Cell), qOut.homogSse[1:n1Cell]/qOut.nobs);
    trueQ = miss(0, 0)*ones(kCell, 1);
endif;
qShipped = 0.1*meanc(diag(H))*ones(kCell, 1);
if mode == 1;
    qRef = trueQ;
else;
    qRef = qShipped;
endif;
