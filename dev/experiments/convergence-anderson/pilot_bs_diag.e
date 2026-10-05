new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;
#include direct_objective.src;
#include pilot_lib.src;

/*
** pilot_bs_diag.e -- post-pilot diagnostic (does not change any pilot
** result). For selected cells where rawStructErr <= 1 but bsErr > 1, it
** splits the bS error by block and compares the symmetry-PROJECTED bS
** (the pilot metric) with an UNPROJECTED recovery of the same raw
** endpoint: homogeneity gives the n-th price row (-sum of the n1 gamma
** rows), adding-up gives the n-th column (-sum of columns, +1 on the
** constant). The unprojected formula is checked against _pilotRecoverBS
** on a cell where quaidsFit's own symmetry fallback fired (symValid=0),
** where the two must coincide.
** Run from this directory: C:\gauss26\tgauss.exe -b -x pilot_bs_diag.e
*/
proc (1) = _diagUnproj(bRaw, nIntFull, n1, linear);
    local n, ng, gam, out, rows1;
    n = n1+1;
    ng = rows(bRaw);
    gam = bRaw[nIntFull+1:nIntFull+n1, 1:n1];
    out = bRaw[1:nIntFull, 1:n1] | gam | (-sumc(gam)') | bRaw[nIntFull+n1+1:ng, 1:n1];
    out = out ~ (-sumc(out'));
    out[1, n] = 1 + out[1, n];
    retp(out);
endp;

proc (1) = _diagBlocks(bS, bT, nIntFull, n, linear);
    local e, a, g, be, la, u, r;
    e = abs(bS - bT);
    r = nIntFull;
    a = maxc(maxc(e[1:r, .]));
    g = maxc(maxc(e[r+1:r+n, .]));
    be = maxc(e[r+n+1, .]');
    if linear;
        la = -1;
        u = maxc(e[r+n+2, .]');
    else;
        la = maxc(e[r+n+2, .]');
        u = maxc(e[r+n+3, .]');
    endif;
    retp(a|g|be|la|u|maxc(maxc(e)));
endp;

/* model, seed, method, start */
cases = { 1 21 2 1, 0 22 3 12, 0 25 3 9, 0 29 3 9, 1 15 3 2, 1 16 3 13, 1 20 3 2, 1 22 3 2, 1 25 3 10, 1 29 3 13, 0 26 3 6, 0 16 1 11 };
struct quaidsControl dCtl;
struct quaidsControl dCtlStone;
struct quaidsOut dOutStone;
cachedModel = -1;
print "cols: model seed method start symValid | projected bS block err [alpha gamma beta lambda u max] | unprojected [same] | raw struct";
ci = 1;
do while ci <= rows(cases);
    cm = cases[ci,1]; cs = cases[ci,2]; cme = cases[ci,3]; cst = cases[ci,4];
    if cm /= cachedModel;
        cells = csvReadM("pilot_out/cells_m" $+ ftos(cm, "%*.*lf", 1, 0) $+ ".csv", 1);
        cachedModel = cm;
    endif;
    idx = selif(seqa(1,1,rows(cells)), cells[.,2].==cs .and cells[.,3].==cme .and cells[.,4].==cst);
    row = cells[idx, .]';
    ng = row[57]; n1 = row[58]; n = n1+1;
    bRaw = reshape(row[59:58+ng*n1], n1, ng)' ~ zeros(ng, 1);
    bSstored = reshape(row[59+ng*n1:58+ng*n1+(ng+1)*n], n, ng+1)';
    dCtl = quaidsControlCreate;
    dCtl.linear = 1-cm; dCtl.maxiter = 100; dCtl.err = .0001; dCtl.homogenous = 1;
    { w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, cs, cm, 1);
    nIntFull = 2;
    if cme == 1;
        /* Anderson cell: stored bS is qOut.bS; recompute (b,b) for the
           comparison of the projection only */
        { bSP, okP, symP, TbP, crP } = _pilotRecoverBS(w, intcpt, prices, totexp, instr, dCtl, bRaw, bRaw);
    else;
        { bSP, okP, symP, TbP, crP } = _pilotRecoverBS(w, intcpt, prices, totexp, instr, dCtl, bRaw, bRaw);
        if maxc(maxc(abs(bSP-bSstored))) > 1e-9;
            print "  WARNING: recomputed bS differs from stored by" maxc(maxc(abs(bSP-bSstored)));
        endif;
    endif;
    bSU = _diagUnproj(bRaw, nIntFull, n1, dCtl.linear);
    pB = _diagBlocks(bSP, trueParams, nIntFull, n, dCtl.linear);
    uB = _diagBlocks(bSU, trueParams, nIntFull, n, dCtl.linear);
    print cm~cs~cme~cst~symP~pB'~uB'~row[11];
    if symP == 0;
        print "  formula check on symmetry-fallback cell: max|unprojected - _pilotRecoverBS|" maxc(maxc(abs(bSU-bSP)));
    endif;
    ci = ci+1;
endo;
print "pilot_bs_diag.e: run complete.";
