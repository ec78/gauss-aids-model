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
** pilot_validate.e -- Stage A of the matched pilot (see PILOT_PROTOCOL.md).
**
** Part 1 (this file): checks that are about the metric and the start
** representation themselves, run before any pilot result is looked at.
**   V1  bS reconstruction: _pilotRecoverBS() reproduces quaidsFit()'s
**       (and _quaidsFitAnderson()'s) own homogB/bS/homogCrit EXACTLY (to
**       rounding) at the same (index, coefficient) pair; plus a shared
**       fixed-point endpoint where the direct representation (b,b) and
**       quaidsFit()'s representation coincide.
**   V2  raw-error layout: _pilotErrs() on synthetic matrices with a
**       distinct known error in every row (u row excluded; lambda
**       included for QUAIDS; column n ignored).
**   V3a start representation: pack -> scale -> unscale -> unpack
**       recovers bStart[.,1:n1] for every start used.
** Part 2: the harness itself on seeds 1, 6, 9, both models, all 16
** starts, all 3 methods -> pilot_validation/cells_m0.csv, _m1.csv.
** Part 3 (pilot_analyze.py --validate): benchmark/smoke-gate
** reproduction, objective consistency, selection-direction unit tests,
** key uniqueness. Part 4 (pilot_resume_test.ps1): resume behavior.
**
** Run from this directory:
**   C:\gauss26\tgauss.exe -b -x pilot_validate.e
*/

nFail = 0;
nCheck = 0;

/* ================= V2: raw-error layout ================= */
print "================================================================";
print "V2: _pilotErrs layout on synthetic matrices";
print "================================================================";
lin = 1;
do while lin >= 0;
    nIntFullT = 2;
    n1T = 4;
    nStructT = nIntFullT+n1T+1+(1-lin);
    bTT = zeros(nStructT+1, n1T+1);
    bHT = bTT;
    bHT[1:nStructT, 1:n1T] = seqa(1, 1, nStructT)/100 .* ones(1, n1T);
    bHT[nStructT+1, 1:n1T] = 100*ones(1, n1T);
    bHT[., n1T+1] = 1000*ones(nStructT+1, 1);
    { errS, errL, errU, errO, errP } = _pilotErrs(bHT, bTT, nIntFullT, n1T, lin);
    if lin;
        lamExp = -1;
    else;
        lamExp = nStructT/100;
    endif;
    ok = abs(errS - nStructT/100) < 1e-15 and abs(errL - lamExp) < 1e-15
        and abs(errU - 100) < 1e-15 and abs(errO - 100) < 1e-15;
    nCheck = nCheck+1;
    if ok;
        print "  PASS V2 linear=" lin " structErr" errS "(= last non-u row," nStructT "rows) lambdaErr" errL "uErr" errU;
    else;
        nFail = nFail+1;
        print "  FAIL V2 linear=" lin " structErr" errS "lambdaErr" errL "uErr" errU "overall" errO;
    endif;
    lin = lin-1;
endo;
print;

/* ================= V1 and V3a on seeds 1, 6, 9 ================= */
struct quaidsControl vCtl;
struct quaidsControl vCtlStone;
struct quaidsControl vCtlX;
struct quaidsOut vOutStone;
struct quaidsOut vOutF;
struct quaidsOut vOutA;
struct quaidsOut vOutAnd;
struct quaidsOut vOutPol;

v1Tol = 1e-8;
vSeeds = 1|6|9;
vModel = 0;
do while vModel <= 1;
    vsi = 1;
    do while vsi <= rows(vSeeds);
        vSeed = vSeeds[vsi];
        print "================================================================";
        print "V1/V3a: model" vModel "seed" vSeed;
        print "================================================================";
        vCtl = quaidsControlCreate;
        vCtl.linear = 1-vModel;
        vCtl.maxiter = 100;
        vCtl.err = .0001;
        vCtl.homogenous = 1;
        { vw, vIntcpt, vPrices, vTotexp, vInstr, vTrue } =
            _quaidsSyntheticDGP(3000, vSeed, vModel, 1);
        vCtlStone = vCtl;
        vCtlStone.maxiter = 1;
        vOutStone = quaidsFit(vw, vIntcpt, vPrices, vTotexp, vInstr, vCtlStone);
        vbStone = vOutStone.homogB;
        vNIntFull = cols(vOutStone.intcptFull);
        vN1 = cols(vPrices)-1;
        vNg = vOutStone.ng;
        vScale = _directObjScale(vbStone, vN1);
        vStarts = _pilotAllStarts(vbStone, vNIntFull, vN1, vModel, vSeed, 16);

        /* V3a: start representation round trip */
        maxRt = 0;
        s = 1;
        do while s <= 16;
            bSt = vStarts[(s-1)*vNg+1:s*vNg, .];
            th = _directObjPack(bSt, vN1)./vScale;
            bBack = _directObjUnpack(th.*vScale, vNg, vN1);
            maxRt = maxc(maxRt | maxc(maxc(abs(bBack[., 1:vN1]-bSt[., 1:vN1]))));
            s = s+1;
        endo;
        nCheck = nCheck+1;
        if maxRt <= 1e-14;
            print "  PASS V3a start round trip, 16 starts: max abs diff" maxRt;
        else;
            nFail = nFail+1;
            print "  FAIL V3a start round trip: max abs diff" maxRt;
        endif;
        vStarts2 = _pilotAllStarts(vbStone, vNIntFull, vN1, vModel, vSeed, 16);
        nCheck = nCheck+1;
        if maxc(maxc(abs(vStarts2-vStarts))) == 0;
            print "  PASS V3a starts are deterministic (redraw identical)";
        else;
            nFail = nFail+1;
            print "  FAIL V3a starts differ on redraw";
        endif;

        /* V1 exact: X = start 2 (not a fixed point) and X = the safe
           Anderson endpoint from the Stone start. */
        vCtlX = vCtl;
        vCtlX.b0 = vbStone;
        vOutAnd = _quaidsFitAnderson(vw, vIntcpt, vPrices, vTotexp, vInstr,
            vCtlX, 0, 8, 1, 1, 1e-10, 10, 10);
        xi = 1;
        do while xi <= 2;
            if xi == 1;
                bX = vStarts[vNg+1:2*vNg, .];
                xName = "start 2";
            else;
                bX = vOutAnd.homogB;
                xName = "Stone-start Anderson endpoint (converged=" $+ ftos(vOutAnd.converged, "%*.*lf", 1, 0) $+ ")";
            endif;
            vCtlX = vCtl;
            vCtlX.maxiter = 2;
            vCtlX.b0 = bX;
            vOutF = quaidsFit(vw, vIntcpt, vPrices, vTotexp, vInstr, vCtlX);
            vOutA = _quaidsFitAnderson(vw, vIntcpt, vPrices, vTotexp, vInstr,
                vCtlX, 0, 0);
            /* index used by quaidsFit's final variance pass */
            if vOutF.iterations == 2;
                { d1, d2, d3, bIdx, d5 } = _pilotRecoverBS(vw, vIntcpt,
                    vPrices, vTotexp, vInstr, vCtl, bX, bX);
            else;
                bIdx = bX;
            endif;
            { bSH, okH, symH, TbH, critH } = _pilotRecoverBS(vw, vIntcpt,
                vPrices, vTotexp, vInstr, vCtl, bIdx, vOutF.homogB);
            dTb = maxc(maxc(abs(TbH-vOutF.homogB)));
            dBS = maxc(maxc(abs(bSH-vOutF.bS)));
            dRec = abs(maxc(maxc(abs(bSH-vTrue))) - maxc(maxc(abs(vOutF.bS-vTrue))));
            dCrit = abs(critH-vOutF.homogCrit);
            dFA = maxc(maxc(abs(vOutA.bS-vOutF.bS)));
            print "  X =" $+ xName;
            print "    quaidsFit iterations" vOutF.iterations "symValid" vOutF.symValid
                "helper ok/symValid" okH symH;
            print "    max|Tb_helper-homogB|" dTb " max|bS_helper-bS|" dBS;
            print "    |recErr_helper-recErr_quaidsFit|" dRec " |crit_helper-homogCrit|" dCrit;
            print "    recErr(quaidsFit bS)" maxc(maxc(abs(vOutF.bS-vTrue)))
                " max|bS_Anderson(mDepth0)-bS_quaidsFit|" dFA;
            nCheck = nCheck+1;
            if okH and symH == vOutF.symValid and dTb <= v1Tol and dBS <= v1Tol
                and dRec <= v1Tol and dCrit <= v1Tol and dFA <= v1Tol;
                print "    PASS V1 exact reproduction (tol" v1Tol ")";
            else;
                nFail = nFail+1;
                print "    FAIL V1 exact reproduction";
            endif;
            xi = xi+1;
        endo;

        /* V1 shared endpoint: polish the Stone-start Anderson endpoint
           with plain fixed-point iterations to err 1e-12, where index
           and coefficient coincide; the direct representation (b,b)
           must then give quaidsFit()'s own bS. */
        vCtlX = vCtl;
        vCtlX.maxiter = 500;
        vCtlX.err = 1e-12;
        vCtlX.b0 = vOutAnd.homogB;
        vOutPol = quaidsFit(vw, vIntcpt, vPrices, vTotexp, vInstr, vCtlX);
        { bSH, okH, symH, TbH, critH } = _pilotRecoverBS(vw, vIntcpt,
            vPrices, vTotexp, vInstr, vCtl, vOutPol.homogB, vOutPol.homogB);
        dBS = maxc(maxc(abs(bSH-vOutPol.bS)));
        dRec = abs(maxc(maxc(abs(bSH-vTrue))) - maxc(maxc(abs(vOutPol.bS-vTrue))));
        print "  shared fixed-point endpoint (polished, converged=" vOutPol.converged
            "finalErr" vOutPol.finalErr "iterations" vOutPol.iterations ")";
        print "    direct representation (b,b) vs quaidsFit: max|dbS|" dBS
            " |dRecErr|" dRec " recErr" maxc(maxc(abs(vOutPol.bS-vTrue)));
        nCheck = nCheck+1;
        if vOutPol.converged and dBS <= 1e-6 and dRec <= 1e-6;
            print "    PASS V1 shared endpoint (tol 1e-6; residual index/coef gap is the 1e-12 stopping rule)";
        elseif not vOutPol.converged;
            print "    INFO V1 shared endpoint: plain iteration did not reach 1e-12 from this endpoint; not scored";
            nCheck = nCheck-1;
        else;
            nFail = nFail+1;
            print "    FAIL V1 shared endpoint";
        endif;
        print;
        vsi = vsi+1;
    endo;
    vModel = vModel+1;
endo;

print "================================================================";
if nFail == 0;
    print "PART 1 (V1/V2/V3a): ALL" nCheck "CHECKS PASSED";
else;
    print "" $+ "PART 1 (V1/V2/V3a): " $+ ftos(nFail, "%*.*lf", 1, 0) $+ " CHECKS FAILED of " $+ ftos(nCheck, "%*.*lf", 1, 0);
    print "STOPPING before the harness run.";
    end;
endif;
print "================================================================";
print;

/* ================= Part 2: harness on the validation seeds ========= */
pilotModels = 0|1;
pilotSeeds = 1|6|9;
pilotStarts = seqa(1, 1, 16);
pilotMethods = 1|2|3;
pilotOutDir = "pilot_validation";
pilotVerbose = 0;
#include pilot_harness_core.src;

print "pilot_validate.e: run complete.";
