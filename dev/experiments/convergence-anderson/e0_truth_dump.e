new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;
#include direct_objective.src;
#include pilot_lib.src;

/*
** e0_truth_dump.e -- experiment E0 of CORRECTNESS_VALIDATION_PROTOCOL.md.
**
** NO-SOLVER truth regeneration for the pilot's seeds 1-30, both models.
** Per (model, seed) it runs only what pilot_harness_core.src's own seed
** setup runs before any method is called:
**   - _quaidsSyntheticDGP(3000, seed, model, 1)        (truth)
**   - quaidsFit(..., maxiter=1)                       (Stone/LA-AIDS: one
**                                                      linear GLS solve,
**                                                      no iteration)
**   - _directObjCoreHomotopy at homotopy 1 for the truth and for bStone
**                                                     (objective
**                                                      evaluations only)
** It calls no LM, continuation, or Anderson routine.
**
** Output e0_out/truth_m<model>.csv, one row per seed (no header):
**   model, seed, ng, n1, critTrue, critStone, vec(trueParams) [(ng+1)*n,
**   column-major, the same layout as the pilot's stored vec(bS)].
**
** Run from this directory: C:\gauss26\tgauss.exe -b -x e0_truth_dump.e
*/

struct quaidsControl eCtl;
struct quaidsControl eCtlStone;
struct quaidsOut eOutStone;

eModel = 0;
do while eModel <= 1;
    eFile = "e0_out/truth_m" $+ ftos(eModel, "%*.*lf", 1, 0) $+ ".csv";
    if filesa(eFile) $/= "";
        print "E0: " $+ eFile $+ " already exists -- STOP (remove it to regenerate).";
        end;
    endif;
    eSeed = 1;
    do while eSeed <= 30;
        eCtl = quaidsControlCreate;
        eCtl.linear = 1-eModel;
        eCtl.maxiter = 100;
        eCtl.err = .0001;
        eCtl.homogenous = 1;
        { ew, eIntcpt, ePrices, eTotexp, eInstr, eTrue } =
            _quaidsSyntheticDGP(3000, eSeed, eModel, 1);
        eCtlStone = eCtl;
        eCtlStone.maxiter = 1;
        eOutStone = quaidsFit(ew, eIntcpt, ePrices, eTotexp, eInstr, eCtlStone);
        ebStone = eOutStone.homogB;
        eIntFull = eOutStone.intcptFull;
        eNIntFull = cols(eIntFull);
        eNGoods = cols(ePrices);
        eN1 = eNGoods-1;
        eNg = eOutStone.ng;
        eTrueRaw = _directTrueRaw(eTrue, eNIntFull, eN1, eNGoods);
        eScale = _directObjScale(ebStone, eN1);

        { efT, efeasT, eresT } = _directObjCoreHomotopy(
            _directObjPack(eTrueRaw, eN1)./eScale, ew, eIntFull, ePrices,
            eTotexp, eOutStone.u, eCtl, eNg, eN1, eScale, 1);
        { efS, efeasS, eresS } = _directObjCoreHomotopy(
            _directObjPack(ebStone, eN1)./eScale, ew, eIntFull, ePrices,
            eTotexp, eOutStone.u, eCtl, eNg, eN1, eScale, 1);
        eCritTrue = miss(0, 0);
        eCritStone = miss(0, 0);
        if efeasT;
            eCritTrue = -efT;
        endif;
        if efeasS;
            eCritStone = -efS;
        endif;

        call _pilotAppendRow(eFile, eModel|eSeed|eNg|eN1|eCritTrue|eCritStone|vec(eTrue));
        print "truth m" eModel "seed" eSeed "critTrue" eCritTrue "critStone" eCritStone;
        eSeed = eSeed+1;
    endo;
    eModel = eModel+1;
endo;
print "e0_truth_dump.e: run complete.";
