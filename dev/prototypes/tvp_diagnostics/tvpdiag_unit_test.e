new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc

ncheck = 0;
nfail = 0;
proc 0 = check(cond, msg);
    ncheck = ncheck + 1;
    if not cond;
        nfail = nfail + 1;
        print "FAIL:" msg;
    endif;
endp;

struct quaidsTVPOut t;
struct _quaidsTVPDiagSingleOut d;
t.mleRetcode = 0;
t.q0 = ones(3, 1);
t.Qfit = 1|1e-4|1e-3;
d = _quaidsTVPDiagSingle(t);
call check(d.status == 3 and d.nScaleCollapse == 1,
    "single: strict < threshold counts one collapse; equality does not count");
call check(d.minQOverQ0 == 1e-4 and d.maxQOverQ0 == 1,
    "single: reports qfit/q0 extrema");

t.mleRetcode = 2;
d = _quaidsTVPDiagSingle(t);
call check(d.status == 1, "single: FIT-FAILED has precedence over scale caution");
t.mleRetcode = 0;
t.Qfit = 1|miss(0, 0)|1;
d = _quaidsTVPDiagSingle(t);
call check(d.status == 2 and scalmiss(d.nScaleCollapse),
    "single: invalid Q is explicit and suppresses ratios");
t.Qfit = 1|2|3;
d = _quaidsTVPDiagSingle(t);
call check(d.status == 0 and d.nScaleCollapse == 0, "single: clean case is OK");

tolInput = 1e-6*200;
{ stable, bestStart, bestLL, maxGap, tol } =
    _tvpDiagMultiSummary(0|0|0, tolInput|0|tolInput/2, 200);
call check(stable == 1 and maxGap == tol,
    "multi: LL tolerance boundary is inclusive");
{ stable, bestStart, bestLL, maxGap, tol } =
    _tvpDiagMultiSummary(0|0|0, tolInput+1e-9|0|tolInput/2, 200);
call check(stable == 0, "multi: gap above tolerance is start-sensitive");
{ stable, bestStart, bestLL, maxGap, tol } =
    _tvpDiagMultiSummary(0|2|0, 10|10|10, 200);
call check(stable == 0, "multi: any nonzero retcode prevents stable label");
call check(bestStart == 1 and bestLL == 10, "multi: exact LL tie selects first start");

print "TVPDIAG_CHECKS" ncheck;
print "TVPDIAG_FAILURES" nfail;
if nfail;
    end;
endif;
print "TVPDIAG_UNIT_PASS";
