new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc;

struct quaidsTVPOut tv;
struct quaidsTVPDiagOut d;
tv.mleRetcode = 0;
tv.q0 = ones(1, 1);
tv.Qfit = ones(1, 1);
d = quaidsTVPDiagnostics(tv, 0);
