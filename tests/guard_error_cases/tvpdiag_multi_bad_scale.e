new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc;

struct quaidsTVPControl ctl;
struct quaidsTVPMultiStartOut d;
d = quaidsTVPMultiStartDiagnostics(0, 0, 0, 0, 0, ctl, 1);
