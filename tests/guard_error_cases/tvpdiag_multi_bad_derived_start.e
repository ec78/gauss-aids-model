new;
library cmlmt, tsmt, sslib;
#include tvpdiag_includes.inc;

struct quaidsTVPControl ctl;
struct quaidsTVPMultiStartOut d;
prices = zeros(2, 3);
d = quaidsTVPMultiStartDiagnostics(zeros(2, 3), prices, zeros(2, 1),
    eye(2), 1e308*ones(7, 1), ctl);
