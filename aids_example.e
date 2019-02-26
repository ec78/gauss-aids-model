new;
#include aid_model.sdf;
#include aidsutil.src
#include aids_rev.src;

seed = 11;
tobs = 1000;
N = 5;

al = round(rndns(1,N-1,seed)*10)/10;
al = al~(1-sumc(al'));
al1 = .5*round(rndns(1,N-1,seed)*10)/10;
al1 = al1~(-sumc(al1'));
ga = round(rndns(N-1,N-1,seed)*10)/10;
ga = xpnd(vech(ga));
ga = ga|(-sumc(ga)');
ga = ga~(-sumc(ga'));
be = .5*round(rndns(1,N-1,seed)*10)/10;
be = be~(-sumc(be'));
la = .01*round(rndns(1,N-1,seed)*10)/10;
la = la~(-sumc(la'));
ro = round(rndns(1,N-1,seed)*10)/10;


prices = 1+rndns(tobs,N,seed);
instr = 5+5*rndns(tobs,1,seed);
intcpt = 2+2*rndns(tobs,1,seed);
u =  .1*rndns(tobs,1,seed);
totexp = .85*instr + u;
e = 2*rndns(tobs,N-1,seed) + u*ro;
e = e~(-sumc(e'));


a_p = sumc( (prices.*(al+intcpt*al1) )') + .5*sumc(((prices*ga).*prices)');
lx = totexp -a_p;
b_p = prices*be';
lx2 = (lx^2)./exp(b_p);

w = al +  prices*ga + lx*be + e +intcpt*al1 + lx2*la ;

xnames = { exog,p1,p2,p3,p4,p5,totexp,instr,w1,w2,w3,w4,w5,instr };

output file=out reset;

// Declare control structure
struct aidsControl aCtl;
aCtl = aidsControlCreate;

// Set parameters
aCtl.linear = 0;

// Maximum iterations
actl.maxiter = 100;

// Homogenous model
actl.homogenous = 1;

// Error tolerance
actl.err = .001;

// Run model
{ b, v, b_s, v_s } =  aids(w, intcpt, prices, totexp, instr, aCtl);

format /ld 5,3;

al|al1|ga|be|la|ro~-sumc(ro');

if aCtl.homogenous;
    b~miss(ones(rows(b),1),1)~b_s;
else;
    b;
endif;
output off;

