new;
/*
** Isolated correctness check for Anderson acceleration (Type-II, depth m)
** on a toy linear fixed point x_{k+1} = A*x_k + c, BEFORE wiring it into
** anything AIDS-related. A is deliberately chosen with spectral radius
** > 1 so plain (undamped) Picard iteration diverges -- exactly the
** "vanishing/exploding" analog discussed: convergence of x_{k+1}=T(x_k)
** is governed by the spectral radius of T's Jacobian (here, just A
** itself, since T is linear).
*/

n = 5;
rndseed 42;
Araw = rndn(n, n)*0.3;
// Scale so the spectral radius is ~1.4 (genuinely divergent under plain
// iteration), by scaling by the actual max abs eigenvalue.
va = eig(Araw);
specrad = maxc(abs(va));
A = Araw * (1.4/specrad);
va2 = eig(A);
print "spectral radius of A (should be ~1.4):" maxc(abs(va2));

c = rndn(n, 1);
xTrue = inv(eye(n) - A)*c;   // (I-A)^-1 * c -- (I-A) is not symmetric/PD, so solpd() doesn't apply
print "true fixed point (first 3 elems):" xTrue[1:3]';

proc (1) = Tmap(x);
    retp(A*x + c);
endp;

/* ---------------------------------------------------------------------
** Plain (undamped) Picard iteration -- expected to diverge.
** --------------------------------------------------------------------- */
x = zeros(n, 1);
k = 1;
diverged = 0;
ok = 1;
do while ok;
    xnew = Tmap(x);
    if maxc(abs(xnew)) > 1e6;
        diverged = 1;
        ok = 0;
    else;
        x = xnew;
        k = k + 1;
        if k > 50;
            ok = 0;
        endif;
    endif;
endo;
print "plain Picard: diverged =" diverged " after" k " iters, maxabs(x) =" maxc(abs(x));

/* ---------------------------------------------------------------------
** Anderson(m) acceleration, Type-II, mixing parameter beta.
**
** History: Xhist (n x depth), Ghist (n x depth) hold the last `depth`
** iterates and their residuals g_i = T(x_i) - x_i, most-recent LAST.
** At each step: build DeltaX = Xhist[.,2:end]-Xhist[.,1:end-1] style
** differences relative to the CURRENT point, solve a small least-squares
** problem for the mixing coefficients, form the new iterate.
** --------------------------------------------------------------------- */
depth = 6;
beta = 1.0;

x = zeros(n, 1);
// GAUSS's zeros(n,0) is out of range -- no true empty-matrix idiom used
// here; track history validity with an explicit counter instead (matches
// this codebase's own "scalar 0 as unset sentinel" convention elsewhere).
Xhist = zeros(n, depth);
Ghist = zeros(n, depth);
histCount = 0;
k = 1;
diverged = 0;
converged = 0;
tol = 1e-10;
ok = 1;
do while ok;
    Tx = Tmap(x);
    g = Tx - x;

    if maxc(abs(g)) < tol;
        converged = 1;
        ok = 0;
    elseif maxc(abs(Tx)) > 1e8;
        diverged = 1;
        ok = 0;
    else;
        if histCount == 0;
            // First step: no history yet, plain step (beta-damped).
            xnew = x + beta*g;
        else;
            // Build difference matrices relative to the CURRENT (x, g),
            // using only the FIRST histCount columns (the valid ones).
            mUse = histCount;
            DX = x - Xhist[., 1:mUse];      // n x mUse, each col = x - x_i
            DG = g - Ghist[., 1:mUse];      // n x mUse, each col = g - g_i
            // Solve min_mixCoef || g - DG*mixCoef ||^2 (least squares,
            // ridge-stabilized since DG can be near-singular). "gamma" is
            // a reserved GAUSS identifier (CLAUDE.md gotcha) -- mixCoef instead.
            mixCoef = invpd(DG'DG + 1e-10*eye(mUse)) * (DG'g);
            xnew = x + beta*g - (DX + beta*DG)*mixCoef;
        endif;

        // Push (x, g) onto history (most-recent in column 1), keep only
        // the last `depth`.
        if histCount < depth;
            Xhist[., histCount+1] = x;
            Ghist[., histCount+1] = g;
            histCount = histCount + 1;
        else;
            Xhist = Xhist[., 2:depth]~x;
            Ghist = Ghist[., 2:depth]~g;
        endif;

        x = xnew;
        k = k + 1;
        if k > 200;
            ok = 0;
        endif;
    endif;
endo;

print "Anderson(" $+ ftocv(depth,1,0) $+ "): converged =" converged " diverged =" diverged " after" k " iters";
print "recovered x (first 3 elems):" x[1:3]';
print "max abs error vs true fixed point:" maxc(abs(x - xTrue));
