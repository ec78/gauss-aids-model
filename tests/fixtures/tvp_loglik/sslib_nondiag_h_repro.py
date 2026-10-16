# sslib_nondiag_h_repro.py
#
# TVP-AIDS promotion P5: hand verification of the one case where GAUSS
# (sslib v2.1.0 kalmanFilterDiffuseTVP) and statsmodels disagree -- case 3
# of tvp_loglik_cases.txt, a NON-diagonal H.
#
# VALIDATION ONLY (numpy + statsmodels_reference.build_design; no runtime
# dependency of the quaids package).
#
# A plain numpy re-implementation of the univariate exact diffuse filter
# (Durbin & Koopman 2012 sec. 5.2/6.4), written from the recursion, with
# two ways to decorrelate a non-diagonal H:
#   "chol" -- sslib's way: rows scaled by Linv = inv(chol(H)), unit
#             observation variances, -0.5*ln|H| Jacobian per period;
#   "ldl"  -- unit-lower-triangular L from H = L D L', variances D, no
#             Jacobian (statsmodels' way).
# Both are exact in real arithmetic. A diffuse ("Case 1") update happens
# when an element's F_inf = Z_j P_inf Z_j' exceeds a threshold; sslib uses
# an absolute 1e-12 and only ends the diffuse phase at the START of a
# period (_kfDiffusePhaseOver).
#
# Result (recorded in README.md): with "chol" and 1e-12, period 4 element 2
# sees a rounding residue F_inf ~ 2e-11 (the last diffuse direction was
# removed by element 1 of the same period, and the 1/sqrt(H) row scaling
# inflates the residue past 1e-12), so it takes a spurious Case-1 update and
# its log-likelihood jumps by about +14. "ldl" (residue ~1e-15), or "chol"
# with a threshold of 1e-10 or 1e-8, matches statsmodels to ~1e-13. The
# exact jump size depends on the residue's rounding, so this replica's total
# (688.26) and sslib's (688.03) differ; the mechanism and the period are the
# same. Diagonal H (cases 1-2) is not row-scaled and agrees for every
# variant.
#
# Run from this directory (tests/fixtures/tvp_loglik/):
#   python sslib_nondiag_h_repro.py

import sys

import numpy as np

sys.dont_write_bytecode = True  # keep the fixture directory clean
from statsmodels_reference import CASES, FIXTURE, REF_TXT, N1, build_design


def univariate_diffuse(y, z, q, h, mode, thr=1e-12, trace_t=None):
    n1, k, nobs = z.shape
    a = np.zeros(k)
    p_inf0 = np.eye(k)
    p_inf = p_inf0.copy()
    p_star = np.zeros((k, k))
    diffuse_on = True
    if abs(h[0, 1]) == 0.0:
        lmat, hrow, logdet = np.eye(n1), np.diag(h).copy(), 0.0
    elif mode == "chol":
        c = np.linalg.cholesky(h)
        lmat, hrow, logdet = np.linalg.inv(c), np.ones(n1), np.log(np.linalg.det(h))
    else:
        c = np.linalg.cholesky(h)
        lmat, hrow, logdet = np.linalg.inv(c / np.diag(c)), np.diag(c) ** 2, 0.0
    ll = np.zeros(nobs)
    for t in range(nobs):
        if diffuse_on and (p_inf * p_inf).sum() <= 1e-10 * (p_inf0 * p_inf0).sum():
            p_inf = np.zeros((k, k))
            diffuse_on = False
        zr = lmat @ z[:, :, t]
        ys = lmat @ y[t]
        s = 0.0
        for j in range(n1):
            zj = zr[j]
            v = ys[j] - zj @ a
            f_inf = zj @ p_inf @ zj
            f_star = zj @ p_star @ zj + hrow[j]
            if trace_t is not None and t + 1 == trace_t:
                print(f"    period {t + 1}, element {j + 1}: F_inf = {f_inf:.3e} ({mode})")
            if f_inf > thr:
                kg = p_inf @ zj / f_inf
                s += np.log(f_inf)
                a = a + kg * v
                p_star = (p_star - np.outer(kg, zj @ p_star) - np.outer(p_star @ zj, kg)
                          + np.outer(kg, kg) * f_star)
                p_inf = p_inf - np.outer(kg, kg) * f_inf
            else:
                kg = p_star @ zj / f_star
                s += np.log(f_star) + v * v / f_star
                a = a + kg * v
                p_star = p_star - np.outer(kg, kg) * f_star
        ll[t] = -0.5 * (n1 * np.log(2 * np.pi) + s) - 0.5 * logdet
        p_star = p_star + np.diag(q)
        p_star = (p_star + p_star.T) / 2
    return ll


def main():
    dat = np.loadtxt(FIXTURE, skiprows=1)
    cases = np.atleast_2d(np.loadtxt(CASES, skiprows=1))
    ref = np.loadtxt(REF_TXT, skiprows=1)[:, 1:]
    y = dat[:, 1:1 + N1]
    z = build_design(dat[:, 4:4 + N1], dat[:, 6], N1)
    for c, row in enumerate(cases):
        q = row[1:8]
        h = np.array([[row[8], row[9]], [row[9], row[10]]])
        variants = [("chol", 1e-12), ("ldl", 1e-12)]
        if h[0, 1] != 0.0:
            variants += [("chol", 1e-10), ("chol", 1e-8)]
        for mode, thr in variants:
            ll = univariate_diffuse(y, z, q, h, mode, thr,
                                    trace_t=4 if (h[0, 1] != 0.0 and thr == 1e-12) else None)
            d = np.abs(ll - ref[:, c])
            print(f"case {int(row[0])} {mode:4s} thr={thr:.0e}: llf = {ll.sum():.10f}, "
                  f"max |diff vs statsmodels| = {d.max():.2e} (first period > 1e-9: "
                  f"{int(np.argmax(d > 1e-9)) + 1 if (d > 1e-9).any() else 'none'})")


if __name__ == "__main__":
    main()
