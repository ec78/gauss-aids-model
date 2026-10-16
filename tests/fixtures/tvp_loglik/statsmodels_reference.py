# statsmodels_reference.py
#
# TVP-AIDS promotion P5: independent fixed-Q log-likelihood reference for
# the TVP-AIDS state space, computed with statsmodels' exact diffuse Kalman
# filter (statsmodels.tsa.statespace.kalman_filter.KalmanFilter with
# initialize_diffuse(); Durbin & Koopman 2012 ch. 5, univariate treatment).
#
# VALIDATION ONLY. Nothing in src/ or the installed quaids package imports
# or needs Python; package users never need it. Requirements: Python 3,
# numpy, statsmodels (the committed reference was produced with the
# versions recorded in tvp_loglik_reference.json).
#
# Why this is independent of the GAUSS path it checks:
#   - The filter is statsmodels' compiled Cython implementation, a separate
#     codebase from sslib's kalmanFilterDiffuseTVP() (gauss-state-space).
#   - The design matrices Z_t are rebuilt HERE from the raw regressors
#     (relative log prices, Stone-deflated log expenditure), from the
#     documented state layout [alpha | gamma unique upper triangle,
#     row-major | beta], not read from GAUSS's _quaidsTVPBuildZ() output.
#   - Q, H and the data come from committed text files that both sides
#     read; neither side's numbers are produced by the other.
#
# Model (identical to the GAUSS side by construction, see the README):
#   w_t[1:n1] = Z_t a_t + eps_t,  eps_t ~ N(0, H)
#   a_(t+1)   = a_t + eta_t,      eta_t ~ N(0, diag(q))
#   every state exactly diffuse (a_1 = 0, P_inf = I, P_star = 0).
# Log-likelihood convention (matches sslib's, checked on a local-level
# toy model while this was written): a period's contribution is
#   -0.5 * sum_j [ ln(2 pi) + ln(F_inf_j) ]          for diffuse elements,
#   -0.5 * sum_j [ ln(2 pi) + ln(F_j) + v_j^2/F_j ]  otherwise,
# summed over ALL periods (loglikelihood_burn = 0 -- statsmodels'
# MLEModel would otherwise burn the diffuse periods by default).
#
# Run from this directory (tests/fixtures/tvp_loglik/):
#   python statsmodels_reference.py            # rewrite the reference
#   python statsmodels_reference.py --check    # recompute and compare to
#                                              # the committed reference

import json
import platform
import sys

import numpy as np
import statsmodels
from statsmodels.tsa.statespace.kalman_filter import KalmanFilter

FIXTURE = "tvp_loglik_fixture.txt"
CASES = "tvp_loglik_cases.txt"
REF_TXT = "tvp_loglik_reference.txt"
REF_JSON = "tvp_loglik_reference.json"

N1 = 2


def build_design(prel, lx, n1):
    """Z_t for every period, shape (n1, k_states, nobs)."""
    nobs = prel.shape[0]
    ngamma = n1 * (n1 + 1) // 2
    k_states = 2 * n1 + ngamma
    gidx = np.zeros((n1, n1), dtype=int)
    cnt = 0
    for i in range(n1):
        for j in range(i, n1):
            gidx[i, j] = gidx[j, i] = cnt
            cnt += 1
    z = np.zeros((n1, k_states, nobs))
    for t in range(nobs):
        for i in range(n1):
            z[i, i, t] = 1.0
            for j in range(n1):
                z[i, n1 + gidx[i, j], t] += prel[t, j]
            z[i, n1 + ngamma + i, t] = lx[t]
    return z


def loglik_obs(y, z, q, h):
    n1, k_states, nobs = z.shape
    kf = KalmanFilter(k_endog=n1, k_states=k_states, k_posdef=k_states)
    kf.bind(np.ascontiguousarray(y, dtype=float))
    kf["design"] = z
    kf["obs_cov"] = h
    kf["transition"] = np.eye(k_states)
    kf["selection"] = np.eye(k_states)
    kf["state_cov"] = np.diag(q)
    kf.initialize_diffuse()
    kf.loglikelihood_burn = 0
    res = kf.filter()
    return np.asarray(res.llf_obs, dtype=float), int(res.nobs_diffuse)


def main():
    dat = np.loadtxt(FIXTURE, skiprows=1)
    cases = np.atleast_2d(np.loadtxt(CASES, skiprows=1))
    y = dat[:, 1:1 + N1]
    prel = dat[:, 4:4 + N1]
    lx = dat[:, 6]
    z = build_design(prel, lx, N1)

    llf = []
    summary = []
    for row in cases:
        cid = int(row[0])
        q = row[1:8]
        h = np.array([[row[8], row[9]], [row[9], row[10]]])
        lo, ndiff = loglik_obs(y, z, q, h)
        llf.append(lo)
        summary.append({"case": cid, "llf_total": float(lo.sum()),
                        "nobs_diffuse": ndiff})
    llf = np.column_stack(llf)

    meta = {
        "description": "statsmodels exact-diffuse fixed-Q TVP-AIDS log-likelihood reference (P5)",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "statsmodels": statsmodels.__version__,
        "platform": platform.platform(),
        "fixture": FIXTURE,
        "cases": CASES,
        "nobs": int(y.shape[0]),
        "k_endog": N1,
        "k_states": int(z.shape[1]),
        "summary": summary,
    }

    if "--check" in sys.argv:
        ref = np.loadtxt(REF_TXT, skiprows=1)[:, 1:]
        d = np.max(np.abs(ref - llf))
        print(f"max |recomputed - committed| per-period llf: {d:.3e}")
        sys.exit(0 if d == 0.0 else 1)

    with open(REF_TXT, "w", newline="\n") as f:
        f.write("t " + " ".join(f"case{s['case']}" for s in summary) + "\n")
        for t in range(llf.shape[0]):
            f.write(" ".join([str(t + 1)] + [repr(float(v)) for v in llf[t]]) + "\n")
    with open(REF_JSON, "w", newline="\n") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")
    for s in summary:
        print(f"case {s['case']}: llf = {s['llf_total']!r}, nobs_diffuse = {s['nobs_diffuse']}")
    print(f"statsmodels {statsmodels.__version__}, numpy {np.__version__}, Python {platform.python_version()}")


if __name__ == "__main__":
    main()
