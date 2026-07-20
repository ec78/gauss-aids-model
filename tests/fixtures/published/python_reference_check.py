# python_reference_check.py
#
# Independent, from-scratch Python cross-check of the same LA-AIDS +
# control-function-IV specification quaidsFit() estimates (see
# quaids_published_validation_test.e), hand-coded from the
# Deaton-Muellbauer (1980) equations rather than using a package (no
# comparably-established Python AIDS package exists the way R has
# micEconAids). Requires numpy, pandas: pip install numpy pandas.
#
# Run from this directory (tests/fixtures/published/):
#   python python_reference_check.py
#
# Status: broadly consistent with both GAUSS and R for most coefficients
# (e.g. equation 1: alpha=-0.342/beta=0.380 here vs GAUSS's
# alpha=-0.345/beta=0.381 -- essentially exact), but shows larger residual
# differences on equation 2 (alpha=-0.032 here vs GAUSS's 0.033 and R's
# 0.033, which agree with each other) than either GAUSS-vs-R comparison
# does. Since R -- an independently authored, published, widely-used
# implementation -- agrees closely with GAUSS on exactly the coefficients
# where this script diverges, that divergence is attributed to this
# from-scratch replica (most likely in how its residual-covariance-weighted
# GLS handles the added IV-residual regressor) rather than to GAUSS. Kept
# here for transparency and as a documented starting point, not as a
# pass/fail reference -- see GOLD_STANDARD_TODO.md Milestone 3.

import numpy as np
import pandas as pd

df = pd.read_csv("blanciforti86_food32.csv")

N = 4
T = len(df)
prices = np.log(df[["pFood1", "pFood2", "pFood3", "pFood4"]].to_numpy())
shares = df[["wFood1", "wFood2", "wFood3", "wFood4"]].to_numpy()
lnX = np.log(df["xFood"].to_numpy())
lnInstr = np.log(df["xAgg"].to_numpy())
lr = prices[:, 0:3] - prices[:, 3:4]  # relative prices, good 4 = reference

# First-stage OLS: totexp on [1, relative prices, absolute reference
# price, instrument] -- matches GAUSS's z = intcpt~prices~instr, where
# `prices` there is the full n-column (3 relative + 1 absolute) matrix.
Z1 = np.column_stack([np.ones(T), lr, prices[:, 3:4], lnInstr])
bIV = np.linalg.lstsq(Z1, lnX, rcond=None)[0]
u = lnX - Z1 @ bIV
print("first-stage R2:", 1 - np.var(u) / np.var(lnX))

# Stone index: fixed mean-share weights over the (correct) absolute prices.
meanw = shares.mean(axis=0)
lx = lnX - prices @ meanw

# Symmetric-GLS system (homogeneity via relative-price reparametrization,
# symmetry via a shared parameter for gamma_ij/gamma_ji), with the
# first-stage residual u as an extra control-function regressor.
X1 = np.column_stack([np.ones(T), lr, lx, u])
resid = np.zeros((T, 3))
for i in range(3):
    b_ols = np.linalg.lstsq(X1, shares[:, i], rcond=None)[0]
    resid[:, i] = shares[:, i] - X1 @ b_ols
Sigma = (resid.T @ resid) / T

nParam = 3 + 6 + 3 + 3  # alpha(3) + gamma-unique(6) + beta(3) + u-coef(3)
Z = np.zeros((3 * T, nParam))
Y = np.zeros(3 * T)
gammaIdx = {(0,0):3,(0,1):4,(0,2):5,(1,0):4,(1,1):6,(1,2):7,(2,0):5,(2,1):7,(2,2):8}
for i in range(3):
    rows = slice(i * T, (i + 1) * T)
    Y[rows] = shares[:, i]
    Z[rows, i] = 1.0
    for j in range(3):
        Z[rows, gammaIdx[(i, j)]] += lr[:, j]
    Z[rows, 9 + i] = lx
    Z[rows, 12 + i] = u

W = np.kron(np.linalg.inv(Sigma), np.eye(T))
theta = np.linalg.solve(Z.T @ W @ Z, Z.T @ W @ Y)

alpha3, beta3, u3 = theta[0:3], theta[9:12], theta[12:15]
gammaSub = np.array([[theta[3], theta[4], theta[5]],
                      [theta[4], theta[6], theta[7]],
                      [theta[5], theta[7], theta[8]]])

alpha = np.append(alpha3, 1.0 - alpha3.sum())
beta = np.append(beta3, -beta3.sum())
ucoef = np.append(u3, -u3.sum())
gamma = np.zeros((4, 4))
gamma[0:3, 0:3] = gammaSub
gamma[0:3, 3] = -gammaSub.sum(axis=1)
gamma[3, 0:3] = -gammaSub.sum(axis=0)
gamma[3, 3] = gammaSub.sum()

np.set_printoptions(precision=6, suppress=True, linewidth=160)
print("alpha:", alpha)
print("beta:", beta)
print("u-coef:", ucoef)
print("gamma:\n", gamma)
