# generate_r_reference.R
#
# Reproduces the R reference numbers hardcoded into
# tests/quaids_published_validation_test.e. Requires the R package
# micEconAids (CRAN, GPL >= 2): install.packages("micEconAids").
#
# Run from this directory (tests/fixtures/published/):
#   Rscript generate_r_reference.R

library(micEconAids)
data("Blanciforti86")

B86 <- Blanciforti86[1:32, ]
priceNames <- c("pFood1", "pFood2", "pFood3", "pFood4")
shareNames <- c("wFood1", "wFood2", "wFood3", "wFood4")

# Same identification strategy as GAUSS's quaidsFit(): instrument
# log(xFood) with log(xAgg) (total aggregate expenditure), treating log
# prices as included/exogenous instruments. R estimates this via 3SLS
# (systemfit); GAUSS uses a control-function/residual-inclusion approach.
# Both are valid, asymptotically consistent IV estimators for the same
# model -- not expected to match bit-for-bit, only closely (observed max
# abs difference ~0.021; see quaids_published_validation_test.e).
for (i in 1:4) {
    B86[[paste0("l", priceNames[i])]] <- log(B86[[priceNames[i]]])
}
B86[["lxAgg"]] <- log(B86[["xAgg"]])
instNames <- c(paste0("l", priceNames), "lxAgg")

ivResult <- aidsEst(priceNames, shareNames, "xFood", data = B86,
    priceIndex = "Ls", instNames = instNames)

cat("alpha:\n"); print(ivResult$coef$alpha)
cat("beta:\n"); print(ivResult$coef$beta)
cat("gamma:\n"); print(ivResult$coef$gamma)
cat("estMethod:", ivResult$est$method, "\n")

# Milestone 9: iterated (nonlinear translog price index) AIDS reference,
# for cross-validating GAUSS's aCtl.linear=1, aCtl.maxiter>1 case --
# method="IL" is micEconAids's Iterated Linear Least Squares Estimator
# (Blundell & Robin 1999): LA-AIDS with priceIndex="Ls" as the starting
# value, then iterates using the (nonlinear) translog price index -- the
# same starting-value/iteration structure quaidsFit() uses. No IV variant
# of method="IL" is offered by micEconAids (instNames is only documented
# for the LA/3SLS path), so this reference is SUR-estimated, not IV --
# GAUSS's iterated fit still instruments log(xFood) (always endogenous in
# this library), so this comparison is across both a different estimation
# algorithm AND an IV-vs-no-IV difference, expect a wider gap than the
# LA-AIDS-vs-3SLS-IV comparison above.
ilResult <- aidsEst(priceNames, shareNames, "xFood", data = B86,
    method = "IL", priceIndex = "Ls")

cat("IL alpha:\n"); print(ilResult$coef$alpha)
cat("IL beta:\n"); print(ilResult$coef$beta)
cat("IL gamma:\n"); print(ilResult$coef$gamma)
cat("IL iterations:", ilResult$ILiter, "\n")
