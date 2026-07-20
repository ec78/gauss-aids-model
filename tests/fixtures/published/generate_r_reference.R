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
