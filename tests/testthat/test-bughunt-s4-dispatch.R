## S4 dispatch ambiguity / shadowing bugs
## DO NOT FIX — these tests document known bugs. See beads issues for fix.

suppressPackageStartupMessages(library(amatrix))
suppressPackageStartupMessages(library(Matrix))

# Shared fixtures
mat   <- matrix(1:4 * 1.0, 2, 2)
A_adge <- adgeMatrix(mat)
A_dge  <- as(as(mat, "dMatrix"), "unpackedMatrix")   # plain dgeMatrix, NOT adgeMatrix
sp_m   <- as(matrix(c(1, 0, 0, 2) * 1.0, 2, 2), "dgCMatrix")
A_adgC <- adgCMatrix(sp_m)

# ── amatrix-7rl ──────────────────────────────────────────────────────────────
# Ops(dgeMatrix, adgeMatrix): Matrix's dMatrix#dMatrix method wins at distance
# (1,1) over amatrix's Ops(ANY,adgeMatrix) at (Inf,0), silently dropping
# backend metadata and returning plain dgeMatrix.
test_that("amatrix-7rl: dgeMatrix + adgeMatrix returns adgeMatrix [BUG: returns dgeMatrix]", {
  result <- A_dge + A_adge
  expect_s4_class(result, "adgeMatrix")
})

# ── amatrix-03p ──────────────────────────────────────────────────────────────
# Same inheritance path in the other direction: adgeMatrix IS-A dMatrix, so
# Ops(dMatrix,dMatrix) at (1,1) beats Ops(adgeMatrix,ANY) at (0,Inf).
test_that("amatrix-03p: adgeMatrix + dgeMatrix returns adgeMatrix [BUG: returns dgeMatrix]", {
  result <- A_adge + A_dge
  expect_s4_class(result, "adgeMatrix")
})

# ── amatrix-scl ──────────────────────────────────────────────────────────────
# Ops(dgCMatrix, adgCMatrix): sparseMatrix#sparseMatrix from Matrix shadows
# amatrix's Ops(ANY,adgCMatrix).  Result is plain dgCMatrix.
test_that("amatrix-scl: dgCMatrix + adgCMatrix returns adgCMatrix [BUG: returns dgCMatrix]", {
  result <- sp_m + A_adgC
  expect_s4_class(result, "adgCMatrix")
})

# ── amatrix-3gv ──────────────────────────────────────────────────────────────
# Mirror of amatrix-scl: adgCMatrix on the left.
test_that("amatrix-3gv: adgCMatrix + dgCMatrix returns adgCMatrix [BUG: returns dgCMatrix]", {
  result <- A_adgC + sp_m
  expect_s4_class(result, "adgCMatrix")
})

# ── amatrix-31q ──────────────────────────────────────────────────────────────
# [<- with missing j: no explicit setReplaceMethod for
# (adgeMatrix, index, missing, *) so the fallback fires through dgeMatrix's
# inherited [<- and the result loses adgeMatrix class.
test_that("amatrix-31q: A[1,] <- 99 preserves adgeMatrix class [BUG: returns dgeMatrix]", {
  tmp <- adgeMatrix(mat)
  tmp[1, ] <- 99
  expect_s4_class(tmp, "adgeMatrix")
})

test_that("amatrix-31q: A[,1] <- 88 preserves adgeMatrix class [BUG: returns dgeMatrix]", {
  tmp <- adgeMatrix(mat)
  tmp[, 1] <- 88
  expect_s4_class(tmp, "adgeMatrix")
})
