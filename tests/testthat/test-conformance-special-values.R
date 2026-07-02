# Special-value conformance suite.
#
# Drives NA / NaN / Inf / -Inf (and mixtures) through the core amatrix ops and
# asserts parity with the identical operation on the plain base-R matrix, for
# every available backend. CPU runs at 1e-10; GPU backends compute in float32
# and are compared at 1e-4.
#
# Documented divergence from base R that this suite accepts:
#   * Row/column reductions (rowSums/colSums/rowMeans/colMeans) return NA where
#     base R returns NaN for a group containing a NaN. Both satisfy is.na(),
#     so the comparison matches NA/NaN positions rather than the exact sentinel.
#   * dist_matrix / kernel_matrix fail fast on non-finite input instead of
#     propagating; that policy is pinned in test-regressions-hunt.R (amatrix-6lg)
#     and those ops are intentionally excluded from this propagation suite.

# -- backend selection --------------------------------------------------------
# Always includes cpu; appends each installed + enabled GPU backend, mirroring
# the skip_if_not_installed + availability-probe idiom used elsewhere. GPU
# probing is env-gated (see tests/testthat/setup.R).
.special_value_backends <- function() {
  backends <- list(list(name = "cpu", tol = 1e-10))
  specs <- list(
    list(pkg = "amatrix.mlx", reg = "amatrix_mlx_register", avail = "amatrix_mlx_is_available", name = "mlx"),
    list(pkg = "amatrix.arrayfire", reg = "amatrix_arrayfire_register", avail = "amatrix_arrayfire_is_available", name = "arrayfire"),
    list(pkg = "amatrix.opencl", reg = "amatrix_opencl_register", avail = "amatrix_opencl_is_available", name = "opencl")
  )
  for (s in specs) {
    if (!requireNamespace(s$pkg, quietly = TRUE)) next
    try(getExportedValue(s$pkg, s$reg)(overwrite = TRUE), silent = TRUE)
    ok <- isTRUE(tryCatch(getExportedValue(s$pkg, s$avail)(), error = function(e) FALSE))
    if (ok) backends[[length(backends) + 1L]] <- list(name = s$name, tol = 1e-4)
  }
  backends
}

# -- NA/NaN/Inf-aware comparison ----------------------------------------------
# Materializes an amatrix result to a plain array/vector and asserts:
#   * identical shape,
#   * identical is.na() pattern (NA and NaN both count as missing), and
#   * value parity (within tol) on all non-missing entries, including +/-Inf.
.materialize <- function(x) {
  if (is.numeric(x) && is.null(dim(x))) {
    return(x)
  }
  as.matrix(x)
}

expect_conforms <- function(actual, reference, tol, label) {
  act <- .materialize(actual)
  testthat::expect_identical(dim(act), dim(reference), label = paste0(label, " [dim]"))
  testthat::expect_identical(
    as.vector(is.na(act)), as.vector(is.na(reference)),
    label = paste0(label, " [NA/NaN pattern]")
  )
  keep <- !is.na(as.vector(reference))
  testthat::expect_equal(
    as.vector(act)[keep], as.vector(reference)[keep],
    tolerance = tol, label = paste0(label, " [values]")
  )
}

# -- shared data with every non-finite flavour --------------------------------
.sv_data <- local({
  set.seed(101)
  A <- matrix(rnorm(12), 3, 4)
  B <- matrix(rnorm(12), 3, 4)
  A[1, 1] <- NA
  A[2, 2] <- NaN
  A[3, 3] <- Inf
  A[1, 4] <- -Inf
  sq <- matrix(rnorm(16), 4, 4)
  sq[1, 2] <- Inf
  sq[3, 4] <- NA
  list(A = A, B = B, sq = sq)
})

test_that("special values: elementwise arithmetic propagates like base R", {
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(.sv_data$A, preferred_backend = bk$name)
    B <- adgeMatrix(.sv_data$B, preferred_backend = bk$name)
    tag <- function(op) paste0("[", bk$name, "] A ", op, " B")
    expect_conforms(A + B, .sv_data$A + .sv_data$B, bk$tol, tag("+"))
    expect_conforms(A - B, .sv_data$A - .sv_data$B, bk$tol, tag("-"))
    expect_conforms(A * B, .sv_data$A * .sv_data$B, bk$tol, tag("*"))
    expect_conforms(A / B, .sv_data$A / .sv_data$B, bk$tol, tag("/"))
  }
})

test_that("special values: matmul/crossprod/tcrossprod propagate like base R", {
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(.sv_data$A, preferred_backend = bk$name) # 3x4
    B <- adgeMatrix(.sv_data$B, preferred_backend = bk$name) # 3x4
    sq <- adgeMatrix(.sv_data$sq, preferred_backend = bk$name) # 4x4
    tag <- function(op) paste0("[", bk$name, "] ", op)
    expect_conforms(A %*% sq, .sv_data$A %*% .sv_data$sq, bk$tol, tag("A %*% sq"))
    expect_conforms(crossprod(A), crossprod(.sv_data$A), bk$tol, tag("crossprod(A)"))
    expect_conforms(crossprod(A, B), crossprod(.sv_data$A, .sv_data$B), bk$tol, tag("crossprod(A,B)"))
    expect_conforms(tcrossprod(A), tcrossprod(.sv_data$A), bk$tol, tag("tcrossprod(A)"))
    expect_conforms(tcrossprod(A, B), tcrossprod(.sv_data$A, .sv_data$B), bk$tol, tag("tcrossprod(A,B)"))
  }
})

test_that("special values: rowSums/colSums/rowMeans/colMeans with and without na.rm", {
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(.sv_data$A, preferred_backend = bk$name)
    Ar <- .sv_data$A
    tag <- function(op) paste0("[", bk$name, "] ", op)
    # without na.rm: amatrix returns NA where base returns NaN (both is.na)
    expect_conforms(rowSums(A), rowSums(Ar), bk$tol, tag("rowSums"))
    expect_conforms(colSums(A), colSums(Ar), bk$tol, tag("colSums"))
    expect_conforms(rowMeans(A), rowMeans(Ar), bk$tol, tag("rowMeans"))
    expect_conforms(colMeans(A), colMeans(Ar), bk$tol, tag("colMeans"))
    # with na.rm=TRUE: missing entries dropped from both sum and mean denominator
    expect_conforms(rowSums(A, na.rm = TRUE), rowSums(Ar, na.rm = TRUE), bk$tol, tag("rowSums na.rm"))
    expect_conforms(colSums(A, na.rm = TRUE), colSums(Ar, na.rm = TRUE), bk$tol, tag("colSums na.rm"))
    expect_conforms(rowMeans(A, na.rm = TRUE), rowMeans(Ar, na.rm = TRUE), bk$tol, tag("rowMeans na.rm"))
    expect_conforms(colMeans(A, na.rm = TRUE), colMeans(Ar, na.rm = TRUE), bk$tol, tag("colMeans na.rm"))
  }
})

test_that("special values: sweep and am_sweep propagate like base R", {
  stat <- c(1, NA, Inf, 2)
  ref <- sweep(.sv_data$A, 2, stat, "-")
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(.sv_data$A, preferred_backend = bk$name)
    tag <- function(op) paste0("[", bk$name, "] ", op)
    expect_conforms(sweep(A, 2, stat, "-"), ref, bk$tol, tag("sweep -"))
    if (exists("am_sweep", mode = "function")) {
      expect_conforms(am_sweep(A, 2, stat, "-"), ref, bk$tol, tag("am_sweep -"))
    }
  }
})

test_that("special values: cbind/rbind preserve non-finite entries", {
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(.sv_data$A, preferred_backend = bk$name)
    B <- adgeMatrix(.sv_data$B, preferred_backend = bk$name)
    tag <- function(op) paste0("[", bk$name, "] ", op)
    expect_conforms(cbind(A, B), cbind(.sv_data$A, .sv_data$B), bk$tol, tag("cbind"))
    expect_conforms(rbind(A, B), rbind(.sv_data$A, .sv_data$B), bk$tol, tag("rbind"))
  }
})

test_that("special values: solve on a matrix containing Inf matches base R", {
  set.seed(202)
  S <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
  S[1, 1] <- Inf
  ref <- solve(S) # base R returns a (NaN-filled) matrix rather than erroring
  for (bk in .special_value_backends()) {
    A <- adgeMatrix(S, preferred_backend = bk$name)
    res <- tryCatch(solve(A), error = function(e) e)
    if (inherits(res, "error")) {
      # A backend may reject a non-finite/singular system; that is acceptable
      # divergence as long as it fails loudly rather than returning wrong values.
      succeed()
    } else {
      expect_conforms(res, ref, bk$tol, paste0("[", bk$name, "] solve(Inf)"))
    }
  }
})
