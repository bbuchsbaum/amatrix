# test-bughunt-coercion.R — coercion round-trip bug hunt tests
# Tagged with beads issue IDs. DO NOT FIX — these are regression anchors.

library(Matrix)

# ── amatrix-chl: complex->adgeMatrix silently discards imaginary parts ────────
# BUG: adgeMatrix() accepts complex input and silently drops imaginary parts,
# emitting only a base-R storage.mode warning. amatrix should stop() with a
# clear message before data is corrupted.
test_that("amatrix-chl: complex input to adgeMatrix errors, not silent loss", {
  m_complex <- matrix(complex(real = 1:4, imaginary = 5:8), 2, 2)
  # Expect an explicit amatrix error — currently only a base warning is emitted
  expect_error(
    adgeMatrix(m_complex),
    regexp = "complex",
    info = "amatrix-chl: complex matrix must be rejected with a clear error"
  )
})

test_that("amatrix-chl: complex->adgeMatrix must not silently discard imaginary parts", {
  m_complex <- matrix(complex(real = 1:4, imaginary = 5:8), 2, 2)
  # If (incorrectly) no error is thrown, imaginary parts must not be silently lost
  result <- tryCatch(adgeMatrix(m_complex), error = function(e) NULL)
  if (!is.null(result)) {
    # Bug is present: coercion succeeded silently — assert imaginary parts survive
    mat_out <- as.matrix(result)
    expect_false(
      isTRUE(all.equal(Im(mat_out), matrix(0, 2, 2))),
      label = "amatrix-chl: imaginary parts must not be silently discarded"
    )
  }
})

# ── amatrix-ymy: setAs(dgCMatrix,adgCMatrix) loses backend/policy metadata ───
# BUG: as(adgCMatrix_obj, "dgCMatrix") strips backend metadata (expected —
# dgCMatrix has no such slots). But the round-trip
#   as(as(x, "dgCMatrix"), "adgCMatrix")
# resets preferred_backend and policy to package defaults instead of preserving
# the originals. Any user who round-trips through dgCMatrix loses their GPU
# backend choice silently.
test_that("amatrix-ymy: backend and policy survive adgCMatrix->dgCMatrix->adgCMatrix round-trip", {
  sp <- Matrix(c(1, 0, 0, 2), 2, 2, sparse = TRUE)
  A_sp <- adgCMatrix(sp, preferred_backend = "mlx", policy = "mlx")

  dg   <- as(A_sp, "dgCMatrix")
  A_sp2 <- as(dg, "adgCMatrix")

  expect_equal(
    A_sp2@preferred_backend, "mlx",
    info = "amatrix-ymy: preferred_backend must survive dgCMatrix round-trip"
  )
  expect_equal(
    A_sp2@policy, "mlx",
    info = "amatrix-ymy: policy must survive dgCMatrix round-trip"
  )
})

# ── amatrix-4sb: no setAs for adgCMatrix<->adgeMatrix cross-class coercion ───
# BUG: as(adgCMatrix_obj, "adgeMatrix") and as(adgeMatrix_obj, "adgCMatrix")
# both throw "no method or default for coercing" errors. These are the natural
# sparse<->dense coercions within the amatrix class family and should work.
test_that("amatrix-4sb: as(adgCMatrix, 'adgeMatrix') converts sparse to dense", {
  sp <- Matrix(c(1, 0, 0, 2, 3, 0), 3, 2, sparse = TRUE,
               dimnames = list(c("r1", "r2", "r3"), c("c1", "c2")))
  A_sp <- adgCMatrix(sp)

  # Currently errors: "no method or default for coercing adgCMatrix to adgeMatrix"
  expect_no_error(
    as(A_sp, "adgeMatrix"),
    message = "amatrix-4sb: as(adgCMatrix, 'adgeMatrix') must not error"
  )

  A_dense <- tryCatch(as(A_sp, "adgeMatrix"), error = function(e) NULL)
  if (!is.null(A_dense)) {
    expect_s4_class(A_dense, "adgeMatrix")
    expect_equal(
      as.matrix(A_dense), as.matrix(sp),
      info = "amatrix-4sb: values must match after sparse->dense coerce"
    )
    expect_equal(
      dimnames(as.matrix(A_dense)), dimnames(sp),
      info = "amatrix-4sb: dimnames must be preserved in sparse->dense coerce"
    )
  }
})

test_that("amatrix-4sb: as(adgeMatrix, 'adgCMatrix') converts dense to sparse", {
  m <- matrix(c(1, 0, 0, 2, 3, 0), 3, 2,
              dimnames = list(c("r1", "r2", "r3"), c("c1", "c2")))
  A <- adgeMatrix(m)

  # Currently errors: "no method or default for coercing adgeMatrix to adgCMatrix"
  expect_no_error(
    as(A, "adgCMatrix"),
    message = "amatrix-4sb: as(adgeMatrix, 'adgCMatrix') must not error"
  )

  A_sp <- tryCatch(as(A, "adgCMatrix"), error = function(e) NULL)
  if (!is.null(A_sp)) {
    expect_s4_class(A_sp, "adgCMatrix")
    expect_equal(
      as.matrix(A_sp), m,
      info = "amatrix-4sb: values must match after dense->sparse coerce"
    )
    expect_equal(
      dimnames(as.matrix(A_sp)), dimnames(m),
      info = "amatrix-4sb: dimnames must be preserved in dense->sparse coerce"
    )
  }
})
