# Regression tests for amatrix-lc1: the deferred-residency path fills the host
# @x slot with a rep(NaN, n) sentinel while the authoritative data lives on the
# device (or in the cached host_x). Genuine user NaN -- and NA / Inf -- payloads
# must never be confused with that sentinel. Detection of "not materialized" is
# gated purely by the host_deferred flag; no consumer may sniff @x values or
# leak the sentinel through a coercion.
#
# On HEAD before the fix, `as(A, "dgeMatrix")` and the CPU host accessor
# (.amatrix_dense_slot_matrix, used by .amatrix_cpu_dense_matrix) read a deferred
# object's @x directly and returned an all-NaN matrix.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# genuine NaN, NA_real_, +/-Inf and finite values mixed together
.nan_sentinel_payload <- function() {
  matrix(
    c(1, NaN, 3, NA_real_, Inf, -Inf, 7, 8, 9),
    nrow = 3, ncol = 3
  )
}

# The internal CPU host accessor stamps dimnames = list(NULL, NULL); drop them so
# value-level (bit-exact NA-vs-NaN) comparisons are not tripped by that attribute.
.strip_dimnames <- function(m) {
  dimnames(m) <- NULL
  m
}

# Build a genuinely deferred adgeMatrix backed by a live (mock or real) resident
# buffer holding `real`. Works on CPU via a recording backend, so no GPU needed.
.build_deferred_mock <- function(real, backend_name) {
  backend <- amatrix:::.amatrix_get_backend(backend_name)
  key <- amatrix:::.amatrix_next_resident_key(backend_name)
  backend$resident_store(key, real)
  A <- amatrix:::new_adgeMatrix_deferred(
    dim = dim(real),
    preferred_backend = backend_name
  )
  amatrix:::.amatrix_bind_resident(A, backend_name, key)
}

# ---------------------------------------------------------------------------
# CPU: non-deferred path (bit-exact, preserves the NA-vs-NaN distinction)
# ---------------------------------------------------------------------------

test_that("non-deferred adgeMatrix preserves genuine NaN/NA/Inf bit-exactly (CPU)", {
  skip_if_not_installed("Matrix")
  real <- .nan_sentinel_payload()
  x <- adgeMatrix(real, preferred_backend = "cpu")

  expect_identical(as.matrix(x), real)
  expect_identical(as.matrix(as(x, "dgeMatrix")), real)
  expect_identical(as.matrix(amatrix_materialize_host(x)), real)
  expect_identical(.strip_dimnames(amatrix:::.amatrix_cpu_dense_matrix(x)), real)
})

test_that("non-deferred adgeMatrix NaN/NA/Inf survive arithmetic and RDS (CPU)", {
  skip_if_not_installed("Matrix")
  real <- .nan_sentinel_payload()
  x <- adgeMatrix(real, preferred_backend = "cpu")

  expect_identical(as.matrix(x + 0), real + 0)
  expect_identical(as.matrix(x * 1), real * 1)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  x2 <- readRDS(tmp)
  expect_identical(as.matrix(x2), real)
})

# ---------------------------------------------------------------------------
# CPU: deferred path via a mock residency backend (the core regression)
# ---------------------------------------------------------------------------

test_that("deferred sentinel never leaks: genuine NaN/NA/Inf survive coercion (CPU mock)", {
  skip_if_not_installed("Matrix")
  counter <- new.env(parent = emptyenv())
  with_registered_backend(
    "nan_sentinel_backend",
    make_recording_backend(counter, resident_supported_ops = "matmul"),
    {
      real <- .nan_sentinel_payload()
      A <- .build_deferred_mock(real, "nan_sentinel_backend")

      expect_true(isTRUE(A@finalizer_env$host_deferred))
      # @x is the non-authoritative sentinel, not the real data.
      expect_true(all(is.nan(A@x)))

      # Every host-facing path must return authoritative data, not the sentinel.
      expect_identical(as.matrix(amatrix_materialize_host(A)), real)
      expect_identical(as.matrix(as(A, "dgeMatrix")), real)
      expect_identical(.strip_dimnames(amatrix:::.amatrix_cpu_dense_matrix(A)), real)
      expect_identical(as.matrix(A), real)

      # Pre-fix, every path above returned an all-NaN matrix.
      expect_false(all(is.nan(as.matrix(A))))
    }
  )
})

test_that("deferred adgeMatrix arithmetic preserves genuine NaN/NA/Inf (CPU mock)", {
  skip_if_not_installed("Matrix")
  counter <- new.env(parent = emptyenv())
  with_registered_backend(
    "nan_sentinel_arith_backend",
    make_recording_backend(counter, resident_supported_ops = "matmul"),
    {
      real <- .nan_sentinel_payload()
      A <- .build_deferred_mock(real, "nan_sentinel_arith_backend")
      expect_identical(as.matrix(A + 0), real + 0)
    }
  )
})

test_that("materialized deferred survives saveRDS/readRDS with genuine NaN (CPU mock)", {
  skip_if_not_installed("Matrix")
  counter <- new.env(parent = emptyenv())
  real <- .nan_sentinel_payload()
  A <- with_registered_backend(
    "nan_sentinel_rds_backend",
    make_recording_backend(counter, resident_supported_ops = "matmul"),
    {
      A <- .build_deferred_mock(real, "nan_sentinel_rds_backend")
      # Materialize so host_x caches the real data (incl. genuine NaN) before
      # the resident buffer/backend go away.
      force(as.matrix(A))
      A
    }
  )

  # Backend unregistered and resident buffer gone, but host_x was cached, so the
  # object is not "dead deferred" and must round-trip its genuine NaN payload.
  expect_false(amatrix:::.amatrix_is_dead_deferred(A))

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(A, tmp)
  A2 <- readRDS(tmp)
  expect_identical(as.matrix(A2), real)
})

# ---------------------------------------------------------------------------
# MLX (GPU): gated. float32 collapses NA_real_ -> NaN, so assert missingness /
# infinity positions and finite values rather than bit-exact NA-vs-NaN.
# ---------------------------------------------------------------------------

test_that("deferred sentinel never leaks on MLX: genuine NaN survives (MLX)", {
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )

  real <- matrix(c(1, NaN, 3, 4, NaN, 6), nrow = 2, ncol = 3)
  h <- resident_handle(real, backend = "mlx")
  A <- amatrix:::as_adgeMatrix.resident_handle(h, defer_host = TRUE)

  expect_true(isTRUE(A@finalizer_env$host_deferred))
  expect_true(all(is.nan(A@x)))

  m_host <- as.matrix(amatrix_materialize_host(A))
  m_dge <- as.matrix(as(A, "dgeMatrix"))

  # No sentinel leak: not all-NaN, coercion paths agree, genuine NaN preserved.
  expect_false(all(is.nan(m_dge)))
  expect_identical(is.nan(m_host), is.nan(real))
  expect_identical(is.nan(m_dge), is.nan(real))
  expect_equal(m_dge, m_host)
  expect_equal(m_host[is.finite(real)], real[is.finite(real)])
})

test_that("deferred sentinel never leaks on MLX: mixed NA/NaN/Inf positions survive (MLX)", {
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )

  real <- matrix(c(1, NA_real_, NaN, Inf, -Inf, 6), nrow = 2, ncol = 3)
  h <- resident_handle(real, backend = "mlx")
  A <- amatrix:::as_adgeMatrix.resident_handle(h, defer_host = TRUE)

  m <- as.matrix(as(A, "dgeMatrix"))

  expect_false(all(is.nan(m)))
  # NA/NaN both collapse to "missing" on a float GPU; positions must survive.
  expect_identical(is.na(m), is.na(real))
  expect_identical(is.infinite(m), is.infinite(real))
  expect_identical(sign(m[is.infinite(real)]), sign(real[is.infinite(real)]))
  expect_equal(m[is.finite(real)], real[is.finite(real)])
})
