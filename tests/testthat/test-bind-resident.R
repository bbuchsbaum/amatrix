test_that("amatrix_bind_resident binds dense matrices to backend residency", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "bind_dense_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = "matmul"
    ),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "bind_dense_backend")
      bound <- amatrix_bind_resident(x, "bind_dense_backend")

      expect_s4_class(bound, "adgeMatrix")
      expect_identical(amatrix:::.amatrix_live_resident_backend(bound), "bind_dense_backend")

      out <- bound %*% diag(2)
      expect_equal(as.matrix(out), as.matrix(bound) %*% diag(2))
      expect_true(counter$resident_store >= 1L)
      expect_true(counter$matmul_resident >= 1L)
    }
  )
})

test_that("amatrix_bind_resident binds sparse matrices for resident SpMM reuse", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "bind_sparse_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      S_host <- Matrix::rsparsematrix(12, 8, density = 0.2)
      B_host <- matrix(rnorm(8 * 3), nrow = 8)

      S <- adgCMatrix(S_host, preferred_backend = "bind_sparse_backend")
      bound <- amatrix_bind_resident(S, "bind_sparse_backend")

      expect_s4_class(bound, "adgCMatrix")
      expect_identical(amatrix:::.amatrix_live_resident_backend(bound), "bind_sparse_backend")

      out <- bound %*% B_host
      expect_equal(as.matrix(out), as.matrix(S_host %*% B_host), tolerance = 1e-10)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 1L)
      expect_false(exists("matmul", envir = counter, inherits = FALSE))
    }
  )
})

test_that("amatrix_bind_resident handles explicit fast-only backend precision", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "bind_fast_only_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE,
      precision_modes = "fast"
    ),
    {
      S_host <- Matrix::rsparsematrix(12, 8, density = 0.2)

      bound <- amatrix_bind_resident(S_host, "bind_fast_only_backend")

      expect_s4_class(bound, "adgCMatrix")
      expect_identical(bound@precision, "fast")
      expect_identical(amatrix:::.amatrix_live_resident_backend(bound), "bind_fast_only_backend")
      expect_true(counter$sparse_resident_store >= 1L)

      expect_error(
        amatrix_bind_resident(
          adgCMatrix(S_host, preferred_backend = "bind_fast_only_backend", precision = "strict"),
          "bind_fast_only_backend"
        ),
        "does not support precision"
      )
    }
  )
})

test_that("auto resident selection can choose and bind a sparse hot-path backend", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "bind_auto_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      S_host <- Matrix::rsparsematrix(20, 10, density = 0.15)
      B_host <- matrix(rnorm(10 * 4), nrow = 10)

      S <- adgCMatrix(
        S_host,
        preferred_backend = "bind_auto_backend",
        policy = "auto",
        precision = "strict"
      )

      expect_identical(
        amatrix_resident_backend_for(S, op = "matmul", y = B_host),
        "bind_auto_backend"
      )

      bound <- amatrix_bind_resident(S, backend = "auto", op = "matmul", y = B_host)

      expect_s4_class(bound, "adgCMatrix")
      expect_identical(amatrix:::.amatrix_live_resident_backend(bound), "bind_auto_backend")
      expect_equal(as.matrix(bound %*% B_host), as.matrix(S_host %*% B_host), tolerance = 1e-10)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 1L)
    }
  )
})

test_that("amatrix_prepare_operands auto-binds both operands for repeated sparse matmul", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "prepare_operands_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      S_host <- Matrix::rsparsematrix(18, 9, density = 0.15)
      B_host <- matrix(rnorm(9 * 5), nrow = 9)

      prep <- amatrix_prepare_operands(
        S_host,
        B_host,
        op = "matmul",
        backend = "auto",
        precision = "strict"
      )

      expect_identical(prep$backend, "prepare_operands_backend")
      expect_s4_class(prep$x, "adgCMatrix")
      expect_s4_class(prep$y, "adgeMatrix")
      expect_identical(amatrix:::.amatrix_live_resident_backend(prep$x), "prepare_operands_backend")
      expect_identical(amatrix:::.amatrix_live_resident_backend(prep$y), "prepare_operands_backend")
      expect_equal(as.matrix(prep$x %*% prep$y), as.matrix(S_host %*% B_host), tolerance = 1e-10)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 1L)
    }
  )
})
