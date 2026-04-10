test_that("amatrix_compile_product reuses a sparse lhs across repeated matmul calls", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "product_plan_backend",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      S_host <- Matrix::rsparsematrix(24, 10, density = 0.15)
      B1 <- matrix(rnorm(10 * 4), nrow = 10)
      B2 <- matrix(rnorm(10 * 4), nrow = 10)

      plan <- amatrix_compile_product(
        adgCMatrix(S_host, preferred_backend = "product_plan_backend", precision = "strict"),
        op = "matmul",
        backend = "auto"
      )

      out1 <- plan(B1)
      out2 <- plan(B2)
      out_direct <- plan(B1, materialize = "matrix")

      expect_s3_class(plan, "am_product_plan")
      expect_true(inherits(plan, "function"))
      expect_true(is.matrix(out_direct))
      expect_equal(as.matrix(out1), as.matrix(S_host %*% B1), tolerance = 1e-10)
      expect_equal(as.matrix(out2), as.matrix(S_host %*% B2), tolerance = 1e-10)
      expect_equal(out_direct, as.matrix(S_host %*% B1), tolerance = 1e-10)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 2L)
      expect_true(counter$spmm_resident >= 1L)
      expect_identical(get0("resident_materialize", envir = counter, inherits = FALSE, ifnotfound = 0L), 0L)
    }
  )
})

test_that("amatrix_compile_product aligns raw inputs to explicit fast-only backends", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "product_plan_fast_only",
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
      S_host <- Matrix::rsparsematrix(24, 10, density = 0.15)
      B <- matrix(rnorm(10 * 4), nrow = 10)

      plan <- amatrix_compile_product(
        S_host,
        op = "matmul",
        backend = "product_plan_fast_only"
      )

      out <- plan(B)

      expect_identical(attr(plan, "amatrix_plan_meta")$precision, "fast")
      expect_equal(as.matrix(out), as.matrix(S_host %*% B), tolerance = 1e-10)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 1L)

      expect_error(
        amatrix_compile_product(
          adgCMatrix(S_host, preferred_backend = "product_plan_fast_only", precision = "strict"),
          op = "matmul",
          backend = "product_plan_fast_only",
          precision = "strict"
        ),
        "does not support precision"
      )
    }
  )
})

test_that("matrix product plans use direct host fallback when resident predicate rejects a call", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(
    counter,
    supported_ops = character(),
    cold_supported_ops = character(),
    resident_supported_ops = character(),
    supports_sparse_ops = "matmul",
    supports_sparse_resident = TRUE
  )
  backend$supports_resident <- function(op, x, y = NULL) {
    !is.null(y) && NCOL(y) > 1L
  }

  with_registered_backend("product_plan_width_gate", backend, {
    S_host <- Matrix::rsparsematrix(24, 10, density = 0.15)
    b1 <- matrix(rnorm(10), nrow = 10, ncol = 1)
    b2 <- matrix(rnorm(10 * 2), nrow = 10, ncol = 2)

    plan <- amatrix_compile_product(
      adgCMatrix(S_host, preferred_backend = "product_plan_width_gate", precision = "strict"),
      op = "matmul",
      backend = "product_plan_width_gate"
    )

    out1 <- plan(b1, materialize = "matrix")
    out2 <- plan(b2, materialize = "matrix")

    expect_equal(out1, as.matrix(S_host %*% b1), tolerance = 1e-10)
    expect_equal(out2, as.matrix(S_host %*% b2), tolerance = 1e-10)
    expect_identical(counter$spmm_resident, 1L)
  })
})

test_that("amatrix_compile_product reuses a dense lhs for repeated sparse-right matmul", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(
    counter,
    supported_ops = "matmul",
    cold_supported_ops = character(),
    resident_supported_ops = character(),
    supports_sparse_ops = "matmul",
    supports_sparse_resident = TRUE
  )
  resident <- new.env(parent = emptyenv())
  sparse_resident <- new.env(parent = emptyenv())

  backend$resident_store <- function(key, x) {
    if (is.null(counter$resident_store)) counter$resident_store <- 0L
    counter$resident_store <- counter$resident_store + 1L
    assign(key, as.matrix(x), envir = resident)
    invisible(key)
  }
  backend$resident_has <- function(key) exists(key, envir = resident, inherits = FALSE)
  backend$resident_drop <- function(key) {
    if (exists(key, envir = resident, inherits = FALSE)) rm(list = key, envir = resident)
    invisible(key)
  }
  backend$resident_materialize <- function(key) get(key, envir = resident, inherits = FALSE)
  backend$sparse_resident_store <- function(key, x_sp) {
    if (is.null(counter$sparse_resident_store)) counter$sparse_resident_store <- 0L
    counter$sparse_resident_store <- counter$sparse_resident_store + 1L
    assign(key, methods::as(x_sp, "dgCMatrix"), envir = sparse_resident)
    invisible(key)
  }
  backend$sparse_resident_has <- function(key) exists(key, envir = sparse_resident, inherits = FALSE)
  backend$sparse_resident_drop <- function(key) {
    if (exists(key, envir = sparse_resident, inherits = FALSE)) rm(list = key, envir = sparse_resident)
    invisible(key)
  }
  backend$dense_sparse_matmul_resident_key <- function(x_key, sp_key, out_key, defer = FALSE) {
    if (is.null(counter$dense_sparse_matmul_resident_key)) counter$dense_sparse_matmul_resident_key <- 0L
    counter$dense_sparse_matmul_resident_key <- counter$dense_sparse_matmul_resident_key + 1L
    value <- as.matrix(
      get(x_key, envir = resident, inherits = FALSE) %*%
        get(sp_key, envir = sparse_resident, inherits = FALSE)
    )
    assign(out_key, value, envir = resident)
    if (isTRUE(defer)) NULL else value
  }

  with_registered_backend("product_plan_dense_sparse", backend, {
    A_host <- matrix(rnorm(6 * 20), nrow = 6)
    S1 <- Matrix::rsparsematrix(20, 8, density = 0.15)
    S2 <- Matrix::rsparsematrix(20, 8, density = 0.10)

    plan <- amatrix_compile_product(
      adgeMatrix(A_host, preferred_backend = "product_plan_dense_sparse", precision = "strict"),
      op = "matmul",
      backend = "auto"
    )

    out1 <- plan(S1)
    out2 <- plan(S2)

    expect_equal(as.matrix(out1), as.matrix(A_host %*% S1), tolerance = 1e-10)
    expect_equal(as.matrix(out2), as.matrix(A_host %*% S2), tolerance = 1e-10)
    expect_true(counter$resident_store >= 1L)
    expect_true(counter$dense_sparse_matmul_resident_key >= 2L)
    expect_identical(get0("resident_materialize", envir = counter, inherits = FALSE, ifnotfound = 0L), 0L)
  })
})
