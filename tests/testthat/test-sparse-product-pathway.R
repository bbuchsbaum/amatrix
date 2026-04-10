test_that("sparse cold products use the shared backend pathway when a backend opts in", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_path",
    make_recording_backend(
      counter,
      supported_ops = c("matmul", "crossprod", "tcrossprod"),
      cold_supported_ops = c("matmul", "crossprod", "tcrossprod"),
      resident_supported_ops = character(),
      supports_sparse_ops = c("matmul", "crossprod", "tcrossprod"),
      supports_sparse_resident = FALSE
    ),
    {
      X_host <- Matrix::rsparsematrix(40, 24, density = 0.08)
      B_host <- matrix(rnorm(24 * 5), nrow = 24, ncol = 5)
      Y_cross <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
      Y_tcross <- matrix(rnorm(7 * 24), nrow = 7, ncol = 24)

      X <- new_adgCMatrix(X_host, preferred_backend = "recording_sparse_path", precision = "strict")

      plan_matmul <- amatrix_backend_plan(X, "matmul", y = B_host)
      plan_crossprod <- amatrix_backend_plan(X, "crossprod", y = Y_cross)
      plan_tcrossprod <- amatrix_backend_plan(X, "tcrossprod", y = Y_tcross)

      expect_identical(plan_matmul$chosen, "recording_sparse_path")
      expect_identical(plan_crossprod$chosen, "recording_sparse_path")
      expect_identical(plan_tcrossprod$chosen, "recording_sparse_path")

      expect_equal(as.matrix(X %*% B_host), as.matrix(X_host %*% B_host), tolerance = 1e-10)
      expect_equal(as.matrix(crossprod(X, Y_cross)), as.matrix(Matrix::crossprod(X_host, Y_cross)), tolerance = 1e-10)
      expect_equal(as.matrix(tcrossprod(X, Y_tcross)), as.matrix(Matrix::tcrossprod(X_host, Y_tcross)), tolerance = 1e-10)

      expect_identical(counter$matmul, 1L)
      expect_identical(counter$crossprod, 1L)
      expect_identical(counter$tcrossprod, 1L)
      expect_false(exists("spmm_resident", envir = counter, inherits = FALSE))
    }
  )
})

test_that("dense-left sparse-right matmul lowers through sparse crossprod", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_rhs_matmul",
    make_recording_backend(
      counter,
      supported_ops = "crossprod",
      cold_supported_ops = "crossprod",
      resident_supported_ops = character(),
      supports_sparse_ops = "crossprod",
      supports_sparse_resident = FALSE
    ),
    {
      A_host <- matrix(rnorm(7 * 24), nrow = 7, ncol = 24)
      S_host <- Matrix::rsparsematrix(24, 11, density = 0.08)

      A <- new_adgeMatrix(A_host, preferred_backend = "cpu", precision = "strict")
      S <- new_adgCMatrix(S_host, preferred_backend = "recording_sparse_rhs_matmul", precision = "strict")

      out <- A %*% S

      expect_true(inherits(out, "adgeMatrix"))
      expect_equal(as.matrix(out), as.matrix(A_host %*% S_host), tolerance = 1e-10)
      expect_identical(counter$crossprod, 1L)
      expect_false(exists("matmul", envir = counter, inherits = FALSE))
    }
  )
})

test_that("dense-left sparse-right matmul can use a direct resident backend path", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(
    counter,
    supported_ops = "matmul",
    cold_supported_ops = "matmul",
    resident_supported_ops = character(),
    supports_sparse_ops = "matmul",
    supports_sparse_resident = TRUE
  )
  resident <- new.env(parent = emptyenv())
  sparse_resident <- new.env(parent = emptyenv())

  backend$resident_store <- function(key, x) {
    if (is.null(counter$resident_store)) {
      counter$resident_store <- 0L
    }
    counter$resident_store <- counter$resident_store + 1L
    assign(key, as.matrix(x), envir = resident)
    invisible(key)
  }
  backend$resident_has <- function(key) exists(key, envir = resident, inherits = FALSE)
  backend$resident_drop <- function(key) {
    if (exists(key, envir = resident, inherits = FALSE)) {
      rm(list = key, envir = resident)
    }
    invisible(key)
  }
  backend$resident_materialize <- function(key) get(key, envir = resident, inherits = FALSE)
  backend$sparse_resident_store <- function(key, x_sp) {
    if (is.null(counter$sparse_resident_store)) {
      counter$sparse_resident_store <- 0L
    }
    counter$sparse_resident_store <- counter$sparse_resident_store + 1L
    assign(key, methods::as(x_sp, "dgCMatrix"), envir = sparse_resident)
    invisible(key)
  }
  backend$sparse_resident_has <- function(key) exists(key, envir = sparse_resident, inherits = FALSE)
  backend$sparse_resident_drop <- function(key) {
    if (exists(key, envir = sparse_resident, inherits = FALSE)) {
      rm(list = key, envir = sparse_resident)
    }
    invisible(key)
  }
  backend$dense_sparse_matmul_resident_key <- function(x_key, sp_key, out_key, defer = FALSE) {
    if (is.null(counter$dense_sparse_matmul_resident_key)) {
      counter$dense_sparse_matmul_resident_key <- 0L
    }
    counter$dense_sparse_matmul_resident_key <- counter$dense_sparse_matmul_resident_key + 1L
    value <- as.matrix(
      get(x_key, envir = resident, inherits = FALSE) %*%
        get(sp_key, envir = sparse_resident, inherits = FALSE)
    )
    assign(out_key, value, envir = resident)
    if (isTRUE(defer)) NULL else value
  }

  with_registered_backend(
    "recording_dense_sparse_direct",
    backend,
    {
      A_host <- matrix(rnorm(8 * 40), nrow = 8, ncol = 40)
      S_host <- Matrix::rsparsematrix(40, 16, density = 0.05)

      A <- new_adgeMatrix(A_host, preferred_backend = "recording_dense_sparse_direct", precision = "strict")
      S <- new_adgCMatrix(S_host, preferred_backend = "recording_dense_sparse_direct", precision = "strict")

      out1 <- A %*% S
      out2 <- A %*% S

      expect_equal(as.matrix(out1), as.matrix(A_host %*% S_host), tolerance = 1e-10)
      expect_equal(as.matrix(out2), as.matrix(A_host %*% S_host), tolerance = 1e-10)
      expect_identical(counter$dense_sparse_matmul_resident_key, 2L)
      expect_false(exists("crossprod", envir = counter, inherits = FALSE))
      expect_false(exists("spmm_resident_key", envir = counter, inherits = FALSE))
    }
  )
})

test_that("dense-left sparse-right crossprod lowers through sparse crossprod", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_rhs_crossprod",
    make_recording_backend(
      counter,
      supported_ops = "crossprod",
      cold_supported_ops = "crossprod",
      resident_supported_ops = character(),
      supports_sparse_ops = "crossprod",
      supports_sparse_resident = FALSE
    ),
    {
      A_host <- matrix(rnorm(24 * 6), nrow = 24, ncol = 6)
      S_host <- Matrix::rsparsematrix(24, 9, density = 0.07)

      A <- new_adgeMatrix(A_host, preferred_backend = "cpu", precision = "strict")
      S <- new_adgCMatrix(S_host, preferred_backend = "recording_sparse_rhs_crossprod", precision = "strict")

      out <- crossprod(A, S)

      expect_true(inherits(out, "adgeMatrix"))
      expect_equal(as.matrix(out), as.matrix(base::crossprod(A_host, S_host)), tolerance = 1e-10)
      expect_identical(counter$crossprod, 1L)
      expect_false(exists("matmul", envir = counter, inherits = FALSE))
    }
  )
})

test_that("dense-left sparse-right tcrossprod lowers through sparse tcrossprod", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_rhs_tcrossprod",
    make_recording_backend(
      counter,
      supported_ops = "tcrossprod",
      cold_supported_ops = "tcrossprod",
      resident_supported_ops = character(),
      supports_sparse_ops = "tcrossprod",
      supports_sparse_resident = FALSE
    ),
    {
      A_host <- matrix(rnorm(5 * 18), nrow = 5, ncol = 18)
      S_host <- Matrix::rsparsematrix(8, 18, density = 0.09)

      A <- new_adgeMatrix(A_host, preferred_backend = "cpu", precision = "strict")
      S <- new_adgCMatrix(S_host, preferred_backend = "recording_sparse_rhs_tcrossprod", precision = "strict")

      out <- tcrossprod(A, S)

      expect_true(inherits(out, "adgeMatrix"))
      expect_equal(as.matrix(out), as.matrix(base::tcrossprod(A_host, S_host)), tolerance = 1e-10)
      expect_identical(counter$tcrossprod, 1L)
      expect_false(exists("matmul", envir = counter, inherits = FALSE))
    }
  )
})

test_that("sparse resident SpMM reuses one stored sparse operand across repeated products", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_resident",
    make_recording_backend(
      counter,
      supported_ops = "matmul",
      cold_supported_ops = "matmul",
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      X_host <- Matrix::rsparsematrix(64, 32, density = 0.05)
      B_host <- matrix(rnorm(32 * 6), nrow = 32, ncol = 6)

      X <- new_adgCMatrix(X_host, preferred_backend = "recording_sparse_resident", precision = "strict")

      out1 <- X %*% B_host
      out2 <- X %*% B_host

      expect_equal(as.matrix(out1), as.matrix(X_host %*% B_host), tolerance = 1e-10)
      expect_equal(as.matrix(out2), as.matrix(X_host %*% B_host), tolerance = 1e-10)

      expect_identical(counter$sparse_resident_store, 1L)
      expect_identical(counter$spmm_resident_key, 2L)
      expect_false(exists("matmul", envir = counter, inherits = FALSE))
    }
  )
})

test_that("sparse resident SpMM can use a resident dense RHS key path", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_resident_key",
    make_recording_backend(
      counter,
      supported_ops = "matmul",
      cold_supported_ops = "matmul",
      resident_supported_ops = character(),
      supports_sparse_ops = "matmul",
      supports_sparse_resident = TRUE
    ),
    {
      X_host <- Matrix::rsparsematrix(64, 32, density = 0.05)
      B_host <- matrix(rnorm(32 * 6), nrow = 32, ncol = 6)

      X <- new_adgCMatrix(X_host, preferred_backend = "recording_sparse_resident_key", precision = "strict")
      B <- new_adgeMatrix(B_host, preferred_backend = "recording_sparse_resident_key", precision = "strict")

      out1 <- X %*% B
      out2 <- X %*% B

      expect_equal(as.matrix(out1), as.matrix(X_host %*% B_host), tolerance = 1e-10)
      expect_equal(as.matrix(out2), as.matrix(X_host %*% B_host), tolerance = 1e-10)

      expect_identical(counter$sparse_resident_store, 1L)
      expect_identical(counter$resident_store, 1L)
      expect_identical(counter$spmm_resident_key, 2L)
      expect_false(exists("spmm_resident", envir = counter, inherits = FALSE))
    }
  )
})

test_that("sparse resident crossprod and tcrossprod route through transposed SpMM", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "recording_sparse_crossprod",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod", "tcrossprod"),
      cold_supported_ops = c("crossprod", "tcrossprod"),
      resident_supported_ops = character(),
      supports_sparse_ops = c("crossprod", "tcrossprod"),
      supports_sparse_resident = TRUE
    ),
    {
      X_host <- Matrix::rsparsematrix(50, 20, density = 0.07)
      Y_cross <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
      Y_tcross <- matrix(rnorm(7 * 20), nrow = 7, ncol = 20)

      X <- new_adgCMatrix(X_host, preferred_backend = "recording_sparse_crossprod", precision = "strict")

      out_cross <- crossprod(X, Y_cross)
      out_tcross <- tcrossprod(X, Y_tcross)

      expect_equal(as.matrix(out_cross), as.matrix(Matrix::crossprod(X_host, Y_cross)), tolerance = 1e-10)
      expect_equal(as.matrix(out_tcross), as.matrix(Matrix::tcrossprod(X_host, Y_tcross)), tolerance = 1e-10)

      expect_identical(counter$sparse_resident_store, 1L)
      expect_identical(counter$spmm_resident_key, 2L)
      expect_false(exists("spmm_resident", envir = counter, inherits = FALSE))
      expect_false(exists("crossprod", envir = counter, inherits = FALSE))
      expect_false(exists("tcrossprod", envir = counter, inherits = FALSE))
    }
  )
})
