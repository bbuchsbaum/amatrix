# Regressions for the fourth batch of pre-CRAN review fixes.

test_that("dist_matrix cosine is a distance, not a similarity", {
  set.seed(1)
  X <- matrix(rnorm(20), 5)
  nx <- sqrt(rowSums(X^2))
  ref <- 1 - tcrossprod(X) / outer(nx, nx)
  diag(ref) <- 0

  D <- as.matrix(am_dist_matrix(X, method = "cosine"))
  expect_equal(D, ref, ignore_attr = TRUE, tolerance = 1e-10)
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0 & D <= 2))

  Dxy <- as.matrix(am_dist_matrix(adgeMatrix(X), X[1:2, ], method = "cosine"))
  expect_equal(Dxy, 1 - tcrossprod(X, X[1:2, ]) / outer(nx, nx[1:2]),
               ignore_attr = TRUE, tolerance = 1e-10)
  # The similarity remains available as a kernel.
  expect_equal(as.matrix(am_kernel_matrix(X, kernel = "cosine")), 1 - ref,
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("products with a plain matrix operand match base above the in-place threshold", {
  set.seed(2)
  n <- 520L
  M <- matrix(rnorm(n * n), n, dimnames = list(paste0("r", seq_len(n)), NULL))
  A <- adgeMatrix(M)
  expect_s4_class(A %*% M, "adgeMatrix")
  expect_equal(as.matrix(A %*% M), M %*% M)
  expect_identical(dimnames(A %*% M), dimnames(M %*% M))
  expect_equal(as.matrix(M %*% A), M %*% M)
  expect_equal(as.matrix(crossprod(A, M)), crossprod(M), ignore_attr = TRUE)
})

test_that("the fallback log is bounded", {
  state <- amatrix:::.amatrix_state
  old_log <- state$fallback_log
  on.exit(state$fallback_log <- old_log, add = TRUE)
  amatrix_fallback_log_reset()
  for (i in seq_len(1100L)) {
    amatrix:::.amatrix_log_fallback(op = "matmul", backend = "fake", reason = as.character(i))
  }
  log <- amatrix_fallback_log()
  expect_lte(nrow(log), 1000L)
  expect_identical(utils::tail(log$reason, 1L), "1100")
})

test_that("amatrix_compile_product leaves the defer_host option as it found it", {
  withr::local_options(amatrix.defer_host = NULL)
  x <- adgeMatrix(matrix(rnorm(20), 5, 4))
  plan <- amatrix_compile_product(x, op = "matmul")
  plan(matrix(rnorm(8), 4, 2))
  expect_null(getOption("amatrix.defer_host"))
})

test_that("sweep() dispatches on amatrix objects and keeps sparse scaling sparse", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  hv <- function(x) if (methods::is(x, "Matrix")) as.matrix(x) else x
  for (margin in 1:2) {
    stats <- if (margin == 1L) c(2, 4) else c(1, 2, 4)
    for (fun in c("-", "+", "*", "/")) {
      ref <- base::sweep(m, margin, stats, fun)
      expect_equal(hv(sweep(A, margin, stats, fun)), ref, ignore_attr = TRUE, info = paste(margin, fun))
      expect_equal(hv(sweep(S, margin, stats, fun)), ref, ignore_attr = TRUE, info = paste(margin, fun))
    }
  }
  # Default operation is "-", as in base, for both sweep() and am_sweep().
  expect_equal(hv(sweep(A, 2, c(1, 2, 4))), base::sweep(m, 2, c(1, 2, 4)), ignore_attr = TRUE)
  expect_equal(hv(am_sweep(A, 2, c(1, 2, 4))), base::sweep(m, 2, c(1, 2, 4)), ignore_attr = TRUE)
  expect_s4_class(sweep(A, 2, c(1, 2, 4)), "adgeMatrix")
  expect_s4_class(sweep(S, 2, c(1, 2, 4), "*"), "adgCMatrix")
  expect_s4_class(am_colscale(S, c(1, 2, 4)), "adgCMatrix")
  expect_s4_class(am_rowscale(S, c(2, 4)), "adgCMatrix")
  expect_equal(hv(am_colscale(S, c(1, 2, 4))), t(t(m) * c(1, 2, 4)), ignore_attr = TRUE)
  expect_equal(hv(am_rowscale(S, c(2, 4))), m * c(2, 4), ignore_attr = TRUE)
  # Non-amatrix inputs are untouched.
  expect_identical(sweep(m, 2, c(1, 2, 4)), base::sweep(m, 2, c(1, 2, 4)))
})

test_that("mean, which and isSymmetric work on amatrix objects", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  expect_identical(mean(A), mean(m))
  expect_equal(mean(S), mean(m))
  expect_equal(mean(A > 2), mean(m > 2))
  expect_identical(as.integer(which(A > 2)), which(m > 2))
  expect_identical(as.integer(which(S > 2)), which(m > 2))
  sq <- matrix(c(2, 1, 1, 3), 2)
  expect_true(isSymmetric(adgeMatrix(sq)))
  expect_true(isSymmetric(adgCMatrix(Matrix::Matrix(sq, sparse = TRUE))))
  expect_false(isSymmetric(adgeMatrix(matrix(1:4 + 0, 2))))
})

test_that("mean() and which() work with only amatrix attached", {
  skip_if_not_installed("callr")
  out <- callr::r(function() {
    library(amatrix)
    A <- adgeMatrix(matrix(c(1, 0, 3, 0, 5, 6), 2, 3))
    list(
      matrix_attached = "package:Matrix" %in% search(),
      mean = mean(A),
      which = as.integer(which(A > 2)),
      sym = isSymmetric(adgeMatrix(diag(2)))
    )
  })
  expect_false(out$matrix_attached)
  expect_identical(out$mean, 2.5)
  expect_identical(out$which, c(3L, 5L, 6L))
  expect_true(out$sym)
})

test_that("dim<- reshapes to a valid object with a new identity", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3, dimnames = list(c("a", "b"), NULL))
  ref <- m
  dim(ref) <- c(3L, 2L)
  A <- adgeMatrix(m)
  d <- A
  dim(d) <- c(3L, 2L)
  expect_true(methods::validObject(d))
  expect_equal(as.matrix(d), ref, ignore_attr = TRUE)
  expect_null(rownames(d))
  expect_false(identical(d@object_id, A@object_id))

  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  ds <- S
  dim(ds) <- c(3L, 2L)
  expect_s4_class(ds, "adgCMatrix")
  expect_equal(as.matrix(ds), ref, ignore_attr = TRUE)
  expect_error({ z <- A; dim(z) <- c(4L, 2L) }, "do not match")
})

test_that("uncalibrated accelerators keep small and transfer-bound cold work on the CPU", {
  state <- amatrix:::.amatrix_state
  old_calibration <- state$calibration
  old_key <- state$calibration_absent_key
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  withr::local_options(amatrix.uncalibrated_min_elements = NULL)
  on.exit({
    state$calibration <- old_calibration
    state$calibration_absent_key <- old_key
  }, add = TRUE)
  state$calibration <- NULL

  ok <- amatrix:::.amatrix_calibration_ok
  small <- adgeMatrix(matrix(0, 16, 16))
  big <- adgeMatrix(matrix(0, 256, 256))
  for (backend in c("mlx", "opencl", "arrayfire")) {
    expect_false(ok(small, "matmul", backend, y = small), info = backend)
    expect_true(ok(big, "matmul", backend, y = big), info = backend)
    expect_true(ok(big, "crossprod", backend), info = backend)
    # Memory-bound ops never pay for a round trip, whatever the size.
    expect_false(ok(big, "ewise", backend), info = backend)
    expect_false(ok(big, "rowSums", backend), info = backend)
  }
  # CPU, custom and test backends are not gated.
  expect_true(ok(small, "matmul", "cpu", y = small))
  expect_true(ok(small, "ewise", "amatrix_test"))
  # The rule can be switched off.
  withr::local_options(amatrix.uncalibrated_min_elements = 0)
  expect_true(ok(small, "matmul", "mlx", y = small))
  expect_true(ok(big, "ewise", "mlx"))
})

test_that("an explicit calibration overrides the uncalibrated rule", {
  state <- amatrix:::.amatrix_state
  old_calibration <- state$calibration
  withr::local_options(amatrix.uncalibrated_min_elements = NULL)
  on.exit(state$calibration <- old_calibration, add = TRUE)
  state$calibration <- list(version = "2", thresholds = list(mlx = list(gemm = 10, ewise = 10)))
  small <- adgeMatrix(matrix(0, 16, 16))
  expect_true(amatrix:::.amatrix_calibration_ok(small, "matmul", "mlx", y = small))
  expect_true(amatrix:::.amatrix_calibration_ok(small, "ewise", "mlx"))
})

test_that("deferred resident results are treated as success, not failure", {
  # A backend whose resident kernels return NULL when defer = TRUE (the result
  # stays on the device). The dense product must use the resident path, give a
  # deferred object, and leave the backend healthy.
  store <- new.env(parent = emptyenv())
  calls <- new.env(parent = emptyenv())
  cpu <- amatrix:::.amatrix_cpu_backend()
  backend <- amatrix_test_backend()
  backend$resident_has <- function(key) exists(key, envir = store, inherits = FALSE)
  backend$resident_store <- function(key, x) { assign(key, as.matrix(x), envir = store); invisible(key) }
  backend$resident_drop <- function(key) { if (exists(key, envir = store, inherits = FALSE)) rm(list = key, envir = store); invisible(TRUE) }
  backend$resident_materialize <- function(key) get(key, envir = store, inherits = FALSE)
  backend$matmul_resident <- function(x_key, y_key, out_key, defer = FALSE) {
    calls$matmul_resident <- (if (is.null(calls$matmul_resident)) 0L else calls$matmul_resident) + 1L
    value <- get(x_key, envir = store) %*% get(y_key, envir = store)
    assign(out_key, value, envir = store)
    if (isTRUE(defer)) NULL else value
  }
  backend$crossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
    xv <- get(x_key, envir = store)
    value <- if (is.null(y_key)) crossprod(xv) else crossprod(xv, get(y_key, envir = store))
    assign(out_key, value, envir = store)
    if (isTRUE(defer)) NULL else value
  }

  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  b <- matrix(c(1, 0, 2, 1, 0, 1, 3, 1, 1), 3, 3)
  withr::local_options(amatrix.defer_host = TRUE)
  out <- amatrix_with_test_backend({
    x <- adgeMatrix(m, backend = "amatrix_test")
    y <- adgeMatrix(b, backend = "amatrix_test")
    z <- x %*% y
    w <- crossprod(z)
    list(
      deferred = isTRUE(z@finalizer_env$host_deferred),
      prod = as.matrix(z),
      cp = as.matrix(w),
      health = amatrix:::.amatrix_backend_health_get("amatrix_test")$status,
      fallbacks = nrow(amatrix_fallback_log()[amatrix_fallback_log()$from_backend == "amatrix_test", ])
    )
  }, backend = backend)

  expect_true(out$deferred)
  expect_equal(out$prod, m %*% b, ignore_attr = TRUE)
  expect_equal(out$cp, crossprod(m %*% b), ignore_attr = TRUE)
  expect_identical(calls$matmul_resident, 1L)
  expect_identical(out$health, "healthy")
  expect_identical(out$fallbacks, 0L)
})

test_that("products keep base dimnames on backend and deferred paths", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  b <- matrix(1:9 + 0, 3, 3, dimnames = list(NULL, c("x", "y", "z")))

  # CPU reference behaviour.
  A <- adgeMatrix(m)
  expect_identical(dimnames(A %*% b), dimnames(m %*% b))
  expect_identical(dimnames(crossprod(A)), dimnames(crossprod(m)))
  expect_identical(dimnames(tcrossprod(A)), dimnames(tcrossprod(m)))
  expect_identical(dimnames(crossprod(A, m)), dimnames(crossprod(m, m)))

  # A backend that returns bare (unnamed) results, as device kernels do.
  backend <- amatrix_test_backend()
  strip <- function(f) { force(f); function(...) unname(as.matrix(f(...))) }
  for (nm in c("matmul", "crossprod", "tcrossprod")) backend[[nm]] <- strip(backend[[nm]])
  out <- amatrix_with_test_backend({
    X <- adgeMatrix(m, backend = "amatrix_test")
    list(mm = dimnames(X %*% b), cp = dimnames(crossprod(X)), tc = dimnames(tcrossprod(X)))
  }, backend = backend)
  expect_identical(out$mm, dimnames(m %*% b))
  expect_identical(out$cp, dimnames(crossprod(m)))
  expect_identical(out$tc, dimnames(tcrossprod(m)))
})
