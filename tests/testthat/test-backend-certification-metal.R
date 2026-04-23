# Metal beta certification smoke
#
# Seed: 20260423
# Dimensions: sparse 48x32 with 160 nonzeros; RHS blocks 32x7, 48x5, 9x32
# Backend / precision / dispatch: metal, fast precision, explicit GPU probe
# R version / platform: captured by the child process sessionInfo()
# Issues: amatrix-x3c.1

test_that("Metal beta certification: explicit sparse product paths", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("amatrix.metal")

  result <- callr::r(
    func = function(repo_dir) {
      Sys.setenv(AMATRIX_METAL_PROBE_GPU = "1")
      pkgload::load_all(repo_dir, quiet = TRUE)

      available <- amatrix.metal::amatrix_metal_enable_probe(register = TRUE)
      if (!isTRUE(available)) {
        return(list(available = FALSE, session = capture.output(sessionInfo())))
      }

      old_options <- options(
        amatrix.metal.spmm_min_nnz = 1L,
        amatrix.metal.spmv_min_nnz = 1L,
        amatrix.metal.resident_spmm_min_nnz = 1L,
        amatrix.metal.resident_spmv_min_nnz = 1L
      )
      on.exit(options(old_options), add = TRUE)

      health <- amatrix_backend_health_probe("metal")

      set.seed(20260423L)
      x_host <- matrix(0, nrow = 48L, ncol = 32L)
      idx <- sample(length(x_host), 160L)
      x_host[idx] <- rnorm(length(idx))
      x_sp <- Matrix::Matrix(x_host, sparse = TRUE)
      x <- adgCMatrix(x_sp, preferred_backend = "metal", precision = "fast")

      y_mm <- matrix(rnorm(32L * 7L), nrow = 32L, ncol = 7L)
      y_cp <- matrix(rnorm(48L * 5L), nrow = 48L, ncol = 5L)
      y_tcp <- matrix(rnorm(9L * 32L), nrow = 9L, ncol = 32L)

      plan_mm <- amatrix_backend_plan(x, "matmul", y = y_mm)
      plan_cp <- amatrix_backend_plan(x, "crossprod", y = y_cp)
      plan_tcp <- amatrix_backend_plan(x, "tcrossprod", y = y_tcp)

      amatrix_fallback_log_reset()
      out_mm <- as.matrix(x %*% y_mm)
      out_cp <- as.matrix(crossprod(x, y_cp))
      out_tcp <- as.matrix(tcrossprod(x, y_tcp))
      fallback <- amatrix_fallback_log()

      list(
        available = TRUE,
        session = capture.output(sessionInfo()),
        bridge = amatrix.metal::amatrix_metal_bridge_info(),
        health = health,
        plan_mm = plan_mm$chosen,
        plan_cp = plan_cp$chosen,
        plan_tcp = plan_tcp$chosen,
        path_mm = plan_mm$chosen_path,
        path_cp = plan_cp$chosen_path,
        path_tcp = plan_tcp$chosen_path,
        err_mm = max(abs(out_mm - as.matrix(x_sp %*% y_mm))),
        err_cp = max(abs(out_cp - as.matrix(Matrix::crossprod(x_sp, y_cp)))),
        err_tcp = max(abs(out_tcp - as.matrix(Matrix::tcrossprod(x_sp, y_tcp)))),
        fallback_rows = nrow(fallback)
      )
    },
    args = list(repo_dir = normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE))
  )

  skip_if_not(isTRUE(result$available), "metal native backend not available")

  expect_identical(result$bridge$available, TRUE)
  expect_identical(result$health$status, "healthy")
  expect_identical(result$plan_mm, "metal")
  expect_identical(result$plan_cp, "metal")
  expect_identical(result$plan_tcp, "metal")
  expect_true(result$path_mm %in% c("cold", "resident"))
  expect_true(result$path_cp %in% c("cold", "resident"))
  expect_true(result$path_tcp %in% c("cold", "resident"))
  expect_lt(result$err_mm, 1e-4)
  expect_lt(result$err_cp, 1e-4)
  expect_lt(result$err_tcp, 1e-4)
  expect_identical(result$fallback_rows, 0L)
})
