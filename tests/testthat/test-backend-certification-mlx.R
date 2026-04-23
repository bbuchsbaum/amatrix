# MLX beta certification smoke
#
# Seed: 20260423
# Dimensions: dense matmul 160x160; Cholesky solve on 96x96 SPD with 8 RHS
# Backend / precision / dispatch: mlx, fast precision, safe fresh Rscript -e worker
# R version / platform: captured by the child process sessionInfo()
# Issues: amatrix-x3c.2, amatrix-goq.2

.mlx_certification_repo_dir <- function() {
  candidates <- unique(c(
    tryCatch(getNamespaceInfo(asNamespace("amatrix"), "path"), error = function(e) NULL),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  ))
  candidates <- Filter(Negate(is.null), candidates)
  matches <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(matches) == 0L) {
    return(NULL)
  }
  normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("MLX beta certification: safe auto startup and claimed fast paths", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )

  repo_dir <- .mlx_certification_repo_dir()
  skip_if(is.null(repo_dir), "source tree not reachable (installed-pkg context)")

  result <- callr::r(
    func = function(repo_dir) {
      Sys.unsetenv("AMATRIX_MLX_PROBE_GPU")
      pkgload::load_all(repo_dir, quiet = TRUE)

      frob <- function(x) sqrt(sum(x * x))
      rel_frob <- function(x, ref) frob(x - ref) / max(frob(ref), .Machine$double.eps)

      status_before <- amatrix_backend_status("mlx")
      health <- amatrix_backend_health_probe("mlx")
      status_after <- amatrix_backend_status("mlx")

      set.seed(20260423L)
      n <- 160L
      x_host <- matrix(rnorm(n * n), n, n)
      y_host <- matrix(rnorm(n * n), n, n)
      x <- adgeMatrix(x_host, mode = "fast")
      y <- adgeMatrix(y_host, mode = "fast")
      matmul_plan <- amatrix_backend_plan(x, "matmul", y = y)

      amatrix_fallback_log_reset()
      matmul <- as.matrix(x %*% y)
      matmul_log <- amatrix_fallback_log()

      z <- matrix(rnorm(320L * 96L), nrow = 320L, ncol = 96L)
      spd <- crossprod(z) + diag(0.75, 96L)
      rhs <- matrix(rnorm(96L * 8L), nrow = 96L, ncol = 8L)
      spd_x <- adgeMatrix(spd, preferred_backend = "mlx", precision = "fast")
      chol_plan <- amatrix_backend_plan(spd_x, "chol")

      amatrix_fallback_log_reset()
      fac <- chol_factor(spd_x)
      sol <- chol_solve(fac, rhs)
      chol_log <- amatrix_fallback_log()
      ref_sol <- solve(spd, rhs)

      list(
        session = capture.output(sessionInfo()),
        startup_env = Sys.getenv("AMATRIX_MLX_PROBE_GPU", unset = ""),
        status_before = status_before,
        health = health,
        status_after = status_after,
        fast_preferred_backend = x@preferred_backend,
        fast_precision = x@precision,
        matmul_chosen = matmul_plan$chosen,
        matmul_rel = rel_frob(matmul, x_host %*% y_host),
        matmul_fallback_rows = nrow(matmul_log),
        chol_chosen = chol_plan$chosen,
        chol_backend = fac@backend,
        chol_precision = fac@precision,
        chol_solve_ref_rel = rel_frob(sol, ref_sol),
        chol_solve_resid_rel = rel_frob(spd %*% sol, rhs),
        chol_fallback_rows = nrow(chol_log)
      )
    },
    args = list(repo_dir = repo_dir)
  )

  expect_identical(result$status_before$available, TRUE)
  expect_identical(result$health$status, "healthy")
  expect_identical(result$status_after$health, "healthy")
  expect_identical(result$fast_preferred_backend, "mlx")
  expect_identical(result$fast_precision, "fast")

  expect_identical(result$matmul_chosen, "mlx")
  expect_lt(result$matmul_rel, 1e-4)
  expect_identical(result$matmul_fallback_rows, 0L)

  expect_identical(result$chol_chosen, "mlx")
  expect_identical(result$chol_backend, "mlx")
  expect_identical(result$chol_precision, "fast")
  expect_lt(result$chol_solve_ref_rel, 5e-6)
  expect_lt(result$chol_solve_resid_rel, 5e-6)
  expect_identical(result$chol_fallback_rows, 0L)
})
