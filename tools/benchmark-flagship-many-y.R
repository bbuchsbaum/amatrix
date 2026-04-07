#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
})

if (!requireNamespace("bench", quietly = TRUE)) {
  message("bench package not installed — skipping benchmark")
  quit(status = 0)
}
library(bench)

run_scenario <- function(n, p, k, label) {
  cat(sprintf("\n=== Scenario: %s  (n=%d, p=%d, k=%d) ===\n", label, n, p, k))

  set.seed(42)
  X_host <- matrix(rnorm(n * p), n, p)
  Y      <- matrix(rnorm(n * k), n, k)

  X_cpu <- adgeMatrix(X_host, preferred_backend = "cpu")

  # Correctness check before timing
  coef_lm  <- do.call(cbind, lapply(seq_len(k), function(j) lm.fit(X_host, Y[, j])$coefficients))
  fit_cpu  <- many_lm(X_cpu, Y, method = "qr", cache = TRUE)
  coef_am  <- as.matrix(fit_cpu$coefficients)
  max_diff <- max(abs(coef_am - coef_lm))
  if (max_diff >= 1e-8) {
    stop(sprintf("Correctness FAILED: max|coef_amatrix - coef_lm| = %g", max_diff))
  }
  cat(sprintf("  Correctness OK: max|coef diff| = %.2e\n", max_diff))

  # Check MLX availability
  mlx_available <- tryCatch({
    X_mlx_test <- adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast")
    many_lm(X_mlx_test, Y, method = "qr", cache = TRUE)
    TRUE
  }, error = function(e) FALSE)

  if (mlx_available) {
    bm <- bench::mark(
      lm_loop = {
        coefs <- vector("list", k)
        for (j in seq_len(k)) coefs[[j]] <- lm.fit(X_host, Y[, j])$coefficients
        invisible(coefs)
      },
      amatrix_cpu = invisible(many_lm(X_cpu, Y, method = "qr", cache = TRUE)),
      amatrix_mlx = invisible(many_lm(
        adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast"),
        Y, method = "qr", cache = TRUE
      )),
      iterations = 3,
      check      = FALSE,
      memory     = FALSE,
      time_unit  = "s"
    )
  } else {
    bm <- bench::mark(
      lm_loop = {
        coefs <- vector("list", k)
        for (j in seq_len(k)) coefs[[j]] <- lm.fit(X_host, Y[, j])$coefficients
        invisible(coefs)
      },
      amatrix_cpu = invisible(many_lm(X_cpu, Y, method = "qr", cache = TRUE)),
      iterations = 3,
      check      = FALSE,
      memory     = FALSE,
      time_unit  = "s"
    )
  }

  expr_labels <- as.character(bm$expression)
  median_lm   <- as.numeric(bm$median[expr_labels == "lm_loop"])
  median_cpu  <- as.numeric(bm$median[expr_labels == "amatrix_cpu"])
  speedup_cpu <- median_lm / median_cpu

  cat(sprintf("  lm_loop      median: %.4f s\n", median_lm))
  cat(sprintf("  amatrix_cpu  median: %.4f s  (%.1fx speedup vs lm_loop)\n",
              median_cpu, speedup_cpu))

  if (mlx_available) {
    median_mlx  <- as.numeric(bm$median[expr_labels == "amatrix_mlx"])
    speedup_mlx <- median_lm / median_mlx
    cat(sprintf("  amatrix_mlx  median: %.4f s  (%.1fx speedup vs lm_loop)\n",
                median_mlx, speedup_mlx))
  } else {
    cat("  amatrix_mlx  [skipped — MLX backend unavailable]\n")
  }
  flush.console()

  bm$scenario <- label
  bm
}

results <- list(
  run_scenario(n = 1000, p =  50, k =  100, label = "small"),
  run_scenario(n = 5000, p = 100, k =  500, label = "large")
)

# Summary table
cat("\n=== Summary ===\n")
summary_rows <- do.call(rbind, lapply(results, function(bm) {
  data.frame(
    scenario   = bm$scenario,
    expression = as.character(bm$expression),  # bch:expr -> character
    median_s   = round(as.numeric(bm$median), 4),
    stringsAsFactors = FALSE
  )
}))
print(summary_rows, row.names = FALSE)
