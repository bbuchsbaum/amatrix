#' Time the same logit_fit() code on CPU and accelerated inputs
#'
#' Generates a synthetic logistic-regression problem, fits it with
#' [logit_fit()] on a plain base matrix (CPU reference), then -- if a healthy
#' GPU backend is installed -- fits the *same* function on an accelerated
#' `adgeMatrix` wrapping of the same data, and reports timings and the
#' maximum absolute coefficient difference.
#'
#' The point of the demo: the algorithm code is untouched; only the input
#' type changes.
#'
#' @param n number of observations.
#' @param p number of predictors (an intercept column is added).
#' @param reps timing repetitions; the median is reported.
#' @param seed RNG seed for the synthetic data.
#'
#' @return invisibly, a data frame with one row per input type
#'   (`backend`, `ms`, `max_coef_diff_vs_cpu`).
#' @export
logit_demo_benchmark <- function(n = 20000L, p = 400L, reps = 3L, seed = 1L) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (p - 1L)), n, p - 1L))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p - 1L)))
  beta_true <- c(0.5, stats::rnorm(p - 1L, sd = 0.2))
  y <- stats::rbinom(n, 1L, stats::plogis(as.numeric(X %*% beta_true)))

  time_fit <- function(input) {
    fit <- NULL
    elapsed <- vapply(seq_len(reps), function(i) {
      t0 <- proc.time()[["elapsed"]]
      fit <<- logit_fit(input, y)
      proc.time()[["elapsed"]] - t0
    }, numeric(1))
    list(ms = stats::median(elapsed) * 1e3, fit = fit)
  }

  cat(sprintf("logit_fit() on n=%d, p=%d (%d reps, median)\n\n", n, p, reps))

  cpu <- time_fit(X)
  cat(sprintf("  %-24s %8.1f ms\n", "base matrix (cpu):", cpu$ms))
  rows <- data.frame(
    backend = "cpu", ms = cpu$ms, max_coef_diff_vs_cpu = 0,
    stringsAsFactors = FALSE
  )

  # Backends register lazily; amatrix_use_gpu() is the documented way to
  # probe and enable whatever is installed. It also flips the session default
  # precision/policy, so restore those afterwards -- this is a benchmark
  # helper, not a session configurator.
  old_prec <- amatrix::amatrix_default_precision()
  old_pol  <- amatrix::amatrix_default_policy()
  on.exit({
    amatrix::amatrix_set_default_precision(old_prec)
    amatrix::amatrix_set_default_policy(old_pol)
  }, add = TRUE)
  amatrix::amatrix_use_gpu(quiet = TRUE)

  status <- amatrix::amatrix_gpu_status()
  avail  <- status[status$available & status$backend != "cpu", , drop = FALSE]
  if (nrow(avail) == 0L) {
    cat("\n  No GPU backend available -- CPU only.",
        "Install one (e.g. amatrix.mlx on Apple Silicon) and rerun.\n")
    return(invisible(rows))
  }

  for (bk in avail$backend) {
    X_acc <- amatrix::adgeMatrix(X, preferred_backend = bk, precision = "fast")
    res <- tryCatch(time_fit(X_acc), error = function(e) e)
    if (inherits(res, "error")) {
      cat(sprintf("  %-24s error: %s\n", paste0(bk, ":"), conditionMessage(res)))
      next
    }
    dcoef <- max(abs(res$fit$coefficients - cpu$fit$coefficients))
    cat(sprintf("  %-24s %8.1f ms   (%.1fx vs cpu, max |coef diff| %.1e)\n",
                paste0("adgeMatrix (", bk, "):"), res$ms, cpu$ms / res$ms, dcoef))
    rows <- rbind(rows, data.frame(
      backend = bk, ms = res$ms, max_coef_diff_vs_cpu = dcoef,
      stringsAsFactors = FALSE
    ))
  }

  cat("\n  Same logit_fit() code in every row -- only the input type changed.\n")
  invisible(rows)
}
