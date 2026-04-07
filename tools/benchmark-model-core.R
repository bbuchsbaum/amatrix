#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
})

bench_expr <- function(fn, reps = 3L) {
  stopifnot(is.function(fn))
  system.time({
    for (i in seq_len(reps)) {
      invisible(fn())
    }
  })[["elapsed"]] / reps
}

shared_response_list <- function(n, p, k, responses, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- matrix(rnorm(p * k), nrow = p, ncol = k)

  Ys <- lapply(seq_len(responses), function(i) {
    X %*% beta + matrix(rnorm(n * k, sd = 0.05), nrow = n, ncol = k)
  })

  list(X = X, Ys = Ys)
}

bench_shared_lm <- function(n = 2000L, p = 32L, k = 8L, responses = 8L, reps = 3L) {
  case <- shared_response_list(n = n, p = p, k = k, responses = responses, seed = 11L)
  x_am <- adgeMatrix(case$X)
  w <- seq_len(n) / n

  data.frame(
    workload = "shared_x_lm",
    mode = c(
      "normal_cache_off",
      "normal_cache_on",
      "qr_cache_off",
      "qr_cache_on",
      "weighted_qr_cache_on"
    ),
    n = n,
    p = p,
    k = k,
    responses = responses,
    lambda = NA_real_,
    elapsed = c(
      bench_expr(function() {
        for (y in case$Ys) {
          many_lm(x_am, y, include_residuals = FALSE, cache = FALSE, method = "normal")
        }
      }, reps = reps),
      bench_expr(function() {
        for (y in case$Ys) {
          many_lm(x_am, y, include_residuals = FALSE, cache = TRUE, method = "normal")
        }
      }, reps = reps),
      bench_expr(function() {
        for (y in case$Ys) {
          many_lm(x_am, y, include_residuals = FALSE, cache = FALSE, method = "qr")
        }
      }, reps = reps),
      bench_expr(function() {
        for (y in case$Ys) {
          many_lm(x_am, y, include_residuals = FALSE, cache = TRUE, method = "qr")
        }
      }, reps = reps),
      bench_expr(function() {
        for (y in case$Ys) {
          many_lm(x_am, y, weights = w, include_residuals = FALSE, cache = TRUE, method = "qr")
        }
      }, reps = reps)
    )
  )
}

bench_shared_ridge <- function(n = 2000L, p = 32L, k = 8L, responses = 8L, reps = 3L, lambda = 0.5) {
  case <- shared_response_list(n = n, p = p, k = k, responses = responses, seed = 17L)
  x_am <- adgeMatrix(case$X)

  data.frame(
    workload = "shared_x_ridge",
    mode = c("cache_off", "cache_on"),
    n = n,
    p = p,
    k = k,
    responses = responses,
    lambda = lambda,
    elapsed = c(
      bench_expr(function() {
        for (y in case$Ys) {
          ridge_fit(
            x_am,
            y,
            lambda = lambda,
            include_fitted = FALSE,
            include_residuals = FALSE,
            cache = FALSE
          )
        }
      }, reps = reps),
      bench_expr(function() {
        for (y in case$Ys) {
          ridge_fit(
            x_am,
            y,
            lambda = lambda,
            include_fitted = FALSE,
            include_residuals = FALSE,
            cache = TRUE
          )
        }
      }, reps = reps)
    )
  )
}

bench_similarity <- function(n = 4000L, p = 64L, reps = 3L) {
  set.seed(29)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x_am <- adgeMatrix(X)
  w <- seq_len(n) / n

  data.frame(
    workload = "similarity",
    mode = c("covariance", "weighted_covariance", "correlation"),
    n = n,
    p = p,
    k = NA_integer_,
    responses = NA_integer_,
    lambda = NA_real_,
    elapsed = c(
      bench_expr(function() {
        covariance(x_am)
      }, reps = reps),
      bench_expr(function() {
        covariance(x_am, weights = w)
      }, reps = reps),
      bench_expr(function() {
        correlation(x_am)
      }, reps = reps)
    )
  )
}

results <- rbind(
  bench_shared_lm(),
  bench_shared_ridge(),
  bench_similarity()
)

print(results)
