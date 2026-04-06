#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
})

parse_csv <- function(x) {
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}

parse_cases <- function(x) {
  parts <- parse_csv(x)
  lapply(parts, function(part) {
    dims <- strsplit(part, "x", fixed = TRUE)[[1L]]
    if (length(dims) != 2L) {
      stop(sprintf("Invalid case '%s'; expected NxM", part), call. = FALSE)
    }
    c(as.integer(dims[[1L]]), as.integer(dims[[2L]]))
  })
}

r_string <- function(x) {
  paste0("\"", gsub("([\"\\\\])", "\\\\\\1", x), "\"")
}

build_command <- function(repo_root, bench_lib, n, p, k, n_oversamples, n_iter) {
  lib_prefix <- if (nzchar(bench_lib)) {
    sprintf(".libPaths(c(%s, .libPaths())); ", r_string(normalizePath(bench_lib)))
  } else {
    ""
  }

  code <- paste0(
    "orig <- .libPaths(); ",
    lib_prefix,
    sprintf("pkgload::load_all(%s, quiet = TRUE); ", r_string(normalizePath(repo_root))),
    "invisible(loadNamespace(\"amatrix.mlx\")); ",
    "options(amatrix.mlx.available = TRUE); ",
    sprintf("n <- %dL; p <- %dL; k <- %dL; n_oversamples <- %dL; n_iter <- %dL; ",
            as.integer(n), as.integer(p), as.integer(k), as.integer(n_oversamples), as.integer(n_iter)),
    "set.seed(20260405L + n + p + k); ",
    "host <- matrix(rnorm(n * p), nrow = n, ncol = p); ",
    "invisible(amatrix.mlx:::amatrix_mlx_rsvd(matrix(rnorm(32L * 16L), nrow = 32L, ncol = 16L), k = 5L, n_oversamples = 4L, n_iter = 1L)); ",
    "t_cpu <- system.time(cpu <- irlba::svdr(host, k = k, extra = n_oversamples, it = n_iter)); ",
    "t_auto <- system.time(fac_auto <- am_svd_factor(adgeMatrix(host, preferred_backend = \"mlx\", precision = \"fast\"), k = k, method = \"auto\", n_oversamples = n_oversamples, n_iter = n_iter)); ",
    "ref <- base::svd(host, nu = k, nv = k)$d[seq_len(k)]; ",
    "plan <- amatrix:::.amatrix_svd_factor_plan(adgeMatrix(host, preferred_backend = \"mlx\", precision = \"fast\"), k, \"auto\", n_oversamples, n_iter); ",
    "err <- max(abs(fac_auto@d - ref) / pmax(abs(ref), 1e-12)); ",
    "cat(sprintf(\"case=%dx%d\\nk=%d\\nn_oversamples=%d\\nn_iter=%d\\nsvdr=%.3f\\nauto=%.3f\\nauto_method=%s\\nauto_backend=%s\\nauto_err=%.4f\\n\", n, p, k, n_oversamples, n_iter, unname(t_cpu[[\\\"elapsed\\\"]]), unname(t_auto[[\\\"elapsed\\\"]]), plan$method, fac_auto@backend, err))"
  )

  paste("Rscript -e", shQuote(code))
}

repo_root <- normalizePath(".")
bench_lib <- Sys.getenv("AMATRIX_BENCH_LIB", "")
k <- as.integer(Sys.getenv("AMATRIX_SVD_K", "20"))
cases <- parse_cases(Sys.getenv("AMATRIX_SVD_CASES", "400x320,500x400,700x560,1000x800"))
oversamples_grid <- as.integer(parse_csv(Sys.getenv("AMATRIX_SVD_OVERSAMPLES_GRID", "5,10,15")))
n_iter_grid <- as.integer(parse_csv(Sys.getenv("AMATRIX_SVD_N_ITER_GRID", "1,2,3")))

grid <- expand.grid(
  case_idx = seq_along(cases),
  n_oversamples = oversamples_grid,
  n_iter = n_iter_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

rows <- lapply(seq_len(nrow(grid)), function(idx) {
  case <- cases[[grid$case_idx[[idx]]]]
  n <- case[[1L]]
  p <- case[[2L]]
  data.frame(
    case = sprintf("%dx%d", n, p),
    k = k,
    n_oversamples = grid$n_oversamples[[idx]],
    n_iter = grid$n_iter[[idx]],
    command = build_command(
      repo_root = repo_root,
      bench_lib = bench_lib,
      n = n,
      p = p,
      k = k,
      n_oversamples = grid$n_oversamples[[idx]],
      n_iter = grid$n_iter[[idx]]
    ),
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)

cat("Notes:\n")
cat("- This script emits direct `Rscript -e` calibration commands because nested/file-entry MLX benchmark harnesses are unstable on this Apple Silicon setup.\n")
cat("- Each command performs a single untimed MLX warm-up call before measuring `am_svd_factor(method = \"auto\")` against CPU `irlba::svdr`.\n")
cat(sprintf("- Current default auto cutoff is min(dim) >= %d.\n", amatrix:::.amatrix_svd_factor_rsvd_min_dim()))
if (nzchar(bench_lib)) {
  cat(sprintf("- Commands prepend AMATRIX_BENCH_LIB=%s inside the R session.\n", bench_lib))
}
cat("\n")

write.table(out, row.names = FALSE, sep = "\t", quote = FALSE)
