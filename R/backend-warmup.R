# amatrix_warm() — pre-compile GPU kernels to eliminate cold-start JIT latency.
#
# GPU frameworks (MLX, ArrayFire) JIT-compile Metal/CUDA kernels on first use.
# That first-call cost can be 200–500 ms per kernel and makes initial benchmarks
# misleading. Call amatrix_warm() once before any timed work to pay that cost
# upfront.
#
# The function runs tiny dummy operations through each requested backend. Errors
# are silently swallowed — warming is always best-effort and never changes
# numerical state.

#' Warm up GPU backends to eliminate cold-start latency
#'
#' Pre-compiles GPU kernels by running tiny dummy operations through
#' each requested backend. Call once before timed work to pay JIT
#' compilation costs upfront. Errors are silently swallowed; warming
#' never alters numerical state.
#'
#' @param backend Character vector of backend names to warm, or
#'   \code{NULL} to warm all non-CPU backends currently registered.
#' @param ops Character vector of operation names to trigger.
#'   Recognised values: \code{"matmul"}, \code{"crossprod"},
#'   \code{"tcrossprod"}, \code{"qr"}, \code{"chol"}, \code{"svd"},
#'   \code{"solve"}.
#' @param size Integer vector of length 2 giving the dimensions
#'   \code{c(nrow, ncol)} of the dummy matrices used during warming.
#' @param quiet Logical; suppress progress messages when \code{TRUE}.
#'
#' @return An invisible named list, one entry per backend, each a list
#'   with elements \code{warmed} (logical) and \code{elapsed_ms}
#'   (numeric milliseconds, or \code{NA} when unavailable).
#'
#' @examples
#' \donttest{
#' results <- amatrix_warm(quiet = TRUE)
#' }
#'
#' @seealso \code{\link{amatrix_backend_names}}
#' @export
amatrix_warm <- function(
  backend = NULL,
  ops     = c("matmul", "crossprod", "qr", "chol"),
  size    = c(64L, 64L),
  quiet   = FALSE
) {
  stopifnot(is.numeric(size), length(size) == 2L, all(size >= 1L))

  if (is.null(backend)) {
    backend <- setdiff(amatrix_backend_names(), "cpu")
  }
  backend <- as.character(backend)

  if (length(backend) == 0L) {
    if (!quiet) message("amatrix_warm: no non-CPU backends available")
    return(invisible(list()))
  }

  results <- vector("list", length(backend))
  names(results) <- backend

  for (be in backend) {
    be_obj <- tryCatch(.amatrix_get_backend(be), error = function(e) NULL)

    if (is.null(be_obj) || !isTRUE(be_obj$available())) {
      results[[be]] <- list(warmed = FALSE, reason = "unavailable", elapsed_ms = NA_real_)
      if (!quiet) message(sprintf("amatrix_warm: '%s' not available, skipping", be))
      next
    }

    # Delegate to backend's own warm method when available (optional contract extension)
    if (is.function(be_obj[["warm"]])) {
      t0 <- proc.time()[["elapsed"]]
      tryCatch(be_obj$warm(ops = ops, size = size), error = function(e) NULL)
      elapsed_ms <- (proc.time()[["elapsed"]] - t0) * 1000
      results[[be]] <- list(warmed = TRUE, elapsed_ms = elapsed_ms)
      if (!quiet) message(sprintf("amatrix_warm: '%s' warmed in %.0f ms", be, elapsed_ms))
      next
    }

    # Generic warm path: run dummy ops to trigger kernel compilation.
    # Use fast precision where available — that is where GPU kernels live.
    precision <- if ("fast" %in% be_obj$precision_modes()) "fast" else "strict"

    nr  <- as.integer(size[[1L]])
    nc  <- as.integer(size[[2L]])

    t0 <- proc.time()[["elapsed"]]
    .amatrix_warm_backend(be_obj, be, ops, nr, nc, precision)
    elapsed_ms <- (proc.time()[["elapsed"]] - t0) * 1000

    results[[be]] <- list(warmed = TRUE, elapsed_ms = elapsed_ms)
    if (!quiet) message(sprintf("amatrix_warm: '%s' warmed in %.0f ms", be, elapsed_ms))
  }

  invisible(results)
}

.amatrix_warm_backend <- function(backend, be_name, ops, nr, nc, precision) {
  caps <- backend$capabilities()

  # Rectangular matrix: used for matmul, crossprod, qr, svd
  X_host <- matrix(seq_len(nr * nc) / (nr * nc + 1L), nr, nc)
  X <- as_adgeMatrix(X_host, preferred_backend = be_name, precision = precision)

  # Square SPD matrix: used for chol, solve
  spd_host <- crossprod(X_host) + diag(nc)
  SPD <- as_adgeMatrix(spd_host, preferred_backend = be_name, precision = precision)

  for (op in ops) {
    tryCatch(
      switch(op,
        matmul     = if ("matmul"    %in% caps) backend$matmul(X, X_host),
        crossprod  = if ("crossprod" %in% caps) backend$crossprod(X),
        tcrossprod = if ("tcrossprod" %in% caps) backend$tcrossprod(X),
        qr         = if ("qr"        %in% caps) backend$qr(X),
        chol       = if ("chol"      %in% caps) backend$chol(SPD),
        svd        = if ("svd"       %in% caps) backend$svd(X, nu = 2L, nv = 2L),
        solve      = if ("solve"     %in% caps) backend$solve(SPD)
      ),
      error = function(e) NULL  # best-effort; never abort
    )
  }

  invisible(NULL)
}
