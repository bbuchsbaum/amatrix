# Testing kit for packages that build on amatrix.
#
# A downstream package's CI usually has no GPU. CPU fallback keeps its tests
# green, but then nothing verifies that the package's code actually reaches a
# backend, or guards GPU-only expectations. These helpers cover both.

#' Is an accelerator backend usable in this session?
#'
#' Reports whether at least one registered non-CPU backend says it is
#' available and has not been marked unhealthy.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @examples
#' amatrix_has_accelerator()
#'
#' @seealso \code{\link{amatrix_skip_if_no_accelerator}},
#'   \code{\link{amatrix_gpu_status}}
#' @export
amatrix_has_accelerator <- function() {
  names <- setdiff(tryCatch(amatrix_backend_names(), error = function(e) "cpu"), "cpu")
  for (name in names) {
    backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
    if (is.null(backend)) next
    available <- isTRUE(tryCatch(backend$available(), error = function(e) FALSE))
    health <- tryCatch(.amatrix_backend_health_get(name)$status, error = function(e) NULL)
    if (available && !identical(health, "unhealthy")) {
      return(TRUE)
    }
  }
  FALSE
}

#' Skip a test when no accelerator backend is usable
#'
#' For use inside \code{testthat} tests of packages that build on amatrix:
#' skips the current test unless \code{\link{amatrix_has_accelerator}()} is
#' \code{TRUE}. Tests of numerical results should normally \emph{not} need
#' this, because amatrix falls back to the CPU; use it for expectations that
#' only make sense on a real device (timings, residency).
#'
#' @return Invisibly \code{TRUE} when the test may proceed; otherwise the test
#'   is skipped.
#'
#' @examples
#' \dontrun{
#' testthat::test_that("device path is faster", {
#'   amatrix_skip_if_no_accelerator()
#'   # ... timing expectations ...
#' })
#' }
#'
#' @seealso \code{\link{amatrix_with_test_backend}}
#' @export
amatrix_skip_if_no_accelerator <- function() {
  if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("amatrix_skip_if_no_accelerator() requires the 'testthat' package", call. = FALSE)
  }
  if (!amatrix_has_accelerator()) {
    testthat::skip("no amatrix accelerator backend is available")
  }
  invisible(TRUE)
}

#' A stand-in backend for tests
#'
#' Builds a backend that satisfies the amatrix backend contract by delegating
#' every operation to the CPU reference implementation while counting calls.
#' Registered under a non-CPU name (see
#' \code{\link{amatrix_with_test_backend}}) it lets a package test, on any
#' machine, that its code really dispatches through a backend, without a GPU.
#' Results are numerically identical to the CPU path.
#'
#' @param ops Character vector of operations the backend claims to support.
#' @param precision_modes Precision modes the backend accepts.
#'
#' @return A backend list suitable for \code{\link{amatrix_register_backend}}.
#'   Its \code{calls} element is an environment holding one integer counter
#'   per operation that has been invoked.
#'
#' @examples
#' backend <- amatrix_test_backend()
#' amatrix_with_test_backend({
#'   x <- adgeMatrix(matrix(1:6, 2, 3), backend = "amatrix_test")
#'   y <- x %*% matrix(1, 3, 2)
#' }, backend = backend)
#' backend$calls$matmul
#'
#' @seealso \code{\link{amatrix_with_test_backend}}
#' @export
amatrix_test_backend <- function(
  ops = c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums"),
  precision_modes = c("strict", "fast")
) {
  stopifnot(is.character(ops), is.character(precision_modes))
  cpu <- .amatrix_cpu_backend()
  calls <- new.env(parent = emptyenv())
  counted <- function(method_name) {
    force(method_name)
    function(...) {
      calls[[method_name]] <- (if (is.null(calls[[method_name]])) 0L else calls[[method_name]]) + 1L
      do.call(cpu[[method_name]], list(...))
    }
  }
  backend <- list(
    capabilities = function() sort(unique(ops)),
    features = function() "dense_f64",
    precision_modes = function() precision_modes,
    available = function() TRUE,
    supports = function(op, x, y = NULL) {
      inherits(x, "adgeMatrix") && x@precision %in% precision_modes && op %in% ops
    },
    calls = calls
  )
  for (method_name in c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums",
                        "solve", "chol", "qr", "svd", "eigen", "diag")) {
    backend[[method_name]] <- counted(method_name)
  }
  backend
}

#' Run code with a test backend registered
#'
#' Registers \code{backend} under \code{name} for the duration of
#' \code{code}, marks it healthy, and restores the previous registry and
#' health state afterwards (also on error). Objects created inside
#' \code{code} with \code{backend = name} dispatch through it.
#'
#' @param code Expression to evaluate.
#' @param backend A backend list; defaults to a fresh
#'   \code{\link{amatrix_test_backend}()}.
#' @param name Name to register the backend under.
#'
#' @return The value of \code{code}.
#'
#' @examples
#' backend <- amatrix_test_backend()
#' result <- amatrix_with_test_backend({
#'   x <- adgeMatrix(matrix(rnorm(12), 4, 3), backend = "amatrix_test")
#'   as.matrix(crossprod(x))
#' }, backend = backend)
#' dim(result)
#' backend$calls$crossprod
#'
#' @seealso \code{\link{amatrix_test_backend}},
#'   \code{\link{amatrix_skip_if_no_accelerator}}
#' @export
amatrix_with_test_backend <- function(code, backend = amatrix_test_backend(), name = "amatrix_test") {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name), !identical(name, "cpu"))
  registry <- .amatrix_state$backends
  existed <- exists(name, envir = registry, inherits = FALSE)
  previous <- if (existed) get(name, envir = registry, inherits = FALSE) else NULL
  .amatrix_backend_health_init()
  health_env <- .amatrix_state$backend_health
  had_health <- exists(name, envir = health_env, inherits = FALSE)
  previous_health <- if (had_health) get(name, envir = health_env, inherits = FALSE) else NULL

  amatrix_register_backend(name, backend, overwrite = TRUE)
  .amatrix_backend_health_mark(name, "healthy", NULL)
  on.exit({
    if (existed) {
      assign(name, previous, envir = registry)
    } else if (exists(name, envir = registry, inherits = FALSE)) {
      rm(list = name, envir = registry)
    }
    if (had_health) {
      assign(name, previous_health, envir = health_env)
    } else if (exists(name, envir = health_env, inherits = FALSE)) {
      rm(list = name, envir = health_env)
    }
  }, add = TRUE)

  force(code)
}
