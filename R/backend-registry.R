.amatrix_optional_backends_enabled <- function() {
  !identical(getOption("amatrix.optional_backends", TRUE), FALSE)
}

.amatrix_optional_backend_specs <- function() {
  list(
    mlx = list(
      package = "amatrix.mlx",
      register_fun = "amatrix_mlx_register",
      enabled = function() TRUE
    ),
    metal = list(
      package = "amatrix.metal",
      register_fun = "amatrix_metal_register",
      enabled = function() isTRUE(getOption("amatrix.enable_metal", FALSE))
    ),
    opencl = list(
      package = "amatrix.opencl",
      register_fun = "amatrix_opencl_register",
      enabled = function() isTRUE(getOption("amatrix.enable_opencl", FALSE))
    ),
    arrayfire = list(
      package = "amatrix.arrayfire",
      register_fun = "amatrix_arrayfire_register",
      enabled = function() isTRUE(getOption("amatrix.enable_arrayfire", FALSE))
    )
  )
}

.amatrix_optional_backend_enabled <- function(spec) {
  enabled <- spec$enabled
  if (is.null(enabled)) {
    return(TRUE)
  }
  isTRUE(tryCatch(enabled(), error = function(e) FALSE))
}

.amatrix_try_register_optional_backend <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))

  if (!.amatrix_optional_backends_enabled()) {
    return(FALSE)
  }

  if (exists(name, envir = .amatrix_state$backends, inherits = FALSE)) {
    return(TRUE)
  }

  spec <- .amatrix_optional_backend_specs()[[name]]
  if (is.null(spec)) {
    return(FALSE)
  }
  if (!.amatrix_optional_backend_enabled(spec)) {
    return(FALSE)
  }

  lib_locs <- unique(c(.libPaths(), .Library, .Library.site))
  ns <- tryCatch(
    loadNamespace(spec$package, lib.loc = lib_locs),
    error = function(e) NULL
  )
  if (is.null(ns)) {
    return(FALSE)
  }

  register_backend <- get0(spec$register_fun, envir = ns, inherits = FALSE)
  if (!is.function(register_backend)) {
    return(FALSE)
  }

  isTRUE(tryCatch({
    register_backend(overwrite = TRUE)
    exists(name, envir = .amatrix_state$backends, inherits = FALSE)
  }, error = function(e) FALSE))
}

#' Register a backend with the amatrix dispatch system
#'
#' Adds a named backend to the session backend registry. The backend
#' must be a named list containing all required callable fields. Once
#' registered, the backend is available for dispatch by any
#' \code{aMatrix} object whose \code{preferred_backend} or \code{policy}
#' slot names it.
#'
#' @param name Character string. Unique identifier for the backend
#'   (e.g. \code{"mlx"}, \code{"opencl"}).
#' @param backend Named list implementing the backend contract. Required
#'   fields: \code{capabilities}, \code{features}, \code{precision_modes}
#'   (each a zero-argument function returning a character vector),
#'   \code{available} (zero-argument logical function), \code{supports},
#'   \code{matmul}, \code{crossprod}, \code{tcrossprod}, \code{ewise},
#'   \code{rowSums}, \code{colSums}.
#' @param overwrite Logical. Allow replacement of an existing registration
#'   with the same \code{name}. Default \code{FALSE}.
#'
#' @return Invisibly, \code{name}.
#'
#' @examples
#' # Minimal no-op backend for illustration only
#' noop <- list(
#'   capabilities   = function() character(),
#'   features       = function() character(),
#'   precision_modes = function() "strict",
#'   available      = function() FALSE,
#'   supports       = function(op, x, y = NULL) FALSE,
#'   matmul         = function(x, y) x,
#'   crossprod      = function(x, y = NULL) x,
#'   tcrossprod     = function(x, y = NULL) x,
#'   ewise          = function(x, y, op) x,
#'   rowSums        = function(x) numeric(nrow(x)),
#'   colSums        = function(x) numeric(ncol(x))
#' )
#' amatrix_register_backend("noop_test", noop, overwrite = TRUE)
#'
#' @seealso \code{\link{amatrix_backend_names}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_register_backend <- function(name, backend, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))

  if (!is.list(backend)) {
    stop("backend must be a named list")
  }

  required_fields <- c(
    "capabilities",
    "features",
    "precision_modes",
    "available",
    "supports",
    "matmul",
    "crossprod",
    "tcrossprod",
    "ewise",
    "rowSums",
    "colSums"
  )
  missing_fields <- setdiff(required_fields, names(backend))
  if (length(missing_fields) > 0L) {
    stop(sprintf("backend is missing required fields: %s", paste(missing_fields, collapse = ", ")))
  }

  if (!is.function(backend$capabilities)) {
    stop("backend$capabilities must be a function")
  }
  if (!is.function(backend$features)) {
    stop("backend$features must be a function")
  }
  if (!is.function(backend$precision_modes)) {
    stop("backend$precision_modes must be a function")
  }

  backend_capabilities <- backend$capabilities()
  backend_features <- backend$features()
  backend_precision_modes <- backend$precision_modes()
  if (!is.character(backend_capabilities)) {
    stop("backend$capabilities() must return a character vector")
  }
  if (!is.character(backend_features)) {
    stop("backend$features() must return a character vector")
  }
  if (!is.character(backend_precision_modes)) {
    stop("backend$precision_modes() must return a character vector")
  }
  if (!all(backend_precision_modes %in% .amatrix_valid_precisions)) {
    stop(sprintf(
      "backend$precision_modes() must be a subset of: %s",
      paste(.amatrix_valid_precisions, collapse = ", ")
    ))
  }

  exists_already <- exists(name, envir = .amatrix_state$backends, inherits = FALSE)
  if (exists_already && !overwrite) {
    stop(sprintf("backend '%s' is already registered", name))
  }

  assign(name, backend, envir = .amatrix_state$backends)
  invisible(name)
}

.amatrix_get_backend <- function(name) {
  backend <- get0(name, envir = .amatrix_state$backends, inherits = FALSE)
  if (is.null(backend) && isTRUE(.amatrix_try_register_optional_backend(name))) {
    backend <- get0(name, envir = .amatrix_state$backends, inherits = FALSE)
  }
  if (is.null(backend)) {
    stop(sprintf("backend '%s' is not registered", name))
  }
  backend
}

#' List names of all registered backends
#'
#' Returns the names of every backend currently in the session registry.
#' When optional backends are enabled (the default), this also attempts
#' to auto-register any installed optional backend packages before
#' returning the list.
#'
#' @return Character vector of registered backend names, sorted
#'   alphabetically. Always includes at least \code{"cpu"}.
#'
#' @examples
#' amatrix_backend_names()
#'
#' @seealso \code{\link{amatrix_backend_status}},
#'   \code{\link{amatrix_register_backend}}
#' @export
amatrix_backend_names <- function() {
  if (.amatrix_optional_backends_enabled()) {
    invisible(lapply(names(.amatrix_optional_backend_specs()), .amatrix_try_register_optional_backend))
  }
  sort(ls(envir = .amatrix_state$backends, all.names = FALSE))
}

#' Query the capabilities of a registered backend
#'
#' Returns the unique capability strings advertised by the named backend,
#' as reported by its \code{capabilities()} function.
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of capability identifiers (e.g.
#'   \code{"matmul"}, \code{"svd"}).
#'
#' @examples
#' amatrix_backend_capabilities("cpu")
#'
#' @seealso \code{\link{amatrix_backend_features}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_capabilities <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$capabilities())
}

#' Query the features of a registered backend
#'
#' Returns the unique feature strings advertised by the named backend,
#' as reported by its \code{features()} function. Features describe
#' optional capabilities such as sparse residency or deferred execution.
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of feature identifiers.
#'
#' @examples
#' amatrix_backend_features("cpu")
#'
#' @seealso \code{\link{amatrix_backend_capabilities}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_features <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$features())
}

#' Query the precision modes supported by a registered backend
#'
#' Returns the precision mode strings advertised by the named backend.
#' Valid values are \code{"strict"} (double precision) and \code{"fast"}
#' (single/mixed precision).
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of precision mode identifiers, a subset of
#'   \code{c("strict", "fast")}.
#'
#' @examples
#' amatrix_backend_precision_modes("cpu")
#'
#' @seealso \code{\link{amatrix_backend_capabilities}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_precision_modes <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$precision_modes())
}

#' Summarise the status of registered backends
#'
#' Returns a data.frame with one row per backend describing its
#' availability, supported precision modes, features, capabilities,
#' and whether it supports GPU residency.
#'
#' @param names Character vector of backend names to query. When
#'   \code{NULL} (default) all registered backends are included,
#'   with optional backends auto-registered first if possible.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{name}{Character. Backend identifier.}
#'     \item{available}{Logical. Whether the backend reports itself
#'       as available on this machine.}
#'     \item{precision_modes}{Character. Comma-separated precision
#'       modes (\code{"strict"}, \code{"fast"}).}
#'     \item{features}{Character. Comma-separated feature strings.}
#'     \item{residency_capable}{Logical. Whether the backend
#'       supports GPU-resident matrix storage.}
#'     \item{capabilities}{Character. Comma-separated operation
#'       capability strings.}
#'   }
#'
#' @examples
#' amatrix_backend_status()
#' amatrix_backend_status("cpu")
#'
#' @seealso \code{\link{amatrix_backend_names}},
#'   \code{\link{amatrix_register_backend}}
#' @export
amatrix_backend_status <- function(names = NULL) {
  if (is.null(names)) {
    if (.amatrix_optional_backends_enabled()) {
      invisible(lapply(names(.amatrix_optional_backend_specs()), .amatrix_try_register_optional_backend))
    }
    names <- amatrix_backend_names()
  }

  stopifnot(is.character(names))

  rows <- lapply(names, function(name) {
    backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
    if (is.null(backend)) {
      stop(sprintf("backend '%s' is not registered", name))
    }
    data.frame(
      name = name,
      available = isTRUE(backend$available()),
      precision_modes = paste(amatrix_backend_precision_modes(name), collapse = ","),
      features = paste(amatrix_backend_features(name), collapse = ","),
      residency_capable = .amatrix_backend_residency_capable(backend),
      capabilities = paste(amatrix_backend_capabilities(name), collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
