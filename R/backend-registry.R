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

amatrix_backend_names <- function() {
  if (.amatrix_optional_backends_enabled()) {
    invisible(lapply(names(.amatrix_optional_backend_specs()), .amatrix_try_register_optional_backend))
  }
  sort(ls(envir = .amatrix_state$backends, all.names = FALSE))
}

amatrix_backend_capabilities <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$capabilities())
}

amatrix_backend_features <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$features())
}

amatrix_backend_precision_modes <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$precision_modes())
}

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
