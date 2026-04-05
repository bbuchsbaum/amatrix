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
  if (is.null(backend)) {
    stop(sprintf("backend '%s' is not registered", name))
  }
  backend
}

amatrix_backend_names <- function() {
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

amatrix_backend_status <- function(names = amatrix_backend_names()) {
  stopifnot(is.character(names))

  rows <- lapply(names, function(name) {
    if (!(name %in% amatrix_backend_names())) {
      stop(sprintf("backend '%s' is not registered", name))
    }

    backend <- .amatrix_get_backend(name)
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
