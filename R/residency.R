.amatrix_next_object_id <- function() {
  .amatrix_state$object_counter <- .amatrix_state$object_counter + 1L
  sprintf("%s:am:%d", .amatrix_state$session_id, .amatrix_state$object_counter)
}

.amatrix_make_finalizer_env <- function(object_id) {
  e <- new.env(parent = emptyenv())
  e$object_id <- object_id
  reg.finalizer(
    e,
    function(env) {
      oid <- get0("object_id", envir = env, inherits = FALSE)
      if (is.null(oid) || !nzchar(oid)) {
        return(invisible(NULL))
      }
      object_key <- paste0("obj:", oid)
      entry <- get0(object_key, envir = .amatrix_state$residency, inherits = FALSE)
      if (is.null(entry)) {
        return(invisible(NULL))
      }
      backend_name <- entry$backend
      resident_key <- entry$resident_key
      rm(list = object_key, envir = .amatrix_state$residency)
      if (identical(backend_name, "mlx")) {
        if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
          return(invisible(NULL))
        }
        drop_fun <- tryCatch(
          get("amatrix_mlx_resident_drop", envir = asNamespace("amatrix.mlx"), inherits = FALSE),
          error = function(e) NULL
        )
        if (!is.null(drop_fun)) {
          try(drop_fun(resident_key), silent = TRUE)
        }
      } else {
        backend <- tryCatch(
          .amatrix_get_backend(backend_name),
          error = function(e) NULL
        )
        if (!is.null(backend) && .amatrix_backend_residency_capable(backend)) {
          if (isTRUE(entry$sparse) && is.function(backend$sparse_resident_drop)) {
            try(backend$sparse_resident_drop(resident_key), silent = TRUE)
          } else {
            try(backend$resident_drop(resident_key), silent = TRUE)
          }
        }
      }
      invisible(NULL)
    },
    onexit = TRUE
  )
  e
}

.amatrix_object_key <- function(x) {
  if (!inherits(x, "aMatrix")) {
    return(NULL)
  }

  paste0("obj:", x@object_id)
}

# Return the residency registry entry for x's source (original matrix when x is
# a transposed view). Returns NULL if x has no src_id or if the source is not
# currently in the registry.
.amatrix_src_resident_entry <- function(x) {
  if (!inherits(x, "aMatrix") || !nzchar(x@src_id)) return(NULL)
  get0(paste0("obj:", x@src_id), envir = .amatrix_state$residency, inherits = FALSE)
}

.amatrix_next_resident_key <- function(backend) {
  .amatrix_state$resident_counter <- .amatrix_state$resident_counter + 1L
  sprintf("%s:%d", backend, .amatrix_state$resident_counter)
}

.amatrix_resident_entry <- function(x) {
  object_key <- .amatrix_object_key(x)
  if (is.null(object_key)) {
    return(NULL)
  }
  get0(object_key, envir = .amatrix_state$residency, inherits = FALSE)
}

.amatrix_bind_resident <- function(x, backend, resident_key, sparse = FALSE) {
  object_key <- .amatrix_object_key(x)
  if (is.null(object_key)) {
    return(x)
  }

  assign(
    object_key,
    list(backend = backend, resident_key = resident_key, sparse = isTRUE(sparse)),
    envir = .amatrix_state$residency
  )

  x
}

.amatrix_drop_resident_binding <- function(x) {
  object_key <- .amatrix_object_key(x)
  if (is.null(object_key)) {
    return(invisible(FALSE))
  }

  if (exists(object_key, envir = .amatrix_state$residency, inherits = FALSE)) {
    rm(list = object_key, envir = .amatrix_state$residency)
    return(invisible(TRUE))
  }

  invisible(FALSE)
}

.amatrix_resident_backend <- function(x) {
  entry <- .amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(NULL)
  }
  entry$backend
}

.amatrix_live_resident_backend <- function(x) {
  entry <- .amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(NULL)
  }

  backend <- .amatrix_get_backend(entry$backend)
  if (!.amatrix_backend_residency_capable(backend) || !isTRUE(backend$resident_has(entry$resident_key))) {
    return(NULL)
  }

  entry$backend
}

.amatrix_resident_key <- function(x, backend = NULL) {
  entry <- .amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(NULL)
  }

  if (!is.null(backend) && !identical(entry$backend, backend)) {
    return(NULL)
  }

  entry$resident_key
}

.amatrix_backend_residency_capable <- function(backend) {
  required <- c("resident_store", "resident_has", "resident_drop")
  all(vapply(required, function(field) is.function(backend[[field]]), logical(1)))
}

.amatrix_backend_supports_resident_op <- function(backend, method) {
  is.function(backend[[paste0(method, "_resident")]])
}

.amatrix_object_is_resident <- function(x, backend_name) {
  if (!inherits(x, "aMatrix")) {
    return(FALSE)
  }

  entry <- .amatrix_resident_entry(x)
  if (is.null(entry) || !identical(entry$backend, backend_name)) {
    return(FALSE)
  }

  backend <- .amatrix_get_backend(backend_name)
  .amatrix_backend_residency_capable(backend) && isTRUE(backend$resident_has(entry$resident_key))
}

amatrix_residency_info <- function(x) {
  stopifnot(inherits(x, "aMatrix"))
  entry <- .amatrix_resident_entry(x)
  pinned_backend <- .amatrix_live_resident_backend(x)

  if (is.null(entry)) {
    return(data.frame(
      backend = NA_character_,
      resident_key = NA_character_,
      pinned_backend = NA_character_,
      live = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  backend <- .amatrix_get_backend(entry$backend)
  data.frame(
    backend = entry$backend,
    resident_key = entry$resident_key,
    pinned_backend = if (is.null(pinned_backend)) NA_character_ else pinned_backend,
    live = isTRUE(backend$resident_has(entry$resident_key)),
    stringsAsFactors = FALSE
  )
}

amatrix_materialize_dense <- function(x) {
  stopifnot(inherits(x, "adgeMatrix"))

  # ── Deferred path: host data not yet downloaded ──────────────────────────
  fenv <- x@finalizer_env
  if (isTRUE(fenv$host_deferred)) {
    if (is.null(fenv$host_x)) {
      # First host access — download from GPU and cache in the shared env
      entry <- .amatrix_resident_entry(x)
      if (!is.null(entry)) {
        backend <- tryCatch(.amatrix_get_backend(entry$backend), error = function(e) NULL)
        if (!is.null(backend) && is.function(backend$resident_materialize) &&
            isTRUE(backend$resident_has(entry$resident_key))) {
          mat <- backend$resident_materialize(entry$resident_key)
          fenv$host_x <- if (is.matrix(mat)) mat else as.matrix(mat)
        }
      }
      if (is.null(fenv$host_x)) {
        stop("deferred adgeMatrix lost its GPU resident data", call. = FALSE)
      }
    }
    return(.amatrix_dense_base(fenv$host_x))
  }

  # ── Eager path (existing logic) ─────────────────────────────────────────
  entry <- .amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(new("dgeMatrix", x = x@x, Dim = x@Dim, Dimnames = x@Dimnames, factors = x@factors))
  }

  backend <- .amatrix_get_backend(entry$backend)
  if (!is.function(backend$resident_materialize) || !isTRUE(backend$resident_has(entry$resident_key))) {
    return(new("dgeMatrix", x = x@x, Dim = x@Dim, Dimnames = x@Dimnames, factors = x@factors))
  }

  materialized <- backend$resident_materialize(entry$resident_key)
  if (inherits(materialized, "dgeMatrix")) {
    return(materialized)
  }
  if (inherits(materialized, "denseMatrix")) {
    return(.amatrix_dense_base(materialized))
  }
  if (is.matrix(materialized)) {
    return(.amatrix_dense_base(materialized))
  }

  stop("resident backend returned an unsupported dense materialization type")
}

amatrix_materialize_host <- function(x) {
  if (inherits(x, "adgeMatrix")) {
    return(amatrix_materialize_dense(x))
  }

  if (inherits(x, "aTransposeView")) {
    return(t(amatrix_materialize_dense(x@source)))
  }

  if (inherits(x, "adgCMatrix")) {
    return(new("dgCMatrix", i = x@i, p = x@p, Dim = x@Dim, Dimnames = x@Dimnames, x = x@x, factors = x@factors))
  }

  x
}
