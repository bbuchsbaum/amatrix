.amatrix_next_object_id <- function() {
  .amatrix_state$object_counter <- .amatrix_state$object_counter + 1L
  sprintf("%s:am:%d", .amatrix_state$session_id, .amatrix_state$object_counter)
}

.amatrix_make_finalizer_env <- function(object_id) {
  e <- new.env(parent = emptyenv())
  e$object_id <- object_id
  # cache_state holds the host_cache_valid flag in a child env that the
  # residency registry can reference WITHOUT creating a cycle back to `e`.
  # Previously the registry stored `finalizer_env = e` directly, which kept
  # `e` alive forever and prevented its reg.finalizer from ever firing —
  # causing residency entries to leak across GC cycles.
  e$cache_state <- new.env(parent = emptyenv())
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
      } else if (identical(backend_name, "mlx") && requireNamespace("amatrix.mlx", quietly = TRUE)) {
        ns <- asNamespace("amatrix.mlx")
        drop_fun <- tryCatch(
          get(
            if (isTRUE(entry$sparse)) "amatrix_mlx_sparse_drop" else "amatrix_mlx_resident_drop",
            envir = ns,
            inherits = FALSE
          ),
          error = function(e) NULL
        )
        if (!is.null(drop_fun)) {
          try(drop_fun(resident_key), silent = TRUE)
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

  # Store cache_state (a forward-only child env), NOT finalizer_env itself —
  # a direct ref to finalizer_env would form a cycle that prevents GC from
  # firing its finalizer, leaking residency entries.
  cache_state <- if (inherits(x, "adgeMatrix")) x@finalizer_env$cache_state else NULL
  assign(
    object_key,
    list(
      backend = backend,
      resident_key = resident_key,
      sparse = isTRUE(sparse),
      cache_state = cache_state
    ),
    envir = .amatrix_state$residency
  )

  if (inherits(x, "adgeMatrix") && !isTRUE(sparse) && !isTRUE(x@finalizer_env$host_deferred)) {
    x@finalizer_env$cache_state$host_cache_valid <- TRUE
  }

  x
}

.amatrix_update_resident_aliases <- function(backend, old_key, new_key = NULL, sparse = FALSE) {
  registry_keys <- ls(envir = .amatrix_state$residency, all.names = TRUE)
  updated <- 0L

  for (registry_key in registry_keys) {
    entry <- get0(registry_key, envir = .amatrix_state$residency, inherits = FALSE)
    if (is.null(entry) ||
        !identical(entry$backend, backend) ||
        !identical(entry$resident_key, old_key) ||
        !identical(isTRUE(entry$sparse), isTRUE(sparse))) {
      next
    }

    if (!is.null(new_key)) {
      entry$resident_key <- new_key
    }

    cs <- entry$cache_state
    if (!is.null(cs) && is.environment(cs)) {
      cs$host_cache_valid <- FALSE
    }

    assign(registry_key, entry, envir = .amatrix_state$residency)
    updated <- updated + 1L
  }

  updated
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

.amatrix_release_resident <- function(x) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  entry <- .amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(.amatrix_get_backend(entry$backend), error = function(e) NULL)
  if (!is.null(backend) && .amatrix_backend_residency_capable(backend)) {
    if (isTRUE(entry$sparse) &&
        is.function(backend$sparse_resident_has) &&
        isTRUE(backend$sparse_resident_has(entry$resident_key)) &&
        is.function(backend$sparse_resident_drop)) {
      try(backend$sparse_resident_drop(entry$resident_key), silent = TRUE)
    } else if (is.function(backend$resident_has) &&
               isTRUE(backend$resident_has(entry$resident_key)) &&
               is.function(backend$resident_drop)) {
      try(backend$resident_drop(entry$resident_key), silent = TRUE)
    }
  }

  .amatrix_drop_resident_binding(x)
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
  if (!.amatrix_backend_residency_capable(backend) ||
      !.amatrix_backend_has_resident_key(backend, entry$resident_key, sparse = isTRUE(entry$sparse))) {
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

.amatrix_backend_has_resident_key <- function(backend, resident_key, sparse = FALSE) {
  if (isTRUE(sparse) && is.function(backend$sparse_resident_has)) {
    return(isTRUE(backend$sparse_resident_has(resident_key)))
  }

  is.function(backend$resident_has) && isTRUE(backend$resident_has(resident_key))
}

.amatrix_sparse_resident_probe_rhs <- function(x, method) {
  if (!inherits(x, "adgCMatrix") || !(method %in% c("matmul", "crossprod", "tcrossprod"))) {
    return(NULL)
  }

  dims <- dim(x)
  if (is.null(dims) || length(dims) != 2L) {
    return(NULL)
  }

  nrow_probe <- switch(
    method,
    matmul = dims[[2L]],
    crossprod = dims[[1L]],
    tcrossprod = dims[[2L]],
    NULL
  )

  if (is.null(nrow_probe) || is.na(nrow_probe) || nrow_probe < 1L) {
    return(NULL)
  }

  matrix(0, nrow = as.integer(nrow_probe), ncol = 1L)
}

.amatrix_backend_supports_resident_op <- function(backend, method, x = NULL, y = NULL) {
  if (inherits(x, "adgCMatrix") && method %in% c("matmul", "crossprod", "tcrossprod")) {
    sparse_fn <- if (identical(method, "matmul") && is.function(backend$spmm_resident_key)) {
      backend$spmm_resident_key
    } else if (is.function(backend$spmm_resident)) {
      backend$spmm_resident
    } else {
      NULL
    }

    if (is.null(sparse_fn)) {
      return(FALSE)
    }

    predicate <- backend[["supports_resident"]]
    if (is.function(predicate) && !is.null(x)) {
      probe_y <- y
      if (is.null(probe_y)) {
        probe_y <- .amatrix_sparse_resident_probe_rhs(x, method)
      }
      return(isTRUE(predicate(method, x, y = probe_y)))
    }

    return(TRUE)
  }

  resident_fn <- backend[[paste0(method, "_resident")]]
  if (!is.function(resident_fn)) {
    return(FALSE)
  }

  predicate <- backend[["supports_resident"]]
  if (is.function(predicate) && !is.null(x)) {
    return(isTRUE(predicate(method, x, y = y)))
  }

  TRUE
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
  .amatrix_backend_residency_capable(backend) &&
    .amatrix_backend_has_resident_key(backend, entry$resident_key, sparse = isTRUE(entry$sparse))
}

#' Query GPU residency state of an aMatrix object
#'
#' Returns a single-row data.frame describing whether \code{x} is
#' currently uploaded to a GPU backend and, if so, which backend holds
#' it and whether that binding is still live (the device buffer still
#' exists).
#'
#' @param x An \code{aMatrix} object.
#'
#' @return A data.frame with one row and columns:
#'   \describe{
#'     \item{backend}{Character. Backend name, or \code{NA} when not
#'       resident.}
#'     \item{resident_key}{Character. Internal device buffer key, or
#'       \code{NA}.}
#'     \item{pinned_backend}{Character. Backend name when the binding
#'       is confirmed live, otherwise \code{NA}.}
#'     \item{live}{Logical. \code{TRUE} when the backend still holds
#'       the buffer identified by \code{resident_key}.}
#'   }
#'
#' @examples
#' m <- adgeMatrix(matrix(1:4, 2, 2))
#' amatrix_residency_info(m)
#'
#' @seealso \code{\link{amatrix_materialize_host}},
#'   \code{\link{amatrix_memory_stats}}
#' @export
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
    live = .amatrix_backend_has_resident_key(backend, entry$resident_key, sparse = isTRUE(entry$sparse)),
    stringsAsFactors = FALSE
  )
}

# ── Residency tripwire ──────────────────────────────────────────────────────
#
# Track 3 mechanism: every real GPU-to-host copy is routed through
# `.amatrix_tripwire_record()`. Tests that enable the tripwire (via the
# `amatrix.residency.tripwire` option or the `AMATRIX_RESIDENCY_TRIPWIRE=1`
# env var) assert that the counter is 0 after an op, which would mean S4
# dispatch silently fell through to a base/Matrix generic and pulled a
# resident tensor home. See planning_docs/quality-tracking.md §7 rule 6.
#
# The tripwire does NOT fire on host-only reads (e.g. an adgeMatrix that
# was never pushed to GPU) — only on actual backend$resident_materialize()
# calls.

.amatrix_tripwire_enabled <- function() {
  opt <- getOption("amatrix.residency.tripwire")
  if (!is.null(opt)) {
    return(isTRUE(opt))
  }
  identical(Sys.getenv("AMATRIX_RESIDENCY_TRIPWIRE", unset = "0"), "1")
}

.amatrix_tripwire_record <- function(op, backend_name, resident_key, dim = NULL) {
  if (!.amatrix_tripwire_enabled()) return(invisible())
  if (is.null(.amatrix_state$tripwire_events)) {
    .amatrix_state$tripwire_count <- 0L
    .amatrix_state$tripwire_events <- list()
  }
  .amatrix_state$tripwire_count <- .amatrix_state$tripwire_count + 1L
  idx <- length(.amatrix_state$tripwire_events) + 1L
  .amatrix_state$tripwire_events[[idx]] <- list(
    op = op,
    backend = backend_name,
    key = resident_key,
    dim = dim,
    timestamp = Sys.time()
  )
  invisible()
}

.amatrix_tripwire_count <- function() {
  .amatrix_state$tripwire_count %||% 0L
}

.amatrix_tripwire_events <- function() {
  .amatrix_state$tripwire_events %||% list()
}

.amatrix_tripwire_reset <- function() {
  .amatrix_state$tripwire_count <- 0L
  .amatrix_state$tripwire_events <- list()
  invisible()
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
          .amatrix_tripwire_record(
            "materialize_dense.deferred",
            entry$backend, entry$resident_key, x@Dim
          )
          mat <- backend$resident_materialize(entry$resident_key)
          host_mat <- if (is.matrix(mat)) mat else as.matrix(mat)
          base::dimnames(host_mat) <- x@Dimnames
          fenv$host_x <- host_mat
        }
      }
      if (is.null(fenv$host_x)) {
        stop("deferred adgeMatrix lost its GPU resident data", call. = FALSE)
      }
    }
    return(.amatrix_dense_base(fenv$host_x))
  }

  if (isTRUE(fenv$cache_state$host_cache_valid)) {
    return(new("dgeMatrix", x = x@x, Dim = x@Dim, Dimnames = x@Dimnames, factors = x@factors))
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

  .amatrix_tripwire_record(
    "materialize_dense.eager",
    entry$backend, entry$resident_key, x@Dim
  )
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

#' Force materialization of an aMatrix to a host Matrix object
#'
#' Downloads any GPU-resident data and returns a standard
#' \code{Matrix}-package object on the host. For \code{adgeMatrix}
#' inputs the result is a \code{dgeMatrix}; for \code{adgCMatrix}
#' inputs the result is a \code{dgCMatrix}; for \code{aTransposeView}
#' the transposed dense host matrix is returned. Host-only objects are
#' returned unchanged.
#'
#' @param x An \code{aMatrix} object (\code{adgeMatrix},
#'   \code{adgCMatrix}, or \code{aTransposeView}).
#'
#' @return A \code{dgeMatrix}, \code{dgCMatrix}, or the original object
#'   if no materialization is needed.
#'
#' @examples
#' m <- adgeMatrix(matrix(1:6, 2, 3))
#' host <- amatrix_materialize_host(m)
#' class(host)
#'
#' @seealso \code{\link{amatrix_residency_info}},
#'   \code{\link{amatrix_gc}}
#' @export
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
