random_dense_case <- function(seed = NULL, nrow_range = 2:5, ncol_range = 2:5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  nr <- sample(nrow_range, 1)
  nc <- sample(ncol_range, 1)

  list(
    nr = nr,
    nc = nc,
    x = matrix(rnorm(nr * nc), nrow = nr, ncol = nc),
    y = matrix(rnorm(nr * nc), nrow = nr, ncol = nc),
    rhs = matrix(rnorm(nc * sample(2:4, 1)), nrow = nc)
  )
}

random_sparse_case <- function(seed = NULL, nrow_range = 3:6, ncol_range = 3:6) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  nr <- sample(nrow_range, 1)
  nc <- sample(ncol_range, 1)
  raw <- matrix(sample(c(0, 0, 0, rnorm(6)), nr * nc, replace = TRUE), nrow = nr, ncol = nc)

  list(
    nr = nr,
    nc = nc,
    raw = raw,
    rhs = matrix(rnorm(nc * sample(2:4, 1)), nrow = nc),
    host = as(Matrix::Matrix(raw, sparse = TRUE), "dgCMatrix")
  )
}

with_registered_backend <- function(name, backend, code) {
  existed <- exists(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
  previous <- if (existed) get(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE) else NULL

  amatrix_register_backend(name, backend, overwrite = TRUE)

  on.exit({
    if (existed) {
      assign(name, previous, envir = amatrix:::.amatrix_state$backends)
    } else if (exists(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
      rm(list = name, envir = amatrix:::.amatrix_state$backends)
    }
  }, add = TRUE)

  force(code)
}

make_recording_backend <- function(
  counter,
  supported_ops = c("matmul", "rowSums", "colSums", "ewise"),
  cold_supported_ops = supported_ops,
  resident_supported_ops = supported_ops,
  precision_modes = c("strict", "fast"),
  features = c("dense_f64"),
  supports_sparse_matmul = FALSE,
  supports_sparse_ops = if (isTRUE(supports_sparse_matmul)) "matmul" else character(),
  supports_sparse_resident = length(supports_sparse_ops) > 0L
) {
  cpu <- amatrix:::.amatrix_cpu_backend()
  resident <- new.env(parent = emptyenv())
  sparse_resident <- new.env(parent = emptyenv())

  wrap <- function(method_name) {
    force(method_name)
    function(...) {
      if (is.null(counter[[method_name]])) {
        counter[[method_name]] <- 0L
      }
      counter[[method_name]] <- counter[[method_name]] + 1L
      do.call(cpu[[method_name]], list(...))
    }
  }

  backend <- list(
    capabilities = function() sort(unique(c(cold_supported_ops, resident_supported_ops))),
    features = function() {
      feats <- features
      if (length(supports_sparse_ops) > 0L && !("sparse_spmm" %in% feats)) {
        feats <- c(feats, "sparse_spmm")
      }
      feats
    },
    precision_modes = function() precision_modes,
    available = function() TRUE,
    supports = function(op, x, y = NULL) {
      precision_ok <- inherits(x, "aMatrix") && (x@precision %in% precision_modes)
      dense_ok <- inherits(x, "adgeMatrix")
      sparse_ok <- inherits(x, "adgCMatrix") && (op %in% supports_sparse_ops)
      if (sparse_ok && op %in% c("crossprod", "tcrossprod") && is.null(y)) {
        sparse_ok <- FALSE
      }

      precision_ok &&
        (dense_ok || sparse_ok) &&
        op %in% cold_supported_ops
    },
    matmul = wrap("matmul"),
    crossprod = wrap("crossprod"),
    tcrossprod = wrap("tcrossprod"),
    ewise = wrap("ewise"),
    rowSums = wrap("rowSums"),
    colSums = wrap("colSums"),
    solve = wrap("solve"),
    chol = wrap("chol"),
    qr = wrap("qr"),
    svd = wrap("svd"),
    eigen = wrap("eigen"),
    diag = wrap("diag"),
    resident_has = function(key) exists(key, envir = resident, inherits = FALSE),
    resident_store = function(key, x) {
      if (is.null(counter$resident_store)) {
        counter$resident_store <- 0L
      }
      counter$resident_store <- counter$resident_store + 1L
      assign(key, as.matrix(x), envir = resident)
      invisible(key)
    },
    resident_drop = function(key) {
      if (is.null(counter$resident_drop)) {
        counter$resident_drop <- 0L
      }
      counter$resident_drop <- counter$resident_drop + 1L
      if (exists(key, envir = resident, inherits = FALSE)) {
        rm(list = key, envir = resident)
      }
      invisible(key)
    },
    resident_materialize = function(key) {
      if (is.null(counter$resident_materialize)) {
        counter$resident_materialize <- 0L
      }
      counter$resident_materialize <- counter$resident_materialize + 1L
      get(key, envir = resident, inherits = FALSE)
    }
  )

  if ("matmul" %in% resident_supported_ops) {
    backend$matmul_resident <- function(x_key, y_key, out_key, defer = FALSE) {
      if (is.null(counter$matmul)) {
        counter$matmul <- 0L
      }
      if (is.null(counter$matmul_resident)) {
        counter$matmul_resident <- 0L
      }
      counter$matmul <- counter$matmul + 1L
      counter$matmul_resident <- counter$matmul_resident + 1L
      value <- resident[[x_key]] %*% resident[[y_key]]
      assign(out_key, value, envir = resident)
      value
    }
  }

  if ("crossprod" %in% resident_supported_ops) {
    backend$crossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
      if (is.null(counter$crossprod)) {
        counter$crossprod <- 0L
      }
      if (is.null(counter$crossprod_resident)) {
        counter$crossprod_resident <- 0L
      }
      counter$crossprod <- counter$crossprod + 1L
      counter$crossprod_resident <- counter$crossprod_resident + 1L
      rhs <- if (is.null(y_key)) resident[[x_key]] else resident[[y_key]]
      value <- crossprod(resident[[x_key]], rhs)
      assign(out_key, value, envir = resident)
      value
    }
  }

  if ("tcrossprod" %in% resident_supported_ops) {
    backend$tcrossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
      if (is.null(counter$tcrossprod)) {
        counter$tcrossprod <- 0L
      }
      if (is.null(counter$tcrossprod_resident)) {
        counter$tcrossprod_resident <- 0L
      }
      counter$tcrossprod <- counter$tcrossprod + 1L
      counter$tcrossprod_resident <- counter$tcrossprod_resident + 1L
      rhs <- if (is.null(y_key)) resident[[x_key]] else resident[[y_key]]
      value <- tcrossprod(resident[[x_key]], rhs)
      assign(out_key, value, envir = resident)
      value
    }
  }

  if ("ewise" %in% resident_supported_ops) {
    backend$ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
      if (is.null(counter$ewise)) {
        counter$ewise <- 0L
      }
      if (is.null(counter$ewise_resident)) {
        counter$ewise_resident <- 0L
      }
      counter$ewise <- counter$ewise + 1L
      counter$ewise_resident <- counter$ewise_resident + 1L
      lhs <- resident[[lhs_key]]
      rhs_value <- if (is.character(rhs)) resident[[rhs]] else rhs
      value <- do.call(op, list(lhs, rhs_value))
      assign(out_key, value, envir = resident)
      value
    }
  }

  if (isTRUE(supports_sparse_resident) && length(supports_sparse_ops) > 0L) {
    backend$sparse_resident_has <- function(key) {
      exists(key, envir = sparse_resident, inherits = FALSE)
    }

    backend$sparse_resident_store <- function(key, x_sp) {
      if (is.null(counter$sparse_resident_store)) {
        counter$sparse_resident_store <- 0L
      }
      counter$sparse_resident_store <- counter$sparse_resident_store + 1L
      assign(key, methods::as(x_sp, "dgCMatrix"), envir = sparse_resident)
      invisible(key)
    }

    backend$sparse_resident_drop <- function(key) {
      if (is.null(counter$sparse_resident_drop)) {
        counter$sparse_resident_drop <- 0L
      }
      counter$sparse_resident_drop <- counter$sparse_resident_drop + 1L
      if (exists(key, envir = sparse_resident, inherits = FALSE)) {
        rm(list = key, envir = sparse_resident)
      }
      invisible(key)
    }

    backend$spmm_resident <- function(sp_key, B, trans_lhs = FALSE) {
      if (is.null(counter$spmm_resident)) {
        counter$spmm_resident <- 0L
      }
      counter$spmm_resident <- counter$spmm_resident + 1L

      if (isTRUE(trans_lhs)) {
        if (is.null(counter$spmm_resident_trans_true)) {
          counter$spmm_resident_trans_true <- 0L
        }
        counter$spmm_resident_trans_true <- counter$spmm_resident_trans_true + 1L
      } else {
        if (is.null(counter$spmm_resident_trans_false)) {
          counter$spmm_resident_trans_false <- 0L
        }
        counter$spmm_resident_trans_false <- counter$spmm_resident_trans_false + 1L
      }

      sp <- get(sp_key, envir = sparse_resident, inherits = FALSE)
      B_mat <- as.matrix(B)
      if (isTRUE(trans_lhs)) {
        as.matrix(Matrix::crossprod(sp, B_mat))
      } else {
        as.matrix(sp %*% B_mat)
      }
    }

    backend$spmm_resident_key <- function(sp_key, y_key, out_key, trans_lhs = FALSE, defer = FALSE) {
      if (is.null(counter$spmm_resident_key)) {
        counter$spmm_resident_key <- 0L
      }
      counter$spmm_resident_key <- counter$spmm_resident_key + 1L

      sp <- get(sp_key, envir = sparse_resident, inherits = FALSE)
      y_mat <- get(y_key, envir = resident, inherits = FALSE)
      value <- if (isTRUE(trans_lhs)) {
        as.matrix(Matrix::crossprod(sp, y_mat))
      } else {
        as.matrix(sp %*% y_mat)
      }
      assign(out_key, value, envir = resident)
      value
    }
  }

  backend
}

optional_backend_specs <- function() {
  list(
    list(
      package = "amatrix.mlx",
      repo_dir = "backends/amatrix.mlx",
      backend = "mlx",
      option = "amatrix.mlx.available",
      disable_option = "amatrix.disable_mlx",
      register_fun = "amatrix_mlx_register",
      capabilities_fun = "amatrix_mlx_capabilities",
      available_fun = "amatrix_mlx_is_available"
    ),
    list(
      package = "amatrix.opencl",
      repo_dir = "backends/amatrix.opencl",
      backend = "opencl",
      option = "amatrix.opencl.available",
      disable_option = "amatrix.disable_opencl",
      register_fun = "amatrix_opencl_register",
      capabilities_fun = "amatrix_opencl_capabilities",
      available_fun = "amatrix_opencl_is_available"
    ),
    list(
      package = "amatrix.arrayfire",
      repo_dir = "backends/amatrix.arrayfire",
      backend = "arrayfire",
      option = "amatrix.arrayfire.available",
      disable_option = "amatrix.disable_arrayfire",
      register_fun = "amatrix_arrayfire_register",
      capabilities_fun = "amatrix_arrayfire_capabilities",
      available_fun = "amatrix_arrayfire_is_available"
    )
  )
}

.optional_backend_repo_dir <- function(repo_dir) {
  if (is.null(repo_dir)) {
    return(NULL)
  }

  candidates <- unique(c(
    repo_dir,
    file.path(getwd(), repo_dir),
    file.path(getwd(), "..", repo_dir),
    file.path(getwd(), "..", "..", repo_dir)
  ))
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  matches <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(matches) == 0L) {
    return(NULL)
  }
  matches[[1L]]
}

skip_if_backend_package_missing <- function(spec) {
  if (is.null(optional_backend_namespace(spec$package))) {
    testthat::skip(sprintf("Optional backend package '%s' is not installed", spec$package))
  }
}

optional_backend_namespace <- function(package) {
  specs <- optional_backend_specs()
  spec_idx <- match(package, vapply(specs, `[[`, character(1), "package"))
  spec <- if (is.na(spec_idx)) NULL else specs[[spec_idx]]
  repo_dir <- if (is.null(spec)) NULL else .optional_backend_repo_dir(spec$repo_dir)
  if (!is.null(repo_dir) &&
      dir.exists(repo_dir) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    if (package %in% loadedNamespaces()) {
      ns <- asNamespace(package)
      ns_path <- tryCatch(getNamespaceInfo(ns, "path"), error = function(e) NULL)
      if (!is.null(ns_path) &&
          identical(
            normalizePath(ns_path, winslash = "/", mustWork = FALSE),
            normalizePath(repo_dir, winslash = "/", mustWork = FALSE)
          )) {
        return(ns)
      }
    }
    # A genuine source checkout can still be unbuildable in this environment
    # (e.g. the backend's native GPU SDK is absent on a CI runner). load_all()
    # then aborts inside compile_dll with "System command 'R' failed"; letting
    # that propagate turns every skip guard into a hard error. Treat an
    # unbuildable/unloadable source tree as "backend unavailable" and fall
    # through to the installed-package lookup below (mirrors the robust
    # ensure_optional_backend_namespace() in tools/benchmark-helpers.R).
    loaded <- tryCatch(
      {
        pkgload::load_all(repo_dir, quiet = TRUE, helpers = FALSE, export_all = FALSE)
        package %in% loadedNamespaces()
      },
      error = function(e) FALSE
    )
    if (isTRUE(loaded)) {
      ns <- asNamespace(package)
      # A source tree can load its R code yet carry a stale or absent compiled
      # bridge; the availability probe then throws '.Call ... not resolved'
      # instead of returning FALSE, and that degraded namespace poisons every
      # consumer for the rest of the run. Validate the probe; if it cannot
      # even run, unload and prefer the installed package below.
      probe_ok <- TRUE
      if (!is.null(spec) && !is.null(spec$available_fun) &&
          exists(spec$available_fun, envir = ns, inherits = FALSE)) {
        probe <- tryCatch(
          get(spec$available_fun, envir = ns, inherits = FALSE)(),
          error = function(e) e
        )
        probe_ok <- !inherits(probe, "error")
      }
      if (probe_ok) {
        return(ns)
      }
      try(pkgload::unload(package), silent = TRUE)
    }
  }

  lib_locs <- unique(c(.libPaths(), .Library, .Library.site))
  tryCatch(
    loadNamespace(package, lib.loc = lib_locs),
    error = function(e) NULL
  )
}

with_optional_backend_available <- function(spec, code) {
  if (!is.null(spec$disable_option)) {
    old_disable <- getOption(spec$disable_option)
    options(structure(list(FALSE), names = spec$disable_option))
    on.exit(options(structure(list(old_disable), names = spec$disable_option)), add = TRUE)
  }

  old <- getOption(spec$option)
  options(structure(list(TRUE), names = spec$option))
  on.exit(options(structure(list(old), names = spec$option)), add = TRUE)
  force(code)
}

backend_package_capabilities <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  if (is.null(ns)) {
    stop(sprintf("backend package '%s' is unavailable", spec$package), call. = FALSE)
  }
  get(spec$capabilities_fun, envir = ns, inherits = FALSE)()
}

backend_package_available <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  if (is.null(ns)) {
    stop(sprintf("backend package '%s' is unavailable", spec$package), call. = FALSE)
  }
  # The availability probe itself can error in degraded harnesses (e.g. the
  # namespace loaded from source without its compiled bridge, so .Call cannot
  # resolve the native symbol). A probe that cannot run means the backend is
  # not available here; report FALSE so guards skip instead of erroring.
  tryCatch(
    isTRUE(get(spec$available_fun, envir = ns, inherits = FALSE)()),
    error = function(e) FALSE
  )
}
