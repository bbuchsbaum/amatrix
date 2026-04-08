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
  features = c("dense_f64")
) {
  cpu <- amatrix:::.amatrix_cpu_backend()
  resident <- new.env(parent = emptyenv())

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
    features = function() features,
    precision_modes = function() precision_modes,
    available = function() TRUE,
    supports = function(op, x, y = NULL) {
      inherits(x, "adgeMatrix") &&
        (x@precision %in% precision_modes) &&
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

  backend
}

optional_backend_specs <- function() {
  list(
    list(
      package = "amatrix.mlx",
      backend = "mlx",
      option = "amatrix.mlx.available",
      enable_option = NULL,
      register_fun = "amatrix_mlx_register",
      capabilities_fun = "amatrix_mlx_capabilities",
      available_fun = "amatrix_mlx_is_available"
    ),
    list(
      package = "amatrix.arrayfire",
      backend = "arrayfire",
      option = "amatrix.arrayfire.available",
      enable_option = "amatrix.enable_arrayfire",
      register_fun = "amatrix_arrayfire_register",
      capabilities_fun = "amatrix_arrayfire_capabilities",
      available_fun = "amatrix_arrayfire_is_available"
    )
  )
}

skip_if_backend_package_missing <- function(spec) {
  if (is.null(optional_backend_namespace(spec$package))) {
    testthat::skip(sprintf("Optional backend package '%s' is not installed", spec$package))
  }
}

optional_backend_namespace <- function(package) {
  lib_locs <- unique(c(.libPaths(), .Library, .Library.site))
  tryCatch(
    loadNamespace(package, lib.loc = lib_locs),
    error = function(e) NULL
  )
}

with_optional_backend_available <- function(spec, code) {
  if (!is.null(spec$enable_option)) {
    old_enable <- getOption(spec$enable_option)
    options(structure(list(TRUE), names = spec$enable_option))
    on.exit(options(structure(list(old_enable), names = spec$enable_option)), add = TRUE)
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
  get(spec$available_fun, envir = ns, inherits = FALSE)()
}
