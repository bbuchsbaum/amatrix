.amatrix_metal_probe_var <- "AMATRIX_METAL_PROBE_GPU"
.amatrix_metal_state <- new.env(parent = emptyenv())
.amatrix_metal_host_resident <- new.env(parent = emptyenv())

.amatrix_metal_probe_enabled <- function() {
  identical(Sys.getenv(.amatrix_metal_probe_var, unset = ""), "1")
}

.amatrix_metal_probe_cache_get <- function() {
  get0("native_available", envir = .amatrix_metal_state, inherits = FALSE)
}

.amatrix_metal_probe_cache_set <- function(value) {
  assign("native_available", isTRUE(value), envir = .amatrix_metal_state)
  invisible(isTRUE(value))
}

.amatrix_metal_probe_cache_clear <- function() {
  if (exists("native_available", envir = .amatrix_metal_state, inherits = FALSE)) {
    rm("native_available", envir = .amatrix_metal_state)
  }
  invisible(NULL)
}

.amatrix_metal_dense_host <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  x_mat
}

.amatrix_metal_sparse_host <- function(x) {
  methods::as(x, "dgCMatrix")
}

.amatrix_metal_rhs_width <- function(y) {
  if (is.null(y)) {
    return(NA_integer_)
  }

  dims <- dim(y)
  if (is.null(dims)) {
    return(1L)
  }

  if (length(dims) == 1L) {
    return(1L)
  }

  as.integer(dims[[2L]])
}

.amatrix_metal_product_width <- function(op, y) {
  if (identical(op, "tcrossprod")) {
    dims <- dim(y)
    if (is.null(dims)) {
      return(1L)
    }
    return(as.integer(dims[[1L]]))
  }

  .amatrix_metal_rhs_width(y)
}

amatrix_metal_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod")
}

amatrix_metal_features <- function() {
  c("sparse_spmm", "resident_sparse", "metal")
}

amatrix_metal_precision_modes <- function() {
  "fast"
}

amatrix_metal_native_available <- function(force = FALSE) {
  if (!force) {
    cached <- .amatrix_metal_probe_cache_get()
    if (!is.null(cached)) {
      return(cached)
    }
  }

  if (!.amatrix_metal_probe_enabled()) {
    return(FALSE)
  }

  available <- isTRUE(.Call("amatrix_metal_native_available_bridge", PACKAGE = "amatrix.metal"))
  .amatrix_metal_probe_cache_set(available)
}

amatrix_metal_enable_probe <- function(register = TRUE) {
  Sys.setenv(AMATRIX_METAL_PROBE_GPU = "1")
  options(amatrix.enable_metal = TRUE)
  .amatrix_metal_probe_cache_clear()

  available <- amatrix_metal_native_available(force = TRUE)
  if (isTRUE(register)) {
    try(amatrix_metal_register(overwrite = TRUE), silent = TRUE)
  }

  invisible(available)
}

amatrix_metal_is_available <- function() {
  isTRUE(getOption("amatrix.metal.available", FALSE)) || isTRUE(amatrix_metal_native_available())
}

amatrix_metal_bridge_info <- function() {
  info <- .Call("amatrix_metal_bridge_info_bridge", PACKAGE = "amatrix.metal")
  info$available <- amatrix_metal_is_available()
  info$capabilities <- amatrix_metal_capabilities()
  info
}

amatrix_metal_profile_enable <- function(enabled = TRUE, reset = FALSE) {
  if (isTRUE(reset)) {
    .Call("amatrix_metal_profile_reset_bridge", PACKAGE = "amatrix.metal")
  }
  .Call(
    "amatrix_metal_profile_set_enabled_bridge",
    isTRUE(enabled),
    PACKAGE = "amatrix.metal"
  )
  invisible(isTRUE(enabled))
}

amatrix_metal_profile_reset <- function() {
  .Call("amatrix_metal_profile_reset_bridge", PACKAGE = "amatrix.metal")
  invisible(TRUE)
}

amatrix_metal_profile <- function(reset = FALSE) {
  out <- .Call("amatrix_metal_profile_bridge", PACKAGE = "amatrix.metal")
  if (isTRUE(reset)) {
    amatrix_metal_profile_reset()
  }
  out
}

amatrix_metal_spmm <- function(x_sp, y, trans_lhs = FALSE) {
  y_mat <- .amatrix_metal_dense_host(y)
  .Call(
    "amatrix_metal_spmm_bridge",
    as.double(x_sp@x),
    as.integer(x_sp@p),
    as.integer(x_sp@i),
    as.integer(x_sp@Dim),
    y_mat,
    as.logical(trans_lhs),
    PACKAGE = "amatrix.metal"
  )
}

amatrix_metal_spmm_resident_key <- function(sp_key, y_key, out_key, trans_lhs = FALSE, defer = FALSE) {
  if (isTRUE(amatrix_metal_native_available())) {
    return(.Call(
      "amatrix_metal_spmm_resident_key_bridge",
      as.character(sp_key),
      as.character(y_key),
      as.character(out_key),
      as.logical(trans_lhs),
      as.logical(defer),
      PACKAGE = "amatrix.metal"
    ))
  }

  y_mat <- amatrix_metal_resident_materialize(y_key)
  value <- .Call(
    "amatrix_metal_spmm_resident_bridge",
    as.character(sp_key),
    y_mat,
    as.logical(trans_lhs),
    PACKAGE = "amatrix.metal"
  )
  amatrix_metal_resident_store(out_key, value)
  if (isTRUE(defer)) {
    return(NULL)
  }
  value
}

amatrix_metal_dense_sparse_matmul_resident_key <- function(x_key, sp_key, out_key, defer = FALSE) {
  if (isTRUE(amatrix_metal_native_available())) {
    return(.Call(
      "amatrix_metal_dense_sparse_matmul_resident_key_bridge",
      as.character(x_key),
      as.character(sp_key),
      as.character(out_key),
      as.logical(defer),
      PACKAGE = "amatrix.metal"
    ))
  }

  x_mat <- amatrix_metal_resident_materialize(x_key)
  value <- t(.Call(
    "amatrix_metal_spmm_resident_bridge",
    as.character(sp_key),
    t(x_mat),
    TRUE,
    PACKAGE = "amatrix.metal"
  ))
  amatrix_metal_resident_store(out_key, value)
  if (isTRUE(defer)) {
    return(NULL)
  }
  value
}

amatrix_metal_sparse_resident_store <- function(key, x_sp) {
  .Call(
    "amatrix_metal_sparse_store_bridge",
    as.character(key),
    as.double(x_sp@x),
    as.integer(x_sp@p),
    as.integer(x_sp@i),
    as.integer(x_sp@Dim),
    PACKAGE = "amatrix.metal"
  )
  invisible(TRUE)
}

amatrix_metal_sparse_resident_has <- function(key) {
  isTRUE(.Call("amatrix_metal_sparse_has_bridge", as.character(key), PACKAGE = "amatrix.metal"))
}

amatrix_metal_sparse_resident_drop <- function(key) {
  .Call("amatrix_metal_sparse_drop_bridge", as.character(key), PACKAGE = "amatrix.metal")
  invisible(TRUE)
}

amatrix_metal_resident_store <- function(key, x) {
  key <- as.character(key)
  x_mat <- .amatrix_metal_dense_host(x)

  if (isTRUE(amatrix_metal_native_available())) {
    .Call("amatrix_metal_dense_store_bridge", key, x_mat, PACKAGE = "amatrix.metal")
  } else {
    assign(key, x_mat, envir = .amatrix_metal_host_resident)
  }
  invisible(TRUE)
}

amatrix_metal_resident_has <- function(key) {
  key <- as.character(key)
  native_has <- if (isTRUE(amatrix_metal_native_available())) {
    isTRUE(.Call("amatrix_metal_dense_has_bridge", key, PACKAGE = "amatrix.metal"))
  } else {
    FALSE
  }

  native_has || exists(key, envir = .amatrix_metal_host_resident, inherits = FALSE)
}

amatrix_metal_resident_drop <- function(key) {
  key <- as.character(key)
  if (isTRUE(amatrix_metal_native_available())) {
    .Call("amatrix_metal_dense_drop_bridge", key, PACKAGE = "amatrix.metal")
  }
  if (exists(key, envir = .amatrix_metal_host_resident, inherits = FALSE)) {
    rm(list = key, envir = .amatrix_metal_host_resident)
  }
  invisible(TRUE)
}

amatrix_metal_resident_materialize <- function(key) {
  key <- as.character(key)
  if (isTRUE(amatrix_metal_native_available()) &&
      isTRUE(.Call("amatrix_metal_dense_has_bridge", key, PACKAGE = "amatrix.metal"))) {
    return(.Call("amatrix_metal_dense_materialize_bridge", key, PACKAGE = "amatrix.metal"))
  }

  if (!exists(key, envir = .amatrix_metal_host_resident, inherits = FALSE)) {
    stop(sprintf("unknown metal resident key '%s'", key), call. = FALSE)
  }
  get(key, envir = .amatrix_metal_host_resident, inherits = FALSE)
}

amatrix_metal_transpose_resident <- function(x_key, out_key) {
  if (isTRUE(amatrix_metal_native_available()) &&
      isTRUE(.Call("amatrix_metal_dense_has_bridge", as.character(x_key), PACKAGE = "amatrix.metal"))) {
    return(invisible(.Call(
      "amatrix_metal_transpose_resident_bridge",
      as.character(x_key),
      as.character(out_key),
      PACKAGE = "amatrix.metal"
    )))
  }

  x_mat <- amatrix_metal_resident_materialize(x_key)
  amatrix_metal_resident_store(out_key, t(x_mat))
  invisible(TRUE)
}

amatrix_metal_backend <- function() {
  cpu <- amatrix:::.amatrix_cpu_backend()

  list(
    capabilities = function() {
      amatrix_metal_capabilities()
    },
    features = function() {
      amatrix_metal_features()
    },
    precision_modes = function() {
      amatrix_metal_precision_modes()
    },
    available = function() {
      amatrix_metal_is_available()
    },
    supports = function(op, x, y = NULL) {
      if (!amatrix_metal_is_available()) {
        return(FALSE)
      }

      if (!inherits(x, "adgCMatrix")) {
        return(FALSE)
      }

      if (!(op %in% c("matmul", "crossprod", "tcrossprod"))) {
        return(FALSE)
      }

      if (is.null(y)) {
        return(FALSE)
      }

      if (!(x@precision %in% amatrix_metal_precision_modes())) {
        return(FALSE)
      }

      nnz <- length(x@x)
      rhs_width <- .amatrix_metal_rhs_width(y)
      min_nnz <- if (!is.na(rhs_width) && rhs_width <= 1L) {
        # Cold sparse MV is rarely worthwhile on Apple GPUs; keep the default
        # threshold high and rely on explicit resident prebinding for hot loops.
        getOption("amatrix.metal.spmv_min_nnz", 2000000L)
      } else {
        # Warm resident SpMM is fast, but first-call upload cost is still high.
        # Default to a conservative cold threshold and let resident reuse bypass
        # this gate once the sparse operand is pinned on device.
        getOption("amatrix.metal.spmm_min_nnz", 1000000L)
      }

      nnz >= as.integer(min_nnz)
    },
    matmul = function(x, y) {
      if (inherits(x, "adgCMatrix")) {
        return(amatrix_metal_spmm(.amatrix_metal_sparse_host(x), y, trans_lhs = FALSE))
      }
      cpu$matmul(x, y)
    },
    crossprod = function(x, y = NULL) {
      if (inherits(x, "adgCMatrix") && !is.null(y)) {
        return(amatrix_metal_spmm(.amatrix_metal_sparse_host(x), y, trans_lhs = TRUE))
      }
      cpu$crossprod(x, y)
    },
    tcrossprod = function(x, y = NULL) {
      if (inherits(x, "adgCMatrix") && !is.null(y)) {
        return(amatrix_metal_spmm(.amatrix_metal_sparse_host(x), t(.amatrix_metal_dense_host(y)), trans_lhs = FALSE))
      }
      cpu$tcrossprod(x, y)
    },
    ewise = function(lhs, rhs = NULL, op) {
      cpu$ewise(lhs, rhs, op)
    },
    rowSums = function(x) {
      cpu$rowSums(x)
    },
    colSums = function(x) {
      cpu$colSums(x)
    },
    resident_store = function(key, x) {
      amatrix_metal_resident_store(key, x)
    },
    resident_has = function(key) {
      isTRUE(amatrix_metal_resident_has(key))
    },
    resident_drop = function(key) {
      amatrix_metal_resident_drop(key)
    },
    resident_materialize = function(key) {
      amatrix_metal_resident_materialize(key)
    },
    supports_resident = function(op, x, y = NULL) {
      if (!amatrix_metal_is_available() ||
          !inherits(x, "adgCMatrix") ||
          !(op %in% c("matmul", "crossprod", "tcrossprod")) ||
          !(x@precision %in% amatrix_metal_precision_modes())) {
        return(FALSE)
      }

      if (is.null(y)) {
        return(TRUE)
      }

      nnz <- length(x@x)
      product_width <- .amatrix_metal_product_width(op, y)
      if (!is.na(product_width) && product_width <= 1L) {
        min_nnz <- getOption(
          "amatrix.metal.resident_spmv_min_nnz",
          getOption("amatrix.metal.spmv_min_nnz", 2000000L)
        )
      } else {
        min_nnz <- getOption("amatrix.metal.resident_spmm_min_nnz", 1L)
      }

      nnz >= as.integer(min_nnz)
    },
    transpose_resident = function(x_key, out_key) {
      amatrix_metal_transpose_resident(x_key, out_key)
    },
    sparse_resident_store = function(key, x_sp) {
      amatrix_metal_sparse_resident_store(key, x_sp)
    },
    sparse_resident_has = function(key) {
      isTRUE(amatrix_metal_sparse_resident_has(key))
    },
    sparse_resident_drop = function(key) {
      amatrix_metal_sparse_resident_drop(key)
    },
    spmm_resident = function(sp_key, B, trans_lhs = FALSE) {
      B_mat <- .amatrix_metal_dense_host(B)
      .Call(
        "amatrix_metal_spmm_resident_bridge",
        as.character(sp_key),
        B_mat,
        as.logical(trans_lhs),
        PACKAGE = "amatrix.metal"
      )
    },
    spmm_resident_key = function(sp_key, y_key, out_key, trans_lhs = FALSE, defer = FALSE) {
      amatrix_metal_spmm_resident_key(
        sp_key = sp_key,
        y_key = y_key,
        out_key = out_key,
        trans_lhs = trans_lhs,
        defer = defer
      )
    },
    dense_sparse_matmul_resident_key = function(x_key, sp_key, out_key, defer = FALSE) {
      amatrix_metal_dense_sparse_matmul_resident_key(
        x_key = x_key,
        sp_key = sp_key,
        out_key = out_key,
        defer = defer
      )
    }
  )
}

amatrix_metal_register <- function(overwrite = TRUE) {
  amatrix_register_backend("metal", amatrix_metal_backend(), overwrite = overwrite)
  invisible("metal")
}
