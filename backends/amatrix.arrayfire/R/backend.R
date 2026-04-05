amatrix_arrayfire_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums", "qr")
}

amatrix_arrayfire_features <- function() {
  c("dense_f32", "qr", "op_resident")
}

amatrix_arrayfire_precision_modes <- function() {
  "fast"
}

amatrix_arrayfire_native_available <- function() {
  .Call("amatrix_arrayfire_native_available_bridge")
}

amatrix_arrayfire_is_available <- function() {
  isTRUE(getOption("amatrix.arrayfire.available", FALSE)) || isTRUE(amatrix_arrayfire_native_available())
}

amatrix_arrayfire_bridge_info <- function() {
  info <- .Call("amatrix_arrayfire_bridge_info_bridge")
  info$available <- amatrix_arrayfire_is_available()
  info$capabilities <- amatrix_arrayfire_capabilities()
  info
}

amatrix_arrayfire_diagnostics <- function() {
  .Call("amatrix_arrayfire_diagnostics_bridge")
}

amatrix_arrayfire_active_backend <- function() {
  amatrix_arrayfire_diagnostics()$active_backend
}

amatrix_arrayfire_set_backend <- function(backend = c("cpu", "opencl", "cuda", "oneapi")) {
  backend <- match.arg(backend)
  backend_id <- switch(
    backend,
    cpu = 1L,
    cuda = 2L,
    opencl = 4L,
    oneapi = 8L
  )
  invisible(.Call("amatrix_arrayfire_set_backend_bridge", as.integer(backend_id)))
}

.amatrix_arrayfire_qr_experimental <- function() {
  isTRUE(getOption("amatrix.arrayfire.experimental_qr", FALSE))
}

.amatrix_arrayfire_qr_safe <- function() {
  identical(amatrix_arrayfire_active_backend(), 1L) || .amatrix_arrayfire_qr_experimental()
}

amatrix_arrayfire_matmul <- function(x, y) {
  x_mat <- as.matrix(x)
  y_mat <- as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_arrayfire_matmul_bridge", x_mat, y_mat)
}

amatrix_arrayfire_crossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) NULL else as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.null(y_mat) && !is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_arrayfire_crossprod_bridge", x_mat, y_mat)
}

amatrix_arrayfire_tcrossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) NULL else as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.null(y_mat) && !is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_arrayfire_tcrossprod_bridge", x_mat, y_mat)
}

amatrix_arrayfire_ewise <- function(lhs, rhs = NULL, op) {
  lhs_mat <- as.matrix(lhs)
  rhs_arg <- rhs

  if (!is.double(lhs_mat)) {
    storage.mode(lhs_mat) <- "double"
  }

  if (is.matrix(rhs_arg)) {
    rhs_arg <- as.matrix(rhs_arg)
    if (!is.double(rhs_arg)) {
      storage.mode(rhs_arg) <- "double"
    }
  } else if (is.numeric(rhs_arg) && length(rhs_arg) == 1L) {
    rhs_arg <- as.double(rhs_arg)
  } else if (!is.null(rhs_arg)) {
    stop("rhs must be NULL, a scalar, or a matrix")
  }

  .Call("amatrix_arrayfire_ewise_bridge", lhs_mat, rhs_arg, as.character(op))
}

amatrix_arrayfire_axis_sums <- function(x, axis) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  .Call("amatrix_arrayfire_sum_axis_bridge", x_mat, as.integer(axis))
}

amatrix_arrayfire_qr <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  .Call("amatrix_arrayfire_qr_bridge", x_mat)
}

amatrix_arrayfire_resident_store <- function(key, x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  invisible(.Call("amatrix_arrayfire_resident_store_bridge", as.character(key), x_mat))
}

amatrix_arrayfire_resident_has <- function(key) {
  isTRUE(.Call("amatrix_arrayfire_resident_has_bridge", as.character(key)))
}

amatrix_arrayfire_resident_drop <- function(key) {
  invisible(.Call("amatrix_arrayfire_resident_drop_bridge", as.character(key)))
}

amatrix_arrayfire_resident_materialize <- function(key) {
  .Call("amatrix_arrayfire_resident_materialize_bridge", as.character(key))
}

amatrix_arrayfire_matmul_resident <- function(x_key, y_key, out_key) {
  .Call("amatrix_arrayfire_matmul_resident_bridge",
        as.character(x_key), as.character(y_key), as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_crossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_crossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_tcrossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_tcrossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_ewise_resident <- function(lhs_key, rhs, op, out_key) {
  rhs_arg <- if (is.character(rhs)) as.character(rhs)
             else if (is.numeric(rhs) && length(rhs) == 1L) as.double(rhs)
             else stop("rhs must be a resident key or numeric scalar")
  .Call("amatrix_arrayfire_ewise_resident_bridge",
        as.character(lhs_key), rhs_arg, as.character(op), as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_rowSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::rowSums(amatrix_arrayfire_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  .Call("amatrix_arrayfire_sum_axis_resident_bridge", as.character(x_key), 0L)
}

amatrix_arrayfire_colSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::colSums(amatrix_arrayfire_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  .Call("amatrix_arrayfire_sum_axis_resident_bridge", as.character(x_key), 1L)
}

amatrix_arrayfire_solve_resident <- function(a_key, b_key = NULL, out_key) {
  a_host <- amatrix_arrayfire_resident_materialize(a_key)
  result <- if (is.null(b_key)) base::solve(a_host)
            else base::solve(a_host, amatrix_arrayfire_resident_materialize(b_key))
  amatrix_arrayfire_resident_store(out_key, result)
  result
}

amatrix_arrayfire_chol_resident <- function(x_key, out_key) {
  result <- base::chol(amatrix_arrayfire_resident_materialize(x_key))
  amatrix_arrayfire_resident_store(out_key, result)
  result
}

.amatrix_arrayfire_forced_available <- function() {
  isTRUE(getOption("amatrix.arrayfire.available", FALSE))
}

.amatrix_arrayfire_product_thresholds <- function() {
  list(
    matmul_min_dim = getOption("amatrix.arrayfire.matmul_min_dim", 512L),
    crossprod_min_dim = getOption("amatrix.arrayfire.crossprod_min_dim", 2048L),
    tcrossprod_min_dim = getOption("amatrix.arrayfire.tcrossprod_min_dim", 2048L),
    ewise_min_dim = getOption("amatrix.arrayfire.ewise_min_dim", 4096L),
    sum_min_dim = getOption("amatrix.arrayfire.sum_min_dim", 4096L),
    qr_min_dim = getOption("amatrix.arrayfire.qr_min_dim", 512L)
  )
}

.amatrix_arrayfire_meets_threshold <- function(x, threshold) {
  dims <- dim(x)
  !is.null(dims) && length(dims) == 2L && max(dims) >= threshold
}

amatrix_arrayfire_backend <- function() {
  cpu <- amatrix:::.amatrix_cpu_backend()
  capabilities <- amatrix_arrayfire_capabilities()
  features <- amatrix_arrayfire_features()
  precision_modes <- amatrix_arrayfire_precision_modes()
  thresholds <- .amatrix_arrayfire_product_thresholds()

  list(
    capabilities = function() {
      capabilities
    },
    features = function() {
      features
    },
    precision_modes = function() {
      precision_modes
    },
    available = function() {
      amatrix_arrayfire_is_available()
    },
    supports = function(op, x, y = NULL) {
      if (!is(x, "adgeMatrix") || !(op %in% capabilities)) {
        return(FALSE)
      }

      if (!(x@precision %in% precision_modes)) {
        return(FALSE)
      }

      if (.amatrix_arrayfire_forced_available()) {
        return(TRUE)
      }

      if (identical(op, "matmul")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$matmul_min_dim))
      }

      if (identical(op, "crossprod")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$crossprod_min_dim))
      }

      if (identical(op, "tcrossprod")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$tcrossprod_min_dim))
      }

      if (identical(op, "ewise")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$ewise_min_dim))
      }

      if (op %in% c("rowSums", "colSums")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$sum_min_dim))
      }

      if (identical(op, "qr")) {
        return(.amatrix_arrayfire_qr_safe() && .amatrix_arrayfire_meets_threshold(x, thresholds$qr_min_dim))
      }

      FALSE
    },
    matmul = function(x, y) {
      amatrix_arrayfire_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      amatrix_arrayfire_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      amatrix_arrayfire_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_arrayfire_ewise(lhs = lhs, rhs = rhs, op = op)
    },
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(cpu$rowSums(x, na.rm = na.rm, dims = dims))
      }
      amatrix_arrayfire_axis_sums(x, axis = 0L)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(cpu$colSums(x, na.rm = na.rm, dims = dims))
      }
      amatrix_arrayfire_axis_sums(x, axis = 1L)
    },
    qr = function(x, ...) {
      if (!.amatrix_arrayfire_qr_safe()) {
        return(cpu$qr(x, ...))
      }
      amatrix_arrayfire_qr(x)
    },
    resident_has = function(key) {
      amatrix_arrayfire_resident_has(key)
    },
    resident_store = function(key, x) {
      amatrix_arrayfire_resident_store(key, x)
    },
    resident_drop = function(key) {
      amatrix_arrayfire_resident_drop(key)
    },
    resident_materialize = function(key) {
      amatrix_arrayfire_resident_materialize(key)
    },
    matmul_resident = function(x_key, y_key, out_key) {
      amatrix_arrayfire_matmul_resident(x_key, y_key, out_key)
    },
    crossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_arrayfire_crossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    tcrossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_arrayfire_tcrossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    ewise_resident = function(lhs_key, rhs, op, out_key) {
      amatrix_arrayfire_ewise_resident(lhs_key, rhs, op, out_key)
    },
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_arrayfire_rowSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_arrayfire_colSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    solve_resident = function(a_key, b_key = NULL, out_key) {
      amatrix_arrayfire_solve_resident(a_key, b_key = b_key, out_key = out_key)
    },
    chol_resident = function(x_key, out_key) {
      amatrix_arrayfire_chol_resident(x_key, out_key)
    }
  )
}

amatrix_arrayfire_register <- function(overwrite = TRUE) {
  amatrix_register_backend("arrayfire", amatrix_arrayfire_backend(), overwrite = overwrite)
  invisible("arrayfire")
}
