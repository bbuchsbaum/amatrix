.amatrix_cpu_dense_matrix <- function(x) {
  if (inherits(x, "adgeMatrix") || inherits(x, "dgeMatrix")) {
    return(.amatrix_dense_slot_matrix(x))
  }
  if (inherits(x, "denseMatrix")) {
    return(.amatrix_dense_slot_matrix(.amatrix_dense_base(x)))
  }
  if (is.matrix(x)) {
    if (!is.double(x)) {
      storage.mode(x) <- "double"
    }
    return(x)
  }
  if (inherits(x, "aMatrix")) {
    return(.amatrix_cpu_dense_matrix(amatrix_materialize_host(x)))
  }
  as.matrix(x)
}

.amatrix_cpu_host_value <- function(x) {
  if (inherits(x, "adgCMatrix")) {
    return(amatrix_materialize_host(x))
  }
  if (inherits(x, "aMatrix")) {
    return(.amatrix_cpu_dense_matrix(x))
  }
  x
}

.amatrix_cpu_backend <- function() {
  capabilities <- c("matmul", "crossprod", "tcrossprod", "ewise", "broadcast_ewise", "argmax", "scatter_mean", "segment_sum", "segment_mean", "rowSums", "colSums", "solve", "chol", "qr", "svd", "eigen", "diag")
  features <- c("dense_f64", "dense_f32", "solve", "chol", "svd", "sparse_spmm")

  list(
    capabilities = function() {
      capabilities
    },
    features = function() {
      features
    },
    precision_modes = function() {
      c("strict", "fast")
    },
    available = function() {
      TRUE
    },
    supports = function(op, x, y = NULL) {
      TRUE
    },
    matmul = function(x, y) {
      x <- .amatrix_cpu_host_value(x)
      y <- .amatrix_cpu_host_value(y)
      if (inherits(x, "sparseMatrix") && inherits(y, "sparseMatrix")) {
        result <- methods::as(x %*% y, "dgCMatrix")
        return(new_adgCMatrix(result))
      }
      if (inherits(x, "sparseMatrix") || inherits(y, "sparseMatrix"))
        return(as.matrix(x %*% y))
      .amatrix_cpu_dense_matrix(x) %*% .amatrix_cpu_dense_matrix(y)
    },
    crossprod = function(x, y = NULL, ...) {
      x <- .amatrix_cpu_host_value(x)
      y <- if (is.null(y)) NULL else .amatrix_cpu_host_value(y)
      if (inherits(x, "sparseMatrix")) {
        result <- if (is.null(y)) Matrix::crossprod(x) else Matrix::crossprod(x, y)
        return(as.matrix(result))
      }
      if (is.null(y)) {
        return(base::crossprod(.amatrix_cpu_dense_matrix(x), ...))
      }
      base::crossprod(.amatrix_cpu_dense_matrix(x), y = .amatrix_cpu_dense_matrix(y), ...)
    },
    tcrossprod = function(x, y = NULL, ...) {
      x <- .amatrix_cpu_host_value(x)
      y <- if (is.null(y)) NULL else .amatrix_cpu_host_value(y)
      if (inherits(x, "sparseMatrix")) {
        result <- if (is.null(y)) Matrix::tcrossprod(x) else Matrix::tcrossprod(x, y)
        return(as.matrix(result))
      }
      if (is.null(y)) {
        return(base::tcrossprod(.amatrix_cpu_dense_matrix(x), ...))
      }
      base::tcrossprod(.amatrix_cpu_dense_matrix(x), y = .amatrix_cpu_dense_matrix(y), ...)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      if (is.null(rhs)) {
        return(do.call(op, c(list(lhs), list(...))))
      }
      do.call(op, c(list(lhs, rhs), list(...)))
    },
    broadcast_ewise = function(x, lhs, v, margin, op, ...) {
      base::sweep(as.matrix(lhs), MARGIN = margin, STATS = v, FUN = op)
    },
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      Matrix::rowSums(x, na.rm = na.rm, dims = dims)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      Matrix::colSums(x, na.rm = na.rm, dims = dims)
    },
    solve = function(x, b = NULL, ...) {
      x <- .amatrix_cpu_host_value(x)
      b <- if (is.null(b)) NULL else .amatrix_cpu_host_value(b)
      if (inherits(x, "sparseMatrix")) {
        if (is.null(b)) return(Matrix::solve(x, ...))
        return(Matrix::solve(x, b, ...))
      }
      if (is.null(b)) {
        return(base::solve(.amatrix_cpu_dense_matrix(x), ...))
      }
      base::solve(.amatrix_cpu_dense_matrix(x), .amatrix_cpu_dense_matrix(b), ...)
    },
    chol = function(x, ...) {
      x <- .amatrix_cpu_host_value(x)
      if (inherits(x, "sparseMatrix")) return(Matrix::chol(x, ...))
      base::chol(.amatrix_cpu_dense_matrix(x), ...)
    },
    qr = function(x, ...) {
      x <- .amatrix_cpu_host_value(x)
      if (inherits(x, "sparseMatrix")) return(Matrix::qr(x, ...))
      base::qr(.amatrix_cpu_dense_matrix(x), ...)
    },
    svd = function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
      x <- .amatrix_cpu_host_value(x)
      if (isTRUE(LINPACK)) {
        stop("LINPACK is not supported", call. = FALSE)
      }
      base::svd(.amatrix_cpu_dense_matrix(x), nu = nu, nv = nv, ...)
    },
    eigen = function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
      base::eigen(.amatrix_cpu_dense_matrix(.amatrix_cpu_host_value(x)), symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
    },
    diag = function(x, nrow, ncol, names = TRUE) {
      x <- .amatrix_cpu_host_value(x)
      if (inherits(x, "sparseMatrix")) return(Matrix::diag(x))
      base::diag(.amatrix_cpu_dense_matrix(x), nrow = nrow, ncol = ncol, names = names)
    }
  )
}
