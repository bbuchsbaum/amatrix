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

# Arithmetic on plain dgeMatrix operands, computed directly on the data slots.
# Matrix's own Arith methods build and validate a new S4 object on every call
# (~0.35 ms), which dominates small elementwise operations. Returns NULL for
# anything outside the simple cases so the caller falls back to Matrix.
.amatrix_cpu_dense_arith <- function(lhs, rhs, op, ...) {
  if (is.null(rhs) || ...length() > 0L || !is.character(op) || length(op) != 1L ||
      !op %in% c("+", "-", "*", "/", "^")) {
    return(NULL)
  }
  plain_dge <- function(v) identical(class(v), structure("dgeMatrix", package = "Matrix"))
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.null(attributes(v))
  if (plain_dge(lhs) && plain_dge(rhs)) {
    if (!identical(lhs@Dim, rhs@Dim) || !identical(lhs@Dimnames, rhs@Dimnames)) {
      return(NULL)
    }
    return(.amatrix_plain_dge(do.call(op, list(lhs@x, rhs@x)), lhs@Dim, lhs@Dimnames))
  }
  if (plain_dge(lhs) && scalar(rhs)) {
    return(.amatrix_plain_dge(as.double(do.call(op, list(lhs@x, rhs))), lhs@Dim, lhs@Dimnames))
  }
  if (scalar(lhs) && plain_dge(rhs)) {
    return(.amatrix_plain_dge(as.double(do.call(op, list(lhs, rhs@x))), rhs@Dim, rhs@Dimnames))
  }
  NULL
}

# Large dense products without operand copies.
#
# The ordinary path turns each adgeMatrix into a base matrix first, which
# duplicates its data (setting `dim` on the shared @x vector), and the base
# result is duplicated again when it is wrapped. For big matrices that is
# 3-4x the memory of the operation itself. Matrix's own dgeMatrix methods
# call BLAS on the @x slots in place and return a dgeMatrix whose @x can be
# adopted as is, so use them once the data is large enough for the copies to
# outweigh Matrix's fixed S4 overhead. Returns NULL when not applicable.
.amatrix_cpu_large_dense_product <- function(kind, x, y = NULL, ...) {
  if (...length() > 0L) return(NULL)
  as_dge <- function(v) {
    if (inherits(v, "adgeMatrix")) return(amatrix_materialize_dense(v))
    if (identical(class(v), structure("dgeMatrix", package = "Matrix"))) return(v)
    NULL
  }
  xh <- as_dge(x)
  if (is.null(xh) || length(xh@x) < 250000L) return(NULL)
  yh <- NULL
  if (!is.null(y)) {
    yh <- as_dge(y)
    if (is.null(yh)) return(NULL)
  }
  res <- switch(kind,
    matmul = if (is.null(yh)) NULL else xh %*% yh,
    crossprod = if (is.null(yh)) Matrix::crossprod(xh) else Matrix::crossprod(xh, yh),
    tcrossprod = if (is.null(yh)) Matrix::tcrossprod(xh) else Matrix::tcrossprod(xh, yh),
    NULL
  )
  if (is.null(res)) return(NULL)
  if (!identical(class(res), structure("dgeMatrix", package = "Matrix"))) {
    # crossprod(x) comes back as a symmetric class; expand to general storage.
    res <- methods::as(methods::as(res, "generalMatrix"), "unpackedMatrix")
    if (!inherits(res, "dgeMatrix")) return(NULL)
  }
  res
}

# A Gram matrix of a sparse operand is returned dense (it usually is, and the
# model code expects a dense p x p matrix) unless that would be both large and
# mostly zeros, in which case densifying it is what runs out of memory.
.amatrix_cpu_sparse_product_value <- function(result) {
  if (inherits(result, "sparseMatrix")) {
    cells <- as.double(nrow(result)) * as.double(ncol(result))
    if (cells > 1e7 && Matrix::nnzero(result) < 0.25 * cells) {
      return(new_adgCMatrix(
        methods::as(methods::as(methods::as(result, "dMatrix"), "generalMatrix"), "CsparseMatrix")
      ))
    }
  }
  as.matrix(result)
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
      big <- .amatrix_cpu_large_dense_product("matmul", x, y)
      if (!is.null(big)) return(big)
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
      big <- .amatrix_cpu_large_dense_product("crossprod", x, y, ...)
      if (!is.null(big)) return(big)
      x <- .amatrix_cpu_host_value(x)
      y <- if (is.null(y)) NULL else .amatrix_cpu_host_value(y)
      if (inherits(x, "sparseMatrix")) {
        result <- if (is.null(y)) Matrix::crossprod(x) else Matrix::crossprod(x, y)
        return(.amatrix_cpu_sparse_product_value(result))
      }
      if (is.null(y)) {
        return(base::crossprod(.amatrix_cpu_dense_matrix(x), ...))
      }
      base::crossprod(.amatrix_cpu_dense_matrix(x), y = .amatrix_cpu_dense_matrix(y), ...)
    },
    tcrossprod = function(x, y = NULL, ...) {
      big <- .amatrix_cpu_large_dense_product("tcrossprod", x, y, ...)
      if (!is.null(big)) return(big)
      x <- .amatrix_cpu_host_value(x)
      y <- if (is.null(y)) NULL else .amatrix_cpu_host_value(y)
      if (inherits(x, "sparseMatrix")) {
        result <- if (is.null(y)) Matrix::tcrossprod(x) else Matrix::tcrossprod(x, y)
        return(.amatrix_cpu_sparse_product_value(result))
      }
      if (is.null(y)) {
        return(base::tcrossprod(.amatrix_cpu_dense_matrix(x), ...))
      }
      base::tcrossprod(.amatrix_cpu_dense_matrix(x), y = .amatrix_cpu_dense_matrix(y), ...)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      fast <- .amatrix_cpu_dense_arith(lhs, rhs, op, ...)
      if (!is.null(fast)) {
        return(fast)
      }
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
