.amatrix_cpu_backend <- function() {
  capabilities <- c("matmul", "crossprod", "tcrossprod", "ewise", "broadcast_ewise", "argmax", "scatter_mean", "rowSums", "colSums", "solve", "chol", "qr", "svd", "eigen", "diag")
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
      x %*% y
    },
    crossprod = function(x, y = NULL, ...) {
      if (is.null(y)) {
        return(base::crossprod(as.matrix(x), ...))
      }
      base::crossprod(as.matrix(x), y = .amatrix_host_arg(y), ...)
    },
    tcrossprod = function(x, y = NULL, ...) {
      if (is.null(y)) {
        return(base::tcrossprod(as.matrix(x), ...))
      }
      base::tcrossprod(as.matrix(x), y = .amatrix_host_arg(y), ...)
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
      if (is.null(b)) {
        return(base::solve(as.matrix(x), ...))
      }
      base::solve(as.matrix(x), as.matrix(.amatrix_host_arg(b)), ...)
    },
    chol = function(x, ...) {
      base::chol(as.matrix(x), ...)
    },
    qr = function(x, ...) {
      base::qr(as.matrix(x), ...)
    },
    svd = function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
      if (isTRUE(LINPACK)) {
        stop("LINPACK is not supported", call. = FALSE)
      }
      base::svd(as.matrix(x), nu = nu, nv = nv, ...)
    },
    eigen = function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
      base::eigen(as.matrix(x), symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
    },
    diag = function(x, nrow, ncol, names = TRUE) {
      base::diag(as.matrix(x), nrow = nrow, ncol = ncol, names = names)
    }
  )
}
