.amatrix_build_dgeMatrix <- function(x, dn = base::dimnames(x)) {
  stopifnot(is.matrix(x))
  storage.mode(x) <- "double"
  if (is.null(dn)) {
    dn <- vector("list", 2L)
  }
  new(
    "dgeMatrix",
    x = as.double(x),
    Dim = as.integer(dim(x)),
    Dimnames = dn,
    factors = list()
  )
}

.amatrix_dense_base <- function(x) {
  if (inherits(x, "dgeMatrix")) {
    return(x)
  }
  if (inherits(x, "denseMatrix")) {
    return(.amatrix_build_dgeMatrix(as.matrix(x), dn = base::dimnames(x)))
  }
  if (is.matrix(x)) {
    return(.amatrix_build_dgeMatrix(x))
  }
  stop("x must be a base matrix or dgeMatrix")
}

.amatrix_sparse_base <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (inherits(x, "sparseMatrix")) {
    base <- x
    if (!inherits(base, "generalMatrix")) {
      base <- as(base, "generalMatrix")
    }
    return(as(base, "dgCMatrix"))
  }
  if (is.matrix(x)) {
    base <- Matrix(x, sparse = TRUE)
    if (!inherits(base, "generalMatrix")) {
      base <- as(base, "generalMatrix")
    }
    return(as(base, "dgCMatrix"))
  }
  stop("x must be a base matrix or dgCMatrix")
}

new_adgeMatrix <- function(
  x,
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  base <- .amatrix_dense_base(x)
  object_id <- .amatrix_next_object_id()
  new(
    "adgeMatrix",
    x = base@x,
    Dim = base@Dim,
    Dimnames = base@Dimnames,
    factors = base@factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    src_id = src_id,
    finalizer_env = .amatrix_make_finalizer_env(object_id)
  )
}

new_adgCMatrix <- function(
  x,
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision()
) {
  base <- .amatrix_sparse_base(x)
  object_id <- .amatrix_next_object_id()
  new(
    "adgCMatrix",
    i = base@i,
    p = base@p,
    Dim = base@Dim,
    Dimnames = base@Dimnames,
    x = base@x,
    factors = base@factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    finalizer_env = .amatrix_make_finalizer_env(object_id)
  )
}

adgeMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgeMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

adgCMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgCMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

as_adgeMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgeMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

as_adgCMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgCMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

setAs("matrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("dgeMatrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("matrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("dgCMatrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("adgeMatrix", "dgeMatrix", function(from) new("dgeMatrix", x = from@x, Dim = from@Dim, Dimnames = from@Dimnames, factors = from@factors))
setAs("adgCMatrix", "dgCMatrix", function(from) new("dgCMatrix", i = from@i, p = from@p, Dim = from@Dim, Dimnames = from@Dimnames, x = from@x, factors = from@factors))
