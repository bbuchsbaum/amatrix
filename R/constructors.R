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

.amatrix_dense_slot_matrix <- function(x) {
  stopifnot(inherits(x, "denseMatrix") || inherits(x, "adgeMatrix"))
  out <- x@x
  dim(out) <- as.integer(x@Dim)
  dimnames(out) <- x@Dimnames
  storage.mode(out) <- "double"
  out
}

.amatrix_new_dense <- function(
  data,
  dim,
  dimnames = NULL,
  factors = list(),
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  if (!is.double(data)) {
    storage.mode(data) <- "double"
  }
  if (is.null(dimnames)) {
    dimnames <- vector("list", 2L)
  }
  object_id <- .amatrix_next_object_id()
  new(
    "adgeMatrix",
    x = as.double(data),
    Dim = as.integer(dim),
    Dimnames = dimnames,
    factors = factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    src_id = src_id,
    finalizer_env = .amatrix_make_finalizer_env(object_id)
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
  if (inherits(x, "dgeMatrix")) {
    return(.amatrix_new_dense(
      data = x@x,
      dim = x@Dim,
      dimnames = x@Dimnames,
      factors = x@factors,
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision,
      src_id = src_id
    ))
  }

  if (is.matrix(x)) {
    return(.amatrix_new_dense(
      data = x,
      dim = dim(x),
      dimnames = base::dimnames(x),
      factors = list(),
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision,
      src_id = src_id
    ))
  }

  base <- .amatrix_dense_base(x)
  .amatrix_new_dense(
    data = base@x,
    dim = base@Dim,
    dimnames = base@Dimnames,
    factors = base@factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    src_id = src_id
  )
}

# Deferred-materialization variant: GPU result exists as a resident key but
# the host copy is not yet downloaded.  @x holds a NaN sentinel (valid
# dgeMatrix structurally) and the actual data lives only on the device until
# the first host access triggers a transparent download.
new_adgeMatrix_deferred <- function(
  dim,
  dimnames = list(NULL, NULL),
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  n <- as.integer(dim[1L]) * as.integer(dim[2L])
  object_id <- .amatrix_next_object_id()
  fenv <- .amatrix_make_finalizer_env(object_id)
  fenv$host_deferred <- TRUE
  fenv$host_x <- NULL

  new(
    "adgeMatrix",
    x = rep(NaN, n),
    Dim = as.integer(dim),
    Dimnames = if (is.null(dimnames)) list(NULL, NULL) else dimnames,
    factors = list(),
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    src_id = src_id,
    finalizer_env = fenv
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
  # resident_handle → adgeMatrix with ownership transfer
  if (inherits(x, "resident_handle")) {
    return(as_adgeMatrix.resident_handle(x))
  }
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

# Internal constructor for zero-copy transpose views.
# Always creates a new object_id distinct from the source so caches and
# the residency registry remain independent.
.new_aTransposeView <- function(source) {
  oid <- .amatrix_next_object_id()
  dn <- source@Dimnames
  new(
    "aTransposeView",
    source           = source,
    Dim              = rev(source@Dim),
    Dimnames         = if (length(dn) == 2L) rev(dn) else list(NULL, NULL),
    preferred_backend = source@preferred_backend,
    policy           = source@policy,
    precision        = source@precision,
    object_id        = oid,
    src_id           = source@object_id,
    finalizer_env    = .amatrix_make_finalizer_env(oid)
  )
}

setAs("matrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("dgeMatrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("matrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("dgCMatrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("adgeMatrix", "dgeMatrix", function(from) new("dgeMatrix", x = from@x, Dim = from@Dim, Dimnames = from@Dimnames, factors = from@factors))
setAs("adgCMatrix", "dgCMatrix", function(from) new("dgCMatrix", i = from@i, p = from@p, Dim = from@Dim, Dimnames = from@Dimnames, x = from@x, factors = from@factors))
