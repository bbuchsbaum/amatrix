.invariant_dense_fixture <- function(
  seed = 20260414L,
  backend = "cpu",
  policy = "auto",
  precision = "strict"
) {
  set.seed(seed)
  host <- matrix(rnorm(12L), nrow = 3L, ncol = 4L)
  dimnames(host) <- list(
    paste0("r", seq_len(nrow(host))),
    paste0("c", seq_len(ncol(host)))
  )

  list(
    host = host,
    dge = Matrix::Matrix(host, sparse = FALSE),
    am = amatrix:::new_adgeMatrix(
      host,
      preferred_backend = backend,
      policy = policy,
      precision = precision
    )
  )
}

.invariant_sparse_fixture <- function(
  backend = "cpu",
  policy = "auto",
  precision = "strict"
) {
  host <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 1L, 3L),
    j = c(1L, 2L, 3L, 4L, 4L),
    x = c(1.25, -2, 3.5, 0.5, -1.75),
    dims = c(3L, 4L),
    dimnames = list(
      paste0("sr", seq_len(3L)),
      paste0("sc", seq_len(4L))
    )
  )
  host <- methods::as(host, "dgCMatrix")

  list(
    host = host,
    matrix = as.matrix(host),
    dgC = host,
    am = amatrix:::new_adgCMatrix(
      host,
      preferred_backend = backend,
      policy = policy,
      precision = precision
    )
  )
}

.invariant_template <- function(e1, e2 = NULL) {
  if (inherits(e1, "aMatrix")) {
    return(e1)
  }
  if (inherits(e2, "aMatrix")) {
    return(e2)
  }
  NULL
}

.invariant_expected_matrix_class <- function(template, host_value) {
  if (is.null(template)) {
    return(NULL)
  }
  if (!(is.matrix(host_value) || inherits(host_value, "Matrix"))) {
    return(NULL)
  }
  if (inherits(template, "adgeMatrix")) {
    return("adgeMatrix")
  }
  if (inherits(template, "adgCMatrix")) {
    if (inherits(host_value, "sparseMatrix")) {
      return("adgCMatrix")
    }
    return("adgeMatrix")
  }
  NULL
}

.invariant_expect_template_metadata <- function(result, template, info = NULL) {
  expect_identical(result@preferred_backend, template@preferred_backend, info = info)
  expect_identical(result@policy, template@policy, info = info)
  expect_identical(result@precision, template@precision, info = info)
}

.invariant_expect_wrapped_result <- function(
  result,
  template,
  host_reference,
  tolerance = 1e-12,
  info = NULL
) {
  expected_class <- .invariant_expected_matrix_class(template, host_reference)
  if (!is.null(expected_class)) {
    expect_true(methods::is(result, expected_class), info = info)
    .invariant_expect_template_metadata(result, template, info = info)
    expect_equal(
      base::dimnames(result),
      base::dimnames(host_reference),
      info = info
    )
    expect_equal(as.matrix(result), as.matrix(host_reference), tolerance = tolerance, info = info)
    return(invisible(NULL))
  }

  if (is.matrix(host_reference) || inherits(host_reference, "Matrix")) {
    expect_equal(as.matrix(result), as.matrix(host_reference), tolerance = tolerance, info = info)
  } else {
    expect_equal(result, host_reference, tolerance = tolerance, info = info)
  }
}

.invariant_mock_resident_backend <- function(counter, inplace = FALSE) {
  backend <- make_recording_backend(counter)

  backend$broadcast_ewise_resident <- function(lhs_key, stats, margin, fun, out_key, defer = FALSE) {
    mat <- backend$resident_materialize(lhs_key)
    value <- sweep(mat, margin, stats, fun)
    backend$resident_store(out_key, value)
    invisible(out_key)
  }

  if (isTRUE(inplace)) {
    backend$broadcast_ewise_resident_inplace <- function(lhs_key, stats, margin, fun) {
      mat <- backend$resident_materialize(lhs_key)
      value <- sweep(mat, margin, stats, fun)
      backend$resident_store(lhs_key, value)
      invisible(lhs_key)
    }
  }

  backend$rowSums_resident <- function(key) {
    base::rowSums(backend$resident_materialize(key))
  }

  backend$colSums_resident <- function(key) {
    base::colSums(backend$resident_materialize(key))
  }

  backend
}
