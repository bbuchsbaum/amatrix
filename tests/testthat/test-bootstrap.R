test_that("dense constructor creates Matrix-compatible subclass", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))

  expect_s4_class(x, "adgeMatrix")
  expect_s4_class(x, "aMatrix")
  expect_equal(dim(x), c(2, 2))
  expect_equal(x@preferred_backend, "cpu")
  expect_equal(x@policy, "auto")
  expect_equal(x@precision, "strict")
})

test_that("sparse constructor creates Matrix-compatible subclass", {
  x <- adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2))

  expect_s4_class(x, "adgCMatrix")
  expect_s4_class(x, "aMatrix")
  expect_equal(dim(x), c(2, 2))
  expect_equal(x@precision, "strict")
  expect_true(inherits(as(x, "dgCMatrix"), "dgCMatrix"))
})

test_that("precision defaults can be configured explicitly", {
  old <- amatrix_default_precision()
  on.exit(amatrix_set_default_precision(old), add = TRUE)

  amatrix_set_default_precision("fast")
  x <- adgeMatrix(matrix(1:4, nrow = 2))

  expect_equal(amatrix_default_precision(), "fast")
  expect_equal(x@precision, "fast")
})

test_that("session defaults can steer backend choice for plain constructors", {
  old_policy <- amatrix_default_policy()
  old_precision <- amatrix_default_precision()
  on.exit({
    amatrix_set_default_policy(old_policy)
    amatrix_set_default_precision(old_precision)
  }, add = TRUE)

  counter <- new.env(parent = emptyenv())

  with_registered_backend("mlx", make_recording_backend(counter, supported_ops = c("matmul"), precision_modes = "fast"), {
    amatrix_set_default_policy("mlx")
    amatrix_set_default_precision("fast")

    x <- adgeMatrix(matrix(1:4, nrow = 2))
    plan <- amatrix_backend_plan(x, "matmul", y = diag(2))

    expect_equal(x@preferred_backend, "mlx")
    expect_equal(x@policy, "mlx")
    expect_equal(x@precision, "fast")
    expect_identical(plan$preferred, c("mlx", "cpu"))
    expect_identical(plan$chosen, "mlx")
  })
})

test_that("fast mode prefers an available accelerator without session setters", {
  old_policy <- amatrix_default_policy()
  old_precision <- amatrix_default_precision()
  on.exit({
    amatrix_set_default_policy(old_policy)
    amatrix_set_default_precision(old_precision)
  }, add = TRUE)

  with_registered_backend("mlx", make_recording_backend(counter = new.env(parent = emptyenv()), supported_ops = c("matmul"), precision_modes = "fast"), {
    amatrix_set_default_policy("auto")
    amatrix_set_default_precision("strict")

    x <- adgeMatrix(matrix(1:4, nrow = 2), mode = "fast")
    plan <- amatrix_backend_plan(x, "matmul", y = diag(2))

    expect_equal(x@preferred_backend, "mlx")
    expect_equal(x@policy, "auto")
    expect_equal(x@precision, "fast")
    # preferred_candidates (R/policy.R) now filters "auto" out of the
    # preferred list so the resolved order is (preferred_backend, cpu).
    expect_identical(plan$preferred, c("mlx", "cpu"))
    expect_identical(plan$chosen, "mlx")
  })
})

test_that("with_amatrix restores defaults on normal exit and error", {
  old_policy <- amatrix_default_policy()
  old_precision <- amatrix_default_precision()
  on.exit({
    amatrix_set_default_policy(old_policy)
    amatrix_set_default_precision(old_precision)
  }, add = TRUE)

  amatrix_set_default_policy("auto")
  amatrix_set_default_precision("strict")

  x <- with_registered_backend("mlx", make_recording_backend(counter = new.env(parent = emptyenv()), supported_ops = c("matmul"), precision_modes = "fast"), {
    with_amatrix(policy = "auto", precision = "fast", {
      expect_equal(amatrix_default_policy(), "auto")
      expect_equal(amatrix_default_precision(), "fast")
      adgeMatrix(matrix(1:4, nrow = 2))
    })
  })

  expect_equal(x@preferred_backend, "mlx")
  expect_equal(x@precision, "fast")
  expect_equal(amatrix_default_policy(), "auto")
  expect_equal(amatrix_default_precision(), "strict")

  expect_error(
    with_amatrix(policy = "cpu", precision = "fast", {
      stop("boom")
    }),
    "boom"
  )

  expect_equal(amatrix_default_policy(), "auto")
  expect_equal(amatrix_default_precision(), "strict")
})

test_that("core operations route through fallback and preserve semantics", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))
  y <- x %*% diag(2)

  expect_s4_class(y, "adgeMatrix")
  expect_equal(as.matrix(y), matrix(1:4, nrow = 2))
  expect_equal(rowSums(x), c(4, 6))
  expect_equal(colSums(x), c(3, 7))
})

test_that("registered backends include cpu", {
  expect_true("cpu" %in% amatrix_backend_names())
})

test_that("subsetting preserves amatrix classes for matrix-like results", {
  x <- adgeMatrix(matrix(1:9, nrow = 3))
  y <- x[1:2, 1:2, drop = FALSE]

  expect_s4_class(y, "adgeMatrix")
  expect_equal(as.matrix(y), matrix(c(1, 2, 4, 5), nrow = 2))
  expect_equal(x[, 1], c(1, 2, 3))
})

test_that("cpu linear algebra methods match host semantics", {
  x <- adgeMatrix(matrix(c(4, 1, 1, 3), nrow = 2))
  chol_x <- chol(x)
  solve_x <- solve(x)
  svd_x <- svd(x)
  eigen_x <- eigen(x, symmetric = TRUE)

  expect_s4_class(chol_x, "adgeMatrix")
  expect_s4_class(solve_x, "adgeMatrix")
  expect_equal(as.matrix(chol_x), chol(matrix(c(4, 1, 1, 3), nrow = 2)))
  expect_equal(as.matrix(solve_x), solve(matrix(c(4, 1, 1, 3), nrow = 2)))
  expect_equal(svd_x$d, base::svd(matrix(c(4, 1, 1, 3), nrow = 2))$d)
  expect_equal(eigen_x$values, base::eigen(matrix(c(4, 1, 1, 3), nrow = 2), symmetric = TRUE)$values)
})

test_that("serialization preserves host semantics without backend state", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))
  roundtrip <- unserialize(serialize(x, NULL))

  expect_s4_class(roundtrip, "adgeMatrix")
  expect_equal(as.matrix(roundtrip), matrix(1:4, nrow = 2))
  expect_equal(roundtrip@preferred_backend, "cpu")
  expect_equal(roundtrip@policy, "auto")
  expect_equal(roundtrip@precision, "strict")
  expect_false("backend_key" %in% slotNames(roundtrip))
  expect_false("state_version" %in% slotNames(roundtrip))
})

test_that("Ops methods preserve numeric amatrix results and comparisons stay logical", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))
  y <- x + 1
  z <- 2 * x
  cmp <- x > 2

  expect_s4_class(y, "adgeMatrix")
  expect_s4_class(z, "adgeMatrix")
  expect_equal(as.matrix(y), matrix(2:5, nrow = 2))
  expect_equal(as.matrix(z), 2 * matrix(1:4, nrow = 2))
  expect_true(is.matrix(cmp) || inherits(cmp, "Matrix"))
  expect_s4_class(cmp, "adlgeMatrix")
  expect_identical(as.matrix(cmp), matrix(c(FALSE, FALSE, TRUE, TRUE), nrow = 2))
})

test_that("dimnames replacement and coercion preserve host semantics", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))
  dimnames(x) <- list(c("r1", "r2"), c("c1", "c2"))

  expect_identical(dimnames(x), list(c("r1", "r2"), c("c1", "c2")))
  expect_identical(rownames(as.matrix(x)), c("r1", "r2"))
  expect_identical(colnames(as.array(x)), c("c1", "c2"))
})

test_that("subassignment preserves amatrix classes and values", {
  x <- adgeMatrix(matrix(1:4, nrow = 2))
  x[1, 2] <- 99

  expect_s4_class(x, "adgeMatrix")
  expect_equal(as.matrix(x), matrix(c(1, 2, 99, 4), nrow = 2))

  s <- adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2))
  s[1, 2] <- 5
  expect_true(inherits(s, "aMatrix"))
  expect_equal(as.matrix(s), matrix(c(1, 0, 5, 1), nrow = 2))
})
