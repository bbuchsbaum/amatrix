# Regression repro metadata
# Seed: none (deterministic dense/sparse logical compare fixtures)
# Dimensions: dense 3 x 4; sparse 3 x 4 with named rows/cols
# Backend / precision / dispatch: cpu / strict / in-session Ops compare path
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-ol8, amatrix-jzk.1

test_that("compare ops preserve an amatrix logical wrapper for dense matrices [amatrix-ol8]", {
  x_host <- matrix(1:12, nrow = 3,
                   dimnames = list(c("r1", "r2", "r3"),
                                   c("c1", "c2", "c3", "c4")))
  x <- new_adgeMatrix(x_host, preferred_backend = "cpu", precision = "strict", policy = "auto")

  gt <- x > 5
  eq <- x == 5
  ne <- x != 5
  lt <- x < 10
  neg <- !(x > 5)
  band <- gt & lt
  bor <- gt | (x == 2)

  expect_true(inherits(gt, "aMatrix"))
  expect_s4_class(gt, "adlgeMatrix")
  expect_true(inherits(eq, "aMatrix"))
  expect_s4_class(eq, "adlgeMatrix")
  expect_s4_class(ne, "adlgeMatrix")
  expect_s4_class(band, "adlgeMatrix")
  expect_s4_class(bor, "adlgeMatrix")
  expect_true(inherits(neg, "aMatrix"))
  expect_s4_class(neg, "adlgeMatrix")

  expect_identical(dim(gt), c(3L, 4L))
  expect_identical(rownames(gt), rownames(x_host))
  expect_identical(colnames(gt), colnames(x_host))
  expect_identical(gt@preferred_backend, "cpu")
  expect_identical(gt@precision, "strict")
  expect_identical(gt@policy, "auto")

  expect_equal(as.matrix(gt), x_host > 5, tolerance = 0)
  expect_equal(as.matrix(eq), x_host == 5, tolerance = 0)
  expect_equal(as.matrix(ne), x_host != 5, tolerance = 0)
  expect_equal(as.matrix(band), (x_host > 5) & (x_host < 10), tolerance = 0)
  expect_equal(as.matrix(bor), (x_host > 5) | (x_host == 2), tolerance = 0)
  expect_equal(as.matrix(neg), !(x_host > 5), tolerance = 0)
})

test_that("compare ops preserve an amatrix logical wrapper for sparse matrices [amatrix-ol8]", {
  x_host <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 3),
    j = c(1, 2, 3, 4, 4),
    x = c(1, 2, 3, 0.5, 1.5),
    dims = c(3, 4),
    dimnames = list(c("r1", "r2", "r3"),
                    c("c1", "c2", "c3", "c4"))
  )
  x <- new_adgCMatrix(as(x_host, "dgCMatrix"), preferred_backend = "cpu", precision = "strict", policy = "auto")

  gt <- x > 1
  ne <- x != 1
  lt <- x < 3
  neg <- !(x > 1)
  band <- gt & lt
  bor <- gt | (x == 0)

  expect_true(inherits(gt, "aMatrix"))
  expect_s4_class(gt, "adlgCMatrix")
  expect_s4_class(ne, "adlgeMatrix")
  expect_s4_class(band, "adlgCMatrix")
  expect_s4_class(bor, "adlgCMatrix")
  expect_true(inherits(neg, "aMatrix"))
  expect_s4_class(neg, "adlgeMatrix")

  expect_identical(dim(gt), c(3L, 4L))
  expect_identical(rownames(gt), rownames(x_host))
  expect_identical(colnames(gt), colnames(x_host))
  expect_identical(gt@preferred_backend, "cpu")
  expect_identical(gt@precision, "strict")
  expect_identical(gt@policy, "auto")

  expect_equal(as.matrix(gt), as.matrix(x_host > 1), tolerance = 0)
  expect_equal(as.matrix(ne), as.matrix(x_host != 1), tolerance = 0)
  expect_equal(as.matrix(band), as.matrix((x_host > 1) & (x_host < 3)), tolerance = 0)
  expect_equal(as.matrix(bor), as.matrix((x_host > 1) | (x_host == 0)), tolerance = 0)
  expect_equal(as.matrix(neg), as.matrix(!(x_host > 1)), tolerance = 0)
})
