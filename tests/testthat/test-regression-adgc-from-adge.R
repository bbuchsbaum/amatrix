# Regression repro for adgCMatrix(adgeMatrix).
# Seed: none; deterministic literal 3 x 4 matrix
# Shape: 3 x 4 dense input coerced to sparse
# Backend: cpu
# Precision mode: strict
# Dispatch path: adgCMatrix()/as_adgCMatrix() constructor
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-dum

test_that("adgCMatrix constructors accept adgeMatrix inputs", {
  host <- matrix(c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0), nrow = 3)
  dense <- adgeMatrix(host, preferred_backend = "cpu", precision = "strict")

  via_constructor <- adgCMatrix(dense)
  via_coercer <- as_adgCMatrix(dense)
  via_s4 <- methods::as(dense, "adgCMatrix")

  expect_s4_class(via_constructor, "adgCMatrix")
  expect_s4_class(via_coercer, "adgCMatrix")
  expect_equal(as.matrix(via_constructor), host)
  expect_equal(as.matrix(via_coercer), host)
  expect_equal(as.matrix(via_s4), host)
  expect_identical(via_constructor@preferred_backend, "cpu")
  expect_identical(via_constructor@precision, "strict")
})
