# Naming policy for the exported API.
#
#   amatrix_*   session, backend and diagnostic functions
#   am_*        everything that computes (kernels, factorizations, spectral
#               methods, model fitters)
#   class constructors / coercions and a few fixed names listed below
#
# The point is twofold: users can discover the API by completing two prefixes,
# and attaching amatrix can never mask a function from another package.

test_that("every exported function follows the naming policy", {
  exports <- getNamespaceExports("amatrix")
  # S4 generics exported for method dispatch on base/Matrix functions.
  generics <- methods::getGenerics(asNamespace("amatrix"))@.Data
  exports <- setdiff(exports, generics)
  exports <- exports[!grepl("^\\.__", exports)]

  allowed_fixed <- c(
    # classes, constructors and coercions
    "adgeMatrix", "adgCMatrix", "as_adgeMatrix", "as_adgCMatrix",
    "as_adgeMatrix.resident_handle", "kron_matrix", "resident_handle",
    # methods exported under their S3 names
    "as.matrix.adgeMatrix", "as.matrix.adgCMatrix", "as.matrix.aTransposeView",
    "as.matrix.KronMatrix", "as.matrix.resident_handle",
    "nrow.resident_handle", "ncol.resident_handle",
    # scoped settings helper
    "with_amatrix"
  )
  offenders <- exports[!grepl("^(amatrix_|am_)", exports) & !exports %in% allowed_fixed]
  expect_identical(
    sort(offenders), character(0),
    info = "new exports must start with amatrix_ or am_ (see the naming policy at the top of this file)"
  )
})

test_that("attaching amatrix masks nothing from base or the recommended packages", {
  exports <- getNamespaceExports("amatrix")
  generics <- methods::getGenerics(asNamespace("amatrix"))@.Data
  plain <- setdiff(exports, generics)
  plain <- plain[!grepl("^\\.", plain)]   # S4 method/class tables
  for (pkg in c("base", "stats", "utils", "methods", "graphics", "grDevices", "Matrix")) {
    clash <- intersect(plain, getNamespaceExports(pkg))
    expect_identical(clash, character(0), info = sprintf("masks %s::%s", pkg, paste(clash, collapse = ", ")))
  }
})

test_that("the downstream testing kit dispatches through a stand-in backend", {
  backend <- amatrix_test_backend()
  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  out <- amatrix_with_test_backend({
    x <- adgeMatrix(m, backend = "amatrix_test")
    list(prod = as.matrix(x %*% t(m)), cp = as.matrix(crossprod(x)), rs = rowSums(x))
  }, backend = backend)
  expect_equal(out$prod, m %*% t(m), ignore_attr = TRUE)
  expect_equal(out$cp, crossprod(m), ignore_attr = TRUE)
  expect_equal(out$rs, rowSums(m))
  expect_identical(backend$calls$matmul, 1L)
  expect_identical(backend$calls$crossprod, 1L)
  expect_identical(backend$calls$rowSums, 1L)
  # Registration and health state are restored afterwards.
  expect_false("amatrix_test" %in% amatrix_backend_names())
  expect_error(amatrix_with_test_backend(stop("boom")), "boom")
  expect_false("amatrix_test" %in% amatrix_backend_names())
  expect_type(amatrix_has_accelerator(), "logical")
})

test_that("mode = 'float32' is an alias for mode = 'fast'", {
  m <- matrix(1:6 + 0, 2, 3)
  a <- adgeMatrix(m, mode = "float32")
  b <- adgeMatrix(m, mode = "fast")
  expect_identical(a@precision, "fast")
  expect_identical(a@precision, b@precision)
  expect_identical(a@preferred_backend, b@preferred_backend)
  expect_identical(adgeMatrix(m, mode = "exact")@precision, "strict")
})
