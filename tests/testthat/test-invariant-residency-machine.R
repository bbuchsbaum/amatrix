# Invariant-driven residency state-machine checks.
#
# Uses a mock residency-capable backend to exercise resident-handle transitions
# without depending on a real accelerator. The checks focus on registry
# rebinding, alias visibility, cache invalidation, and ownership transfer.

suppressPackageStartupMessages(library(amatrix))

.with_resident_case <- function(name, backend, code) {
  with_registered_backend(name, backend, force(code))
}

.expect_resident_alias_state <- function(x, h, expected, key = NULL, info = NULL) {
  if (!is.null(key)) {
    expect_identical(h$resident_key, key, info = info)
    expect_identical(amatrix:::.amatrix_resident_key(x), key, info = info)
  }
  expect_equal(amatrix:::as.matrix.resident_handle(h), expected, tolerance = 1e-12, info = info)
  expect_equal(as.matrix(x), expected, tolerance = 1e-12, info = info)
}

test_that("resident handle replacement-path sweep rebinds aliases to the new key", {
  counter <- new.env(parent = emptyenv())
  backend <- .invariant_mock_resident_backend(counter, inplace = FALSE)

  .with_resident_case("mock-machine-rebind", backend, {
    host <- matrix(1:9, nrow = 3L, ncol = 3L)
    x <- adgeMatrix(host, preferred_backend = "mock-machine-rebind")
    x <- amatrix_bind_resident(x, backend = "mock-machine-rebind")
    h <- resident_handle(x, backend = "mock-machine-rebind")

    key0 <- h$resident_key
    expected <- sweep(host, 1L, c(1, 2, 3), "*")
    am_sweep_inplace(h, 1L, c(1, 2, 3), "*")
    key1 <- h$resident_key

    expect_false(identical(key1, key0), info = "replacement-path sweep must allocate a new key")
    expect_false(backend$resident_has(key0), info = "old resident key must be dropped after rebind")
    .expect_resident_alias_state(x, h, expected, key = key1, info = "post-sweep alias state")
    expect_false(isTRUE(h$owns_key), info = "handle should not claim ownership while aliases share the key")

    expected <- expected + 5
    am_ewise_inplace(h, 5, "+")
    key2 <- h$resident_key

    expect_false(identical(key2, key1), info = "ewise replacement path must rotate the key")
    expect_false(backend$resident_has(key1), info = "prior resident key must be dropped after ewise rebind")
    .expect_resident_alias_state(x, h, expected, key = key2, info = "post-ewise alias state")
  })
})

test_that("resident handle in-place sweep invalidates alias host cache without rebinding", {
  counter <- new.env(parent = emptyenv())
  backend <- .invariant_mock_resident_backend(counter, inplace = TRUE)

  .with_resident_case("mock-machine-inplace", backend, {
    host <- matrix(seq_len(12L), nrow = 3L, ncol = 4L)
    x <- adgeMatrix(host, preferred_backend = "mock-machine-inplace")
    x <- amatrix_bind_resident(x, backend = "mock-machine-inplace")
    h <- resident_handle(x, backend = "mock-machine-inplace")

    key0 <- h$resident_key
    expect_equal(as.matrix(x), host, tolerance = 1e-12, info = "baseline host materialization")

    expected <- sweep(host, 2L, c(10, 20, 30, 40), "+")
    am_sweep_inplace(h, 2L, c(10, 20, 30, 40), "+")

    expect_identical(h$resident_key, key0, info = "in-place backend should keep the same resident key")
    expect_identical(amatrix:::.amatrix_resident_key(x), key0, info = "alias should keep the shared key")
    .expect_resident_alias_state(x, h, expected, key = key0, info = "in-place alias state")
    expect_false(isTRUE(h$owns_key), info = "shared in-place key should remain alias-owned")
  })
})

test_that("resident handle ownership transfers cleanly to adgeMatrix", {
  counter <- new.env(parent = emptyenv())
  backend <- .invariant_mock_resident_backend(counter, inplace = FALSE)

  .with_resident_case("mock-machine-transfer", backend, {
    host <- matrix(rnorm(9L), nrow = 3L, ncol = 3L)
    x <- adgeMatrix(host, preferred_backend = "mock-machine-transfer")
    x <- amatrix_bind_resident(x, backend = "mock-machine-transfer")
    h <- resident_handle(x, backend = "mock-machine-transfer")
    key0 <- h$resident_key

    y <- amatrix:::as_adgeMatrix.resident_handle(h, defer_host = TRUE)
    expect_false(isTRUE(h$active), info = "handle becomes inert after ownership transfer")
    expect_false(isTRUE(h$owns_key), info = "handle must relinquish key ownership after transfer")
    expect_null(h$resident_key, info = "handle no longer points at the resident key")
    expect_identical(amatrix:::.amatrix_resident_key(y), key0, info = "transferred object inherits resident key")
    expect_true(backend$resident_has(key0), info = "resident key stays live after transfer")
    expect_equal(as.matrix(y), host, tolerance = 1e-12, info = "transferred object materializes correctly")

    drops_before <- if (is.null(counter$resident_drop)) 0L else counter$resident_drop
    rm(h)
    gc(verbose = FALSE)
    drops_after <- if (is.null(counter$resident_drop)) 0L else counter$resident_drop
    expect_identical(drops_after, drops_before, info = "garbage collecting inert handle must not drop transferred key")
    expect_true(backend$resident_has(key0), info = "key remains live after handle GC")
  })
})
