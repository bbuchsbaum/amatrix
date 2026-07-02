## Regression tests for dispatch/registry bugs found by dispatch-hunter.
## Each test is tagged with its amatrix issue ID.
## These tests are EXPECTED TO FAIL until the bugs are fixed.

.make_minimal_backend <- function(precision_modes = c("strict", "fast"),
                                  available = TRUE) {
  list(
    capabilities    = function() c("matmul"),
    features        = function() character(),
    precision_modes = function() precision_modes,
    available       = function() available,
    supports        = function(op, x, y = NULL) TRUE,
    matmul          = function(x, y) x,
    crossprod       = function(x, y = NULL) x,
    tcrossprod      = function(x, y = NULL) x,
    ewise           = function(x, y, op) x,
    rowSums         = function(x) numeric(nrow(x)),
    colSums         = function(x) numeric(ncol(x))
  )
}

# amatrix-or5: amatrix_register_backend accepts empty precision_modes silently,
# then backend is never dispatched because precision_compatible is always FALSE.
test_that("amatrix-or5: register_backend rejects empty precision_modes()", {
  be <- .make_minimal_backend(precision_modes = character())
  expect_error(
    amatrix_register_backend("test_empty_prec", be, overwrite = TRUE),
    regexp = "precision_modes"
  )
})

# amatrix-or5 corollary: a backend with empty precision_modes that somehow gets
# registered is never chosen as preferred backend, even when available.
test_that("amatrix-or5: backend with empty precision_modes is never dispatched", {
  be <- .make_minimal_backend(precision_modes = character())
  # Force registration by bypassing the public API (simulate the bug)
  assign("test_empty_prec2", be, envir = .amatrix_state$backends)
  on.exit(
    if (exists("test_empty_prec2", envir = .amatrix_state$backends, inherits = FALSE))
      rm("test_empty_prec2", envir = .amatrix_state$backends),
    add = TRUE
  )

  m <- new_adgeMatrix(matrix(1:4, 2, 2),
                      preferred_backend = "test_empty_prec2",
                      policy = "auto",
                      precision = "strict")
  plan <- amatrix_backend_plan(m, "matmul")
  # Backends with no declared precision support are invalid and must never win
  # dispatch, even if they are inserted into the registry out-of-band.
  expect_equal(plan$chosen, "cpu")
})

# amatrix-y85: 'opencl' is in .amatrix_auto_fast_backend_order() but not in
# .amatrix_valid_policies, so amatrix_set_default_policy(default_fast_backend())
# throws an error when opencl is the fastest registered backend.
test_that("amatrix-y85: opencl backend can be used as default policy", {
  opencl_be <- .make_minimal_backend()
  amatrix_register_backend("opencl", opencl_be, overwrite = TRUE)
  on.exit(
    if (exists("opencl", envir = .amatrix_state$backends, inherits = FALSE))
      rm("opencl", envir = .amatrix_state$backends),
    add = TRUE
  )

  fb <- .amatrix_default_fast_backend()
  if (identical(fb, "opencl")) {
    # Regression amatrix-y85: opencl must be accepted as a policy value.
    expect_no_error(amatrix_set_default_policy(fb))
  } else {
    skip("opencl was not selected as fastest backend in this session")
  }
})

# amatrix-y85 corollary: .amatrix_valid_policies must include every backend
# name that .amatrix_auto_fast_backend_order() can return.
test_that("amatrix-y85: all fast backend order names are valid policies", {
  fast_order <- .amatrix_auto_fast_backend_order()
  valid_pol  <- .amatrix_valid_policies
  missing    <- setdiff(fast_order, c(valid_pol, "auto", "cpu"))
  expect_equal(
    missing, character(0),
    label = paste("backends in fast order but not valid policies:", paste(missing, collapse = ", "))
  )
})
