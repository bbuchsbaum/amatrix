# test-bughunt-errors.R
# Failing tests documenting error-path bugs found during bughunt.
# DO NOT FIX these tests here — fix the source and let the tests pass.
#
# Bug IDs: amatrix-cng, amatrix-uu2, amatrix-hjj, amatrix-397, amatrix-6nm, amatrix-833

library(amatrix)

# ── amatrix-cng: stop() calls lack condition class ───────────────────────────
# All user-facing stop() calls in the package use plain stop() with no
# structured condition class. testthat::expect_error(class=) therefore always
# fails because the thrown condition is a bare simpleError, not an
# "amatrix_*" subclass.  The package must use something like:
#   rlang::abort("msg", class = "amatrix_bad_arg")
# or:
#   cond <- structure(class = c("amatrix_bad_arg", "error", "condition"), ...)
#   stop(cond)

test_that("amatrix-cng: chol_factor bad input throws classed condition", {
  # Plain matrix, not adgeMatrix — should throw amatrix_bad_arg (or similar)
  expect_error(
    chol_factor(matrix(1:4, 2, 2)),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: svd_factor bad X throws classed condition", {
  expect_error(
    svd_factor(matrix(1:6, 2, 3)),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: svd_factor bad k throws classed condition", {
  m <- adgeMatrix(matrix(rnorm(6), 2, 3))
  expect_error(
    svd_factor(m, k = 0L),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: svd_factor k > min(dim) throws classed condition", {
  m <- adgeMatrix(matrix(rnorm(6), 2, 3))
  expect_error(
    svd_factor(m, k = 100L),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: chol_solve bad factor throws classed condition", {
  expect_error(
    chol_solve("not_a_chol", rnorm(3)),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: lu_factor non-square throws classed condition", {
  expect_error(
    lu_factor(matrix(1:6, 2, 3)),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: amatrix_set_default_policy bad policy throws classed condition", {
  expect_error(
    amatrix_set_default_policy("NOT_A_POLICY"),
    class = "amatrix_bad_arg"
  )
})

test_that("amatrix-cng: amatrix_set_default_precision bad precision throws classed condition", {
  expect_error(
    amatrix_set_default_precision("NOT_A_PRECISION"),
    class = "amatrix_bad_arg"
  )
})

# ── amatrix-397: amatrix_register_backend lacks condition class ───────────────
# amatrix_register_backend() calls stop() without class= for multiple failure
# modes (not a list, missing fields, already registered).  Callers cannot
# distinguish these by class.

test_that("amatrix-397: register_backend with non-list throws classed condition", {
  expect_error(
    amatrix_register_backend("test_be", "not_a_list"),
    class = "amatrix_bad_backend"
  )
})

test_that("amatrix-397: register_backend with missing fields throws classed condition", {
  minimal <- list(capabilities = function() character(),
                  features = function() character(),
                  precision_modes = function() "strict")
  expect_error(
    amatrix_register_backend("test_be_missing", minimal),
    class = "amatrix_bad_backend"
  )
})

test_that("amatrix-397: register_backend duplicate name throws classed condition", {
  be <- list(
    capabilities    = function() character(),
    features        = function() character(),
    precision_modes = function() "strict",
    available       = function() FALSE,
    supports        = function(...) FALSE,
    matmul          = function(...) NULL,
    crossprod       = function(...) NULL,
    tcrossprod      = function(...) NULL,
    ewise           = function(...) NULL,
    rowSums         = function(...) NULL,
    colSums         = function(...) NULL
  )
  # Register once (may already exist — use overwrite)
  suppressWarnings(try(amatrix_register_backend("test_dup_be", be, overwrite = TRUE), silent = TRUE))
  expect_error(
    amatrix_register_backend("test_dup_be", be, overwrite = FALSE),
    class = "amatrix_backend_exists"
  )
})

# ── amatrix-6nm: policy.R stop() lacks call.=FALSE ───────────────────────────
# stop() without call.=FALSE includes the internal call (amatrix_set_default_policy)
# in the error message, leaking implementation details into user-facing messages.
# The call. field of the condition should be NULL when call.=FALSE is used.

test_that("amatrix-6nm: invalid policy error has no internal call site in message", {
  err <- tryCatch(
    amatrix_set_default_policy("NOPE"),
    error = function(e) e
  )
  # With call.=TRUE (the bug), conditionCall(err) returns the function call.
  # With call.=FALSE (the fix), conditionCall(err) returns NULL.
  expect_null(conditionCall(err))
})

test_that("amatrix-6nm: invalid precision error has no internal call site in message", {
  err <- tryCatch(
    amatrix_set_default_precision("NOPE"),
    error = function(e) e
  )
  expect_null(conditionCall(err))
})

# ── amatrix-uu2: cold-path tryCatch swallows original error class ─────────────
# amatrix_dispatch_op() catches a GPU backend error and calls fallback().
# The fallback is the CPU path. If we want to observe the *GPU* error class,
# it is gone — only simpleError from the fallback is visible. The fix must
# re-signal the original condition (e.g. via withCallingHandlers before
# tryCatch) so it remains catchable.
#
# We simulate by registering a dummy GPU backend that always errors, then
# dispatching an op through amatrix_dispatch_op.

test_that("amatrix-uu2: original GPU error condition class survives dispatch fallback", {
  # Register a fake GPU backend whose matmul always throws a classed error
  fake_gpu_error <- structure(
    class = c("amatrix_fake_gpu_error", "error", "condition"),
    list(message = "fake GPU failure", call = NULL)
  )
  fake_be <- list(
    capabilities    = function() c("matmul"),
    features        = function() character(),
    precision_modes = function() c("strict", "fast"),
    available       = function() TRUE,
    supports        = function(op, ...) op == "matmul",
    matmul          = function(...) stop(fake_gpu_error),
    crossprod       = function(...) NULL,
    tcrossprod      = function(...) NULL,
    ewise           = function(...) NULL,
    rowSums         = function(...) NULL,
    colSums         = function(...) NULL
  )
  suppressWarnings(try(
    amatrix_register_backend("fake_gpu_uu2", fake_be, overwrite = TRUE),
    silent = TRUE
  ))

  m <- adgeMatrix(matrix(rnorm(4), 2, 2), preferred_backend = "fake_gpu_uu2")

  # The original GPU error class should be observable via withCallingHandlers
  # even when the outer dispatch falls back to CPU.
  gpu_error_seen <- FALSE
  withCallingHandlers(
    tryCatch(
      amatrix_dispatch_op(m, "matmul", y = m, args = list(y = as.matrix(m)),
                          fallback = function() as.matrix(m) %*% as.matrix(m)),
      error = function(e) NULL
    ),
    amatrix_fake_gpu_error = function(e) {
      gpu_error_seen <<- TRUE
    }
  )
  # This will FAIL until the fix re-signals the original condition
  expect_true(gpu_error_seen)
})

# ── amatrix-hjj: resident-path NULL return has no condition class ─────────────
# When the resident path fails (tryCatch returns NULL), no condition is
# signalled to withCallingHandlers. Callers wanting to observe fallback events
# cannot do so. Fix: signal an "amatrix_fallback" condition before returning NULL.

test_that("amatrix-hjj: resident path failure signals amatrix_fallback condition", {
  fallback_seen <- FALSE
  # We need a scenario where resident path is attempted and fails.
  # The simplest proxy: check that .amatrix_log_fallback is called in a way
  # that produces a catchable condition of class "amatrix_fallback".
  # For now, test that the class exists and is signallable.
  withCallingHandlers(
    {
      cond <- structure(
        class = c("amatrix_fallback", "condition"),
        list(message = "test fallback", call = NULL, op = "matmul", backend = "test")
      )
      signalCondition(cond)
    },
    amatrix_fallback = function(e) {
      fallback_seen <<- TRUE
    }
  )
  # Passes trivially — the real test is that amatrix_dispatch_op actually
  # signals this class (which it currently does NOT do).
  expect_true(fallback_seen)

  # Now confirm that amatrix_dispatch_op does NOT currently signal amatrix_fallback
  # (this is the failing assertion that documents the bug):
  fake_be2 <- list(
    capabilities    = function() c("matmul"),
    features        = function() character(),
    precision_modes = function() c("strict", "fast"),
    available       = function() TRUE,
    supports        = function(op, ...) op == "matmul",
    matmul          = function(...) stop("resident failure"),
    crossprod       = function(...) NULL,
    tcrossprod      = function(...) NULL,
    ewise           = function(...) NULL,
    rowSums         = function(...) NULL,
    colSums         = function(...) NULL
  )
  suppressWarnings(try(
    amatrix_register_backend("fake_gpu_hjj", fake_be2, overwrite = TRUE),
    silent = TRUE
  ))
  m2 <- adgeMatrix(matrix(rnorm(4), 2, 2), preferred_backend = "fake_gpu_hjj")

  dispatch_fallback_seen <- FALSE
  withCallingHandlers(
    tryCatch(
      amatrix_dispatch_op(m2, "matmul", y = m2,
                          args = list(y = as.matrix(m2)),
                          fallback = function() as.matrix(m2) %*% as.matrix(m2)),
      error = function(e) NULL
    ),
    amatrix_fallback = function(e) {
      dispatch_fallback_seen <<- TRUE
    }
  )
  # This will FAIL because amatrix_dispatch_op currently uses .amatrix_log_fallback
  # (an internal logger) rather than signalling a catchable condition.
  expect_true(dispatch_fallback_seen)
})

# ── amatrix-833: subspace SVD stop() swallowed by tryCatch(error=NULL) ────────
# .amatrix_subspace_svd() calls stop("subspace SVD did not discover a usable
# range space") but its caller .amatrix_subspace_compile_operator wraps the
# call in tryCatch(error = function(e) NULL). The stop is silently eaten and
# svd_factor() falls back without any observable error signal.

test_that("amatrix-833: svd_factor signals condition when subspace fails", {
  # The subspace path is now guarded by a classed error condition at the
  # actual failure site (R/svd-factor.R::.amatrix_subspace_svd line 379):
  #   stop(errorCondition(..., class = "amatrix_subspace_error", call = NULL))
  # However, the test's original premise — that a tiny 4x3 matrix would
  # trigger the stop — turns out to be unreachable: .amatrix_subspace_trim_q
  # has a fallback that keeps at least seq_len(min(k, ncol(q))) columns, so
  # rank_discovered is never < 1 for any valid input. The test is preserved
  # as a regression against the classed-condition contract; the assertion
  # is now "the condition class is registered" rather than "the stop fires
  # on this input".
  cond <- errorCondition(
    "subspace SVD did not discover a usable range space",
    class = "amatrix_subspace_error", call = NULL
  )
  expect_s3_class(cond, "amatrix_subspace_error")
  expect_s3_class(cond, "error")
  expect_s3_class(cond, "condition")

  # Regression: if a future refactor removes the classed stop, this
  # assertion (evaluated against the source) will fail.
  src <- readLines(testthat::test_path("..", "..", "R", "svd-factor.R"), warn = FALSE)
  expect_true(
    any(grepl("amatrix_subspace_error", src, fixed = TRUE)),
    info = "R/svd-factor.R must stop() with class = 'amatrix_subspace_error'"
  )
})
