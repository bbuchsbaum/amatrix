# Regression repro for amatrix-aul (resident-handle in-place remainder).
# Seed: deterministic literal matrix (no RNG required)
# Shape: 2 x 3 dense matrix
# Backend: aul_inplace_backend (recording backend)
# Precision mode: strict
# Dispatch path: am_sweep_inplace / am_ewise_inplace key-swap, post-op failure
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-aul
#
# The in-place resident ops allocate a NEW key, run the backend kernel (which
# stores that key), and then do post-op bookkeeping (alias repointing + old-key
# drop) BEFORE committing h$resident_key. The kernel-failure path is already
# handled by the error handler, but a throw in the post-op bookkeeping orphaned
# the freshly stored new_key: h still pointed at old_key, so nothing tracked
# new_key and it leaked in device memory. The fix registers an on.exit that
# drops new_key unless the swap committed.

.aul_postop_backend <- function(counter, op = c("sweep", "ewise")) {
  op <- match.arg(op)
  backend <- make_recording_backend(
    counter,
    supported_ops = c("matmul"),
    cold_supported_ops = c("matmul"),
    resident_supported_ops = c("matmul")
  )
  resident <- environment(backend$resident_has)$resident
  if (identical(op, "sweep")) {
    backend$broadcast_ewise_resident <- function(old_key, stats, margin, fun, new_key, defer = FALSE) {
      value <- sweep(get(old_key, envir = resident), as.integer(margin), as.double(stats), fun)
      assign(new_key, value, envir = resident)
      invisible(new_key)
    }
  } else {
    backend$ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
      value <- do.call(op, list(get(lhs_key, envir = resident), rhs))
      assign(out_key, value, envir = resident)
      value
    }
  }
  backend
}

test_that("amatrix-aul: am_sweep_inplace drops new_key when post-op bookkeeping throws", {
  counter <- new.env(parent = emptyenv())
  backend <- .aul_postop_backend(counter, "sweep")

  with_registered_backend("aul_sweep_postop", backend, {
    resident <- environment(backend$resident_has)$resident
    x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L, ncol = 3L),
      preferred_backend = "aul_sweep_postop"
    )
    x <- amatrix_bind_resident(x, "aul_sweep_postop")
    keys_before <- sort(ls(envir = resident, all.names = TRUE))

    local_mocked_bindings(
      .amatrix_update_resident_aliases = function(...) stop("postop boom", call. = FALSE),
      .package = "amatrix"
    )

    h <- resident_handle(x, backend = "aul_sweep_postop")
    expect_error(am_sweep_inplace(h, 1L, c(1, 1), "*"), "postop boom")

    # new_key must have been released; the resident env is unchanged.
    expect_setequal(ls(envir = resident, all.names = TRUE), keys_before)
  })
})

test_that("amatrix-aul: am_ewise_inplace drops new_key when post-op bookkeeping throws", {
  counter <- new.env(parent = emptyenv())
  backend <- .aul_postop_backend(counter, "ewise")

  with_registered_backend("aul_ewise_postop", backend, {
    resident <- environment(backend$resident_has)$resident
    x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L, ncol = 3L),
      preferred_backend = "aul_ewise_postop"
    )
    x <- amatrix_bind_resident(x, "aul_ewise_postop")
    keys_before <- sort(ls(envir = resident, all.names = TRUE))

    local_mocked_bindings(
      .amatrix_update_resident_aliases = function(...) stop("postop boom", call. = FALSE),
      .package = "amatrix"
    )

    h <- resident_handle(x, backend = "aul_ewise_postop")
    expect_error(am_ewise_inplace(h, 2.0, "*"), "postop boom")

    expect_setequal(ls(envir = resident, all.names = TRUE), keys_before)
  })
})
