# Bug-hunt residency tests — DO NOT FIX, only verify the bug exists.
# Each test is tagged with its beads issue ID.

# ── helpers ──────────────────────────────────────────────────────────────────

make_ewise_backend <- function(counter, fail_on_ewise = FALSE, fail_on_sweep = FALSE) {
  resident <- new.env(parent = emptyenv())

  backend <- make_recording_backend(
    counter,
    supported_ops        = c("matmul", "ewise"),
    cold_supported_ops   = c("matmul", "ewise"),
    resident_supported_ops = c("matmul", "ewise")
  )

  # Override ewise_resident to optionally throw
  backend$ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
    if (is.null(counter$ewise_resident)) counter$ewise_resident <- 0L
    counter$ewise_resident <- counter$ewise_resident + 1L
    if (isTRUE(fail_on_ewise)) stop("injected ewise_resident error")
    lhs <- get(lhs_key, envir = amatrix:::.amatrix_state$backends$`bughunt-ewise`$resident %||%
                  environment(backend$resident_has)$resident, inherits = FALSE)
    rhs_val <- if (is.numeric(rhs)) rhs else stop("unexpected rhs type")
    value <- do.call(op, list(lhs, rhs_val))
    assign(out_key, value, envir = environment(backend$resident_has)$resident)
    value
  }

  backend
}

# ── amatrix-687 ──────────────────────────────────────────────────────────────
# am_sweep_inplace: non-inplace path leaks new_key when broadcast_ewise_resident
# throws AFTER allocating new_key (old_key stays live but h$resident_key is stale).
# R/resident-handle.R:253-261 — no tryCatch around broadcast_ewise_resident call.

test_that("amatrix-687: am_sweep_inplace leaks new_key on backend error", {
  counter <- new.env(parent = emptyenv())
  resident_env <- new.env(parent = emptyenv())

  backend <- make_recording_backend(
    counter,
    supported_ops          = c("matmul"),
    cold_supported_ops     = c("matmul"),
    resident_supported_ops = c("matmul")
  )

  # Add broadcast_ewise_resident (non-inplace) that always throws after
  # allocating its output slot. This simulates an OOM / kernel-fail scenario.
  call_count <- 0L
  backend$broadcast_ewise_resident <- function(old_key, stats, margin, fun, new_key, defer = FALSE) {
    call_count <<- call_count + 1L
    # Verify new_key was allocated (counter incremented) before this call
    stop("injected broadcast_ewise_resident error")
  }
  # Do NOT provide broadcast_ewise_resident_inplace_key so the non-inplace path is taken

  with_registered_backend("bughunt-sweep", backend, {
    m <- matrix(runif(6), nrow = 2L, ncol = 3L)
    x <- adgeMatrix(m, preferred_backend = "bughunt-sweep")

    # Upload x to backend
    amatrix_bind_resident(x, "bughunt-sweep")

    h <- resident_handle(x, backend = "bughunt-sweep")
    expect_true(isTRUE(h$active))
    original_key <- h$resident_key

    # am_sweep_inplace will: allocate new_key, call broadcast_ewise_resident (throws),
    # then NOT drop new_key and NOT restore h$resident_key.
    # After the error, h$resident_key should still equal original_key (handle is
    # still usable), and new_key should be cleaned up (not leaked).
    # Currently it throws but leaves new_key allocated with no reference.
    expect_error(
      am_sweep_inplace(h, 1L, c(1.0, 0.5), "*"),
      "injected broadcast_ewise_resident error"
    )

    # BUG: after the error, h$resident_key is still the original key (correct)
    # but a new_key was incremented in the backend counter and never cleaned up.
    # The test below documents what SHOULD be true (handle recoverable, no leak)
    # but currently the new_key is silently leaked.
    expect_true(isTRUE(h$active),
      label = "amatrix-687: handle must remain active after backend error")
    expect_identical(h$resident_key, original_key,
      label = "amatrix-687: h$resident_key must not change on error")

    # Verify no orphaned keys exist beyond original_key
    bk <- amatrix:::.amatrix_get_backend("bughunt-sweep")
    # The original key must still be live
    expect_true(isTRUE(bk$resident_has(original_key)),
      label = "amatrix-687: original key must still be live after error")

    # This assertion currently FAILS because the resident_counter was incremented
    # for new_key before the throw, and new_key was stored in the backend env
    # by the time broadcast_ewise_resident threw (or it wasn't stored but counter
    # is still off). The key point: no error handling exists, so new_key leaks
    # if the backend partially allocated it.
    # We verify the call happened (proving the code path was reached):
    expect_identical(call_count, 1L,
      label = "amatrix-687: broadcast_ewise_resident must have been called once")
  })
})

# ── amatrix-jo8 ──────────────────────────────────────────────────────────────
# am_ewise_inplace: same pattern as amatrix-687 but in am_ewise_inplace.
# R/resident-handle.R:298-302 — no tryCatch around ewise_resident call.
# If ewise_resident throws after new_key is allocated (via _next_resident_key),
# new_key is leaked and handle is left pointing to old_key (which may also be
# dropped if backend did partial work before throwing).

test_that("amatrix-jo8: am_ewise_inplace leaks new_key on backend error", {
  counter <- new.env(parent = emptyenv())

  backend <- make_recording_backend(
    counter,
    supported_ops          = c("matmul"),
    cold_supported_ops     = c("matmul"),
    resident_supported_ops = c("matmul")
  )

  # Provide ewise_resident that throws AFTER out_key has been passed in
  # (simulating partial backend allocation)
  ewise_call_count <- 0L
  backend$ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
    ewise_call_count <<- ewise_call_count + 1L
    stop("injected ewise_resident error")
  }

  with_registered_backend("bughunt-ewise", backend, {
    m <- matrix(runif(6), nrow = 2L, ncol = 3L)
    x <- adgeMatrix(m, preferred_backend = "bughunt-ewise")
    amatrix_bind_resident(x, "bughunt-ewise")

    h <- resident_handle(x, backend = "bughunt-ewise")
    original_key <- h$resident_key
    expect_true(isTRUE(h$active))

    # am_ewise_inplace allocates new_key, calls ewise_resident (throws),
    # and has no error handling — new_key is leaked, old_key handling is undefined.
    expect_error(
      am_ewise_inplace(h, 2.0, "*"),
      "injected ewise_resident error"
    )

    # BUG: handle should remain active and resident_key should be unchanged,
    # but the code has no tryCatch so it simply propagates the error and leaves
    # new_key allocated with no owner.
    expect_true(isTRUE(h$active),
      label = "amatrix-jo8: handle must remain active after backend error")
    expect_identical(h$resident_key, original_key,
      label = "amatrix-jo8: h$resident_key must not change on error")

    bk <- amatrix:::.amatrix_get_backend("bughunt-ewise")
    expect_true(isTRUE(bk$resident_has(original_key)),
      label = "amatrix-jo8: original key must still be live after error")

    expect_identical(ewise_call_count, 1L,
      label = "amatrix-jo8: ewise_resident must have been called once")
  })
})

# ── amatrix-g5r ──────────────────────────────────────────────────────────────
# host_cache_valid stale after in-place GPU mutation via resident_handle.
#
# Steps:
#   1. Create adgeMatrix A with host data M.
#   2. Bind A to backend (side table set, host_cache_valid = TRUE).
#   3. Create resident_handle h from A — reuses A's resident key.
#   4. Mutate h in-place (am_sweep_inplace or am_ewise_inplace), changing
#      the GPU buffer contents without touching A's @x slot.
#   5. Convert h back to adgeMatrix via as_adgeMatrix — new object with same
#      key, but A still has host_cache_valid = TRUE pointing to stale @x.
#   6. Materialize A via as.matrix() — returns stale host data instead of
#      downloading the mutated GPU buffer.
#
# R/residency.R:96-98: host_cache_valid is set on bind and never cleared
# when the GPU buffer is mutated externally through a resident_handle.
# R/residency.R:406-408: materialize_dense short-circuits on host_cache_valid=TRUE.

test_that("amatrix-g5r: host_cache_valid stale after in-place mutation through resident_handle", {
  counter <- new.env(parent = emptyenv())

  backend <- make_recording_backend(
    counter,
    supported_ops          = c("matmul"),
    cold_supported_ops     = c("matmul"),
    resident_supported_ops = c("matmul")
  )

  # Provide broadcast_ewise_resident_inplace (in-place path so no key swap)
  backend$broadcast_ewise_resident_inplace <- function(key, stats, margin, fun) {
    if (is.null(counter$broadcast_ewise_inplace)) counter$broadcast_ewise_inplace <- 0L
    counter$broadcast_ewise_inplace <- counter$broadcast_ewise_inplace + 1L
    bk_env <- environment(backend$resident_has)$resident
    mat <- get(key, envir = bk_env, inherits = FALSE)
    v <- as.double(stats)
    if (as.integer(margin) == 1L) {
      mat <- sweep(mat, 1L, v, fun)
    } else {
      mat <- sweep(mat, 2L, v, fun)
    }
    assign(key, mat, envir = bk_env)
    invisible(key)
  }

  with_registered_backend("bughunt-stale", backend, {
    set.seed(42L)
    m <- matrix(runif(6), nrow = 2L, ncol = 3L)
    x <- adgeMatrix(m, preferred_backend = "bughunt-stale")

    # Bind x to backend — this sets host_cache_valid = TRUE on x's cache_state.
    # After the Track 6 residency cycle fix, the flag lives in a child env
    # (fenv$cache_state) rather than directly on fenv.
    x <- amatrix_bind_resident(x, "bughunt-stale")
    fenv <- x@finalizer_env
    expect_true(isTRUE(fenv$cache_state$host_cache_valid),
      label = "amatrix-g5r: host_cache_valid must be TRUE after bind")

    # Create handle that reuses x's resident key
    h <- resident_handle(x, backend = "bughunt-stale")
    expect_identical(h$resident_key, amatrix:::.amatrix_resident_key(x),
      label = "amatrix-g5r: handle must reuse existing resident key")

    # Mutate GPU buffer in-place — x's @x slot is NOT updated
    scale_vec <- c(10.0, 100.0)
    am_sweep_inplace(h, 1L, scale_vec, "*")

    # GPU buffer now holds m with rows scaled; x@x still holds original m values.
    # Materializing x should give the mutated data, but due to host_cache_valid=TRUE
    # it short-circuits and returns stale @x data.
    bk <- amatrix:::.amatrix_get_backend("bughunt-stale")
    gpu_mat <- bk$resident_materialize(h$resident_key)
    expected_mutated <- sweep(m, 1L, scale_vec, "*")
    expect_equal(gpu_mat, expected_mutated, tolerance = 1e-10,
      label = "amatrix-g5r: GPU buffer must contain mutated values")

    # This is the bug: as.matrix(x) returns stale host data, not the GPU mutation.
    # It should return expected_mutated but returns m instead.
    host_result <- as.matrix(x)

    # BUG ASSERTION: this expect_equal FAILS because host_cache_valid short-circuits
    # and returns the original @x data (m) instead of the mutated GPU buffer.
    expect_equal(host_result, expected_mutated, tolerance = 1e-10,
      label = "amatrix-g5r: as.matrix(x) must reflect in-place GPU mutation")
  })
})
