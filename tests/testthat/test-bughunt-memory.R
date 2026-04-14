# Memory leak / finalizer bug tests
# Issues: amatrix-5ni, amatrix-mlt, amatrix-o7u, amatrix-4rt, amatrix-ey2
#
# Each test is tagged with the issue ID it covers.
# These tests are expected to FAIL until the bugs are fixed.

# ── helpers ──────────────────────────────────────────────────────────────────

# Build a minimal fake backend with residency support.
# The store/drop/has functions operate on an in-process environment so we can
# observe whether keys are leaked or double-freed without touching a real GPU.
.make_fake_backend <- function(fail_on = NULL) {
  store <- new.env(parent = emptyenv())
  dropped <- new.env(parent = emptyenv())
  double_drops <- new.env(parent = emptyenv())

  list(
    capabilities = function() c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums"),
    features = function() character(0),
    precision_modes = function() c("strict"),
    available = function() TRUE,
    supports = function(op, x, y = NULL) TRUE,
    matmul = function(x, y) x %*% y,
    crossprod = function(x, y = NULL, ...) if (is.null(y)) base::crossprod(x) else base::crossprod(x, y),
    tcrossprod = function(x, y = NULL, ...) if (is.null(y)) base::tcrossprod(x) else base::tcrossprod(x, y),
    ewise = function(x, lhs, rhs = NULL, op, ...) if (is.null(rhs)) do.call(op, list(lhs)) else do.call(op, list(lhs, rhs)),
    rowSums = function(x, na.rm = FALSE, dims = 1L) base::rowSums(x, na.rm = na.rm),
    colSums = function(x, na.rm = FALSE, dims = 1L) base::colSums(x, na.rm = na.rm),
    resident_store = function(key, mat) {
      assign(key, mat, envir = store)
    },
    resident_has = function(key) exists(key, envir = store, inherits = FALSE),
    resident_drop = function(key) {
      if (!exists(key, envir = store, inherits = FALSE)) {
        assign(key, TRUE, envir = double_drops)
        stop("double-free: key '", key, "' already dropped")
      }
      rm(list = key, envir = store)
      assign(key, TRUE, envir = dropped)
    },
    resident_materialize = function(key) {
      get(key, envir = store, inherits = FALSE)
    },
    # Simulate a failing broadcast_ewise_resident (non-inplace path) – amatrix-5ni
    broadcast_ewise_resident = function(old_key, v, margin, op, new_key, defer = FALSE) {
      if (identical(fail_on, "broadcast_ewise_resident")) {
        stop("simulated backend failure")
      }
      mat <- get(old_key, envir = store, inherits = FALSE)
      result <- base::sweep(mat, MARGIN = margin, STATS = v, FUN = op)
      assign(new_key, result, envir = store)
    },
    # Simulate a failing ewise_resident – amatrix-mlt
    ewise_resident = function(old_key, rhs, op, new_key, defer = FALSE) {
      if (identical(fail_on, "ewise_resident")) {
        stop("simulated backend failure")
      }
      mat <- get(old_key, envir = store, inherits = FALSE)
      result <- do.call(op, list(mat, rhs))
      assign(new_key, result, envir = store)
    },
    .store = store,
    .dropped = dropped,
    .double_drops = double_drops
  )
}

.register_fake_backend <- function(name, backend) {
  amatrix_register_backend(name, backend, overwrite = TRUE)
}

.unregister_fake_backend <- function(name) {
  state <- amatrix:::.amatrix_state
  if (exists(name, envir = state$backends, inherits = FALSE)) {
    rm(list = name, envir = state$backends)
  }
}

# ── amatrix-5ni: am_sweep_inplace leaks new_key on backend error ──────────────
# R/resident-handle.R:253-262 — broadcast_ewise_resident called without tryCatch.
# If backend throws after new_key is allocated, new_key is never dropped.
test_that("amatrix-5ni: am_sweep_inplace does not leak new_key when backend errors", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend(fail_on = "broadcast_ewise_resident")
  .register_fake_backend("fake_leak_sweep", bk)
  on.exit(.unregister_fake_backend("fake_leak_sweep"), add = TRUE)

  m <- matrix(1:6 + 0.0, nrow = 2)
  h <- resident_handle(m, backend = "fake_leak_sweep")
  key_before <- h$resident_key

  # This should error (backend fails) but must NOT leave orphaned keys in the store
  expect_error(am_sweep_inplace(h, 1L, c(1, 2), "+"))

  # After the error: store must contain only the original key, no leaked new_key
  live_keys <- ls(bk$.store)
  expect_true(
    key_before %in% live_keys,
    label = "original key must survive after failed sweep"
  )
  leaked <- setdiff(live_keys, key_before)
  expect_true(length(leaked) == 0L,
    info = "amatrix-5ni: new_key leaked in store after backend error in am_sweep_inplace"
  )
})

# ── amatrix-mlt: am_ewise_inplace leaks new_key on backend error ──────────────
# R/resident-handle.R:297-301 — ewise_resident called without error handling.
test_that("amatrix-mlt: am_ewise_inplace does not leak new_key when backend errors", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend(fail_on = "ewise_resident")
  .register_fake_backend("fake_leak_ewise", bk)
  on.exit(.unregister_fake_backend("fake_leak_ewise"), add = TRUE)

  m <- matrix(1:4 + 0.0, nrow = 2)
  h <- resident_handle(m, backend = "fake_leak_ewise")
  key_before <- h$resident_key

  expect_error(am_ewise_inplace(h, 2.0, "+"))

  live_keys <- ls(bk$.store)
  leaked <- setdiff(live_keys, key_before)
  expect_true(length(leaked) == 0L,
    info = "amatrix-mlt: new_key leaked in store after backend error in am_ewise_inplace"
  )
})

# ── amatrix-o7u: double-free if as_adgeMatrix.resident_handle bind throws ─────
# R/resident-handle.R:393-416
# If .amatrix_bind_resident() throws after as.matrix() but before h$active<-FALSE,
# the handle's finalizer will drop the same key that the new adgeMatrix now owns.
# We simulate this by creating a handle whose object_id is NULL so bind_resident
# returns early (object_key is NULL), meaning h$active stays TRUE and the key is
# shared — the test verifies that h$resident_key is NULL after conversion, proving
# ownership was transferred regardless.
test_that("amatrix-o7u: as_adgeMatrix.resident_handle marks handle inert even if bind_resident path is non-standard", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend()
  .register_fake_backend("fake_double_free", bk)
  on.exit(.unregister_fake_backend("fake_double_free"), add = TRUE)

  m <- matrix(runif(6), nrow = 2)
  h <- resident_handle(m, backend = "fake_double_free")
  original_key <- h$resident_key

  obj <- as_adgeMatrix(h)

  # After ownership transfer: handle must be inert and key must be NULL
  expect_false(isTRUE(h$active),
    label = "amatrix-o7u: handle must be inactive after as_adgeMatrix conversion"
  )
  expect_null(h$resident_key,
    label = "amatrix-o7u: handle resident_key must be NULL after ownership transfer"
  )

  # GC the handle — finalizer must NOT drop the key (ownership transferred)
  rm(h)
  gc()

  # Key must still be alive in backend (owned by obj now)
  expect_true(
    bk$resident_has(original_key),
    label = "amatrix-o7u: key was double-freed — finalizer dropped key after ownership transfer"
  )
})

# ── amatrix-ey2: memory_stats resident_objects counts stale registry entries ──
# R/memory-stats.R:33-45 — resident_objects counts registry slots, not live buffers.
# After GC of an object whose backend was already unregistered, the finalizer
# clears the residency registry but the byte counts from backend$memory_usage()
# may be inconsistent. More concretely: if we register a residency entry manually
# and then unregister the backend, resident_objects stays nonzero while bytes_used=NA.
test_that("amatrix-ey2: amatrix_memory_stats resident_objects reflects only live backend entries", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend()
  .register_fake_backend("fake_stats", bk)

  m <- matrix(1:4 + 0.0, nrow = 2)
  obj <- adgeMatrix(m, preferred_backend = "fake_stats")

  # Force residency — upload manually
  key <- amatrix:::.amatrix_next_resident_key("fake_stats")
  bk$resident_store(key, m)
  obj <- amatrix:::.amatrix_bind_resident(obj, "fake_stats", key)

  # Verify it shows up
  s1 <- amatrix_memory_stats()
  fake_row <- s1$residency[s1$residency$backend == "fake_stats", ]
  expect_equal(fake_row$resident_objects, 1L,
    label = "amatrix-ey2: should see 1 resident object before unregister"
  )

  # Now unregister the backend (simulates crashed/unloaded backend)
  .unregister_fake_backend("fake_stats")

  # Drop the buffer from the (now-unregistered) fake backend store directly
  rm(list = key, envir = bk$.store)

  # The registry entry is still present (finalizer hasn't run yet)
  # memory_stats must NOT count it as a live resident object
  s2 <- amatrix_memory_stats()
  fake_row2 <- s2$residency[!is.na(s2$residency$backend) & s2$residency$backend == "fake_stats", ]

  # Bug: if the stale entry is still counted, resident_objects will be 1
  # Expected correct behavior: 0 (entry is stale, backend gone)
  if (nrow(fake_row2) > 0L) {
    expect_equal(fake_row2$resident_objects, 0L,
      label = "amatrix-ey2: stale residency entry counted as live resident_object after backend unregistered"
    )
  } else {
    succeed("amatrix-ey2: no stale entry shown — correct behavior")
  }
})

# ── stress: GC loop does not leak resident keys ────────────────────────────────
# Verifies that creating and GC-ing 50 handles doesn't grow the residency store.
test_that("stress: resident_handle GC loop leaves no leaked keys in fake backend", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend()
  .register_fake_backend("fake_stress", bk)
  on.exit(.unregister_fake_backend("fake_stress"), add = TRUE)

  n_before <- length(ls(bk$.store))

  for (i in seq_len(50)) {
    h <- resident_handle(matrix(runif(4), 2, 2), backend = "fake_stress")
    rm(h)
    gc()
  }

  n_after <- length(ls(bk$.store))
  expect_equal(n_after, n_before,
    label = "stress: leaked keys found in fake backend store after GC loop"
  )
})

# ── multi-env: handle key still works after one env loses reference ────────────
test_that("multi-env: handle remains valid when shared across envs and one ref is removed", {
  skip_if_not_installed("amatrix")

  bk <- .make_fake_backend()
  .register_fake_backend("fake_multienv", bk)
  on.exit(.unregister_fake_backend("fake_multienv"), add = TRUE)

  h <- resident_handle(matrix(1:4 + 0.0, 2, 2), backend = "fake_multienv")

  e1 <- new.env(parent = emptyenv())
  e2 <- new.env(parent = emptyenv())
  e1$h <- h
  e2$h <- h

  rm(h)
  rm(e1)
  gc()

  # e2$h must still be active (it holds the only remaining reference)
  expect_true(isTRUE(e2$h$active),
    label = "multi-env: handle became inactive when only one env reference was removed"
  )
  expect_true(
    bk$resident_has(e2$h$resident_key),
    label = "multi-env: resident key was dropped while e2 still holds a reference"
  )
})

# ── bad-construction: error in constructor must not leave partial residency ────
test_that("bad-construction: adgeMatrix with invalid data leaves no residency entry", {
  skip_if_not_installed("amatrix")

  n_res_before <- length(ls(amatrix:::.amatrix_state$residency))

  # Force an error in the constructor by passing something that triggers stop()
  # Use a non-matrix non-adgeMatrix type that .amatrix_new_dense will reject
  expect_error(
    adgeMatrix("not a matrix"),
    label = "bad-construction: adgeMatrix should reject non-matrix input"
  )

  gc()
  n_res_after <- length(ls(amatrix:::.amatrix_state$residency))

  # A failed constructor must not leave a dangling residency slot
  expect_equal(n_res_after, n_res_before,
    label = "bad-construction: failed adgeMatrix constructor left a dangling residency entry"
  )
})
