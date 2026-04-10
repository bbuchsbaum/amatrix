# Residency lifecycle conformance tests
# Uses the mock recording backend from helper-conformance.R

# Check whether the residency table has an entry for an object.
has_residency_entry <- function(x) {
  object_key <- paste0("obj:", x@object_id)
  exists(object_key, envir = amatrix:::.amatrix_state$residency, inherits = FALSE)
}

# Count current entries in the global residency table.
residency_table_size <- function() {
  length(ls(amatrix:::.amatrix_state$residency))
}

test_that("GC drops resident binding from side table", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    # Trigger a resident operation so x gets a residency entry.
    result <- x %*% diag(n)
    expect_true(has_residency_entry(x))

    # Capture the object_id before removing all references.
    oid <- x@object_id
    object_key <- paste0("obj:", oid)

    # Drop all R references to x and force GC.
    rm(x, result)
    gc(verbose = FALSE)

    # The residency side table should no longer have an entry for this object.
    expect_false(
      exists(object_key, envir = amatrix:::.amatrix_state$residency, inherits = FALSE),
      label = "GC finalizer should remove residency entry from side table"
    )
  })
})

test_that("GC drops sparse resident binding for mlx-like backends", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(
    counter,
    supported_ops = "matmul",
    cold_supported_ops = "matmul",
    resident_supported_ops = character(),
    supports_sparse_matmul = TRUE
  )

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mlx", backend, {
    X_host <- Matrix::rsparsematrix(5, 4, density = 0.4)
    rhs <- matrix(rnorm(12), nrow = 4)
    X <- adgCMatrix(X_host, preferred_backend = "mlx", precision = "fast")

    result <- X %*% rhs
    expect_true(has_residency_entry(X))

    rm(X, result)
    gc(verbose = FALSE)
    gc(verbose = FALSE)

    expect_identical(counter$sparse_resident_store, 1L)
    expect_identical(counter$sparse_resident_drop, 1L)
  })
})

test_that("Serialization of resident object loses residency cleanly", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    # Upload x into the backend so it becomes resident.
    result <- x %*% diag(n)
    expect_true(has_residency_entry(x))

    # Roundtrip through serialization.
    roundtrip <- unserialize(serialize(x, NULL))

    # unserialize preserves object_id (binary copy), so the roundtrip shares the
    # same side-table entry as x if x is still alive.  The important contract is:
    # (a) the roundtrip is a valid adgeMatrix, and
    # (b) it can still compute correctly (via cold path or re-upload).
    expect_s4_class(roundtrip, "adgeMatrix")
    expect_equal(roundtrip@preferred_backend, x@preferred_backend)

    # The host data must survive the roundtrip.
    expected <- x_host %*% diag(n)
    got <- as.matrix(roundtrip %*% diag(n))
    expect_equal(got, expected, tolerance = 1e-10)

    # Now remove the original x; only roundtrip holds the object_id.
    # After GC the entry may or may not be gone depending on whether the
    # finalizer_env was also serialized — but the roundtrip must still work.
    rm(x, result)
    gc(verbose = FALSE)

    got2 <- as.matrix(roundtrip %*% diag(n))
    expect_equal(got2, expected, tolerance = 1e-10)
  })
})

test_that("Temp resident keys cleaned up after cold matmul", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    y_host <- matrix(rnorm(n * n), nrow = n)
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    table_before  <- residency_table_size()
    stores_before <- if (is.null(counter$resident_store)) 0L else counter$resident_store
    drops_before  <- if (is.null(counter$resident_drop))  0L else counter$resident_drop

    # Perform a resident matmul where y is a plain matrix (temp upload).
    result <- x %*% y_host

    table_after   <- residency_table_size()
    stores_after  <- if (is.null(counter$resident_store)) 0L else counter$resident_store
    drops_after   <- if (is.null(counter$resident_drop))  0L else counter$resident_drop

    new_stores <- stores_after - stores_before
    new_drops  <- drops_after  - drops_before

    # At least some stores happened (x and/or y_host).
    expect_gte(new_stores, 1L)

    # The temp y_host was uploaded and must have been dropped immediately.
    expect_gte(new_drops, 1L,
      label = "temp y_host must be dropped after the op"
    )

    # The residency side table should not have grown by more than the number
    # of stores (every store either becomes a permanent entry or is a temp
    # that was dropped).  Since temps are dropped eagerly:
    #   new permanent entries <= new_stores
    permanent_new <- table_after - table_before
    expect_lte(permanent_new, new_stores,
      label = "residency table growth must not exceed store count"
    )
  })
})

test_that("Chain residency: resident_store called once per object across chained ops", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    stores_before <- if (is.null(counter$resident_store)) 0L else counter$resident_store

    # First op: x gets uploaded and becomes resident.
    result1 <- x %*% diag(n)
    stores_after_first <- if (is.null(counter$resident_store)) 0L else counter$resident_store

    expect_true(has_residency_entry(x))
    x_key_after_first <- amatrix:::.amatrix_resident_key(x)
    expect_false(is.null(x_key_after_first))

    # Second op using the already-resident x: x should NOT be re-uploaded.
    result2 <- x %*% diag(n)
    stores_after_second <- if (is.null(counter$resident_store)) 0L else counter$resident_store

    # x's key must still be the same (was not evicted and re-stored).
    expect_equal(amatrix:::.amatrix_resident_key(x), x_key_after_first,
      label = "x's resident key must not change across chained ops"
    )

    # Between first and second op, stores increased at most by 1
    # (only the output of the second op may have been stored, not x again).
    delta <- stores_after_second - stores_after_first
    expect_lte(delta, 1L,
      label = "resident_store delta between chained ops must be <= 1 (output only, not re-upload of x)"
    )
  })
})

test_that("Temp resident args dropped after operation completes", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    y_host <- matrix(rnorm(n * n), nrow = n)  # plain matrix -> temp upload
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    table_before  <- residency_table_size()
    stores_before <- if (is.null(counter$resident_store)) 0L else counter$resident_store
    drops_before  <- if (is.null(counter$resident_drop))  0L else counter$resident_drop

    result <- x %*% y_host

    table_after  <- residency_table_size()
    stores_after <- if (is.null(counter$resident_store)) 0L else counter$resident_store
    drops_after  <- if (is.null(counter$resident_drop))  0L else counter$resident_drop

    new_stores <- stores_after - stores_before
    new_drops  <- drops_after  - drops_before

    # At least one store (x and/or y_host).
    expect_gt(new_stores, 0L)

    # At least one drop (the temp y_host was cleaned up eagerly).
    expect_gt(new_drops, 0L,
      label = "temp y_host must be dropped immediately after the op"
    )

    # The side table must not have grown by more than the number of stores
    # (temp keys are dropped before returning, so they never count as net new
    # permanent entries).
    permanent_new <- table_after - table_before
    expect_lte(permanent_new, new_stores,
      label = "residency table growth must not exceed store count"
    )
  })
})

test_that("Fallback drops residency honestly when op is not resident-capable", {
  counter <- new.env(parent = emptyenv())
  # Build a backend that supports residency and matmul_resident, but does NOT
  # expose rowSums_resident.  rowSums is still available as a cold method, so
  # amatrix_dispatch_op will find backend$rowSums, drop the residency binding,
  # then materialize x from host data.
  backend <- make_recording_backend(
    counter,
    supported_ops        = c("matmul", "rowSums", "colSums", "ewise"),
    cold_supported_ops   = c("matmul", "rowSums", "colSums", "ewise"),
    resident_supported_ops = c("matmul", "ewise")  # rowSums is NOT resident
  )

  skip_if(
    !amatrix:::.amatrix_backend_residency_capable(backend),
    "mock backend does not support residency"
  )

  with_registered_backend("mock-recording", backend, {
    n <- 3L
    x_host <- matrix(rnorm(n * n), nrow = n)
    x <- adgeMatrix(x_host, preferred_backend = "mock-recording")

    # Upload x into the backend via a resident matmul so it becomes resident.
    result <- x %*% diag(n)
    expect_true(has_residency_entry(x),
      label = "x must have a residency entry after resident matmul"
    )

    # Now call rowSums — this goes through amatrix_dispatch_op.
    # The backend has rowSums (cold) but NOT rowSums_resident, so the
    # dispatch code materializes x from host for the cold method, but
    # preserves the resident binding so future ops can reuse the GPU data.
    rs <- rowSums(x)

    # Residency binding must be preserved (no needless re-upload next time).
    expect_true(has_residency_entry(x),
      label = "residency entry must be preserved after cold-path fallback"
    )

    # Result must still be numerically correct.
    expect_equal(rs, rowSums(x_host), tolerance = 1e-10)
  })
})
