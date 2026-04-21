# Track 5 — backend health, fallback telemetry, and balanced-mode deprecation.
#
# These tests verify:
#   - .amatrix_backend_health_{mark,get} state machine
#   - amatrix_backend_health_probe() runs a canary and marks health
#   - amatrix_fallback_log() / amatrix_fallback_log_reset() API
#   - amatrix_backend_status() surfaces health columns
#   - Balanced mode emits a one-shot deprecation warning and maps to "exact"
#   - Dispatch auto-fallback: a failing fake backend logs an event and the
#     CPU fallback returns the correct result

test_that("backend health defaults to 'unprobed' and can be marked", {
  amatrix:::.amatrix_backend_health_init()
  # Reset state for known-fixed backend ids.
  amatrix:::.amatrix_backend_health_mark("cpu", "unprobed", NA_character_)

  rec <- amatrix:::.amatrix_backend_health_get("cpu")
  expect_identical(rec$status, "unprobed")

  amatrix:::.amatrix_backend_health_mark("cpu", "healthy", NULL)
  rec <- amatrix:::.amatrix_backend_health_get("cpu")
  expect_identical(rec$status, "healthy")
  expect_true(inherits(rec$timestamp, "POSIXct"))

  amatrix:::.amatrix_backend_health_mark("cpu", "unhealthy", "test")
  rec <- amatrix:::.amatrix_backend_health_get("cpu")
  expect_identical(rec$status, "unhealthy")
  expect_identical(rec$reason, "test")
})

test_that("amatrix_backend_health_probe runs canary on cpu and marks healthy", {
  # Reset before probe.
  amatrix:::.amatrix_backend_health_mark("cpu", "unprobed", NA_character_)

  result <- amatrix_backend_health_probe("cpu")
  expect_identical(result$status, "healthy")
  expect_null(result$reason)

  # Status reflects in amatrix_backend_status().
  status <- amatrix_backend_status("cpu")
  expect_identical(status$health[status$name == "cpu"], "healthy")
})

test_that("health probe marks unregistered backends unhealthy", {
  result <- amatrix_backend_health_probe("nonexistent-fake-backend-xyz")
  expect_identical(result$status, "unhealthy")
  expect_true(grepl("not registered", result$reason))
})

test_that("unhealthy backend is excluded from automatic planning", {
  fake <- list(
    capabilities    = function() "matmul",
    features        = function() "dense_f64",
    precision_modes = function() "strict",
    available       = function() TRUE,
    supports        = function(op, x, y = NULL) identical(op, "matmul"),
    matmul          = function(x, y) x %*% y,
    crossprod       = function(x, y = NULL, ...) base::crossprod(x, y),
    tcrossprod      = function(x, y = NULL, ...) base::tcrossprod(x, y),
    ewise           = function(x, lhs, rhs = NULL, op, ...) if (is.null(rhs)) do.call(op, list(lhs)) else do.call(op, list(lhs, rhs)),
    rowSums         = function(x, ...) base::rowSums(x, ...),
    colSums         = function(x, ...) base::colSums(x, ...)
  )
  amatrix_register_backend("fake_unhealthy", fake, overwrite = TRUE)
  on.exit({
    if (exists("fake_unhealthy", envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
      rm(list = "fake_unhealthy", envir = amatrix:::.amatrix_state$backends)
    }
  }, add = TRUE)

  amatrix:::.amatrix_backend_health_mark("fake_unhealthy", "unhealthy", "forced test failure")

  x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "fake_unhealthy", policy = "auto", precision = "strict")
  plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
  fake_entry <- Filter(function(entry) identical(entry$name, "fake_unhealthy"), plan$candidates)[[1]]

  expect_identical(plan$chosen, "cpu")
  expect_identical(fake_entry$health, "unhealthy")
  expect_identical(fake_entry$health_reason, "forced test failure")
  expect_false(fake_entry$health_eligible)
  expect_false(fake_entry$supported)
})

test_that("amatrix_fallback_log() starts empty and accumulates events", {
  amatrix_fallback_log_reset()
  empty <- amatrix_fallback_log()
  expect_s3_class(empty, "data.frame")
  expect_identical(nrow(empty), 0L)
  expect_identical(
    names(empty),
    c("timestamp", "op", "from_backend", "to_backend", "reason")
  )

  amatrix:::.amatrix_log_fallback(
    op = "matmul",
    backend = "mlx",
    reason = "test event",
    from_backend = "mlx",
    to_backend = "cpu"
  )
  amatrix:::.amatrix_log_fallback(
    op = "svd",
    backend = "arrayfire",
    reason = "another event",
    from_backend = "arrayfire",
    to_backend = "cpu"
  )

  log <- amatrix_fallback_log()
  expect_identical(nrow(log), 2L)
  expect_identical(log$op, c("matmul", "svd"))
  expect_identical(log$from_backend, c("mlx", "arrayfire"))
  expect_identical(log$to_backend, c("cpu", "cpu"))

  amatrix_fallback_log_reset()
  expect_identical(nrow(amatrix_fallback_log()), 0L)
})

test_that("amatrix_backend_status() includes health and health_reason columns", {
  status <- amatrix_backend_status("cpu")
  expect_true(all(c("name", "available", "health", "health_reason") %in% names(status)))
  expect_identical(status$name, "cpu")
})

test_that("balanced mode emits a one-shot deprecation warning and maps to exact", {
  reset_balanced_flag <- function() {
    st <- amatrix:::.amatrix_state
    st$balanced_deprecation_warned <- FALSE
  }
  reset_balanced_flag()

  # The first call should warn AND return the same triple as mode="exact".
  expect_warning(
    amatrix:::.amatrix_resolve_mode("balanced", backend = NULL, preferred_backend = NULL,
                                    policy = NULL, precision = NULL),
    regexp = "balanced.*deprecated"
  )

  # Reset, then capture the value (suppressing the warning this time) and
  # assert it matches the exact-mode triple.
  reset_balanced_flag()
  triple_balanced <- suppressWarnings(
    amatrix:::.amatrix_resolve_mode("balanced", backend = NULL, preferred_backend = NULL,
                                    policy = NULL, precision = NULL)
  )
  triple_exact <- amatrix:::.amatrix_resolve_mode(
    "exact", backend = NULL, preferred_backend = NULL,
    policy = NULL, precision = NULL
  )
  expect_identical(triple_balanced, triple_exact)
  expect_identical(triple_balanced$preferred_backend, "cpu")
  expect_identical(triple_balanced$precision, "strict")

  # Second "balanced" call in the same session (flag already set) is silent.
  # expect_silent would spuriously fail on benign stray messages; use
  # warnings-collection instead.
  warnings_seen <- character(0)
  withCallingHandlers(
    amatrix:::.amatrix_resolve_mode("balanced", backend = NULL, preferred_backend = NULL,
                                    policy = NULL, precision = NULL),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(
    grep("balanced.*deprecated", warnings_seen, value = TRUE),
    0L
  )
})

test_that("auto-fallback dispatch logs when a registered backend raises", {
  # Register a deliberately-broken fake backend that errors on matmul. Verify
  # that amatrix_dispatch_op catches the error, logs a fallback, and the
  # caller-provided fallback produces the correct result.
  fake <- list(
    capabilities    = function() c("matmul", "crossprod"),
    features        = function() "dense_f64",
    precision_modes = function() "strict",
    available       = function() TRUE,
    supports        = function(op, x, y = NULL) TRUE,
    matmul          = function(x, y) stop("fake backend always errors", call. = FALSE),
    crossprod       = function(x, y = NULL, ...) stop("fake backend always errors", call. = FALSE),
    tcrossprod      = function(x, y = NULL, ...) stop("fake backend always errors", call. = FALSE),
    ewise           = function(x, y, op) x,
    rowSums         = function(x, ...) numeric(nrow(x)),
    colSums         = function(x, ...) numeric(ncol(x))
  )
  amatrix_register_backend("fake_broken", fake, overwrite = TRUE)
  on.exit({
    if (exists("fake_broken", envir = amatrix:::.amatrix_state$backends)) {
      rm("fake_broken", envir = amatrix:::.amatrix_state$backends)
    }
  }, add = TRUE)

  amatrix_fallback_log_reset()
  amatrix:::.amatrix_backend_health_mark("fake_broken", "unprobed", NA_character_)

  x_host <- matrix(rnorm(12), 3, 4)
  y_host <- matrix(rnorm(12), 4, 3)
  # policy="auto" so the dispatcher honours preferred_backend; with an
  # explicit non-auto policy it would win over preferred_backend per
  # R/policy.R::.amatrix_backend_preference.
  x <- adgeMatrix(x_host, preferred_backend = "fake_broken",
                  policy = "auto", precision = "strict")
  y <- adgeMatrix(y_host, preferred_backend = "fake_broken",
                  policy = "auto", precision = "strict")

  result <- amatrix_dispatch_op(
    x, "matmul", method = "matmul", y = y,
    args = list(y = amatrix_materialize_host(y)),
    fallback = function() x_host %*% y_host
  )
  expect_equal(result, x_host %*% y_host, tolerance = 1e-10)

  log <- amatrix_fallback_log()
  expect_true(nrow(log) >= 1L)
  expect_true(any(grepl("fake_broken", log$from_backend)))
  expect_true(any(grepl("runtime error", log$reason)))

  # Backend should have been marked unhealthy.
  h <- amatrix:::.amatrix_backend_health_get("fake_broken")
  expect_identical(h$status, "unhealthy")

  amatrix_fallback_log_reset()
})
