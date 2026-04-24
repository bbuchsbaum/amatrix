# Regression repro for .onLoad CPU backend overwrite state invalidation.
# Seed: none; deterministic literal calibration state
# Shape: no numeric payload; backend registry startup path only
# Backend: cpu
# Precision mode: strict
# Dispatch path: .onLoad -> cpu backend registration
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-lei

test_that(".onLoad preserves existing CPU calibration and health state", {
  state <- amatrix:::.amatrix_state
  old_calibration <- state$calibration
  old_backend_health <- state$backend_health
  old_session_id <- state$session_id
  old_cpu <- get("cpu", envir = state$backends, inherits = FALSE)

  on.exit({
    state$calibration <- old_calibration
    state$backend_health <- old_backend_health
    state$session_id <- old_session_id
    assign("cpu", old_cpu, envir = state$backends)
  }, add = TRUE)

  state$backend_health <- new.env(parent = emptyenv())
  state$calibration <- list(
    version = "2",
    sys_hash = amatrix:::.amatrix_sys_hash(),
    calibrated_at = Sys.time(),
    thresholds = list(cpu = list(gemm = 123L)),
    results = data.frame(
      backend = "cpu",
      op = "gemm",
      elements = 123L,
      gpu_wins = FALSE,
      stringsAsFactors = FALSE
    )
  )
  amatrix:::.amatrix_backend_health_mark("cpu", "healthy", NULL)

  amatrix:::.onLoad("", "amatrix")

  expect_identical(state$calibration$thresholds$cpu$gemm, 123L)
  expect_true(any(state$calibration$results$backend == "cpu"))
  expect_identical(amatrix:::.amatrix_backend_health_get("cpu")$status, "healthy")
})
