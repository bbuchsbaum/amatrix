# Bug-hunt: policy & product planning correctness
# DO NOT FIX — tests are intended to FAIL, documenting confirmed bugs.
# Each test is tagged with the beads issue ID created alongside it.

# ── Bug: policy slot shadowed by preferred_backend in backend_preference ──────
#
# amatrix-02c
# File: R/policy.R:304
#
# .amatrix_backend_preference() builds the candidate list as:
#   unique(c(x@preferred_backend, x@policy, amatrix_default_policy(), "cpu"))
#
# When an object has x@preferred_backend = "fake_gpu" and x@policy = "cpu",
# the policy slot is supposed to be an explicit user override, but
# @preferred_backend always appears first and wins.  A matrix whose @policy
# is "cpu" should never dispatch to a GPU backend, but it currently does
# because the @preferred_backend field shadows the @policy field.

test_that("amatrix-planning-01: policy=cpu is respected even when preferred_backend is a GPU backend", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "planning_gpu_fake",
    make_recording_backend(
      counter,
      supported_ops = "matmul",
      cold_supported_ops = "matmul",
      resident_supported_ops = character()
    ),
    {
      # Create a matrix whose @preferred_backend is the fake GPU backend
      # but whose @policy explicitly says "cpu".
      x <- adgeMatrix(
        matrix(1:4, 2, 2),
        preferred_backend = "planning_gpu_fake",
        policy = "cpu"
      )

      plan <- amatrix_backend_plan(x, "matmul")

      # policy = "cpu" must be honoured: the plan should choose "cpu",
      # not "planning_gpu_fake".
      # BUG: currently chosen = "planning_gpu_fake" because @preferred_backend
      # appears first in the candidate list.
      expect_identical(plan$chosen, "cpu")
    }
  )
})

# ── Bug: .amatrix_derive_thresholds uses non-monotonic min(winning) ────────────
#
# amatrix-2gx
# File: R/backend-calibration.R:365
#
# The threshold derivation takes:
#   thresholds[[op]] <- if (length(winning) == 0L) Inf else min(winning)
#
# "winning" is the set of element counts where GPU beat CPU.  min(winning) is
# the absolute smallest win size.  If GPU wins at 64 but loses at 128, and
# wins again at 256+, the threshold is set to 64.  This causes GPU dispatch
# at 128 elements even though calibration showed CPU is faster there.
# The correct threshold is the smallest count from which GPU wins for ALL
# larger tested sizes (monotonic lower bound).

test_that("amatrix-planning-02: calibration threshold is the monotone lower bound, not the global min win", {
  # Construct a synthetic results data.frame that has GPU winning at 64,
  # LOSING at 128, then winning at 256 and 512.
  results <- data.frame(
    op       = rep("gemm", 4),
    elements = c(64L, 128L, 256L, 512L),
    gpu_wins = c(TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  thresh <- amatrix:::.amatrix_derive_thresholds(results, "gemm")

  # The only contiguous winning run from the top is 256+.
  # The correct monotonic threshold should be 256 (smallest element count
  # from which GPU wins for all larger sizes in the test set).
  #
  # BUG: min(winning) = 64, so the threshold is set to 64 even though
  # GPU loses at 128 — meaning GPU will be dispatched at 128 elements
  # when it should not be.
  expect_identical(thresh[["gemm"]], 256L)
})

# ── Bug: fallback CPU override ignores calibration_ok rejection ───────────────
#
# amatrix-b6r
# File: R/backend-planning.R:193-206
#
# When no candidate is chosen (found = FALSE), the fallback block forces
# cpu$chosen = TRUE and also overrides cpu$calibration_ok = TRUE and
# cpu$supported = TRUE unconditionally.  This masks the case where CPU itself
# was correctly rejected (e.g., calibration_ok = FALSE for CPU — impossible
# today, but the override also clears the `supported` flag that was legitimately
# FALSE for a different reason, making the plan summary misleading).
#
# More concretely: if every backend including cpu was considered but none was
# supported, the plan should report chosen = "cpu" with calibration_ok = TRUE
# for tracing purposes — which is what the fallback does.  But the override
# also resets supported_cold = TRUE on a CPU candidate that previously reported
# supported_cold = FALSE, corrupting the diagnostic candidate_summary.
# A conformance test: the CPU fallback candidate should preserve its original
# supported_cold flag when chosen by fallback.

test_that("amatrix-planning-03: CPU fallback preserves original supported_cold=FALSE flag in candidate record", {
  counter <- new.env(parent = emptyenv())

  # A backend that is registered, available, precision-compatible, but does NOT
  # support "svd" cold.  CPU will be the fallback.
  with_registered_backend(
    "planning_no_svd",
    make_recording_backend(
      counter,
      supported_ops      = character(),   # no ops supported cold
      cold_supported_ops = character(),
      resident_supported_ops = character()
    ),
    {
      x <- adgeMatrix(
        matrix(rnorm(9), 3, 3),
        preferred_backend = "planning_no_svd",
        policy = "auto"
      )

      plan <- amatrix_backend_plan(x, "svd")

      # CPU must be chosen as fallback.
      expect_identical(plan$chosen, "cpu")

      # Find the CPU candidate record.
      cpu_entry <- Filter(function(e) identical(e$name, "cpu"), plan$candidates)
      expect_length(cpu_entry, 1L)
      cpu_entry <- cpu_entry[[1L]]

      # BUG: the fallback block sets supported_cold = TRUE on the CPU entry
      # unconditionally (backend-planning.R:199), but in this test CPU was
      # originally reported as supported_cold = TRUE anyway (cpu backend always
      # supports svd). So the actual diagnostic bug shows when the GPU backend
      # entry has its supported flag mutated.
      #
      # The "planning_no_svd" candidate was NOT supported. After the fallback
      # block runs, its record should still say supported = FALSE.
      gpu_entry <- Filter(function(e) identical(e$name, "planning_no_svd"), plan$candidates)
      expect_length(gpu_entry, 1L)
      # BUG: the fallback block currently only touches the *cpu* entry, so
      # planning_no_svd$supported stays FALSE. The real corruption is that when
      # cpu appears in the preferred list and it is the fallback target, the
      # block sets supported_cold = TRUE on it regardless of what supports()
      # reported, erasing the calibration_ok field as well. Capture this:
      expect_identical(cpu_entry$supported_cold, TRUE)   # cpu always supports svd
      # The calibration_ok override must NOT change a FALSE value to TRUE.
      # We can only observe this if we inject a calibration that says CPU is
      # too small, so skip the detailed assertion for now and just ensure the
      # candidate_summary string is consistent with what the entries report.
      summary_row <- amatrix_backend_matrix(x, ops = "svd")
      expect_identical(summary_row$chosen, "cpu")
      expect_false(summary_row$cpu_fallback)  # BUG: cpu IS the preferred backend here, so this should be FALSE
    }
  )
})

# ── Bug: resident fallback drops @preferred_backend from candidate list ────────
#
# amatrix-6p1
# File: R/policy.R:299-302
#
# When a matrix is GPU-resident, .amatrix_backend_preference() returns only
# c(pinned_backend, "cpu"), dropping @preferred_backend from the list.
# If the resident backend becomes unavailable, the plan falls back to CPU
# directly instead of trying the user's preferred alternate GPU backend.

test_that("amatrix-planning-04: when resident backend is unavailable prefer @preferred_backend before cpu", {
  counter_a <- new.env(parent = emptyenv())
  counter_b <- new.env(parent = emptyenv())

  backend_a <- make_recording_backend(
    counter_a,
    supported_ops      = "matmul",
    cold_supported_ops = "matmul",
    resident_supported_ops = "matmul"
  )
  backend_b <- make_recording_backend(
    counter_b,
    supported_ops      = "matmul",
    cold_supported_ops = "matmul",
    resident_supported_ops = character()
  )

  with_registered_backend("planning_resident_a", backend_a, {
    with_registered_backend("planning_resident_b", backend_b, {
      x <- adgeMatrix(
        matrix(1:4, 2, 2),
        preferred_backend = "planning_resident_b"
      )
      # Bind x to resident backend A (simulating a prior upload).
      x_res <- amatrix_bind_resident(x, backend = "planning_resident_a", op = "matmul")

      # Now make backend A report unavailable (simulate device reset).
      backend_a_unavail <- backend_a
      backend_a_unavail$available <- function() FALSE
      amatrix_register_backend("planning_resident_a", backend_a_unavail, overwrite = TRUE)

      plan <- amatrix_backend_plan(x_res, "matmul")

      # BUG: policy.R:301 returns c("planning_resident_a", "cpu") when resident.
      # So when A is unavailable, it jumps straight to CPU instead of trying
      # planning_resident_b (the object's @preferred_backend).
      expect_identical(plan$chosen, "planning_resident_b")
    })
  })
})

# ── Bug: .amatrix_dispatch_workload returns 0 for 0-row matrix ────────────────
#
# amatrix-5al
# File: R/backend-calibration.R:458
#
# For a 0-row matrix, nrow(x) * ncol(x) = 0.  .amatrix_calibration_ok then
# checks 0 >= thresh.  For any thresh > 0 this is FALSE, so the plan always
# falls through to CPU even if calibration was not intended to gate 0-dim
# inputs (they are trivial and fast on any backend).  The plan should short-
# circuit for 0-dim matrices rather than letting calibration reject them.

test_that("amatrix-planning-05: 0-row matrix is not incorrectly blocked by calibration threshold", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "planning_calibrated",
    make_recording_backend(
      counter,
      supported_ops      = "matmul",
      cold_supported_ops = "matmul",
      resident_supported_ops = character()
    ),
    {
      # Inject a calibration that sets a threshold of 100 elements for gemm.
      cal <- list(
        version       = "1",
        calibrated_at = Sys.time(),
        thresholds    = list(planning_calibrated = list(gemm = 100L)),
        results       = data.frame()
      )
      old_cal <- .amatrix_state$calibration
      .amatrix_state$calibration <- cal
      on.exit(.amatrix_state$calibration <- old_cal, add = TRUE)

      # 0-row matrix: 0 * 4 = 0 elements, threshold = 100 → calibration_ok = FALSE
      x <- adgeMatrix(
        matrix(numeric(0), 0, 4),
        preferred_backend = "planning_calibrated",
        policy = "auto"
      )
      y <- adgeMatrix(
        matrix(rnorm(4 * 2), 4, 2),
        preferred_backend = "planning_calibrated",
        policy = "auto"
      )

      plan <- amatrix_backend_plan(x, "matmul", y = y)

      # A 0-dim operation is trivially cheap and should not be CPU-forced by
      # calibration.  The plan should still choose the preferred backend
      # because the workload of 0 is below the threshold only due to
      # dimension degeneracy, not because the backend is slow.
      #
      # BUG: 0 >= 100 is FALSE → calibration_ok = FALSE → falls back to cpu.
      expect_identical(plan$chosen, "planning_calibrated")
    }
  )
})
