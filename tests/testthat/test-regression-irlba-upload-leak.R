# Regression repro for amatrix-4rt.
# Seed: deterministic literal matrix (no RNG required for the leak assertion)
# Shape: 4 x 4 dense matrix
# Backend: lanczos_leak_backend (recording backend)
# Precision mode: strict
# Dispatch path: block_lanczos source/right operator compilation
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-4rt
#
# block_lanczos() uploads the source operator (A_left) to a resident backend
# BEFORE it registers the on.exit that drops it. If building the right operator
# (A_right) throws, the source operator's device buffer leaks because the
# cleanup on.exit was never reached. The fix registers the source-drop on.exit
# immediately after the upload, before the right-operator build.

test_that("amatrix-4rt: source operator key is dropped when right operator build fails", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "lanczos_leak_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      backend <- get(
        "lanczos_leak_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident

      local_mocked_bindings(
        .amatrix_block_lanczos_source_operator = function(A) {
          key <- amatrix:::.amatrix_next_resident_key("lanczos_leak_backend")
          backend$resident_store(key, matrix(0, nrow = 4L, ncol = 4L))
          list(backend = "lanczos_leak_backend", resident_key = key, temporary = TRUE)
        },
        .amatrix_block_lanczos_right_operator = function(A, source_operator = NULL) {
          stop("right operator boom", call. = FALSE)
        },
        .package = "amatrix"
      )

      A <- adgeMatrix(matrix(c(
        4, 1, 0, 0,
        1, 4, 1, 0,
        0, 1, 4, 1,
        0, 0, 1, 4
      ), nrow = 4L), preferred_backend = "lanczos_leak_backend")

      expect_error(block_lanczos(A, nv = 1L), "right operator boom")

      # The source operator upload must not leak: after the failed build the
      # recording backend's resident environment must hold no buffers.
      expect_length(ls(envir = resident, all.names = TRUE), 0L)
      expect_true(isTRUE(counter$resident_drop >= 1L))
    }
  )
})
