# Regression repros for amatrix-8kj secondary-error cleanup.
# Seed: deterministic literal matrix/vector (no RNG required)
# Shape: 3 x 2 dense matrix
# Backend: weighted_cleanup_backend
# Precision mode: strict
# Dispatch path: resident weighted crossprod / tcrossprod
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-8kj

test_that("amatrix-8kj: crossprod_weighted drops resident temp keys after wrap failure", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "weighted_cleanup_backend",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod", "broadcast_ewise"),
      cold_supported_ops = c("crossprod", "broadcast_ewise"),
      resident_supported_ops = c("crossprod", "broadcast_ewise")
    ),
    {
      backend <- get(
        "weighted_cleanup_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident
      backend$broadcast_ewise_resident <- function(lhs_key, stats, margin, fun, out_key, defer = FALSE) {
        if (is.null(counter$broadcast_ewise_resident)) {
          counter$broadcast_ewise_resident <- 0L
        }
        counter$broadcast_ewise_resident <- counter$broadcast_ewise_resident + 1L
        value <- sweep(resident[[lhs_key]], margin, stats, FUN = fun)
        assign(out_key, value, envir = resident)
        value
      }
      assign("weighted_cleanup_backend", backend, envir = amatrix:::.amatrix_state$backends)

      x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L),
        preferred_backend = "weighted_cleanup_backend"
      )
      x <- amatrix_bind_resident(x, "weighted_cleanup_backend")
      x_key <- amatrix:::.amatrix_resident_key(x, backend = "weighted_cleanup_backend")

      local_mocked_bindings(
        .amatrix_rewrap_value = function(...) {
          stop("weighted wrap boom", call. = FALSE)
        },
        .package = "amatrix"
      )

      expect_error(crossprod_weighted(x, c(1, 4, 9)), "weighted wrap boom")
      expect_identical(sort(ls(envir = resident, all.names = FALSE)), x_key)
      expect_true(counter$crossprod_resident >= 1L)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})

test_that("amatrix-8kj: tcrossprod_weighted drops resident temp keys after wrap failure", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "weighted_cleanup_backend_t",
    make_recording_backend(
      counter,
      supported_ops = c("tcrossprod", "broadcast_ewise"),
      cold_supported_ops = c("tcrossprod", "broadcast_ewise"),
      resident_supported_ops = c("tcrossprod", "broadcast_ewise")
    ),
    {
      backend <- get(
        "weighted_cleanup_backend_t",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident
      backend$broadcast_ewise_resident <- function(lhs_key, stats, margin, fun, out_key, defer = FALSE) {
        if (is.null(counter$broadcast_ewise_resident)) {
          counter$broadcast_ewise_resident <- 0L
        }
        counter$broadcast_ewise_resident <- counter$broadcast_ewise_resident + 1L
        value <- sweep(resident[[lhs_key]], margin, stats, FUN = fun)
        assign(out_key, value, envir = resident)
        value
      }
      assign("weighted_cleanup_backend_t", backend, envir = amatrix:::.amatrix_state$backends)

      x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L),
        preferred_backend = "weighted_cleanup_backend_t"
      )
      x <- amatrix_bind_resident(x, "weighted_cleanup_backend_t")
      x_key <- amatrix:::.amatrix_resident_key(x, backend = "weighted_cleanup_backend_t")

      local_mocked_bindings(
        .amatrix_rewrap_value = function(...) {
          stop("weighted wrap boom", call. = FALSE)
        },
        .package = "amatrix"
      )

      expect_error(tcrossprod_weighted(x, c(1, 4, 9)), "weighted wrap boom")
      expect_identical(sort(ls(envir = resident, all.names = FALSE)), x_key)
      expect_true(counter$tcrossprod_resident >= 1L)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})
