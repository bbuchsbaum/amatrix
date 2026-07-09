# Regression guard for amatrix-36q.
# Seed: deterministic literal matrix/vector (no RNG required)
# Shape: 3 x 2 dense matrix
# Backend: row_weight_backend (recording backend)
# Precision mode: strict
# Dispatch path: .amatrix_apply_row_weights resident scaling, op failure
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-36q
#
# amatrix-36q reported two "net-new double-drop" sites (backend-planning.R:395,
# models-lm.R:620) where an error-handler drop of out_key was followed by an
# unconditional second drop of the same key. Both were converted to the
# on.exit + release-flag shape (commits 2de7b852 / f3e9ef37), so out_key is now
# released exactly once on every path. This guard pins that invariant on the
# models-lm.R site: the resident op stores out_key and then throws; the call
# must fall back to CPU and release out_key exactly once (never twice).

test_that("amatrix-36q: weighted row scaling releases out_key exactly once when the resident op fails", {
  counter <- new.env(parent = emptyenv())
  drop_counts <- new.env(parent = emptyenv())

  with_registered_backend(
    "row_weight_backend",
    make_recording_backend(
      counter,
      supported_ops = c("broadcast_ewise"),
      cold_supported_ops = c("broadcast_ewise"),
      resident_supported_ops = c("broadcast_ewise")
    ),
    {
      backend <- get(
        "row_weight_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident

      # Count every raw resident_drop call (whether or not the key still exists),
      # so a redundant second drop of the same key is visible.
      orig_drop <- backend$resident_drop
      backend$resident_drop <- function(key) {
        prev <- drop_counts[[key]]
        drop_counts[[key]] <- (if (is.null(prev)) 0L else prev) + 1L
        orig_drop(key)
      }
      # Allocate out_key on device, then fail — this is the exact window the
      # 36q double-drop lived in.
      backend$broadcast_ewise_resident <- function(lhs_key, stats, margin, fun, out_key, defer = FALSE) {
        assign(out_key, sweep(resident[[lhs_key]], margin, stats, FUN = fun), envir = resident)
        stop("resident op boom", call. = FALSE)
      }
      assign("row_weight_backend", backend, envir = amatrix:::.amatrix_state$backends)

      x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L),
        preferred_backend = "row_weight_backend"
      )
      x <- amatrix_bind_resident(x, "row_weight_backend")
      x_key <- amatrix:::.amatrix_resident_key(x, "row_weight_backend")

      # The op failure is caught internally and the call falls back to CPU.
      weighted <- amatrix:::.amatrix_apply_row_weights(x, c(1, 4, 9))
      expect_true(inherits(weighted, "adgeMatrix"))
      expect_equal(
        as.matrix(amatrix_materialize_host(weighted)),
        matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L, ncol = 2L) * sqrt(c(1, 4, 9)),
        ignore_attr = TRUE
      )

      # out_key was allocated then dropped exactly once; never twice.
      out_keys <- setdiff(ls(envir = drop_counts, all.names = TRUE), x_key)
      expect_true(length(out_keys) >= 1L)
      for (k in out_keys) {
        expect_lte(drop_counts[[k]], 1L)
      }
      # The input binding must survive the op failure.
      expect_true(backend$resident_has(x_key))
    }
  )
})
