# Regression repros for amatrix-15n.
# Seed: deterministic literal matrices (no RNG required)
# Shape: 2 x 2 dense matrices
# Backend: pref_backend / policy_backend / summary_backend
# Precision mode: strict
# Dispatch path: cold planner selection and backend summary inspection
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-15n

test_that("amatrix-15n: planner prefers preferred_backend before policy fallback", {
  preferred_counter <- new.env(parent = emptyenv())
  policy_counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "mlx",
    make_recording_backend(policy_counter, supported_ops = "matmul"),
    {
      with_registered_backend(
        "arrayfire",
        make_recording_backend(preferred_counter, supported_ops = "matmul"),
        {
          x <- adgeMatrix(
            matrix(1:4, nrow = 2),
            preferred_backend = "arrayfire",
            policy = "mlx"
          )

          plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
          result <- x %*% diag(2)

          expect_identical(plan$preferred, c("arrayfire", "mlx", "cpu"))
          expect_identical(plan$chosen, "arrayfire")
          expect_s4_class(result, "adgeMatrix")
          expect_equal(as.matrix(result), matrix(1:4, nrow = 2))
          expect_equal(preferred_counter$matmul, 1L)
          expect_equal(if (is.null(policy_counter$matmul)) 0L else policy_counter$matmul, 0L)
        }
      )
    }
  )
})

test_that("amatrix-15n: backend matrix reports cpu fallback after non-cpu first preference", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "summary_backend",
    make_recording_backend(counter, supported_ops = c("matmul", "ewise")),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "summary_backend")
      summary <- amatrix_backend_matrix(
        x,
        ops = c("matmul", "solve"),
        y_map = list(matmul = diag(2))
      )

      expect_identical(summary$chosen, c("summary_backend", "cpu"))
      expect_identical(summary$cpu_fallback, c(FALSE, TRUE))
    }
  )
})
