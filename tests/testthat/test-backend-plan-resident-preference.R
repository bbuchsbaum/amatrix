test_that("backend plan prefers resident path when both cold and resident support are available", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "resident_preferred",
    make_recording_backend(
      counter,
      supported_ops = c("matmul", "ewise"),
      cold_supported_ops = c("matmul", "ewise"),
      resident_supported_ops = c("matmul", "ewise")
    ),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "resident_preferred")

      cold_plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
      resident_x <- amatrix_bind_resident(x * 2, backend = "resident_preferred")
      resident_plan <- amatrix_backend_plan(resident_x, "matmul", y = diag(2))
      summary <- amatrix_backend_matrix(resident_x, ops = "matmul", y_map = list(matmul = diag(2)))

      expect_identical(cold_plan$chosen, "resident_preferred")
      expect_identical(cold_plan$chosen_path, "cold")

      expect_identical(resident_plan$chosen, "resident_preferred")
      expect_identical(resident_plan$chosen_path, "resident")
      expect_true(resident_plan$candidates[[1]]$supported_cold)
      expect_true(resident_plan$candidates[[1]]$supported_resident)
      expect_true(resident_plan$candidates[[1]]$resident_active)

      expect_identical(summary$chosen, "resident_preferred")
      expect_identical(summary$chosen_path, "resident")
      expect_true(summary$resident_reuse)
    }
  )
})
