test_that("backend_for short-circuits live resident backends", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "resident_choice_fastpath",
    make_recording_backend(
      counter,
      supported_ops = "ewise",
      cold_supported_ops = "ewise",
      resident_supported_ops = c("matmul", "ewise")
    ),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "resident_choice_fastpath")
      backend <- get(
        "resident_choice_fastpath",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident_key <- amatrix:::.amatrix_next_resident_key("resident_choice_fastpath")
      backend$resident_store(resident_key, as.matrix(x))
      resident_x <- amatrix:::.amatrix_bind_resident(x, "resident_choice_fastpath", resident_key)

      expect_identical(
        amatrix:::.amatrix_live_resident_backend(resident_x),
        "resident_choice_fastpath"
      )

      fast_backend <- backend
      fast_backend$capabilities <- function() {
        stop("planner fast path regressed: capabilities() should not be queried")
      }
      fast_backend$features <- function() {
        stop("planner fast path regressed: features() should not be queried")
      }
      assign(
        "resident_choice_fastpath",
        fast_backend,
        envir = amatrix:::.amatrix_state$backends
      )

      choice <- amatrix:::.amatrix_backend_for(resident_x, "matmul", y = diag(2))

      expect_identical(choice$name, "resident_choice_fastpath")
      expect_equal(as.matrix(resident_x %*% diag(2)), as.matrix(resident_x) %*% diag(2))
      expect_true(counter$matmul_resident >= 1L)
    }
  )
})
