# amatrix_use_gpu() / amatrix_gpu_status() — the user-facing GPU
# enablement surface. Tests are environment-tolerant: they must pass on
# machines with any subset of backend packages installed (including
# none), so structural assertions dominate and backend-specific ones
# are gated.

test_that("amatrix_gpu_status returns one diagnosable row per GPU backend", {
  st <- amatrix_gpu_status()

  expect_s3_class(st, "data.frame")
  expect_setequal(
    colnames(st),
    c("backend", "package", "installed", "registered", "available",
      "health", "reason")
  )
  expect_setequal(st$backend, c("mlx", "metal", "arrayfire", "opencl"))
  expect_type(st$installed, "logical")
  # every non-fully-active backend has a non-empty human-readable reason
  needs_reason <- !(st$registered & st$available & st$health == "healthy")
  expect_true(all(nzchar(st$reason[needs_reason]) & !is.na(st$reason[needs_reason])))
})

test_that("amatrix_use_gpu validates the backend argument", {
  expect_error(amatrix_use_gpu(backend = "nonsense"))
})

test_that("amatrix_use_gpu is quiet when asked and returns name or FALSE", {
  expect_silent(res <- amatrix_use_gpu(quiet = TRUE))
  expect_true(identical(res, FALSE) || (is.character(res) && nzchar(res)))
  if (is.character(res)) {
    st <- amatrix_gpu_status()
    row <- st[st$backend == res, ]
    expect_true(row$registered)
    expect_identical(row$health, "healthy")
  }
})

test_that(".onAttach startup note is silent without backend packages and one line with them", {
  specs <- amatrix:::.amatrix_optional_backend_specs()
  any_installed <- any(vapply(
    specs, function(s) nzchar(system.file(package = s$package)), logical(1)
  ))

  msgs <- character()
  withCallingHandlers(
    amatrix:::.onAttach("", "amatrix"),
    packageStartupMessage = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  if (any_installed) {
    expect_length(msgs, 1L)
    expect_match(msgs, "amatrix GPU backends", fixed = TRUE)
    expect_match(msgs, "amatrix_gpu_status()", fixed = TRUE)
  } else {
    expect_length(msgs, 0L)
  }
})
