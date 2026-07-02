# Regression repro metadata
# Seed: none (deterministic startup-policy probes)
# Dimensions: none; fresh-process runtime-context checks only
# Backend / precision / dispatch: mlx/opencl startup policy in source-tree subprocesses
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-goq.2, amatrix-goq.3
#
# CONTRACT CHANGE 2026-07-01: MLX probing is now default-ON in every launch
# context (including direct `Rscript file.R`), with AMATRIX_MLX_PROBE_GPU=0 /
# options(amatrix.auto_probe = FALSE) as opt-outs and a contained
# child-process first probe as the crash safety net. The old file-entry
# probe skip was retired after the upstream NSRangeException stopped
# reproducing on mlx >= 0.31 (26/26 certified runs; see
# planning_docs/mlx-spectral-benchmark-instability.md and
# tools/certify-mlx-file-entry.R).

.run_backend_policy_probe <- function(mode = c("expr", "file"), body_lines, env = character()) {
  mode <- match.arg(mode)
  script_path <- tempfile(fileext = ".R")
  on.exit(unlink(script_path), add = TRUE)

  writeLines(body_lines, script_path)

  args <- if (identical(mode, "expr")) {
    expr <- sprintf("source(%s, local = TRUE)", dQuote(script_path))
    c("--vanilla", "-e", shQuote(expr))
  } else {
    c("--vanilla", script_path)
  }

  system2(
    R.home("bin/Rscript"),
    args = args,
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
}

.backend_policy_repo_dir <- function() {
  candidates <- unique(c(
    tryCatch(getNamespaceInfo(asNamespace("amatrix"), "path"), error = function(e) NULL),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  ))
  candidates <- Filter(Negate(is.null), candidates)
  matches <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(matches) == 0L) {
    return(NULL)
  }
  normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("mlx startup policy allows probing in direct file-entry context", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    "Sys.unsetenv('AMATRIX_MLX_PROBE_GPU')",
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "policy <- amatrix:::.amatrix_optional_backend_probe_policy('mlx')",
    "cat('allowed=', policy$allowed, '\\n', sep='')"
  )

  output <- .run_backend_policy_probe("file", lines)

  expect_true(any(grepl("^allowed=TRUE$", output)))
})

test_that("mlx probing honors the AMATRIX_MLX_PROBE_GPU=0 opt-out", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "status <- amatrix_backend_status('mlx')",
    "policy <- amatrix:::.amatrix_optional_backend_probe_policy('mlx')",
    "cat('allowed=', policy$allowed, '\\n', sep='')",
    "cat('available=', status$available[[1]], '\\n', sep='')",
    "cat('reason=', status$health_reason[[1]], '\\n', sep='')"
  )

  output <- .run_backend_policy_probe(
    "file", lines,
    env = c("AMATRIX_MLX_PROBE_GPU=0")
  )

  expect_true(any(grepl("^allowed=FALSE$", output)))
  expect_true(any(grepl("^available=FALSE$", output)))
  expect_true(any(grepl("AMATRIX_MLX_PROBE_GPU=0", output, fixed = TRUE)))
})

test_that("mlx probing honors options(amatrix.auto_probe = FALSE)", {
  old <- options(amatrix.auto_probe = FALSE)
  on.exit(options(old), add = TRUE)
  env_backup <- Sys.getenv("AMATRIX_MLX_PROBE_GPU", unset = NA)
  Sys.unsetenv("AMATRIX_MLX_PROBE_GPU")
  on.exit({
    if (!is.na(env_backup)) Sys.setenv(AMATRIX_MLX_PROBE_GPU = env_backup)
  }, add = TRUE)

  policy <- amatrix:::.amatrix_optional_backend_probe_policy("mlx")

  expect_false(policy$allowed)
  expect_match(policy$reason, "amatrix.auto_probe", fixed = TRUE)
})

test_that("contained gpu probe survives a crashing child and blocks the probe env", {
  env_backup <- Sys.getenv("AMATRIX_TEST_FAKE_PROBE_GPU", unset = NA)
  Sys.unsetenv("AMATRIX_TEST_FAKE_PROBE_GPU")
  on.exit({
    Sys.unsetenv("AMATRIX_TEST_FAKE_PROBE_GPU")
    if (!is.na(env_backup)) Sys.setenv(AMATRIX_TEST_FAKE_PROBE_GPU = env_backup)
    rm(list = intersect("fakecrash", ls(amatrix:::.amatrix_contained_probe_cache)),
       envir = amatrix:::.amatrix_contained_probe_cache)
  }, add = TRUE)

  spec <- list(
    package = "amatrix.fake",
    probe_env = "AMATRIX_TEST_FAKE_PROBE_GPU",
    available_bridge = "none",
    probe_expr = "q(status = 139)"
  )
  verdict <- amatrix:::.amatrix_contained_gpu_probe("fakecrash", spec)

  expect_false(verdict)
  expect_identical(Sys.getenv("AMATRIX_TEST_FAKE_PROBE_GPU"), "0")
})

test_that("contained gpu probe passes a healthy child and exports the probe env", {
  env_backup <- Sys.getenv("AMATRIX_TEST_FAKE_PROBE_GPU", unset = NA)
  Sys.unsetenv("AMATRIX_TEST_FAKE_PROBE_GPU")
  on.exit({
    Sys.unsetenv("AMATRIX_TEST_FAKE_PROBE_GPU")
    if (!is.na(env_backup)) Sys.setenv(AMATRIX_TEST_FAKE_PROBE_GPU = env_backup)
    rm(list = intersect("fakeok", ls(amatrix:::.amatrix_contained_probe_cache)),
       envir = amatrix:::.amatrix_contained_probe_cache)
  }, add = TRUE)

  spec <- list(
    package = "amatrix.fake",
    probe_env = "AMATRIX_TEST_FAKE_PROBE_GPU",
    available_bridge = "none",
    probe_expr = "cat('AMATRIX-PROBE-OK')"
  )
  verdict <- amatrix:::.amatrix_contained_gpu_probe("fakeok", spec)

  expect_true(verdict)
  expect_identical(Sys.getenv("AMATRIX_TEST_FAKE_PROBE_GPU"), "1")
})

test_that("mlx startup policy allows probe in -e context", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "policy <- amatrix:::.amatrix_optional_backend_probe_policy('mlx')",
    "cat('allowed=', policy$allowed, '\\n', sep='')"
  )

  output <- .run_backend_policy_probe("expr", lines)

  expect_true(any(grepl("^allowed=TRUE$", output)))
})

test_that("mlx optional backend registration activates the safe gpu probe hook", {
  fake_ns <- new.env(parent = emptyenv())
  fake_ns$called <- FALSE
  fake_ns$amatrix_mlx_enable_gpu_probe <- function() {
    fake_ns$called <- TRUE
    TRUE
  }

  activated <- amatrix:::.amatrix_activate_optional_backend_probe(
    "mlx",
    fake_ns,
    list(
      probe_env = "AMATRIX_MLX_PROBE_GPU",
      enable_probe_fun = "amatrix_mlx_enable_gpu_probe"
    )
  )

  expect_true(activated)
  expect_true(fake_ns$called)
})

test_that("opencl startup policy requires explicit probe opt-in", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    "Sys.unsetenv('AMATRIX_OPENCL_PROBE_GPU')",
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "status <- amatrix_backend_status('opencl')",
    "policy <- amatrix:::.amatrix_optional_backend_probe_policy('opencl')",
    "cat('allowed=', policy$allowed, '\\n', sep='')",
    "cat('available=', status$available[[1]], '\\n', sep='')",
    "cat('reason=', status$health_reason[[1]], '\\n', sep='')"
  )

  output <- .run_backend_policy_probe("expr", lines)

  expect_true(any(grepl("^allowed=FALSE$", output)))
  expect_true(any(grepl("^available=FALSE$", output)))
  expect_true(any(grepl("AMATRIX_OPENCL_PROBE_GPU=1", output, fixed = TRUE)))
})

test_that("opencl startup policy allows explicit probe runs", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "policy <- amatrix:::.amatrix_optional_backend_probe_policy('opencl')",
    "cat('allowed=', policy$allowed, '\\n', sep='')"
  )

  output <- .run_backend_policy_probe(
    "expr",
    lines,
    env = c("AMATRIX_OPENCL_PROBE_GPU=1")
  )

  expect_true(any(grepl("^allowed=TRUE$", output)))
})
