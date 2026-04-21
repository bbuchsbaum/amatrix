# Regression repro metadata
# Seed: none (deterministic startup-policy probes)
# Dimensions: none; fresh-process runtime-context checks only
# Backend / precision / dispatch: mlx/opencl startup policy in source-tree subprocesses
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-goq.2, amatrix-goq.3

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

test_that("mlx startup policy skips probe in direct file-entry context", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  repo_dir <- .backend_policy_repo_dir()
  skip_if(is.null(repo_dir),
          "source tree not reachable (installed-pkg context)")

  lines <- c(
    sprintf("pkgload::load_all(%s, quiet = TRUE)", shQuote(repo_dir)),
    "status <- amatrix_backend_status('mlx')",
    "cat('available=', status$available[[1]], '\\n', sep='')",
    "cat('health=', status$health[[1]], '\\n', sep='')",
    "cat('reason=', status$health_reason[[1]], '\\n', sep='')"
  )

  output <- .run_backend_policy_probe("file", lines)

  expect_true(any(grepl("^available=FALSE$", output)))
  expect_true(any(grepl("^health=unprobed$", output)))
  expect_true(any(grepl("unsafe file-entry Rscript context", output, fixed = TRUE)))
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
