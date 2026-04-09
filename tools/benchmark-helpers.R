local_backend_libpaths <- function() {
  candidates <- c(".tmp/opencl-lib", ".tmp/lib", ".tmp/backends-lib")
  Filter(dir.exists, candidates)
}

load_benchmark_amatrix <- function() {
  prepare_benchmark_libpaths()

  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }

  suppressPackageStartupMessages(library(amatrix))
  invisible(TRUE)
}

prepare_benchmark_libpaths <- function() {
  lib_candidates <- c(local_backend_libpaths(), .libPaths())
  lib_candidates <- unique(normalizePath(lib_candidates, winslash = "/", mustWork = FALSE))

  if (length(lib_candidates) > 0L) {
    .libPaths(c(lib_candidates, .libPaths()))
  }

  invisible(.libPaths())
}

ensure_optional_backend_namespace <- function(package, repo_dir = NULL) {
  prepare_benchmark_libpaths()

  if (!is.null(repo_dir) && dir.exists(repo_dir) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_dir, quiet = TRUE, helpers = FALSE, export_all = FALSE)
    if (package %in% loadedNamespaces()) {
      return(asNamespace(package))
    }
  }

  if (requireNamespace(package, quietly = TRUE)) {
    return(loadNamespace(package))
  }

  NULL
}

.benchmark_arrayfire_requested <- function() {
  identical(Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE", unset = ""), "1") ||
    identical(Sys.getenv("AMATRIX_ARRAYFIRE_PROBE_GPU", unset = ""), "1")
}

.benchmark_optional_backend_specs <- function(include_arrayfire = .benchmark_arrayfire_requested()) {
  specs <- list(
    mlx = list(
      package = "amatrix.mlx",
      repo_dir = "backends/amatrix.mlx",
      name = "mlx",
      precision = "fast",
      register_fun = "amatrix_mlx_register",
      available_fun = "amatrix_mlx_is_available",
      options = c(amatrix.mlx.available = TRUE),
      env = NULL,
      available_args = list()
    ),
    opencl = list(
      package = "amatrix.opencl",
      repo_dir = "backends/amatrix.opencl",
      name = "opencl",
      precision = "fast",
      register_fun = "amatrix_opencl_register",
      available_fun = "amatrix_opencl_native_available",
      options = c(amatrix.enable_opencl = TRUE),
      env = c(AMATRIX_OPENCL_PROBE_GPU = "1"),
      available_args = list(force = TRUE)
    )
  )

  if (isTRUE(include_arrayfire)) {
    specs$arrayfire <- list(
      package = "amatrix.arrayfire",
      repo_dir = "backends/amatrix.arrayfire",
      name = "arrayfire",
      precision = "fast",
      register_fun = "amatrix_arrayfire_register",
      available_fun = "amatrix_arrayfire_is_available",
      options = c(amatrix.enable_arrayfire = TRUE, amatrix.arrayfire.available = TRUE),
      env = c(AMATRIX_ARRAYFIRE_PROBE_GPU = "1"),
      available_args = list()
    )
  }

  specs
}

.benchmark_enable_backend <- function(spec) {
  ns <- ensure_optional_backend_namespace(spec$package, repo_dir = spec$repo_dir)
  if (is.null(ns)) {
    return(FALSE)
  }

  if (!is.null(spec$env)) {
    do.call(Sys.setenv, as.list(spec$env))
  }
  if (!is.null(spec$options)) {
    options(as.list(spec$options))
  }

  try(get(spec$register_fun, envir = ns)(overwrite = TRUE), silent = TRUE)
  available <- try(do.call(get(spec$available_fun, envir = ns), spec$available_args), silent = TRUE)
  isTRUE(available)
}

available_benchmark_backends <- function(
  include_cpu = TRUE,
  include_mlx = TRUE,
  include_opencl = TRUE,
  include_arrayfire = .benchmark_arrayfire_requested()
) {
  backends <- list()

  if (isTRUE(include_cpu)) {
    backends$cpu <- list(name = "cpu", precision = "strict")
  }

  specs <- .benchmark_optional_backend_specs(include_arrayfire = include_arrayfire)
  wanted <- c(
    if (isTRUE(include_mlx)) "mlx",
    if (isTRUE(include_opencl)) "opencl",
    if (isTRUE(include_arrayfire)) "arrayfire"
  )

  for (name in wanted) {
    spec <- specs[[name]]
    if (is.null(spec)) {
      next
    }
    if (.benchmark_enable_backend(spec)) {
      backends[[name]] <- list(name = spec$name, precision = spec$precision)
    }
  }

  backends
}

benchmark_backend_names <- function(...) {
  backends <- available_benchmark_backends(...)
  vapply(backends, `[[`, character(1), "name")
}
