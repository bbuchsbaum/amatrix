local_backend_libpaths <- function() {
  candidates <- c(".tmp/opencl-lib", ".tmp/lib", ".tmp/backends-lib", ".tmp/metal-lib")
  Filter(dir.exists, candidates)
}

.benchmark_debug <- function(...) {
  path <- Sys.getenv("AMATRIX_BENCHMARK_LAUNCH_DEBUG", unset = "")
  if (!nzchar(path)) {
    return(invisible(FALSE))
  }

  line <- paste(..., collapse = "")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), line), file = path, append = TRUE)
  invisible(TRUE)
}

.benchmark_debug_state <- function(tag) {
  env_names <- c(
    "AMATRIX_OPENCL_PROBE_GPU",
    "AMATRIX_METAL_PROBE_GPU",
    "R_LIBS",
    "R_LIBS_USER",
    "DYLD_LIBRARY_PATH",
    "LD_LIBRARY_PATH",
    "PATH"
  )
  env_values <- vapply(
    env_names,
    function(name) sprintf("%s=%s", name, Sys.getenv(name, unset = "<unset>")),
    character(1)
  )

  .benchmark_debug(
    tag,
    " ; wd=", getwd(),
    " ; libpaths=", paste(.libPaths(), collapse = " | "),
    " ; env=", paste(env_values, collapse = " ; ")
  )
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
  .benchmark_debug_state(sprintf("ensure_optional_backend_namespace start package=%s repo_dir=%s", package, repo_dir %||% "<null>"))

  if (!is.null(repo_dir) && dir.exists(repo_dir) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload_error <- NULL
    source_ns <- tryCatch(
      {
        pkgload::load_all(repo_dir, quiet = TRUE, helpers = FALSE, export_all = FALSE)
        if (package %in% loadedNamespaces()) {
          asNamespace(package)
        } else {
          NULL
        }
      },
      error = function(e) {
        pkgload_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(source_ns)) {
      .benchmark_debug(
        "ensure_optional_backend_namespace loaded via pkgload package=", package,
        " ; ns_path=", getNamespaceInfo(source_ns, "path")
      )
      return(source_ns)
    }
    if (!is.null(pkgload_error)) {
      warning(
        sprintf("pkgload::load_all('%s') failed: %s — falling back to installed package", repo_dir, pkgload_error),
        call. = FALSE, immediate. = TRUE
      )
      .benchmark_debug(
        "ensure_optional_backend_namespace pkgload FAILED package=", package,
        " ; error=", pkgload_error
      )
    }
  }

  if (requireNamespace(package, quietly = TRUE)) {
    ns <- loadNamespace(package)
    ns_path <- getNamespaceInfo(ns, "path")
    .benchmark_debug(
      "ensure_optional_backend_namespace loaded installed package=", package,
      " ; ns_path=", ns_path
    )
    if (!is.null(repo_dir) && dir.exists(repo_dir)) {
      warning(
        sprintf(
          "Backend '%s' loaded from installed path (%s) instead of source (%s) — installed .so may be stale",
          package, ns_path, repo_dir
        ),
        call. = FALSE, immediate. = TRUE
      )
    }
    return(ns)
  }

  .benchmark_debug("ensure_optional_backend_namespace failed package=", package)
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
      options = c(amatrix.enable_opencl = TRUE, amatrix.opencl.factor_gpu = TRUE),
      env = c(AMATRIX_OPENCL_PROBE_GPU = "1"),
      available_args = list(force = TRUE)
    ),
    metal = list(
      package = "amatrix.metal",
      repo_dir = "backends/amatrix.metal",
      name = "metal",
      precision = "fast",
      register_fun = "amatrix_metal_register",
      available_fun = "amatrix_metal_native_available",
      options = c(amatrix.enable_metal = TRUE, amatrix.metal.available = TRUE),
      env = c(AMATRIX_METAL_PROBE_GPU = "1"),
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
  if (!is.null(spec$env)) {
    do.call(Sys.setenv, as.list(spec$env))
  }
  if (!is.null(spec$options)) {
    options(as.list(spec$options))
  }
  .benchmark_debug_state(sprintf("enable_backend start %s", spec$name))

  ns <- ensure_optional_backend_namespace(spec$package, repo_dir = spec$repo_dir)
  if (is.null(ns)) {
    .benchmark_debug("enable_backend ", spec$name, ": namespace unavailable")
    return(FALSE)
  }

  try(get(spec$register_fun, envir = ns)(overwrite = TRUE), silent = TRUE)
  available <- try(do.call(get(spec$available_fun, envir = ns), spec$available_args), silent = TRUE)
  ns_path <- getNamespaceInfo(ns, "path")
  so_path <- file.path(ns_path, "libs", paste0(spec$package, .Platform$dynlib.ext))
  if (!file.exists(so_path)) {
    so_path <- file.path(ns_path, "src", paste0(spec$package, .Platform$dynlib.ext))
  }
  .benchmark_debug(
    "enable_backend ", spec$name,
    " ; ns_path=", ns_path,
    " ; so_path=", so_path,
    " ; so_exists=", file.exists(so_path),
    " ; so_mtime=", if (file.exists(so_path)) format(file.mtime(so_path), "%Y-%m-%d %H:%M:%S") else "<missing>",
    ": env=", paste(names(spec$env %||% character()), unlist(spec$env %||% character()), collapse = ","),
    " ; options=", paste(names(spec$options %||% character()), unlist(spec$options %||% character()), collapse = ","),
    " ; available_result=", if (inherits(available, "try-error")) as.character(available) else as.character(available)
  )

  if (identical(spec$name, "opencl")) {
    diag <- try(get("amatrix_opencl_diagnostics", envir = ns, inherits = FALSE)(), silent = TRUE)
    if (!inherits(diag, "try-error") && is.list(diag)) {
      diag_parts <- vapply(
        names(diag),
        function(name) sprintf("%s=%s", name, as.character(diag[[name]])),
        character(1)
      )
      .benchmark_debug("enable_backend opencl diagnostics ; ", paste(diag_parts, collapse = " ; "))
    }
  }

  isTRUE(available)
}

available_benchmark_backends <- function(
  include_cpu = TRUE,
  include_mlx = TRUE,
  include_metal = TRUE,
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
    if (isTRUE(include_metal)) "metal",
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
      .benchmark_debug("available_backends admitted ", spec$name)
    } else {
      .benchmark_debug("available_backends rejected ", spec$name)
    }
  }

  backends
}

benchmark_backend_names <- function(...) {
  backends <- available_benchmark_backends(...)
  vapply(backends, `[[`, character(1), "name")
}

r_string_literal <- function(x) {
  encodeString(x, quote = "\"")
}

benchmark_rscript_source_args <- function(script_path, args = character(), working_dir = getwd(), main_call = NULL) {
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  working_dir <- normalizePath(working_dir, winslash = "/", mustWork = TRUE)
  expr <- sprintf("setwd(%s); source(%s, local = globalenv())",
    r_string_literal(working_dir),
    r_string_literal(script_path)
  )
  if (!is.null(main_call) && nzchar(main_call)) {
    expr <- sprintf("%s; %s", expr, main_call)
  }
  c("-e", expr, "--args", args)
}

benchmark_system2_capture <- function(command, args) {
  warned_status <- NULL
  quoted_args <- vapply(args, shQuote, character(1), USE.NAMES = FALSE)
  output <- withCallingHandlers(
    system2(command, quoted_args, stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )

  list(
    output = output,
    status = attr(output, "status") %||% warned_status %||% 0L
  )
}

# ---------------------------------------------------------------------------
# prime_backend: force host->device upload before timing begins
# ---------------------------------------------------------------------------

prime_backend <- function(obj, backend) {
  if (identical(backend, "cpu") || is.null(backend)) {
    return(invisible(obj))
  }
  tryCatch(
    {
      primed <- amatrix::amatrix_bind_resident(obj, backend = backend)
      invisible(primed)
    },
    error = function(e) invisible(obj)
  )
}

# ---------------------------------------------------------------------------
# append_benchmark_history: record every run's summary rows to an append-only
# history CSV so baseline.csv updates never lose historical measurements.
# ---------------------------------------------------------------------------

.benchmark_history_columns <- function() {
  c(
    "timestamp", "git_sha", "op", "size_label", "backend", "variant",
    "median_ms", "sd_ms", "n_reps", "host_os", "host_cpu"
  )
}

.benchmark_git_sha <- function() {
  sha <- tryCatch(
    suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (length(sha) == 0L || !nzchar(sha[[1L]])) "unknown" else sha[[1L]]
}

.benchmark_host_cpu <- function() {
  info <- tryCatch(benchmarkme::get_cpu(), error = function(e) NULL)
  if (is.list(info) && !is.null(info$model_name)) {
    return(info$model_name)
  }
  sys <- Sys.info()
  sys[["machine"]] %||% "unknown"
}

append_benchmark_history <- function(df, path) {
  cols <- .benchmark_history_columns()

  if (is.null(df) || nrow(df) == 0L) {
    return(invisible(character()))
  }

  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  git_sha <- .benchmark_git_sha()
  host_os <- paste(Sys.info()[["sysname"]], Sys.info()[["release"]], sep = "/")
  host_cpu <- .benchmark_host_cpu()

  backend_col <- if ("requested_backend" %in% names(df)) {
    df$requested_backend
  } else {
    rep(NA_character_, nrow(df))
  }

  row <- data.frame(
    timestamp  = rep(timestamp, nrow(df)),
    git_sha    = rep(git_sha, nrow(df)),
    op         = df$op %||% rep(NA_character_, nrow(df)),
    size_label = df$size_label %||% rep(NA_character_, nrow(df)),
    backend    = backend_col,
    variant    = df$variant %||% rep(NA_character_, nrow(df)),
    median_ms  = df$median_ms %||% rep(NA_real_, nrow(df)),
    sd_ms      = df$sd_ms %||% rep(NA_real_, nrow(df)),
    n_reps     = df$n_reps %||% rep(NA_integer_, nrow(df)),
    host_os    = rep(host_os, nrow(df)),
    host_cpu   = rep(host_cpu, nrow(df)),
    stringsAsFactors = FALSE
  )
  row <- row[, cols, drop = FALSE]

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  need_header <- !file.exists(path) || file.info(path)$size == 0
  utils::write.table(
    row,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = need_header,
    append = !need_header,
    qmethod = "double"
  )
  invisible(path)
}

# ---------------------------------------------------------------------------
# assert_backend_accuracy: compare GPU result against CPU reference
# ---------------------------------------------------------------------------

TOLERANCES <- list(
  svd             = 1e-4,
  rsvd            = 1e-4,
  svd_factor_subspace = 1e-4,
  chol            = 1e-6,
  matmul          = 1e-5,
  crossprod       = 1e-5,
  tcrossprod      = 1e-5,
  dist            = 1e-5,
  many_lm         = 1e-5
)

assert_backend_accuracy <- function(ref, gpu, op, tol = NULL) {
  if (is.null(tol)) {
    tol <- TOLERANCES[[op]]
    if (is.null(tol)) {
      stop(sprintf("assert_backend_accuracy: no tolerance registered for op '%s'", op), call. = FALSE)
    }
  }
  ref_mat <- as.matrix(ref)
  gpu_mat <- as.matrix(gpu)
  rel_err <- norm(ref_mat - gpu_mat, type = "F") /
    max(norm(ref_mat, type = "F"), .Machine$double.eps)
  if (rel_err > tol) {
    stop(sprintf(
      "accuracy regression: op=%s tol=%.2e rel_err=%.2e",
      op, tol, rel_err
    ), call. = FALSE)
  }
  rel_err
}
