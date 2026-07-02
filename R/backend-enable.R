# User-facing GPU enablement and diagnostics.
#
# amatrix_use_gpu() is the documented one-liner that turns on GPU
# acceleration: it walks the fast-backend preference order, enables and
# health-probes each installed backend (explicit user consent, so
# opt-in backends like opencl/arrayfire get their probe env set), and
# adopts the first healthy one as the session default.
#
# amatrix_gpu_status() is the one-call answer to "why am I not on the
# GPU?": a per-backend table of installed/registered/available/health
# with the registry's recorded reasons.

.amatrix_backend_enable_wrappers <- c(
  mlx = "amatrix_mlx_enable_gpu_probe",
  metal = "amatrix_metal_enable_probe",
  arrayfire = "amatrix_arrayfire_enable_probe",
  opencl = "amatrix_opencl_enable_probe"
)

.amatrix_gpu_install_hint <- function() {
  repos_hint <- 'install.packages("amatrix.<backend>", repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))'
  backend_hint <- if (identical(Sys.info()[["sysname"]], "Darwin") &&
                      identical(Sys.info()[["machine"]], "arm64")) {
    "amatrix.mlx (Apple Silicon GPU via MLX)"
  } else {
    "amatrix.opencl (portable OpenCL) or amatrix.arrayfire (ArrayFire)"
  }
  sprintf("no GPU backend packages installed; install %s, e.g. %s",
          backend_hint, repos_hint)
}

#' Enable GPU acceleration for this session
#'
#' Finds, enables, and health-checks an installed GPU backend, then
#' adopts it as the session default for \code{"fast"}-precision work.
#' On Apple Silicon with \pkg{amatrix.mlx} installed this is usually
#' unnecessary: MLX probing is on by default and activates on first
#' use. Call this for the opt-in backends (\pkg{amatrix.opencl},
#' \pkg{amatrix.arrayfire}, \pkg{amatrix.metal}), to force a specific
#' backend, or to get an explicit confirmation line.
#'
#' GPU backends compute in float32 (\code{"fast"} precision, conformance
#' tolerance ~1e-4); \code{"strict"} float64 work always stays on the
#' CPU reference backend regardless of this setting.
#'
#' Side effect: on success this sets the session default precision to
#' \code{"fast"} (and, when \code{backend} is given explicitly, the
#' session default policy to that backend) so subsequent matrices route
#' to the GPU without per-object arguments. Undo with
#' \code{amatrix_set_default_precision("strict")} /
#' \code{amatrix_set_default_policy("auto")}.
#'
#' @param backend Optional backend name (\code{"mlx"}, \code{"metal"},
#'   \code{"arrayfire"}, \code{"opencl"}). Default \code{NULL} tries
#'   the automatic preference order and adopts the first healthy one.
#' @param quiet Logical; suppress the status messages. Default
#'   \code{FALSE}.
#'
#' @return Invisibly, the name of the enabled backend, or \code{FALSE}
#'   if no GPU backend could be enabled.
#'
#' @examples
#' status <- amatrix_gpu_status()
#' if (interactive()) amatrix_use_gpu()
#'
#' @seealso \code{\link{amatrix_gpu_status}},
#'   \code{\link{amatrix_backend_status}},
#'   \code{\link{amatrix_set_default_precision}}
#' @export
amatrix_use_gpu <- function(backend = NULL, quiet = FALSE) {
  specs <- .amatrix_optional_backend_specs()
  order <- intersect(.amatrix_auto_fast_backend_order(), names(specs))
  if (!is.null(backend)) {
    backend <- match.arg(backend, order)
    order <- backend
  }

  say <- function(...) if (!isTRUE(quiet)) message(...)
  failures <- character()

  for (name in order) {
    spec <- specs[[name]]

    if (!nzchar(system.file(package = spec$package))) {
      failures[name] <- sprintf("package %s not installed", spec$package)
      next
    }

    ns <- tryCatch(loadNamespace(spec$package), error = function(e) NULL)
    if (is.null(ns)) {
      failures[name] <- sprintf("package %s failed to load", spec$package)
      next
    }

    # Explicit user consent: activate the backend's probe (sets its
    # AMATRIX_*_PROBE_GPU env), routed through containment for backends
    # that declare it.
    if (isTRUE(spec$contained_probe) &&
        !isTRUE(.amatrix_contained_gpu_probe(name, spec))) {
      health <- .amatrix_backend_health_get(name)
      failures[name] <- health$reason %||% "isolated GPU probe failed"
      next
    }
    wrapper <- .amatrix_backend_enable_wrappers[[name]]
    enable_fun <- get0(wrapper, envir = ns, inherits = FALSE)
    if (is.function(enable_fun)) {
      tryCatch(enable_fun(), error = function(e) NULL)
    }

    if (!isTRUE(.amatrix_try_register_optional_backend(name))) {
      health <- .amatrix_backend_health_get(name)
      failures[name] <- health$reason %||% "registration failed"
      next
    }

    health <- amatrix_backend_health_probe(name)
    if (!identical(health$status, "healthy")) {
      failures[name] <- health$reason %||% "health probe failed"
      next
    }

    amatrix_set_default_precision("fast")
    if (!is.null(backend)) {
      amatrix_set_default_policy(backend)
    }

    say(sprintf(
      "amatrix: GPU enabled - %s backend (float32 'fast' precision, ~1e-4 vs float64; 'strict' float64 stays on CPU). amatrix_gpu_status() for details.",
      name
    ))
    return(invisible(name))
  }

  if (length(failures) == 0L) {
    say(sprintf("amatrix: %s", .amatrix_gpu_install_hint()))
  } else {
    say("amatrix: no GPU backend could be enabled:")
    for (name in names(failures)) {
      say(sprintf("  - %s: %s", name, failures[name]))
    }
    installed_any <- any(vapply(
      specs, function(s) nzchar(system.file(package = s$package)), logical(1)
    ))
    if (!installed_any) {
      say(sprintf("  %s", .amatrix_gpu_install_hint()))
    }
  }
  invisible(FALSE)
}

#' GPU backend status: why am I (not) on the GPU?
#'
#' One row per known GPU backend with the state of every gate between
#' "installed" and "computing on the GPU": package installed, backend
#' registered, device available, health, and the registry's recorded
#' reason when something is off.
#'
#' @return A data frame with columns \code{backend}, \code{package},
#'   \code{installed}, \code{registered}, \code{available},
#'   \code{health}, and \code{reason}.
#'
#' @examples
#' amatrix_gpu_status()
#'
#' @seealso \code{\link{amatrix_use_gpu}},
#'   \code{\link{amatrix_backend_status}}, \code{\link{amatrix_explain}}
#' @export
amatrix_gpu_status <- function() {
  specs <- .amatrix_optional_backend_specs()
  order <- intersect(.amatrix_auto_fast_backend_order(), names(specs))

  rows <- lapply(order, function(name) {
    spec <- specs[[name]]
    installed <- nzchar(system.file(package = spec$package))
    registered <- exists(name, envir = .amatrix_state$backends, inherits = FALSE)
    available <- FALSE
    if (registered) {
      backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
      available <- isTRUE(tryCatch(backend$available(), error = function(e) FALSE))
    }
    health <- .amatrix_backend_health_get(name)
    reason <- health$reason %||% NA_character_
    if (is.na(reason) && installed && !registered) {
      policy <- .amatrix_optional_backend_probe_policy(name, spec)
      reason <- if (!isTRUE(policy$allowed)) {
        policy$reason
      } else if (isTRUE(spec$auto_probe)) {
        "not yet registered (activates on first use)"
      } else {
        "opt-in backend; enable with amatrix_use_gpu()"
      }
    }
    if (is.na(reason) && registered && !identical(health$status, "healthy")) {
      reason <- if (available) {
        "registered; health canary not yet run (runs on first dispatch or amatrix_use_gpu())"
      } else {
        "registered but device unavailable; see amatrix_backend_status()"
      }
    }
    if (!installed) {
      reason <- sprintf("package %s not installed", spec$package)
    }
    data.frame(
      backend = name,
      package = spec$package,
      installed = installed,
      registered = registered,
      available = available,
      health = health$status %||% "unprobed",
      reason = reason,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
