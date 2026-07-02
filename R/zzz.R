.amatrix_state <- new.env(parent = emptyenv())
.amatrix_state$backends <- new.env(parent = emptyenv())
.amatrix_state$default_policy <- "auto"
.amatrix_state$default_precision <- "strict"
.amatrix_state$residency <- new.env(parent = emptyenv())
.amatrix_state$model_cache <- new.env(parent = emptyenv())
.amatrix_state$resident_counter <- 0L
.amatrix_state$object_counter <- 0L
.amatrix_state$session_id <- ""

.amatrix_backend_registration_valid <- function(backend) {
  required_fields <- c(
    "capabilities",
    "features",
    "precision_modes",
    "available",
    "supports",
    "matmul",
    "crossprod",
    "tcrossprod",
    "ewise",
    "rowSums",
    "colSums"
  )

  if (!is.list(backend)) return(FALSE)
  if (length(setdiff(required_fields, names(backend))) > 0L) return(FALSE)
  if (!all(vapply(backend[required_fields], is.function, logical(1)))) return(FALSE)

  capabilities <- tryCatch(backend$capabilities(), error = function(e) NULL)
  features <- tryCatch(backend$features(), error = function(e) NULL)
  precision_modes <- tryCatch(backend$precision_modes(), error = function(e) NULL)

  is.character(capabilities) &&
    is.character(features) &&
    is.character(precision_modes) &&
    length(precision_modes) > 0L &&
    all(precision_modes %in% .amatrix_valid_precisions)
}

.amatrix_register_cpu_backend_on_load <- function() {
  existing <- get0("cpu", envir = .amatrix_state$backends, inherits = FALSE)
  if (is.null(existing)) {
    amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite = FALSE)
  } else if (!isTRUE(.amatrix_backend_registration_valid(existing))) {
    amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite = TRUE)
  }
  invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
  .amatrix_state$session_id <- paste0(
    format(Sys.time(), "%Y%m%d%H%M%OS6"), "-",
    as.hexmode(sample.int(2^31 - 1L, 1L))
  )
  ns <- asNamespace(pkgname)
  registerS3method("as.matrix", "KronMatrix", get("as.matrix.KronMatrix", envir = ns), envir = ns)
  registerS3method("as.matrix", "resident_handle", get("as.matrix.resident_handle", envir = ns), envir = ns)
  .amatrix_cache_init()
  .amatrix_register_cpu_backend_on_load()
}

# TRUE when the user has asked for a silent attach, via either
# options(amatrix.quiet_startup = TRUE) or AMATRIX_QUIET ("1"/"true").
.amatrix_quiet_startup <- function() {
  if (isTRUE(getOption("amatrix.quiet_startup", FALSE))) {
    return(TRUE)
  }
  tolower(Sys.getenv("AMATRIX_QUIET", unset = "")) %in% c("1", "true")
}

.onAttach <- function(libname, pkgname) {
  # One-line GPU visibility note. Cheap checks only: installed-package
  # lookups (no namespace loads) and pure policy evaluation — no
  # probing, registration, or subprocess work at attach time. Silent
  # for pure-CPU users (no backend packages installed), and fully
  # suppressible via options(amatrix.quiet_startup = TRUE) or the
  # AMATRIX_QUIET environment variable.
  if (.amatrix_quiet_startup()) {
    return(invisible())
  }
  # optional_backends disabled globally: no optional backend will be
  # used this session, so there is nothing worth announcing.
  if (!.amatrix_optional_backends_enabled()) {
    return(invisible())
  }
  notes <- tryCatch({
    specs <- .amatrix_optional_backend_specs()
    parts <- character()
    for (name in intersect(.amatrix_auto_fast_backend_order(), names(specs))) {
      spec <- specs[[name]]
      if (!nzchar(system.file(package = spec$package))) next
      # Per-backend opt-out (e.g. options(amatrix.disable_mlx = TRUE)):
      # omit a disabled backend from the note entirely. If every
      # installed backend is disabled, `parts` stays empty and nothing
      # is printed.
      if (!.amatrix_optional_backend_enabled(spec)) next
      policy <- .amatrix_optional_backend_probe_policy(name, spec)
      parts[name] <- if (isTRUE(policy$allowed) && isTRUE(spec$auto_probe)) {
        sprintf("%s ready (activates on first use)", name)
      } else if (isTRUE(policy$allowed)) {
        sprintf("%s installed (enable with amatrix_use_gpu())", name)
      } else {
        sprintf("%s installed but disabled (%s)", name, policy$reason)
      }
    }
    parts
  }, error = function(e) character())

  if (length(notes) > 0L) {
    packageStartupMessage(
      "amatrix GPU backends: ", paste(notes, collapse = "; "),
      ". See amatrix_gpu_status()."
    )
  }
}
