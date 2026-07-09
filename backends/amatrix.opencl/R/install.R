# ---------------------------------------------------------------------------
# CLBlast runtime installation.
#
# amatrix.opencl never links CLBlast: it is loaded at run time from the system
# search path, from AMATRIX_CLBLAST_LIB, or from the package's user data
# directory. amatrix_install_clblast() populates that data directory with the
# official CLBlast release binary so a user can enable GPU BLAS/linear algebra
# with a single call, most importantly on Windows where CLBlast is rarely
# preinstalled.
# ---------------------------------------------------------------------------

.amatrix_opencl_data_dir <- function() {
  tools::R_user_dir("amatrix.opencl", "data")
}

# Tell the native loader (src/amatrix_cl_loader.c) about the user data
# directory so a CLBlast library placed there is found at probe time. This is
# pure string bookkeeping and touches no device; safe to call at load time.
.amatrix_opencl_register_clblast_dir <- function(dir = .amatrix_opencl_data_dir()) {
  try(
    .Call("amatrix_opencl_set_clblast_path_bridge", as.character(dir), PACKAGE = "amatrix.opencl"),
    silent = TRUE
  )
  invisible(dir)
}

# Human-readable runtime availability reasons from the loader.
amatrix_opencl_availability_reason <- function() {
  .Call("amatrix_opencl_availability_reason_bridge", PACKAGE = "amatrix.opencl")
}

.amatrix_opencl_installed_clblast <- function(dir = .amatrix_opencl_data_dir()) {
  if (!dir.exists(dir)) {
    return(NULL)
  }
  pattern <- if (.Platform$OS.type == "windows") {
    "clblast\\.dll$"
  } else {
    "libclblast.*\\.(so|dylib).*$"
  }
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

# Resolve the official CLBlast release asset for this platform. Returns a list
# with `supported = FALSE` and a `hint` for platform/arch combinations that
# have no official release binary (notably macOS arm64 and Linux aarch64).
.amatrix_opencl_clblast_asset <- function(version) {
  sysname <- Sys.info()[["sysname"]]
  arch <- R.version$arch
  base <- sprintf("https://github.com/CNugteren/CLBlast/releases/download/%s/", version)

  if (identical(sysname, "Windows")) {
    return(list(
      supported = TRUE,
      url = paste0(base, sprintf("CLBlast-%s-windows-x64.7z", version)),
      archive = "7z",
      lib_pattern = "clblast\\.dll$"
    ))
  }
  if (identical(sysname, "Linux") && grepl("x86_64|amd64", arch)) {
    return(list(
      supported = TRUE,
      url = paste0(base, sprintf("CLBlast-%s-linux-x86_64.tar.gz", version)),
      archive = "tar.gz",
      lib_pattern = "libclblast\\.so.*$"
    ))
  }
  if (identical(sysname, "Darwin") && grepl("x86_64", arch)) {
    return(list(
      supported = TRUE,
      url = paste0(base, sprintf("CLBlast-%s-macos-x86_64.tar.gz", version)),
      archive = "tar.gz",
      lib_pattern = "libclblast.*\\.dylib$"
    ))
  }

  hint <- if (identical(sysname, "Darwin")) {
    paste(
      "no official CLBlast release binary is published for macOS arm64.",
      "Install it with Homebrew ('brew install clblast') - amatrix.opencl",
      "finds the Homebrew library automatically - or set AMATRIX_CLBLAST_LIB",
      "to a libclblast.dylib you built for this architecture."
    )
  } else if (identical(sysname, "Linux")) {
    paste(
      "no official CLBlast release binary is published for this architecture.",
      "Install it with your package manager (e.g. `apt-get install libclblast-dev`)",
      "or set AMATRIX_CLBLAST_LIB to a libclblast.so you built."
    )
  } else {
    paste(
      "no official CLBlast release binary is available for this platform.",
      "Set AMATRIX_CLBLAST_LIB to a CLBlast shared library."
    )
  }
  list(supported = FALSE, hint = hint)
}

.amatrix_opencl_extract_clblast <- function(archive_path, asset, dest_dir, quiet) {
  exdir <- tempfile("amatrix-clblast-")
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(exdir, recursive = TRUE, force = TRUE), add = TRUE)

  if (identical(asset$archive, "tar.gz")) {
    ok <- tryCatch({
      utils::untar(archive_path, exdir = exdir)
      TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) {
      return(NULL)
    }
  } else if (identical(asset$archive, "7z")) {
    seven <- Sys.which(c("7z", "7za", "7zr"))
    seven <- seven[nzchar(seven)]
    if (length(seven) == 0L) {
      return(NULL)
    }
    ok <- tryCatch(
      identical(as.integer(system2(
        seven[[1L]],
        c("x", "-y", paste0("-o", exdir), archive_path),
        stdout = if (isTRUE(quiet)) FALSE else "",
        stderr = if (isTRUE(quiet)) FALSE else ""
      )), 0L),
      error = function(e) FALSE
    )
    if (!isTRUE(ok)) {
      return(NULL)
    }
  } else {
    return(NULL)
  }

  hits <- list.files(exdir, pattern = asset$lib_pattern, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!is.na(file.info(hits)$size)]
  if (length(hits) == 0L) {
    return(NULL)
  }
  # Prefer the real (largest) file over version symlinks.
  hits <- hits[order(file.info(hits)$size, decreasing = TRUE)]
  primary <- hits[[1L]]
  dest <- file.path(dest_dir, basename(primary))
  if (!file.copy(primary, dest, overwrite = TRUE)) {
    return(NULL)
  }
  # Copy versioned siblings (e.g. libclblast.so.1) so the loader's bare-name
  # candidates resolve too.
  siblings <- list.files(dirname(primary), pattern = asset$lib_pattern,
                         full.names = TRUE, ignore.case = TRUE)
  for (sib in setdiff(siblings, primary)) {
    file.copy(sib, file.path(dest_dir, basename(sib)), overwrite = TRUE)
  }
  dest
}

#' Download and install the CLBlast runtime library
#'
#' \code{amatrix.opencl} loads CLBlast (the OpenCL BLAS library that powers
#' GPU matrix multiply, \code{crossprod}, triangular solves, and the dense
#' factorizations) at run time; it is never linked at build time. On most
#' systems CLBlast is not preinstalled. This helper downloads the official
#' CLBlast release binary for the current platform from GitHub and places its
#' shared library in the package's user data directory
#' (\code{tools::R_user_dir("amatrix.opencl", "data")}), where the runtime
#' loader searches for it. It is most useful on Windows, where CLBlast is
#' rarely available otherwise.
#'
#' The download requires consent: in an interactive session you are prompted
#' before anything is fetched; in a non-interactive session you must pass
#' \code{force = TRUE} to proceed. Nothing about this function runs at package
#' load or during \code{R CMD check}; the OpenCL device itself is only touched
#' later, by the gated probe (see \code{amatrix_opencl_enable_probe}).
#'
#' No official CLBlast release binary is published for macOS arm64 (Apple
#' Silicon) or Linux aarch64. On those platforms the function stops with a
#' pointer to the system package manager (e.g. \code{brew install clblast});
#' the Homebrew library is discovered automatically. You can always bypass this
#' helper entirely by setting the \code{AMATRIX_CLBLAST_LIB} environment
#' variable to the full path of a CLBlast shared library.
#'
#' @param force Logical. Reinstall even if CLBlast is already present in the
#'   data directory, and consent to downloading without an interactive prompt.
#'   Required in non-interactive sessions. Default \code{FALSE}.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#' @param version Character. CLBlast release version (git tag) to install.
#'   Default \code{"1.7.0"}.
#'
#' @return Invisibly, the path to the installed CLBlast shared library, or
#'   \code{NULL} if the user declined the download.
#'
#' @seealso \code{amatrix_opencl_enable_probe},
#'   \code{amatrix_opencl_diagnostics}
#' @export
amatrix_install_clblast <- function(force = FALSE, quiet = FALSE, version = "1.7.0") {
  stopifnot(
    is.logical(force), length(force) == 1L,
    is.logical(quiet), length(quiet) == 1L,
    is.character(version), length(version) == 1L, nzchar(version)
  )
  say <- function(...) if (!isTRUE(quiet)) message(...)
  dest_dir <- .amatrix_opencl_data_dir()

  existing <- .amatrix_opencl_installed_clblast(dest_dir)
  if (!is.null(existing) && !isTRUE(force)) {
    say("CLBlast is already installed at ", existing,
        ". Use force = TRUE to reinstall.")
    .amatrix_opencl_register_clblast_dir(dest_dir)
    return(invisible(existing))
  }

  asset <- .amatrix_opencl_clblast_asset(version)
  if (!isTRUE(asset$supported)) {
    stop("amatrix_install_clblast(): ", asset$hint, call. = FALSE)
  }

  if (!isTRUE(force)) {
    if (!interactive()) {
      stop(
        "amatrix_install_clblast(): downloading CLBlast requires consent. ",
        "Call it interactively, or pass force = TRUE to consent in a ",
        "non-interactive session.",
        call. = FALSE
      )
    }
    answer <- tolower(trimws(readline(sprintf(
      "Download the official CLBlast %s binary from\n  %s\ninto %s ? [y/N] ",
      version, asset$url, dest_dir
    ))))
    if (!answer %in% c("y", "yes")) {
      say("Aborted; nothing was downloaded.")
      return(invisible(NULL))
    }
  }

  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  archive_path <- tempfile(fileext = paste0(".", asset$archive))
  on.exit(unlink(archive_path, force = TRUE), add = TRUE)

  say("Downloading ", asset$url, " ...")
  dl <- tryCatch(
    utils::download.file(asset$url, archive_path, mode = "wb", quiet = isTRUE(quiet)),
    error = function(e) e
  )
  if (inherits(dl, "error") || !file.exists(archive_path) ||
      is.na(file.info(archive_path)$size) || file.info(archive_path)$size == 0) {
    stop(
      "amatrix_install_clblast(): download failed from ", asset$url,
      if (inherits(dl, "error")) paste0(" (", conditionMessage(dl), ")") else "",
      call. = FALSE
    )
  }

  lib_path <- .amatrix_opencl_extract_clblast(archive_path, asset, dest_dir, quiet = quiet)
  if (is.null(lib_path)) {
    detail <- if (identical(asset$archive, "7z")) {
      paste0(
        "The Windows CLBlast release is a .7z archive; install 7-Zip so that ",
        "'7z' is on the PATH and retry, or extract clblast.dll manually into "
      )
    } else {
      "Extract the CLBlast shared library manually into "
    }
    stop(
      "amatrix_install_clblast(): could not extract the CLBlast library from ",
      "the downloaded archive. ", detail, dest_dir, ".",
      call. = FALSE
    )
  }

  .amatrix_opencl_register_clblast_dir(dest_dir)
  .amatrix_opencl_probe_cache_clear()
  say(
    "CLBlast installed: ", lib_path, "\n",
    "Enable the GPU backend with amatrix::amatrix_use_gpu() ",
    "(or set AMATRIX_OPENCL_PROBE_GPU=1)."
  )
  invisible(lib_path)
}
