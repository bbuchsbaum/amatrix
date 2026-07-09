# Regression: OpenCL/CLBlast must be resolved at RUNTIME, never linked.
#
# Bug: `R CMD check` on core amatrix aborted the R process during "checking
# dependencies in R code" with an uncaught C++ exception:
#   libc++abi: terminating due to uncaught exception of type cl::Error: clGetDeviceIDs
# The amatrix.opencl shared object linked OpenCL and CLBlast at load time, so a
# C++ cl::Error thrown inside the OpenCL/CLBlast stack (as happens on a headless
# builder where device enumeration misbehaves) propagated through our C frame to
# std::terminate. A process abort during check is fatal on CRAN and R-universe.
#
# Fix: amatrix.opencl links against neither library. Every OpenCL and CLBlast
# symbol is resolved at run time via dlopen()/LoadLibrary() (src/amatrix_cl_loader.c),
# so loading the namespace maps nothing OpenCL, and the gated probe degrades to
# "unavailable" instead of aborting.
#
# Reproduction metadata:
#   - Trigger: loadNamespace("amatrix.opencl") + gated availability probe.
#   - Dispatch path: cold namespace load; probe both disabled and force-enabled.
#   - Platform of original report: aarch64-apple-darwin (R 4.5.1); the abort
#     only manifests where device enumeration errors (headless builder), so the
#     subprocess guard passes on a GPU-equipped host and protects CI.
#   - Issues: amatrix-tw4, bd-01KX33B03AHS382F70FNAE07F2 (mote).

test_that("loading amatrix.opencl and probing never aborts the R process", {
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript) || file.exists(paste0(rscript, ".exe")),
              "Rscript not available")

  child <- '
    if (!requireNamespace("amatrix.opencl", quietly = TRUE)) {
      cat("AMATRIX-OPENCL-NOT-INSTALLED\\n"); quit(save = "no", status = 0)
    }
    loadNamespace("amatrix.opencl")
    if (!exists("amatrix_opencl_availability_reason",
                envir = asNamespace("amatrix.opencl"), inherits = FALSE)) {
      # The library visible to this child holds a build that predates the
      # runtime loader; there is nothing meaningful to probe.
      cat("AMATRIX-OPENCL-STALE-INSTALL\\n"); quit(save = "no", status = 0)
    }
    # Probe disabled: must report unavailable, touch no device.
    Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "0")
    invisible(amatrix.opencl::amatrix_opencl_native_available(force = TRUE))
    # Probe force-enabled: may be TRUE or FALSE, must NEVER abort or error out.
    Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
    invisible(tryCatch(
      amatrix.opencl::amatrix_opencl_native_available(force = TRUE),
      error = function(e) FALSE
    ))
    reason <- amatrix.opencl:::amatrix_opencl_availability_reason()
    stopifnot(is.list(reason),
              is.logical(reason$opencl_loaded),
              is.logical(reason$clblast_loaded),
              is.character(reason$opencl_reason))
    cat("AMATRIX-OPENCL-LOAD-OK\\n")
  '

  old_rlibs <- Sys.getenv("R_LIBS", unset = NA)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  on.exit({
    if (is.na(old_rlibs)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old_rlibs)
  }, add = TRUE)

  out <- tryCatch(
    suppressWarnings(system2(
      rscript, c("--vanilla", "-e", shQuote(child)),
      stdout = TRUE, stderr = TRUE
    )),
    error = function(e) structure(character(), status = -1L)
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  text <- paste(out, collapse = "\n")

  if (grepl("AMATRIX-OPENCL-NOT-INSTALLED", text, fixed = TRUE)) {
    skip("amatrix.opencl is not installed in a library the subprocess can load")
  }
  if (grepl("AMATRIX-OPENCL-STALE-INSTALL", text, fixed = TRUE)) {
    skip("installed amatrix.opencl predates the runtime loader; reinstall to test")
  }

  # A C++ abort leaves these fingerprints and/or a non-zero exit status.
  expect_false(grepl("libc\\+\\+abi|terminating due to uncaught|cl::Error", text),
               info = text)
  expect_identical(as.integer(status), 0L, info = text)
  expect_true(grepl("AMATRIX-OPENCL-LOAD-OK", text, fixed = TRUE), info = text)
})

test_that("the compiled shared object has no OpenCL or CLBlast link-time dependency", {
  dll <- getLoadedDLLs()[["amatrix.opencl"]]
  skip_if(is.null(dll), "amatrix.opencl DLL not loaded in this session")
  so_path <- dll[["path"]]
  skip_if(is.null(so_path) || !nzchar(so_path) || !file.exists(so_path),
          "shared object path unavailable")

  sysname <- Sys.info()[["sysname"]]
  deps <- NULL
  if (identical(sysname, "Darwin") && nzchar(Sys.which("otool"))) {
    deps <- suppressWarnings(system2("otool", c("-L", so_path), stdout = TRUE, stderr = TRUE))
  } else if (identical(sysname, "Linux")) {
    tool <- Sys.which(c("objdump", "readelf"))
    tool <- tool[nzchar(tool)]
    if (length(tool) > 0L) {
      flag <- if (grepl("objdump", tool[[1L]])) "-p" else "-d"
      deps <- suppressWarnings(system2(tool[[1L]], c(flag, so_path), stdout = TRUE, stderr = TRUE))
    }
  }
  skip_if(is.null(deps), "no tool available to inspect shared-object dependencies")

  linked <- paste(deps, collapse = "\n")
  # After the runtime-loading refactor the object must not name any OpenCL ICD
  # loader or CLBlast library among its load-time dependencies.
  expect_false(grepl("OpenCL\\.framework|libOpenCL|OpenCL\\.dll", linked, ignore.case = TRUE),
               info = linked)
  expect_false(grepl("clblast", linked, ignore.case = TRUE), info = linked)
})
