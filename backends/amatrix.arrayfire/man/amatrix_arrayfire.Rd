\name{amatrix_arrayfire_register}
\alias{amatrix_arrayfire_register}
\alias{amatrix_arrayfire_is_available}
\alias{amatrix_arrayfire_native_available}
\alias{amatrix_arrayfire_capabilities}
\alias{amatrix_arrayfire_diagnostics}
\alias{amatrix_arrayfire_enable_probe}
\alias{amatrix_arrayfire_bridge_info}
\alias{amatrix_arrayfire_backend}
\title{ArrayFire Backend for amatrix}
\description{
  Functions for registering and inspecting the ArrayFire GPU backend.
  The backend is optional: if ArrayFire is not installed the package
  loads but all GPU operations fall back to CPU.
}
\usage{
amatrix_arrayfire_register(overwrite = TRUE)
amatrix_arrayfire_enable_probe(register = TRUE)
amatrix_arrayfire_is_available()
amatrix_arrayfire_native_available()
amatrix_arrayfire_capabilities()
amatrix_arrayfire_diagnostics()
amatrix_arrayfire_bridge_info()
amatrix_arrayfire_backend()
}
\arguments{
  \item{overwrite}{logical; whether to overwrite a previously registered
    backend of the same name (default \code{TRUE}).}
  \item{register}{logical; whether to register the backend with
    \pkg{amatrix} after enabling the runtime probe (default \code{TRUE}).}
}
\value{
  \code{amatrix_arrayfire_register} returns the backend name invisibly.

  \code{amatrix_arrayfire_enable_probe} returns the probed availability
  invisibly.

  \code{amatrix_arrayfire_is_available} and
  \code{amatrix_arrayfire_native_available} return a scalar logical.

  \code{amatrix_arrayfire_capabilities} returns a character vector of
  supported operation names.

  \code{amatrix_arrayfire_diagnostics} and \code{amatrix_arrayfire_bridge_info}
  return named lists with diagnostic information.

  \code{amatrix_arrayfire_backend} returns the backend contract list used
  internally by \pkg{amatrix}.
}
\details{
  Call \code{amatrix_arrayfire_register()} once per session (or place it in
  \file{.Rprofile}) to make the ArrayFire backend available for dispatch.
  ArrayFire registration is opt-in; loading \pkg{amatrix.arrayfire} no longer
  auto-registers the backend unless \code{options(amatrix.enable_arrayfire = TRUE)}
  is set before package load.

  \code{amatrix_arrayfire_enable_probe()} activates the runtime probe for the
  current session by setting \env{AMATRIX_ARRAYFIRE_PROBE_GPU=1}, enables
  optional registration, and returns the probed availability invisibly.

  \code{amatrix_arrayfire_is_available()} returns \code{TRUE} when ArrayFire
  was found at build time and initialises successfully at runtime. Availability
  checks are probe-gated so ordinary package load and backend listing stay safe.

  \code{amatrix_arrayfire_diagnostics()} returns details about the active
  compute backend, device count, and LAPACK availability.

  The \code{ARRAYFIRE_PREFIX} environment variable can be set before
  installation to point to a non-standard ArrayFire location.
}
\seealso{
  \code{amatrix::amatrix_backend_names()},
  \code{amatrix::amatrix_backend_status()}
}
