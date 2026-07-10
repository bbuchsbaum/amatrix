\name{amatrix.opencl-api}
\alias{amatrix_opencl_backend}
\alias{amatrix_opencl_bridge_info}
\alias{amatrix_opencl_capabilities}
\alias{amatrix_opencl_diagnostics}
\alias{amatrix_opencl_enable_probe}
\alias{amatrix_opencl_init_status}
\alias{amatrix_opencl_is_available}
\alias{amatrix_opencl_native_available}
\alias{amatrix_opencl_register}
\title{OpenCL backend API surface}
\description{
  Registration, availability probing, and diagnostics for the OpenCL/CLBlast
  backend. These functions implement the \pkg{amatrix} backend contract and
  are primarily consumed by the core dispatch layer. OpenCL and CLBlast are
  located at run time; \code{amatrix_opencl_diagnostics()} and
  \code{amatrix_opencl_init_status()} report which libraries were found and
  why the backend is or is not available. End users normally interact with
  the backend only through \code{\link[amatrix]{amatrix_use_gpu}},
  \code{\link[amatrix]{amatrix_gpu_status}}, and
  \code{\link{amatrix_install_clblast}}.
}
\keyword{internal}
