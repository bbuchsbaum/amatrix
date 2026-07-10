\name{amatrix.metal-api}
\alias{amatrix_metal_backend}
\alias{amatrix_metal_bridge_info}
\alias{amatrix_metal_capabilities}
\alias{amatrix_metal_enable_probe}
\alias{amatrix_metal_is_available}
\alias{amatrix_metal_native_available}
\alias{amatrix_metal_register}
\title{Metal backend API surface}
\description{
  Registration, availability probing, and capability introspection for the
  experimental Metal sparse backend. These functions implement the
  \pkg{amatrix} backend contract and are primarily consumed by the core
  dispatch layer; \code{amatrix_metal_enable_probe()} activates the explicit
  device probe. End users normally interact with the backend only through
  \code{\link[amatrix]{amatrix_use_gpu}} and
  \code{\link[amatrix]{amatrix_gpu_status}}.
}
\keyword{internal}
