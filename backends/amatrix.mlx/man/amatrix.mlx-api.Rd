\name{amatrix.mlx-api}
\alias{amatrix_mlx_backend}
\alias{amatrix_mlx_bridge_info}
\alias{amatrix_mlx_capabilities}
\alias{amatrix_mlx_is_available}
\alias{amatrix_mlx_native_available}
\alias{amatrix_mlx_register}
\alias{amatrix_mlx_tcrossprod_host_bridge}
\title{MLX backend API surface}
\description{
  Registration, availability probing, capability introspection, and host
  bridge entry points for the MLX backend. These functions implement the
  \pkg{amatrix} backend contract and are primarily consumed by the core
  dispatch layer: \code{amatrix_mlx_register()} registers the backend,
  \code{amatrix_mlx_is_available()} / \code{amatrix_mlx_native_available()}
  report device availability (the first probe runs in a contained child
  process), \code{amatrix_mlx_backend()} returns the backend implementation
  list, \code{amatrix_mlx_capabilities()} the supported operation set, and
  \code{amatrix_mlx_bridge_info()} build/runtime details. End users normally
  interact with the backend only through
  \code{\link[amatrix]{amatrix_use_gpu}} and
  \code{\link[amatrix]{amatrix_gpu_status}}.
}
\keyword{internal}
