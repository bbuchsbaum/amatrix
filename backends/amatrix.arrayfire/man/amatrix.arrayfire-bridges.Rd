\name{amatrix.arrayfire-bridges}
\alias{amatrix_arrayfire_dist_sq_host_bridge}
\alias{amatrix_arrayfire_kernel_host_bridge}
\alias{amatrix_arrayfire_lanczos_drop_bridge}
\alias{amatrix_arrayfire_lanczos_upload_bridge}
\alias{amatrix_arrayfire_lanczos_warm_bridge}
\alias{amatrix_arrayfire_tcrossprod_host_bridge}
\title{ArrayFire host bridge entry points}
\description{
  Low-level host bridges used by the core \pkg{amatrix} dispatch layer for
  pairwise distance, kernel matrix, Lanczos residency, and transposed
  cross-product paths on the ArrayFire backend. Not intended for direct use;
  see \code{\link{amatrix_arrayfire_register}} and
  \code{\link[amatrix]{amatrix_use_gpu}}.
}
\keyword{internal}
