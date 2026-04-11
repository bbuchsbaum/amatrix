\name{amatrix_resident_backend_for}
\alias{amatrix_resident_backend_for}
\title{Choose a Residency-Capable Accelerator Backend}
\usage{
amatrix_resident_backend_for(x, op = NULL, y = NULL)
}
\arguments{
  \item{x}{An \code{aMatrix}.}

  \item{op}{Optional operation name such as \code{"matmul"}.}

  \item{y}{Optional rhs object used when checking resident-op support.}
}
\value{
A backend name, or \code{NULL} when no residency-capable accelerator backend is
available for the requested operation.
}
\description{
Returns the first available non-CPU backend that can keep \code{x} resident for
the requested operation. This is intended for package authors who want repeated
work to stay on the fastest available accelerator without hardcoding backend
names such as \code{"metal"} or \code{"mlx"}.
}
