\name{amatrix_bind_resident}
\alias{amatrix_bind_resident}
\title{Bind an amatrix Object to Resident Backend Storage}
\usage{
amatrix_bind_resident(x, backend = NULL, op = NULL, y = NULL)
}
\arguments{
  \item{x}{An \code{adgeMatrix}, \code{adgCMatrix}, base matrix, or sparse
  \code{Matrix} object.}

  \item{backend}{Backend name, \code{"auto"}, or \code{NULL}. When left
  \code{NULL} or set to \code{"auto"}, \code{amatrix} picks the first
  residency-capable accelerator backend that supports the requested resident
  operation.}

  \item{op}{Optional operation name such as \code{"matmul"} used when
  selecting an automatic resident backend.}

  \item{y}{Optional rhs object used when checking resident-op support for
  automatic backend selection.}
}
\value{
An \code{adgeMatrix} or \code{adgCMatrix} with a live resident binding on
\code{backend}. When no suitable accelerator backend is available in automatic
mode, returns \code{x} unchanged.
}
\description{
Uploads a dense or sparse matrix to a residency-capable backend and returns the
corresponding \code{aMatrix} object with a live resident binding. This is most
useful for repeated GPU work where paying the upload cost once is preferable to
relying on cold-path dispatch.
}
