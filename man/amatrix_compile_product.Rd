\name{amatrix_compile_product}
\alias{amatrix_compile_product}
\title{Compile a Reusable Matrix-Product Plan}
\usage{
amatrix_compile_product(
  x,
  op = c("matmul", "crossprod", "tcrossprod"),
  backend = "auto",
  precision = amatrix_default_precision(),
  policy = amatrix_default_policy()
)
}
\arguments{
  \item{x}{Fixed left operand.}

  \item{op}{Product primitive: \code{"matmul"}, \code{"crossprod"}, or
  \code{"tcrossprod"}.}

  \item{backend}{Backend name or \code{"auto"}.}

  \item{precision}{Precision to use when wrapping base matrices.}

  \item{policy}{Policy to use when wrapping base matrices.}
}
\value{
A callable \code{am_product_plan}.
}
\description{
Prepares a fixed left operand for repeated products, choosing and binding a
resident accelerator backend when beneficial. The returned object is a callable
function, so repeated right-hand sides can be applied without rethinking
backend selection each time.
}
