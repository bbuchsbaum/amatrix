#' Internal S3 methods for amatrix helper classes
#'
#' These S3 methods implement standard base generics (\code{as.matrix},
#' \code{dim}, \code{nrow}, \code{ncol}) for internal amatrix helper
#' classes (\code{KronMatrix}, \code{resident_handle}). They are not
#' part of the public user-facing API — use the generics directly. This
#' help page exists only to satisfy R CMD check.
#'
#' @param x A \code{KronMatrix} or \code{resident_handle} object.
#' @param ... Additional arguments passed to base methods.
#'
#' @return For \code{as.matrix} methods: a plain R \code{matrix}. For
#'   \code{dim}, \code{nrow}, \code{ncol}: an integer (or length-2
#'   integer vector for \code{dim}).
#'
#' @name amatrix-s3-methods
#' @keywords internal
NULL
