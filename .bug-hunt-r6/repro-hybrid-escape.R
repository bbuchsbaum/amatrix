suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
})

mk_dense <- function() {
  x <- adgeMatrix(matrix(c(1, 4, 2, 5, 3, 6), nrow = 3, byrow = FALSE))
  dimnames(x) <- list(c("r1", "r2", "r3"), c("c1", "c2"))
  x
}

mk_sparse <- function() {
  x <- adgCMatrix(Matrix::Matrix(matrix(c(1, 0, 0, 2, 3, 0), nrow = 3), sparse = TRUE))
  dimnames(x) <- list(c("r1", "r2", "r3"), c("c1", "c2"))
  x
}

capture_probe <- function(label, expr) {
  out <- tryCatch(
    {
      value <- force(expr)
      list(
        ok = TRUE,
        class = paste(class(value), collapse = ","),
        typeof = typeof(value),
        dim = paste(dim(value) %||% integer(), collapse = "x"),
        rownames = paste(rownames(value) %||% character(), collapse = ","),
        colnames = paste(colnames(value) %||% character(), collapse = ","),
        value = utils::capture.output(str(value))
      )
    },
    error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    }
  )
  cat("\n== ", label, " ==\n", sep = "")
  if (!isTRUE(out$ok)) {
    cat("ERROR: ", out$error, "\n", sep = "")
    return(invisible(NULL))
  }
  cat("class: ", out$class, "\n", sep = "")
  cat("typeof: ", out$typeof, "\n", sep = "")
  cat("dim: ", out$dim, "\n", sep = "")
  cat("rownames: ", out$rownames, "\n", sep = "")
  cat("colnames: ", out$colnames, "\n", sep = "")
  cat(paste(out$value, collapse = "\n"), "\n", sep = "")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

x_dense <- mk_dense()
x_sparse <- mk_sparse()

capture_probe("dense::as.matrix", as.matrix(x_dense))
capture_probe("dense::methods::as('matrix')", methods::as(x_dense, "matrix"))
capture_probe("dense::as.array", as.array(x_dense))
capture_probe("dense::as.numeric", as.numeric(x_dense))
capture_probe("dense::as.vector", as.vector(x_dense))
capture_probe("dense::Matrix::Matrix", Matrix::Matrix(x_dense))
capture_probe("dense::data.matrix", data.matrix(x_dense))
capture_probe("dense::escape_matmul", as.matrix(x_dense) %*% c(1, 10))

capture_probe("sparse::as.matrix", as.matrix(x_sparse))
capture_probe("sparse::methods::as('matrix')", methods::as(x_sparse, "matrix"))
capture_probe("sparse::Matrix::Matrix", Matrix::Matrix(x_sparse))
capture_probe("sparse::data.matrix", data.matrix(x_sparse))

cat("\n== sessionInfo ==\n")
print(sessionInfo())
