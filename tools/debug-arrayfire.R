#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(amatrix.arrayfire))

cat("bridge info:\n")
print(amatrix_arrayfire_bridge_info())

cat("\ndiagnostics:\n")
print(amatrix_arrayfire_diagnostics())

if (is.loaded("amatrix_arrayfire_native_available_bridge")) {
  cat("native available bridge loaded\n")
}

x <- matrix(c(1, 2, 3, 4), nrow = 2)

cat("backend available:", amatrix_arrayfire_is_available(), "\n")
cat("calling backend$matmul() through current safe path\n")
print(amatrix_arrayfire_backend()$matmul(x, diag(2)))
