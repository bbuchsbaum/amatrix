suppressPackageStartupMessages({
  library(amatrix)
})

x <- amatrix:::new_adgeMatrix_deferred(
  dim = c(2L, 3L),
  preferred_backend = "cpu"
)

run_probe <- function(label, fn) {
  cat("\n== ", label, " ==\n", sep = "")
  out <- tryCatch(fn(x), error = function(e) e)
  if (inherits(out, "error")) {
    cat("ERROR: ", conditionMessage(out), "\n", sep = "")
    return(invisible(NULL))
  }
  print(out)
  cat("class: ", paste(class(out), collapse = ","), "\n", sep = "")
  cat("anyNaN: ", any(is.nan(out)), "\n", sep = "")
}

run_probe("as.matrix", as.matrix)
run_probe("as.array", as.array)
run_probe("as.numeric", as.numeric)
run_probe("as.vector", as.vector)

cat("\n== sessionInfo ==\n")
print(sessionInfo())
