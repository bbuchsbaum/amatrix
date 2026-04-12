.onLoad <- function(libname, pkgname) {
  try(amatrix_opencl_register(overwrite = TRUE), silent = TRUE)
}
