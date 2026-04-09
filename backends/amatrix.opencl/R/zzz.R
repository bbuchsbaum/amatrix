.onLoad <- function(libname, pkgname) {
  if (isTRUE(getOption("amatrix.enable_opencl", FALSE))) {
    amatrix_opencl_register(overwrite = TRUE)
  }
}
