.onLoad <- function(libname, pkgname) {
  if (isTRUE(getOption("amatrix.enable_arrayfire", FALSE))) {
    amatrix_arrayfire_register(overwrite = TRUE)
  }
}
