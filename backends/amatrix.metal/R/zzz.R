.onLoad <- function(libname, pkgname) {
  if (isTRUE(getOption("amatrix.enable_metal", FALSE))) {
    amatrix_metal_register(overwrite = TRUE)
  }
}
