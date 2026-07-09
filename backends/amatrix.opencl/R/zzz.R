.onLoad <- function(libname, pkgname) {
  # Point the native runtime loader at the package's user data directory so a
  # CLBlast library installed by amatrix_install_clblast() is found at probe
  # time. Pure string bookkeeping; touches no device and never probes.
  try(.amatrix_opencl_register_clblast_dir(), silent = TRUE)
  try(amatrix_opencl_register(overwrite = TRUE), silent = TRUE)
}
