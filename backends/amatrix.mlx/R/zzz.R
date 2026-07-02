.onLoad <- function(libname, pkgname) {
  # Register unconditionally (cheap, metadata-only); device probing is
  # separately guarded by the probe policy in the core registry.
  try(amatrix_mlx_register(overwrite = TRUE), silent = TRUE)
}
