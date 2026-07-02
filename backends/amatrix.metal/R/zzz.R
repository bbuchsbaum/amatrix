.onLoad <- function(libname, pkgname) {
  # Register unconditionally (cheap, metadata-only); device probing is
  # separately guarded (AMATRIX_METAL_PROBE_GPU / amatrix_use_gpu()).
  # The former options(amatrix.enable_metal) load gate made the backend
  # invisible after install; core honors amatrix.disable_metal instead.
  try(amatrix_metal_register(overwrite = TRUE), silent = TRUE)
}
