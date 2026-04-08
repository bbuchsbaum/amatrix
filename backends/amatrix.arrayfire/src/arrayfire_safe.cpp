#include <cstdlib>
#include <cstring>

#ifdef HAVE_ARRAYFIRE
#include <arrayfire.h>
#include <exception>
#endif

static bool amatrix_arrayfire_probe_enabled_cpp() {
  const char* probe = std::getenv("AMATRIX_ARRAYFIRE_PROBE_GPU");
  return probe != nullptr && std::strcmp(probe, "1") == 0;
}

static void amatrix_arrayfire_zero_diag_cpp(int* init_ok, int* backends,
                                            int* devices, int* active_backend,
                                            int* lapack_available) {
  if (init_ok != nullptr) {
    *init_ok = 0;
  }
  if (backends != nullptr) {
    *backends = 0;
  }
  if (devices != nullptr) {
    *devices = 0;
  }
  if (active_backend != nullptr) {
    *active_backend = 0;
  }
  if (lapack_available != nullptr) {
    *lapack_available = 0;
  }
}

extern "C" int amatrix_arrayfire_safe_native_available_cpp(void) {
#ifdef HAVE_ARRAYFIRE
  if (!amatrix_arrayfire_probe_enabled_cpp()) {
    return 0;
  }

  try {
    return af::getDeviceCount() > 0 ? 1 : 0;
  } catch (const af::exception&) {
    return 0;
  } catch (const std::exception&) {
    return 0;
  } catch (...) {
    return 0;
  }
#else
  return 0;
#endif
}

extern "C" int amatrix_arrayfire_safe_diagnostics_cpp(int* init_ok, int* backends,
                                                      int* devices, int* active_backend,
                                                      int* lapack_available) {
  amatrix_arrayfire_zero_diag_cpp(init_ok, backends, devices, active_backend,
                                  lapack_available);

#ifdef HAVE_ARRAYFIRE
  if (!amatrix_arrayfire_probe_enabled_cpp()) {
    return 0;
  }

  try {
    if (backends != nullptr) {
      *backends = af::getAvailableBackends();
    }
    if (devices != nullptr) {
      *devices = af::getDeviceCount();
    }
    if (active_backend != nullptr) {
      *active_backend = static_cast<int>(af::getActiveBackend());
    }
    if (lapack_available != nullptr) {
      *lapack_available = af::isLAPACKAvailable() ? 1 : 0;
    }
    if (init_ok != nullptr) {
      *init_ok = 1;
    }
    return 1;
  } catch (const af::exception&) {
    return 0;
  } catch (const std::exception&) {
    return 0;
  } catch (...) {
    return 0;
  }
#else
  return 0;
#endif
}

extern "C" int amatrix_arrayfire_safe_set_backend_cpp(int backend_id) {
#ifdef HAVE_ARRAYFIRE
  try {
    af::setBackend(static_cast<af::Backend>(backend_id));
    return 1;
  } catch (const af::exception&) {
    return 0;
  } catch (const std::exception&) {
    return 0;
  } catch (...) {
    return 0;
  }
#else
  (void) backend_id;
  return 0;
#endif
}
