#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "arrayfire_bridges.h"

/* Each raw bridge is fronted by a C++ firewall wrapper `<name>_guarded`
 * defined in arrayfire_guard.cpp. We register the guarded wrapper under the
 * bridge's original R name, so every .Call entry into ArrayFire code passes
 * through try/catch(...) and a crashing OpenCL driver surfaces as an R error
 * instead of terminating the process. The shared bridge list keeps the
 * registered arity in lockstep with the wrapper definitions. */
#define AMATRIX_AF_DECLARE_GUARD(name, arity) \
  extern SEXP name##_guarded(AMATRIX_AF_PROTO_##arity);
AMATRIX_AF_BRIDGES(AMATRIX_AF_DECLARE_GUARD)
#undef AMATRIX_AF_DECLARE_GUARD

static const R_CallMethodDef call_methods[] = {
#define AMATRIX_AF_REGISTER_GUARD(name, arity) \
  {#name, (DL_FUNC) &name##_guarded, arity},
    AMATRIX_AF_BRIDGES(AMATRIX_AF_REGISTER_GUARD)
#undef AMATRIX_AF_REGISTER_GUARD
    {NULL, NULL, 0}
};

void R_init_amatrix_arrayfire(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
