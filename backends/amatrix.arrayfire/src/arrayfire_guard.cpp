/* arrayfire_guard.cpp — C++ firewall around every ArrayFire .Call bridge.
 *
 * The raw bridges in arrayfire_bridge.c call the ArrayFire C API (af_*)
 * directly from C. On some hosts ArrayFire's OpenCL backend throws a C++
 * exception (cl::Error from clGetDeviceIDs during lazy device enumeration)
 * that is NOT converted to an af_err at the C boundary. With no C++ catch
 * frame between R and the throw, the exception unwinds through the C bridge
 * and reaches std::terminate -> SIGABRT, killing the process.
 *
 * For every registered bridge we generate an `extern "C"` wrapper that runs
 * the raw bridge inside try/catch(...). cl::Error and af::exception both
 * derive from std::exception, and catch(...) is the final net; any escaping
 * C++ exception is converted into an ordinary R error (Rf_error) so the
 * process survives. The raw bridges' own error() calls (Rf_error/longjmp)
 * pass straight through the try block unaffected, so the success path and
 * every existing error path are preserved exactly.
 *
 * The bridge list lives in arrayfire_bridges.h so that the wrapper arity and
 * the arity registered in init.c can never drift apart.
 */

#include <cstdio>
#include <exception>

#include "arrayfire_bridges.h"

/* Scratch buffer for the caught exception message. The wrapper writes it in
 * the catch block and immediately calls Rf_error afterwards, with no R
 * evaluation in between, so a single shared buffer is safe even under nested
 * .Call invocations. */
static char amatrix_af_guard_msg[512];

#define AMATRIX_AF_MAKE_GUARD(name, arity)                                    \
  extern "C" SEXP name(AMATRIX_AF_PROTO_##arity);                             \
  extern "C" SEXP name##_guarded(AMATRIX_AF_PARAMS_##arity) {                 \
    try {                                                                     \
      return name(AMATRIX_AF_FWD_##arity);                                    \
    } catch (const std::exception& e) {                                       \
      std::snprintf(amatrix_af_guard_msg, sizeof(amatrix_af_guard_msg),       \
                    "amatrix.arrayfire: ArrayFire raised a C++ exception in "  \
                    "%s(): %s",                                               \
                    #name, e.what());                                         \
    } catch (...) {                                                           \
      std::snprintf(amatrix_af_guard_msg, sizeof(amatrix_af_guard_msg),       \
                    "amatrix.arrayfire: ArrayFire raised an unknown C++ "      \
                    "exception in %s()",                                      \
                    #name);                                                   \
    }                                                                         \
    Rf_error("%s", amatrix_af_guard_msg);                                     \
    return R_NilValue; /* unreachable: Rf_error does not return */            \
  }

AMATRIX_AF_BRIDGES(AMATRIX_AF_MAKE_GUARD)

#undef AMATRIX_AF_MAKE_GUARD
