#ifndef AMATRIX_ARRAYFIRE_BRIDGES_H
#define AMATRIX_ARRAYFIRE_BRIDGES_H

/* Single source of truth for every .Call bridge exported by this package.
 *
 * AMATRIX_AF_BRIDGES(X) expands X(name, arity) once per registered bridge,
 * where `name` is the R-visible .Call symbol / raw C implementation and
 * `arity` is its number of SEXP arguments.
 *
 * Two translation units consume this list:
 *   - arrayfire_guard.cpp generates, for every bridge, an `extern "C"`
 *     C++ wrapper `<name>_guarded` that runs the raw C bridge inside a
 *     try/catch(...). ArrayFire's OpenCL backend can throw a C++ exception
 *     (cl::Error from clGetDeviceIDs during device enumeration) that would
 *     otherwise unwind through the C bridge and hit std::terminate/SIGABRT.
 *     The wrapper converts any escaping C++ exception into an ordinary R
 *     error instead of crashing the process.
 *   - init.c registers each `<name>_guarded` under the original R name, so
 *     every entry point from R into ArrayFire code passes through the
 *     firewall. Keeping the list here guarantees the wrapper's arity and
 *     the registered arity can never drift apart.
 */
#include <Rinternals.h>

/* Prototype argument lists (types only) for extern declarations. */
#define AMATRIX_AF_PROTO_0 void
#define AMATRIX_AF_PROTO_1 SEXP
#define AMATRIX_AF_PROTO_2 SEXP, SEXP
#define AMATRIX_AF_PROTO_3 SEXP, SEXP, SEXP
#define AMATRIX_AF_PROTO_4 SEXP, SEXP, SEXP, SEXP
#define AMATRIX_AF_PROTO_5 SEXP, SEXP, SEXP, SEXP, SEXP
#define AMATRIX_AF_PROTO_6 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
#define AMATRIX_AF_PROTO_7 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
#define AMATRIX_AF_PROTO_8 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP

/* Named parameter lists for wrapper definitions (arrayfire_guard.cpp). */
#define AMATRIX_AF_PARAMS_0 void
#define AMATRIX_AF_PARAMS_1 SEXP a0
#define AMATRIX_AF_PARAMS_2 SEXP a0, SEXP a1
#define AMATRIX_AF_PARAMS_3 SEXP a0, SEXP a1, SEXP a2
#define AMATRIX_AF_PARAMS_4 SEXP a0, SEXP a1, SEXP a2, SEXP a3
#define AMATRIX_AF_PARAMS_5 SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4
#define AMATRIX_AF_PARAMS_6 SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5
#define AMATRIX_AF_PARAMS_7 SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5, SEXP a6
#define AMATRIX_AF_PARAMS_8 SEXP a0, SEXP a1, SEXP a2, SEXP a3, SEXP a4, SEXP a5, SEXP a6, SEXP a7

/* Forwarded argument lists for wrapper definitions (arrayfire_guard.cpp). */
#define AMATRIX_AF_FWD_0
#define AMATRIX_AF_FWD_1 a0
#define AMATRIX_AF_FWD_2 a0, a1
#define AMATRIX_AF_FWD_3 a0, a1, a2
#define AMATRIX_AF_FWD_4 a0, a1, a2, a3
#define AMATRIX_AF_FWD_5 a0, a1, a2, a3, a4
#define AMATRIX_AF_FWD_6 a0, a1, a2, a3, a4, a5
#define AMATRIX_AF_FWD_7 a0, a1, a2, a3, a4, a5, a6
#define AMATRIX_AF_FWD_8 a0, a1, a2, a3, a4, a5, a6, a7

#define AMATRIX_AF_BRIDGES(X)                                       \
  X(amatrix_arrayfire_native_available_bridge,        0)            \
  X(amatrix_arrayfire_bridge_info_bridge,             0)            \
  X(amatrix_arrayfire_diagnostics_bridge,             0)            \
  X(amatrix_arrayfire_set_backend_bridge,             1)            \
  X(amatrix_arrayfire_matmul_bridge,                  2)            \
  X(amatrix_arrayfire_crossprod_bridge,               2)            \
  X(amatrix_arrayfire_tcrossprod_bridge,              2)            \
  X(amatrix_arrayfire_ewise_bridge,                   3)            \
  X(amatrix_arrayfire_sum_axis_bridge,                2)            \
  X(amatrix_arrayfire_qr_bridge,                      1)            \
  X(amatrix_arrayfire_resident_store_bridge,          2)            \
  X(amatrix_arrayfire_resident_has_bridge,            1)            \
  X(amatrix_arrayfire_resident_drop_bridge,           1)            \
  X(amatrix_arrayfire_resident_materialize_bridge,    1)            \
  X(amatrix_arrayfire_matmul_resident_bridge,         3)            \
  X(amatrix_arrayfire_crossprod_resident_bridge,      3)            \
  X(amatrix_arrayfire_tcrossprod_resident_bridge,     3)            \
  X(amatrix_arrayfire_ewise_resident_bridge,          4)            \
  X(amatrix_arrayfire_broadcast_ewise_resident_bridge, 5)           \
  X(amatrix_arrayfire_argreduce_bridge,               3)            \
  X(amatrix_arrayfire_scatter_mean_bridge,            3)            \
  X(amatrix_arrayfire_segment_sum_bridge,             4)            \
  X(amatrix_arrayfire_segment_mean_bridge,            4)            \
  X(amatrix_arrayfire_sum_axis_resident_bridge,       2)            \
  X(am_af_lanczos_bidiag_bridge,                      3)            \
  X(am_af_lbz_upload_A_bridge,                        1)            \
  X(am_af_lbz_drop_A_bridge,                          0)            \
  X(am_af_lanczos_warm_bridge,                        4)            \
  X(am_af_dist_sq_bridge,                             2)            \
  X(am_af_kernel_bridge,                              6)            \
  X(am_af_kernel_resident_bridge,                     8)            \
  X(amatrix_arrayfire_matmul_correct_bridge,          2)            \
  X(amatrix_arrayfire_crossprod_correct_bridge,       2)            \
  X(amatrix_arrayfire_tcrossprod_correct_bridge,      2)            \
  X(amatrix_arrayfire_qr_q_correct_bridge,            1)            \
  X(amatrix_arrayfire_chol_bridge,                    1)            \
  X(amatrix_arrayfire_chol_resident_bridge,           2)            \
  X(amatrix_arrayfire_solve_bridge,                   2)            \
  X(amatrix_arrayfire_solve_resident_bridge,          3)            \
  X(amatrix_arrayfire_qr_Q_resident_bridge,           2)            \
  X(amatrix_arrayfire_svd_bridge,                     3)            \
  X(amatrix_arrayfire_svd_safe_bridge,                0)            \
  X(amatrix_arrayfire_bdc_bidiag_bridge,              1)            \
  X(amatrix_arrayfire_bdc_orgbr_bridge,               6)            \
  X(amatrix_arrayfire_bdc_dbdsdc_bridge,              3)            \
  X(amatrix_arrayfire_spmm_bridge,                    6)            \
  X(amatrix_arrayfire_sparse_store_bridge,            6)            \
  X(amatrix_arrayfire_sparse_has_bridge,              1)            \
  X(amatrix_arrayfire_sparse_drop_bridge,             1)            \
  X(amatrix_arrayfire_spmm_resident_bridge,           3)

#endif /* AMATRIX_ARRAYFIRE_BRIDGES_H */
