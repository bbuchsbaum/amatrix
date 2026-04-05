#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP amatrix_arrayfire_native_available_bridge(void);
extern SEXP amatrix_arrayfire_bridge_info_bridge(void);
extern SEXP amatrix_arrayfire_diagnostics_bridge(void);
extern SEXP amatrix_arrayfire_set_backend_bridge(SEXP backend);
extern SEXP amatrix_arrayfire_matmul_bridge(SEXP x, SEXP y);
extern SEXP amatrix_arrayfire_crossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_arrayfire_tcrossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_arrayfire_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op);
extern SEXP amatrix_arrayfire_sum_axis_bridge(SEXP x, SEXP axis);
extern SEXP amatrix_arrayfire_qr_bridge(SEXP x);

static const R_CallMethodDef call_methods[] = {
    {"amatrix_arrayfire_native_available_bridge", (DL_FUNC) &amatrix_arrayfire_native_available_bridge, 0},
    {"amatrix_arrayfire_bridge_info_bridge", (DL_FUNC) &amatrix_arrayfire_bridge_info_bridge, 0},
    {"amatrix_arrayfire_diagnostics_bridge", (DL_FUNC) &amatrix_arrayfire_diagnostics_bridge, 0},
    {"amatrix_arrayfire_set_backend_bridge", (DL_FUNC) &amatrix_arrayfire_set_backend_bridge, 1},
    {"amatrix_arrayfire_matmul_bridge", (DL_FUNC) &amatrix_arrayfire_matmul_bridge, 2},
    {"amatrix_arrayfire_crossprod_bridge", (DL_FUNC) &amatrix_arrayfire_crossprod_bridge, 2},
    {"amatrix_arrayfire_tcrossprod_bridge", (DL_FUNC) &amatrix_arrayfire_tcrossprod_bridge, 2},
    {"amatrix_arrayfire_ewise_bridge", (DL_FUNC) &amatrix_arrayfire_ewise_bridge, 3},
    {"amatrix_arrayfire_sum_axis_bridge", (DL_FUNC) &amatrix_arrayfire_sum_axis_bridge, 2},
    {"amatrix_arrayfire_qr_bridge", (DL_FUNC) &amatrix_arrayfire_qr_bridge, 1},
    {NULL, NULL, 0}
};

void R_init_amatrix_arrayfire(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
