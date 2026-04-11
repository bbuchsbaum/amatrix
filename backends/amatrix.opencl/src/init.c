#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP amatrix_opencl_native_available_bridge(void);
extern SEXP amatrix_opencl_bridge_info_bridge(void);
extern SEXP amatrix_opencl_diagnostics_bridge(void);
extern SEXP amatrix_opencl_matmul_bridge(SEXP x, SEXP y);
extern SEXP amatrix_opencl_crossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_opencl_tcrossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_opencl_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op);
extern SEXP amatrix_opencl_broadcast_ewise_bridge(SEXP lhs, SEXP v, SEXP margin, SEXP op);
extern SEXP amatrix_opencl_sum_axis_bridge(SEXP x, SEXP axis);
extern SEXP amatrix_opencl_sparse_store_bridge(SEXP key, SEXP values, SEXP p, SEXP i, SEXP dim);
extern SEXP amatrix_opencl_sparse_has_bridge(SEXP key);
extern SEXP amatrix_opencl_sparse_drop_bridge(SEXP key);
extern SEXP amatrix_opencl_spmm_bridge(SEXP values, SEXP p, SEXP i, SEXP dim, SEXP b, SEXP trans_lhs);
extern SEXP amatrix_opencl_resident_store_bridge(SEXP key, SEXP x);
extern SEXP amatrix_opencl_resident_has_bridge(SEXP key);
extern SEXP amatrix_opencl_resident_dim_bridge(SEXP key);
extern SEXP amatrix_opencl_resident_drop_bridge(SEXP key);
extern SEXP amatrix_opencl_resident_materialize_bridge(SEXP key);
extern SEXP amatrix_opencl_spmm_resident_bridge(SEXP key, SEXP b, SEXP trans_lhs);
extern SEXP amatrix_opencl_spmm_resident_key_bridge(SEXP sp_key, SEXP y_key, SEXP out_key, SEXP trans_lhs, SEXP defer);
extern SEXP amatrix_opencl_chol_resident_bridge(SEXP x_key, SEXP out_key);
extern SEXP amatrix_opencl_solve_resident_bridge(SEXP a_key, SEXP b_key, SEXP out_key);
extern SEXP amatrix_opencl_solve_triangular_resident_bridge(SEXP factor_key, SEXP rhs_key, SEXP out_key, SEXP lower, SEXP transpose);
extern SEXP amatrix_opencl_chol_solve_resident_bridge(SEXP factor_key, SEXP rhs_key, SEXP out_key);
extern SEXP amatrix_opencl_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_opencl_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_opencl_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_opencl_matmul_resident_host_bridge(SEXP x_key, SEXP y);
extern SEXP amatrix_opencl_crossprod_resident_host_bridge(SEXP x_key, SEXP y);
extern SEXP amatrix_opencl_matmul_resident_host_into_bridge(SEXP x_key, SEXP y, SEXP out_key);
extern SEXP amatrix_opencl_crossprod_resident_host_into_bridge(SEXP x_key, SEXP y, SEXP out_key);
extern SEXP amatrix_opencl_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key);
extern SEXP amatrix_opencl_broadcast_ewise_resident_bridge(SEXP lhs_key, SEXP v, SEXP margin, SEXP op, SEXP out_key);
extern SEXP amatrix_opencl_broadcast_ewise_resident_inplace_bridge(SEXP lhs_key, SEXP v, SEXP margin, SEXP op);
extern SEXP amatrix_opencl_sum_axis_resident_bridge(SEXP x_key, SEXP axis);
extern SEXP amatrix_opencl_sum_axis_resident_key_bridge(SEXP x_key, SEXP axis, SEXP out_key);

static const R_CallMethodDef call_methods[] = {
    {"amatrix_opencl_native_available_bridge", (DL_FUNC) &amatrix_opencl_native_available_bridge, 0},
    {"amatrix_opencl_bridge_info_bridge", (DL_FUNC) &amatrix_opencl_bridge_info_bridge, 0},
    {"amatrix_opencl_diagnostics_bridge", (DL_FUNC) &amatrix_opencl_diagnostics_bridge, 0},
    {"amatrix_opencl_matmul_bridge", (DL_FUNC) &amatrix_opencl_matmul_bridge, 2},
    {"amatrix_opencl_crossprod_bridge", (DL_FUNC) &amatrix_opencl_crossprod_bridge, 2},
    {"amatrix_opencl_tcrossprod_bridge", (DL_FUNC) &amatrix_opencl_tcrossprod_bridge, 2},
    {"amatrix_opencl_ewise_bridge", (DL_FUNC) &amatrix_opencl_ewise_bridge, 3},
    {"amatrix_opencl_broadcast_ewise_bridge", (DL_FUNC) &amatrix_opencl_broadcast_ewise_bridge, 4},
    {"amatrix_opencl_sum_axis_bridge", (DL_FUNC) &amatrix_opencl_sum_axis_bridge, 2},
    {"amatrix_opencl_sparse_store_bridge", (DL_FUNC) &amatrix_opencl_sparse_store_bridge, 5},
    {"amatrix_opencl_sparse_has_bridge", (DL_FUNC) &amatrix_opencl_sparse_has_bridge, 1},
    {"amatrix_opencl_sparse_drop_bridge", (DL_FUNC) &amatrix_opencl_sparse_drop_bridge, 1},
    {"amatrix_opencl_spmm_bridge", (DL_FUNC) &amatrix_opencl_spmm_bridge, 6},
    {"amatrix_opencl_resident_store_bridge", (DL_FUNC) &amatrix_opencl_resident_store_bridge, 2},
    {"amatrix_opencl_resident_has_bridge", (DL_FUNC) &amatrix_opencl_resident_has_bridge, 1},
    {"amatrix_opencl_resident_dim_bridge", (DL_FUNC) &amatrix_opencl_resident_dim_bridge, 1},
    {"amatrix_opencl_resident_drop_bridge", (DL_FUNC) &amatrix_opencl_resident_drop_bridge, 1},
    {"amatrix_opencl_resident_materialize_bridge", (DL_FUNC) &amatrix_opencl_resident_materialize_bridge, 1},
    {"amatrix_opencl_spmm_resident_bridge", (DL_FUNC) &amatrix_opencl_spmm_resident_bridge, 3},
    {"amatrix_opencl_spmm_resident_key_bridge", (DL_FUNC) &amatrix_opencl_spmm_resident_key_bridge, 5},
    {"amatrix_opencl_chol_resident_bridge", (DL_FUNC) &amatrix_opencl_chol_resident_bridge, 2},
    {"amatrix_opencl_solve_resident_bridge", (DL_FUNC) &amatrix_opencl_solve_resident_bridge, 3},
    {"amatrix_opencl_solve_triangular_resident_bridge", (DL_FUNC) &amatrix_opencl_solve_triangular_resident_bridge, 5},
    {"amatrix_opencl_chol_solve_resident_bridge", (DL_FUNC) &amatrix_opencl_chol_solve_resident_bridge, 3},
    {"amatrix_opencl_matmul_resident_bridge", (DL_FUNC) &amatrix_opencl_matmul_resident_bridge, 3},
    {"amatrix_opencl_crossprod_resident_bridge", (DL_FUNC) &amatrix_opencl_crossprod_resident_bridge, 3},
    {"amatrix_opencl_tcrossprod_resident_bridge", (DL_FUNC) &amatrix_opencl_tcrossprod_resident_bridge, 3},
    {"amatrix_opencl_matmul_resident_host_bridge", (DL_FUNC) &amatrix_opencl_matmul_resident_host_bridge, 2},
    {"amatrix_opencl_crossprod_resident_host_bridge", (DL_FUNC) &amatrix_opencl_crossprod_resident_host_bridge, 2},
    {"amatrix_opencl_matmul_resident_host_into_bridge", (DL_FUNC) &amatrix_opencl_matmul_resident_host_into_bridge, 3},
    {"amatrix_opencl_crossprod_resident_host_into_bridge", (DL_FUNC) &amatrix_opencl_crossprod_resident_host_into_bridge, 3},
    {"amatrix_opencl_ewise_resident_bridge", (DL_FUNC) &amatrix_opencl_ewise_resident_bridge, 4},
    {"amatrix_opencl_broadcast_ewise_resident_bridge", (DL_FUNC) &amatrix_opencl_broadcast_ewise_resident_bridge, 5},
    {"amatrix_opencl_broadcast_ewise_resident_inplace_bridge", (DL_FUNC) &amatrix_opencl_broadcast_ewise_resident_inplace_bridge, 4},
    {"amatrix_opencl_sum_axis_resident_bridge", (DL_FUNC) &amatrix_opencl_sum_axis_resident_bridge, 2},
    {"amatrix_opencl_sum_axis_resident_key_bridge", (DL_FUNC) &amatrix_opencl_sum_axis_resident_key_bridge, 3},
    {NULL, NULL, 0}
};

void R_init_amatrix_opencl(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
