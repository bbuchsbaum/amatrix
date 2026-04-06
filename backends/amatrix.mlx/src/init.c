#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP amatrix_mlx_matmul_bridge(SEXP x, SEXP y);
extern SEXP amatrix_mlx_crossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_mlx_tcrossprod_bridge(SEXP x, SEXP y);
extern SEXP amatrix_mlx_solve_triangular_bridge(SEXP a, SEXP b, SEXP upper);
extern SEXP amatrix_mlx_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op);
extern SEXP amatrix_mlx_sum_axis_bridge(SEXP x, SEXP axis);
extern SEXP amatrix_mlx_qr_bridge(SEXP x, SEXP q_key);
extern SEXP amatrix_mlx_tsqr_build_bridge(SEXP x, SEXP block_rows, SEXP q_keys, SEXP top_q_key, SEXP top_r_key, SEXP r_stack_key);
extern SEXP amatrix_mlx_qr_qty_key_bridge(SEXP q_key, SEXP y);
extern SEXP amatrix_mlx_qr_qy_key_bridge(SEXP q_key, SEXP y);
extern SEXP amatrix_mlx_qr_coef_key_bridge(SEXP q_key, SEXP r, SEXP y);
extern SEXP amatrix_mlx_tsqr_coef_resident_bridge(SEXP q_keys, SEXP block_rows, SEXP top_q_key, SEXP top_r_key, SEXP y);
extern SEXP amatrix_mlx_tsqr_coef_key_bridge(SEXP q_keys, SEXP block_rows, SEXP top_q_key, SEXP r, SEXP y);
extern SEXP amatrix_mlx_native_available_bridge(void);
extern SEXP amatrix_mlx_bridge_info_bridge(void);
extern SEXP amatrix_mlx_resident_has_bridge(SEXP key);
extern SEXP amatrix_mlx_resident_store_bridge(SEXP key, SEXP x);
extern SEXP amatrix_mlx_resident_drop_bridge(SEXP key);
extern SEXP amatrix_mlx_resident_materialize_bridge(SEXP key);
extern SEXP amatrix_mlx_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_mlx_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_mlx_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_mlx_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key);
extern SEXP amatrix_mlx_rsvd_bridge(SEXP x_r, SEXP k_r, SEXP n_oversamples_r, SEXP n_iter_r);
extern SEXP amatrix_mlx_chol_solve_bridge(SEXP A_r, SEXP B_r);
extern SEXP amatrix_mlx_chol_factor_bridge(SEXP X_r);
extern SEXP amatrix_mlx_eigh_bridge(SEXP A_r);

static const R_CallMethodDef call_methods[] = {
    {"amatrix_mlx_native_available_bridge", (DL_FUNC) &amatrix_mlx_native_available_bridge, 0},
    {"amatrix_mlx_bridge_info_bridge", (DL_FUNC) &amatrix_mlx_bridge_info_bridge, 0},
    {"amatrix_mlx_matmul_bridge", (DL_FUNC) &amatrix_mlx_matmul_bridge, 2},
    {"amatrix_mlx_crossprod_bridge", (DL_FUNC) &amatrix_mlx_crossprod_bridge, 2},
    {"amatrix_mlx_tcrossprod_bridge", (DL_FUNC) &amatrix_mlx_tcrossprod_bridge, 2},
    {"amatrix_mlx_solve_triangular_bridge", (DL_FUNC) &amatrix_mlx_solve_triangular_bridge, 3},
    {"amatrix_mlx_ewise_bridge", (DL_FUNC) &amatrix_mlx_ewise_bridge, 3},
    {"amatrix_mlx_sum_axis_bridge", (DL_FUNC) &amatrix_mlx_sum_axis_bridge, 2},
    {"amatrix_mlx_qr_bridge", (DL_FUNC) &amatrix_mlx_qr_bridge, 2},
    {"amatrix_mlx_tsqr_build_bridge", (DL_FUNC) &amatrix_mlx_tsqr_build_bridge, 6},
    {"amatrix_mlx_qr_qty_key_bridge", (DL_FUNC) &amatrix_mlx_qr_qty_key_bridge, 2},
    {"amatrix_mlx_qr_qy_key_bridge", (DL_FUNC) &amatrix_mlx_qr_qy_key_bridge, 2},
    {"amatrix_mlx_qr_coef_key_bridge", (DL_FUNC) &amatrix_mlx_qr_coef_key_bridge, 3},
    {"amatrix_mlx_tsqr_coef_resident_bridge", (DL_FUNC) &amatrix_mlx_tsqr_coef_resident_bridge, 5},
    {"amatrix_mlx_tsqr_coef_key_bridge", (DL_FUNC) &amatrix_mlx_tsqr_coef_key_bridge, 5},
    {"amatrix_mlx_resident_has_bridge", (DL_FUNC) &amatrix_mlx_resident_has_bridge, 1},
    {"amatrix_mlx_resident_store_bridge", (DL_FUNC) &amatrix_mlx_resident_store_bridge, 2},
    {"amatrix_mlx_resident_drop_bridge", (DL_FUNC) &amatrix_mlx_resident_drop_bridge, 1},
    {"amatrix_mlx_resident_materialize_bridge", (DL_FUNC) &amatrix_mlx_resident_materialize_bridge, 1},
    {"amatrix_mlx_matmul_resident_bridge", (DL_FUNC) &amatrix_mlx_matmul_resident_bridge, 3},
    {"amatrix_mlx_crossprod_resident_bridge", (DL_FUNC) &amatrix_mlx_crossprod_resident_bridge, 3},
    {"amatrix_mlx_tcrossprod_resident_bridge", (DL_FUNC) &amatrix_mlx_tcrossprod_resident_bridge, 3},
    {"amatrix_mlx_ewise_resident_bridge", (DL_FUNC) &amatrix_mlx_ewise_resident_bridge, 4},
    {"amatrix_mlx_rsvd_bridge",           (DL_FUNC) &amatrix_mlx_rsvd_bridge,           4},
    {"amatrix_mlx_chol_solve_bridge",     (DL_FUNC) &amatrix_mlx_chol_solve_bridge,     2},
    {"amatrix_mlx_chol_factor_bridge",    (DL_FUNC) &amatrix_mlx_chol_factor_bridge,    1},
    {"amatrix_mlx_eigh_bridge",           (DL_FUNC) &amatrix_mlx_eigh_bridge,           1},
    {NULL, NULL, 0}
};

void R_init_amatrix_mlx(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
