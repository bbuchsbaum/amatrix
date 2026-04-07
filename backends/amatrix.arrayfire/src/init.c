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
extern SEXP amatrix_arrayfire_resident_store_bridge(SEXP key, SEXP x);
extern SEXP amatrix_arrayfire_resident_has_bridge(SEXP key);
extern SEXP amatrix_arrayfire_resident_drop_bridge(SEXP key);
extern SEXP amatrix_arrayfire_resident_materialize_bridge(SEXP key);
extern SEXP amatrix_arrayfire_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_arrayfire_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_arrayfire_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key);
extern SEXP amatrix_arrayfire_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key);
extern SEXP amatrix_arrayfire_broadcast_ewise_resident_bridge(SEXP lhs_key, SEXP v, SEXP margin_r, SEXP op, SEXP out_key);
extern SEXP amatrix_arrayfire_argreduce_bridge(SEXP lhs_key, SEXP axis_r, SEXP is_max_r);
extern SEXP amatrix_arrayfire_scatter_mean_bridge(SEXP lhs_key, SEXP labels_r, SEXP K_r);
extern SEXP amatrix_arrayfire_segment_sum_bridge(SEXP x_key, SEXP labels_r, SEXP K_r, SEXP out_key);
extern SEXP amatrix_arrayfire_segment_mean_bridge(SEXP x_key, SEXP labels_r, SEXP K_r, SEXP out_key);
extern SEXP amatrix_arrayfire_sum_axis_resident_bridge(SEXP x_key, SEXP axis);
extern SEXP am_af_lanczos_bidiag_bridge(SEXP A_r, SEXP v0_r, SEXP k_r);
extern SEXP am_af_lbz_upload_A_bridge(SEXP A_r);
extern SEXP am_af_lbz_drop_A_bridge(void);
extern SEXP am_af_lanczos_warm_bridge(SEXP V_warm_r, SEXP U_warm_r, SEXP p0_r, SEXP k_r);
extern SEXP am_af_dist_sq_bridge(SEXP X_r, SEXP Y_r);
extern SEXP am_af_kernel_bridge(SEXP X_r, SEXP Y_r, SEXP kernel_r, SEXP sigma_r, SEXP degree_r, SEXP coef_r);
extern SEXP amatrix_arrayfire_matmul_correct_bridge(SEXP A_r, SEXP B_r);
extern SEXP amatrix_arrayfire_crossprod_correct_bridge(SEXP A_r, SEXP B_r);
extern SEXP amatrix_arrayfire_tcrossprod_correct_bridge(SEXP A_r, SEXP B_r);
extern SEXP amatrix_arrayfire_qr_q_correct_bridge(SEXP A_r);
extern SEXP amatrix_arrayfire_chol_bridge(SEXP x);
extern SEXP amatrix_arrayfire_chol_resident_bridge(SEXP x_key, SEXP out_key);
extern SEXP amatrix_arrayfire_solve_bridge(SEXP a, SEXP b);
extern SEXP amatrix_arrayfire_solve_resident_bridge(SEXP a_key, SEXP b_key, SEXP out_key);
extern SEXP amatrix_arrayfire_qr_Q_resident_bridge(SEXP x_key, SEXP q_out_key);
extern SEXP amatrix_arrayfire_svd_bridge(SEXP x, SEXP nu_r, SEXP nv_r);
extern SEXP amatrix_arrayfire_svd_safe_bridge(void);
extern SEXP amatrix_arrayfire_bdc_bidiag_bridge(SEXP A_r);
extern SEXP amatrix_arrayfire_bdc_orgbr_bridge(SEXP vect_r, SEXP A_r, SEXP tau_r, SEXP M_r, SEXP N_r, SEXP K_r);
extern SEXP amatrix_arrayfire_bdc_dbdsdc_bridge(SEXP d_r, SEXP e_r, SEXP uplo_r);

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
    {"amatrix_arrayfire_resident_store_bridge",        (DL_FUNC) &amatrix_arrayfire_resident_store_bridge,        2},
    {"amatrix_arrayfire_resident_has_bridge",          (DL_FUNC) &amatrix_arrayfire_resident_has_bridge,          1},
    {"amatrix_arrayfire_resident_drop_bridge",         (DL_FUNC) &amatrix_arrayfire_resident_drop_bridge,         1},
    {"amatrix_arrayfire_resident_materialize_bridge",  (DL_FUNC) &amatrix_arrayfire_resident_materialize_bridge,  1},
    {"amatrix_arrayfire_matmul_resident_bridge",       (DL_FUNC) &amatrix_arrayfire_matmul_resident_bridge,       3},
    {"amatrix_arrayfire_crossprod_resident_bridge",    (DL_FUNC) &amatrix_arrayfire_crossprod_resident_bridge,    3},
    {"amatrix_arrayfire_tcrossprod_resident_bridge",   (DL_FUNC) &amatrix_arrayfire_tcrossprod_resident_bridge,   3},
    {"amatrix_arrayfire_ewise_resident_bridge",                  (DL_FUNC) &amatrix_arrayfire_ewise_resident_bridge,                  4},
    {"amatrix_arrayfire_broadcast_ewise_resident_bridge",        (DL_FUNC) &amatrix_arrayfire_broadcast_ewise_resident_bridge,        5},
    {"amatrix_arrayfire_argreduce_bridge",                       (DL_FUNC) &amatrix_arrayfire_argreduce_bridge,                       3},
    {"amatrix_arrayfire_scatter_mean_bridge",                    (DL_FUNC) &amatrix_arrayfire_scatter_mean_bridge,                    3},
    {"amatrix_arrayfire_segment_sum_bridge",                     (DL_FUNC) &amatrix_arrayfire_segment_sum_bridge,                     4},
    {"amatrix_arrayfire_segment_mean_bridge",                    (DL_FUNC) &amatrix_arrayfire_segment_mean_bridge,                    4},
    {"amatrix_arrayfire_sum_axis_resident_bridge",     (DL_FUNC) &amatrix_arrayfire_sum_axis_resident_bridge,     2},
    {"am_af_lanczos_bidiag_bridge",                    (DL_FUNC) &am_af_lanczos_bidiag_bridge,                    3},
    {"am_af_lbz_upload_A_bridge",                      (DL_FUNC) &am_af_lbz_upload_A_bridge,                      1},
    {"am_af_lbz_drop_A_bridge",                        (DL_FUNC) &am_af_lbz_drop_A_bridge,                        0},
    {"am_af_lanczos_warm_bridge",                      (DL_FUNC) &am_af_lanczos_warm_bridge,                      4},
    {"am_af_dist_sq_bridge",                           (DL_FUNC) &am_af_dist_sq_bridge,                           2},
    {"am_af_kernel_bridge",                            (DL_FUNC) &am_af_kernel_bridge,                            6},
    {"amatrix_arrayfire_matmul_correct_bridge",        (DL_FUNC) &amatrix_arrayfire_matmul_correct_bridge,        2},
    {"amatrix_arrayfire_crossprod_correct_bridge",     (DL_FUNC) &amatrix_arrayfire_crossprod_correct_bridge,     2},
    {"amatrix_arrayfire_tcrossprod_correct_bridge",    (DL_FUNC) &amatrix_arrayfire_tcrossprod_correct_bridge,    2},
    {"amatrix_arrayfire_qr_q_correct_bridge",          (DL_FUNC) &amatrix_arrayfire_qr_q_correct_bridge,          1},
    {"amatrix_arrayfire_chol_bridge",                  (DL_FUNC) &amatrix_arrayfire_chol_bridge,                  1},
    {"amatrix_arrayfire_chol_resident_bridge",         (DL_FUNC) &amatrix_arrayfire_chol_resident_bridge,         2},
    {"amatrix_arrayfire_solve_bridge",                 (DL_FUNC) &amatrix_arrayfire_solve_bridge,                 2},
    {"amatrix_arrayfire_solve_resident_bridge",        (DL_FUNC) &amatrix_arrayfire_solve_resident_bridge,        3},
    {"amatrix_arrayfire_qr_Q_resident_bridge",         (DL_FUNC) &amatrix_arrayfire_qr_Q_resident_bridge,         2},
    {"amatrix_arrayfire_svd_bridge",                   (DL_FUNC) &amatrix_arrayfire_svd_bridge,                   3},
    {"amatrix_arrayfire_svd_safe_bridge",              (DL_FUNC) &amatrix_arrayfire_svd_safe_bridge,              0},
    {"amatrix_arrayfire_bdc_bidiag_bridge",            (DL_FUNC) &amatrix_arrayfire_bdc_bidiag_bridge,            1},
    {"amatrix_arrayfire_bdc_orgbr_bridge",             (DL_FUNC) &amatrix_arrayfire_bdc_orgbr_bridge,             6},
    {"amatrix_arrayfire_bdc_dbdsdc_bridge",            (DL_FUNC) &amatrix_arrayfire_bdc_dbdsdc_bridge,            3},
    {NULL, NULL, 0}
};

void R_init_amatrix_arrayfire(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
