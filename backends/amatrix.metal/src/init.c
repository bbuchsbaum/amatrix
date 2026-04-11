#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP amatrix_metal_native_available_bridge(void);
extern SEXP amatrix_metal_bridge_info_bridge(void);
extern SEXP amatrix_metal_profile_set_enabled_bridge(SEXP enabled);
extern SEXP amatrix_metal_profile_reset_bridge(void);
extern SEXP amatrix_metal_profile_bridge(void);
extern SEXP amatrix_metal_spmm_bridge(SEXP values, SEXP p, SEXP i, SEXP dim, SEXP b, SEXP trans_lhs);
extern SEXP amatrix_metal_sparse_store_bridge(SEXP key, SEXP values, SEXP p, SEXP i, SEXP dim);
extern SEXP amatrix_metal_sparse_has_bridge(SEXP key);
extern SEXP amatrix_metal_sparse_drop_bridge(SEXP key);
extern SEXP amatrix_metal_dense_store_bridge(SEXP key, SEXP x);
extern SEXP amatrix_metal_dense_has_bridge(SEXP key);
extern SEXP amatrix_metal_dense_drop_bridge(SEXP key);
extern SEXP amatrix_metal_dense_materialize_bridge(SEXP key);
extern SEXP amatrix_metal_spmm_resident_bridge(SEXP key, SEXP b, SEXP trans_lhs);
extern SEXP amatrix_metal_spmm_resident_key_bridge(SEXP sp_key, SEXP y_key, SEXP out_key, SEXP trans_lhs, SEXP defer);
extern SEXP amatrix_metal_dense_sparse_matmul_resident_key_bridge(SEXP x_key, SEXP sp_key, SEXP out_key, SEXP defer);
extern SEXP amatrix_metal_transpose_resident_bridge(SEXP x_key, SEXP out_key);

static const R_CallMethodDef call_methods[] = {
    {"amatrix_metal_native_available_bridge", (DL_FUNC) &amatrix_metal_native_available_bridge, 0},
    {"amatrix_metal_bridge_info_bridge", (DL_FUNC) &amatrix_metal_bridge_info_bridge, 0},
    {"amatrix_metal_profile_set_enabled_bridge", (DL_FUNC) &amatrix_metal_profile_set_enabled_bridge, 1},
    {"amatrix_metal_profile_reset_bridge", (DL_FUNC) &amatrix_metal_profile_reset_bridge, 0},
    {"amatrix_metal_profile_bridge", (DL_FUNC) &amatrix_metal_profile_bridge, 0},
    {"amatrix_metal_spmm_bridge", (DL_FUNC) &amatrix_metal_spmm_bridge, 6},
    {"amatrix_metal_sparse_store_bridge", (DL_FUNC) &amatrix_metal_sparse_store_bridge, 5},
    {"amatrix_metal_sparse_has_bridge", (DL_FUNC) &amatrix_metal_sparse_has_bridge, 1},
    {"amatrix_metal_sparse_drop_bridge", (DL_FUNC) &amatrix_metal_sparse_drop_bridge, 1},
    {"amatrix_metal_dense_store_bridge", (DL_FUNC) &amatrix_metal_dense_store_bridge, 2},
    {"amatrix_metal_dense_has_bridge", (DL_FUNC) &amatrix_metal_dense_has_bridge, 1},
    {"amatrix_metal_dense_drop_bridge", (DL_FUNC) &amatrix_metal_dense_drop_bridge, 1},
    {"amatrix_metal_dense_materialize_bridge", (DL_FUNC) &amatrix_metal_dense_materialize_bridge, 1},
    {"amatrix_metal_spmm_resident_bridge", (DL_FUNC) &amatrix_metal_spmm_resident_bridge, 3},
    {"amatrix_metal_spmm_resident_key_bridge", (DL_FUNC) &amatrix_metal_spmm_resident_key_bridge, 5},
    {"amatrix_metal_dense_sparse_matmul_resident_key_bridge", (DL_FUNC) &amatrix_metal_dense_sparse_matmul_resident_key_bridge, 4},
    {"amatrix_metal_transpose_resident_bridge", (DL_FUNC) &amatrix_metal_transpose_resident_bridge, 2},
    {NULL, NULL, 0}
};

void R_init_amatrix_metal(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
