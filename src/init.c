#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP amatrix_block_reorth_bridge(SEXP z, SEXP basis, SEXP return_projection);
extern SEXP amatrix_block_reorth_prefix_bridge(SEXP z, SEXP basis, SEXP basis_cols, SEXP return_projection);
extern SEXP amatrix_block_thin_qr_bridge(SEXP z);

static const R_CallMethodDef call_methods[] = {
  {"amatrix_block_reorth_bridge", (DL_FUNC) &amatrix_block_reorth_bridge, 3},
  {"amatrix_block_reorth_prefix_bridge", (DL_FUNC) &amatrix_block_reorth_prefix_bridge, 4},
  {"amatrix_block_thin_qr_bridge", (DL_FUNC) &amatrix_block_thin_qr_bridge, 1},
  {NULL, NULL, 0}
};

void R_init_amatrix(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
