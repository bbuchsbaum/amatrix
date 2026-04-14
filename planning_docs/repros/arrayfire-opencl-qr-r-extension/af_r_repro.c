#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <arrayfire.h>
#include <stdlib.h>

static void fail_af(const char* stage, af_err err) {
  error("%s failed with error code %d", stage, (int)err);
}

static void copy_r_to_row_major_f32(float* out, const double* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[i * ncol + j] = (float)in[i + nrow * j];
    }
  }
}

static void copy_row_major_f32_to_r(double* out, const float* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[i + nrow * j] = in[i * ncol + j];
    }
  }
}

static af_array matrix_from_r(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  dim_t dims[2] = {nrow, ncol};
  size_t size = (size_t)nrow * (size_t)ncol;
  float* buffer = (float*)malloc(size * sizeof(float));
  af_array out = 0;
  af_err err;

  if (buffer == NULL) error("malloc failed");
  copy_r_to_row_major_f32(buffer, REAL(x), nrow, ncol);
  err = af_create_array(&out, buffer, 2, dims, f32);
  free(buffer);
  if (err != AF_SUCCESS) fail_af("af_create_array", err);
  return out;
}

static SEXP result_to_r_matrix(const af_array arr) {
  dim_t d0 = 0, d1 = 0, d2 = 0, d3 = 0;
  af_err err = af_get_dims(&d0, &d1, &d2, &d3, arr);
  int nrow = (int)d0;
  int ncol = (int)d1;
  size_t size = (size_t)nrow * (size_t)ncol;
  float* buffer;
  SEXP out;

  if (err != AF_SUCCESS) fail_af("af_get_dims", err);
  buffer = (float*)malloc(size * sizeof(float));
  if (buffer == NULL) error("malloc failed");

  err = af_get_data_ptr(buffer, arr);
  if (err != AF_SUCCESS) {
    free(buffer);
    fail_af("af_get_data_ptr", err);
  }

  out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  copy_row_major_f32_to_r(REAL(out), buffer, nrow, ncol);
  free(buffer);
  UNPROTECT(1);
  return out;
}

SEXP af_r_set_backend(SEXP backend_) {
  af_err err;
  if (!isInteger(backend_) || XLENGTH(backend_) != 1) error("backend must be scalar integer");
  err = af_set_backend((af_backend)INTEGER(backend_)[0]);
  if (err != AF_SUCCESS) fail_af("af_set_backend", err);
  return ScalarLogical(1);
}

SEXP af_r_diag(void) {
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  af_backend active = AF_BACKEND_DEFAULT;
  bool lapack = false;
  int devices = 0;
  af_err err;

  err = af_init();
  if (err != AF_SUCCESS) fail_af("af_init", err);
  err = af_get_active_backend(&active);
  if (err != AF_SUCCESS) fail_af("af_get_active_backend", err);
  err = af_is_lapack_available(&lapack);
  if (err != AF_SUCCESS) fail_af("af_is_lapack_available", err);
  err = af_get_device_count(&devices);
  if (err != AF_SUCCESS) fail_af("af_get_device_count", err);

  SET_STRING_ELT(names, 0, mkChar("active_backend"));
  SET_STRING_ELT(names, 1, mkChar("lapack_available"));
  SET_STRING_ELT(names, 2, mkChar("device_count"));
  SET_VECTOR_ELT(out, 0, ScalarInteger((int)active));
  SET_VECTOR_ELT(out, 1, ScalarLogical(lapack));
  SET_VECTOR_ELT(out, 2, ScalarInteger(devices));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP af_r_qr(SEXP x) {
  af_array ax = 0, ax_t = 0, q = 0, r = 0, q_t = 0, r_t = 0, tau = 0;
  af_err err;
  SEXP q_r = R_NilValue, r_r = R_NilValue, out = R_NilValue, names = R_NilValue;

  if (!isReal(x) || !isMatrix(x)) error("x must be numeric matrix");

  Rprintf("[af_r_qr] create_input\n"); R_FlushConsole();
  ax = matrix_from_r(x);

  Rprintf("[af_r_qr] transpose_input\n"); R_FlushConsole();
  err = af_transpose(&ax_t, ax, false);
  if (err != AF_SUCCESS) fail_af("af_transpose(input)", err);

  Rprintf("[af_r_qr] af_qr\n"); R_FlushConsole();
  err = af_qr(&q, &r, &tau, ax_t);
  if (err != AF_SUCCESS) fail_af("af_qr", err);

  Rprintf("[af_r_qr] transpose_q\n"); R_FlushConsole();
  err = af_transpose(&q_t, q, false);
  if (err != AF_SUCCESS) fail_af("af_transpose(q)", err);

  Rprintf("[af_r_qr] transpose_r\n"); R_FlushConsole();
  err = af_transpose(&r_t, r, false);
  if (err != AF_SUCCESS) fail_af("af_transpose(r)", err);

  Rprintf("[af_r_qr] materialize_q\n"); R_FlushConsole();
  PROTECT(q_r = result_to_r_matrix(q_t));
  Rprintf("[af_r_qr] materialize_r\n"); R_FlushConsole();
  PROTECT(r_r = result_to_r_matrix(r_t));

  PROTECT(out = allocVector(VECSXP, 2));
  PROTECT(names = allocVector(STRSXP, 2));
  SET_VECTOR_ELT(out, 0, q_r);
  SET_VECTOR_ELT(out, 1, r_r);
  SET_STRING_ELT(names, 0, mkChar("q"));
  SET_STRING_ELT(names, 1, mkChar("r"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(4);

  if (tau) af_release_array(tau);
  if (r_t) af_release_array(r_t);
  if (q_t) af_release_array(q_t);
  if (r) af_release_array(r);
  if (q) af_release_array(q);
  if (ax_t) af_release_array(ax_t);
  if (ax) af_release_array(ax);
  return out;
}

static const R_CallMethodDef call_methods[] = {
  {"af_r_set_backend", (DL_FUNC)&af_r_set_backend, 1},
  {"af_r_diag", (DL_FUNC)&af_r_diag, 0},
  {"af_r_qr", (DL_FUNC)&af_r_qr, 1},
  {NULL, NULL, 0}
};

void R_init_af_r_repro(DllInfo* dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
