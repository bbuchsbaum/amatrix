#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_ARRAYFIRE
#include <arrayfire.h>
#endif

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

static SEXP arrayfire_named_list2(const char* name1, SEXP value1, const char* name2, SEXP value2) {
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SEXP names = PROTECT(allocVector(STRSXP, 2));

  SET_VECTOR_ELT(out, 0, value1);
  SET_VECTOR_ELT(out, 1, value2);
  SET_STRING_ELT(names, 0, mkChar(name1));
  SET_STRING_ELT(names, 1, mkChar(name2));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

static void arrayfire_debug_stage(const char* stage) {
  const char* enabled = getenv("AMATRIX_AF_DEBUG");
  if (enabled == NULL || strcmp(enabled, "1") != 0) {
    return;
  }

  Rprintf("[arrayfire-qr] %s\n", stage);
  R_FlushConsole();
}

static void* arrayfire_xmalloc(size_t bytes) {
  void* ptr = malloc(bytes);
  if (ptr == NULL) {
    error("arrayfire bridge failed to allocate host buffer");
  }
  return ptr;
}

#ifdef HAVE_ARRAYFIRE

static SEXP arrayfire_result_to_r_matrix(const af_array arr) {
  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], arr);

  int nrow = (int)dims[0];
  int ncol = (int)dims[1];
  size_t size = (size_t)nrow * (size_t)ncol;
  float* buffer = (float*) arrayfire_xmalloc(size * sizeof(float));

  af_get_data_ptr(buffer, arr);

  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  copy_row_major_f32_to_r(REAL(out), buffer, nrow, ncol);
  free(buffer);
  UNPROTECT(1);
  return out;
}

static af_array arrayfire_matrix_from_r(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  dim_t dims[2] = {nrow, ncol};
  size_t size = (size_t)nrow * (size_t)ncol;
  float* buffer = (float*) arrayfire_xmalloc(size * sizeof(float));
  af_array out = 0;

  copy_r_to_row_major_f32(buffer, REAL(x), nrow, ncol);
  af_create_array(&out, buffer, 2, dims, f32);
  free(buffer);
  return out;
}

static af_err arrayfire_scalar_like(af_array *out, double value, const af_array like) {
  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], like);
  return af_constant(out, value, 2, dims, f32);
}

#endif

SEXP amatrix_arrayfire_native_available_bridge(void) {
#ifdef HAVE_ARRAYFIRE
  af_err err_init = af_init();
  int devices = 0;

  if (err_init != AF_SUCCESS) {
    return ScalarLogical(0);
  }

  if (af_get_device_count(&devices) != AF_SUCCESS) {
    return ScalarLogical(0);
  }

  return ScalarLogical(devices > 0);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_arrayfire_bridge_info_bridge(void) {
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, mkChar("compiled"));
  SET_STRING_ELT(names, 1, mkChar("native"));
  SET_STRING_ELT(names, 2, mkChar("engine"));

  SET_VECTOR_ELT(out, 0, ScalarLogical(1));
#ifdef HAVE_ARRAYFIRE
  SET_VECTOR_ELT(out, 1, amatrix_arrayfire_native_available_bridge());
  SET_VECTOR_ELT(out, 2, mkString("arrayfire-c"));
#else
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
  SET_VECTOR_ELT(out, 2, mkString("mock-c-bridge"));
#endif

  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP amatrix_arrayfire_diagnostics_bridge(void) {
  SEXP out = PROTECT(allocVector(VECSXP, 6));
  SEXP names = PROTECT(allocVector(STRSXP, 6));

  SET_STRING_ELT(names, 0, mkChar("compiled"));
  SET_STRING_ELT(names, 1, mkChar("init_ok"));
  SET_STRING_ELT(names, 2, mkChar("available_backends"));
  SET_STRING_ELT(names, 3, mkChar("device_count"));
  SET_STRING_ELT(names, 4, mkChar("active_backend"));
  SET_STRING_ELT(names, 5, mkChar("lapack_available"));

  SET_VECTOR_ELT(out, 0, ScalarLogical(1));

#ifdef HAVE_ARRAYFIRE
  af_err err_init = af_init();
  int backends = 0;
  int devices = 0;
  af_backend active = AF_BACKEND_DEFAULT;
  bool lapack = false;

  af_get_available_backends(&backends);
  af_get_device_count(&devices);
  af_get_active_backend(&active);
  af_is_lapack_available(&lapack);

  SET_VECTOR_ELT(out, 1, ScalarLogical(err_init == AF_SUCCESS));
  SET_VECTOR_ELT(out, 2, ScalarInteger(backends));
  SET_VECTOR_ELT(out, 3, ScalarInteger(devices));
  SET_VECTOR_ELT(out, 4, ScalarInteger((int)active));
  SET_VECTOR_ELT(out, 5, ScalarLogical(lapack));
#else
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
  SET_VECTOR_ELT(out, 2, ScalarInteger(0));
  SET_VECTOR_ELT(out, 3, ScalarInteger(0));
  SET_VECTOR_ELT(out, 4, ScalarInteger(0));
  SET_VECTOR_ELT(out, 5, ScalarLogical(0));
#endif

  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP amatrix_arrayfire_set_backend_bridge(SEXP backend) {
  if (!isInteger(backend) || XLENGTH(backend) != 1) {
    error("backend must be a scalar integer");
  }

#ifdef HAVE_ARRAYFIRE
  af_err err = af_set_backend((af_backend) INTEGER(backend)[0]);
  if (err != AF_SUCCESS) {
    error("af_set_backend failed");
  }
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_arrayfire_matmul_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x) || !isReal(y) || !isMatrix(y)) {
    error("x and y must be numeric matrices");
  }

#ifdef HAVE_ARRAYFIRE
  af_array ax = 0, ay = 0, out = 0;
  SEXP result = R_NilValue;
  af_err err;

  ax = arrayfire_matrix_from_r(x);
  ay = arrayfire_matrix_from_r(y);
  err = af_matmul(&out, ax, ay, AF_MAT_NONE, AF_MAT_NONE);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ay) af_release_array(ay);
    if (out) af_release_array(out);
    error("af_matmul failed");
  }

  result = arrayfire_result_to_r_matrix(out);

  af_release_array(ax);
  af_release_array(ay);
  af_release_array(out);
  return result;
#else
  SEXP x_dim = getAttrib(x, R_DimSymbol);
  SEXP y_dim = getAttrib(y, R_DimSymbol);
  int x_nrow = INTEGER(x_dim)[0];
  int x_ncol = INTEGER(x_dim)[1];
  int y_nrow = INTEGER(y_dim)[0];
  int y_ncol = INTEGER(y_dim)[1];

  if (x_ncol != y_nrow) {
    error("non-conformable arguments");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, x_nrow, y_ncol));
  double *x_ptr = REAL(x);
  double *y_ptr = REAL(y);
  double *out_ptr = REAL(out);

  for (int j = 0; j < y_ncol; ++j) {
    for (int i = 0; i < x_nrow; ++i) {
      double acc = 0.0;
      for (int k = 0; k < x_ncol; ++k) {
        acc += x_ptr[i + x_nrow * k] * y_ptr[k + y_nrow * j];
      }
      out_ptr[i + x_nrow * j] = acc;
    }
  }

  UNPROTECT(1);
  return out;
#endif
}

SEXP amatrix_arrayfire_crossprod_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }
  if (!isNull(y) && (!isReal(y) || !isMatrix(y))) {
    error("y must be NULL or a numeric matrix");
  }

#ifdef HAVE_ARRAYFIRE
  af_array ax = 0, ay = 0, out = 0;
  SEXP result = R_NilValue;
  af_err err;

  ax = arrayfire_matrix_from_r(x);
  if (!isNull(y)) {
    ay = arrayfire_matrix_from_r(y);
  } else {
    ay = arrayfire_matrix_from_r(x);
  }

  err = af_matmul(&out, ax, ay, AF_MAT_NONE, AF_MAT_TRANS);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ay) af_release_array(ay);
    if (out) af_release_array(out);
    error("af_matmul failed for crossprod");
  }

  result = arrayfire_result_to_r_matrix(out);

  af_release_array(ax);
  af_release_array(ay);
  af_release_array(out);
  return result;
#else
  if (isNull(y)) {
    return amatrix_arrayfire_matmul_bridge(x, x);
  }

  SEXP call = PROTECT(lang3(install("crossprod"), x, y));
  SEXP out = PROTECT(eval(call, R_BaseEnv));
  UNPROTECT(2);
  return out;
#endif
}

SEXP amatrix_arrayfire_tcrossprod_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }
  if (!isNull(y) && (!isReal(y) || !isMatrix(y))) {
    error("y must be NULL or a numeric matrix");
  }

#ifdef HAVE_ARRAYFIRE
  af_array ax = 0, ay = 0, out = 0;
  SEXP result = R_NilValue;
  af_err err;

  ax = arrayfire_matrix_from_r(x);
  if (!isNull(y)) {
    ay = arrayfire_matrix_from_r(y);
  } else {
    ay = arrayfire_matrix_from_r(x);
  }

  err = af_matmul(&out, ax, ay, AF_MAT_TRANS, AF_MAT_NONE);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ay) af_release_array(ay);
    if (out) af_release_array(out);
    error("af_matmul failed for tcrossprod");
  }

  result = arrayfire_result_to_r_matrix(out);

  af_release_array(ax);
  af_release_array(ay);
  af_release_array(out);
  return result;
#else
  if (isNull(y)) {
    SEXP call = PROTECT(lang2(install("tcrossprod"), x));
    SEXP out = PROTECT(eval(call, R_BaseEnv));
    UNPROTECT(2);
    return out;
  }

  SEXP call = PROTECT(lang3(install("tcrossprod"), x, y));
  SEXP out = PROTECT(eval(call, R_BaseEnv));
  UNPROTECT(2);
  return out;
#endif
}

SEXP amatrix_arrayfire_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op) {
  if (!isReal(lhs) || !isMatrix(lhs) || !isString(op) || LENGTH(op) != 1) {
    error("lhs must be a numeric matrix and op must be a string");
  }
  if (!isNull(rhs) && !isReal(rhs)) {
    error("rhs must be NULL, a numeric scalar, or a numeric matrix");
  }

#ifdef HAVE_ARRAYFIRE
  af_array alhs = 0, arhs = 0, out = 0;
  af_err err = AF_SUCCESS;
  const char *op_name = CHAR(STRING_ELT(op, 0));
  SEXP result = R_NilValue;

  alhs = arrayfire_matrix_from_r(lhs);

  if (!isNull(rhs)) {
    if (isMatrix(rhs)) {
      arhs = arrayfire_matrix_from_r(rhs);
    } else if (XLENGTH(rhs) == 1) {
      err = arrayfire_scalar_like(&arhs, REAL(rhs)[0], alhs);
      if (err != AF_SUCCESS) {
        if (alhs) af_release_array(alhs);
        error("af_constant failed for scalar rhs");
      }
    } else {
      if (alhs) af_release_array(alhs);
      error("rhs must be a scalar or matrix");
    }
  } else {
    if (alhs) af_release_array(alhs);
    error("rhs must not be NULL for arrayfire ewise");
  }

  if (strcmp(op_name, "+") == 0) {
    err = af_add(&out, alhs, arhs, false);
  } else if (strcmp(op_name, "-") == 0) {
    err = af_sub(&out, alhs, arhs, false);
  } else if (strcmp(op_name, "*") == 0) {
    err = af_mul(&out, alhs, arhs, false);
  } else {
    if (alhs) af_release_array(alhs);
    if (arhs) af_release_array(arhs);
    error("unsupported arrayfire ewise operator");
  }

  if (err != AF_SUCCESS) {
    if (alhs) af_release_array(alhs);
    if (arhs) af_release_array(arhs);
    if (out) af_release_array(out);
    error("arrayfire ewise op failed");
  }

  result = arrayfire_result_to_r_matrix(out);
  af_release_array(alhs);
  af_release_array(arhs);
  af_release_array(out);
  return result;
#else
  SEXP call = PROTECT(lang3(install(CHAR(STRING_ELT(op, 0))), lhs, rhs));
  SEXP out = PROTECT(eval(call, R_BaseEnv));
  UNPROTECT(2);
  return out;
#endif
}

SEXP amatrix_arrayfire_sum_axis_bridge(SEXP x, SEXP axis) {
  if (!isReal(x) || !isMatrix(x) || !isInteger(axis) || LENGTH(axis) != 1) {
    error("x must be a numeric matrix and axis must be a scalar integer");
  }

#ifdef HAVE_ARRAYFIRE
  af_array ax = 0, out = 0;
  af_err err;
  int axis_value = INTEGER(axis)[0];
  SEXP result;

  ax = arrayfire_matrix_from_r(x);
  err = af_sum(&out, ax, axis_value);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (out) af_release_array(out);
    error("arrayfire sum failed");
  }

  result = arrayfire_result_to_r_matrix(out);
  af_release_array(ax);
  af_release_array(out);

  if (axis_value == 0) {
    SEXP vec = PROTECT(allocVector(REALSXP, XLENGTH(result)));
    for (R_xlen_t i = 0; i < XLENGTH(result); ++i) {
      REAL(vec)[i] = REAL(result)[i];
    }
    UNPROTECT(1);
    return vec;
  }

  if (axis_value == 1) {
    SEXP vec = PROTECT(allocVector(REALSXP, XLENGTH(result)));
    for (R_xlen_t i = 0; i < XLENGTH(result); ++i) {
      REAL(vec)[i] = REAL(result)[i];
    }
    UNPROTECT(1);
    return vec;
  }

  return result;
#else
  SEXP call = PROTECT(lang2(install(INTEGER(axis)[0] == 0 ? "colSums" : "rowSums"), x));
  SEXP out = PROTECT(eval(call, R_BaseEnv));
  UNPROTECT(2);
  return out;
#endif
}

SEXP amatrix_arrayfire_qr_bridge(SEXP x) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

#ifdef HAVE_ARRAYFIRE
  af_array ax = 0, ax_t = 0, q = 0, r = 0, q_t = 0, r_t = 0, tau = 0;
  SEXP q_r = R_NilValue;
  SEXP r_r = R_NilValue;
  SEXP result = R_NilValue;
  af_err err;

  arrayfire_debug_stage("create_input");
  ax = arrayfire_matrix_from_r(x);
  arrayfire_debug_stage("transpose_input");
  err = af_transpose(&ax_t, ax, false);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ax_t) af_release_array(ax_t);
    error("af_transpose failed before af_qr");
  }

  arrayfire_debug_stage("af_qr");
  err = af_qr(&q, &r, &tau, ax_t);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ax_t) af_release_array(ax_t);
    if (q) af_release_array(q);
    if (r) af_release_array(r);
    if (tau) af_release_array(tau);
    error("af_qr failed");
  }

  arrayfire_debug_stage("transpose_q");
  err = af_transpose(&q_t, q, false);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ax_t) af_release_array(ax_t);
    if (q) af_release_array(q);
    if (r) af_release_array(r);
    if (q_t) af_release_array(q_t);
    if (r_t) af_release_array(r_t);
    if (tau) af_release_array(tau);
    error("af_transpose failed for Q");
  }

  arrayfire_debug_stage("transpose_r");
  err = af_transpose(&r_t, r, false);
  if (err != AF_SUCCESS) {
    if (ax) af_release_array(ax);
    if (ax_t) af_release_array(ax_t);
    if (q) af_release_array(q);
    if (r) af_release_array(r);
    if (q_t) af_release_array(q_t);
    if (r_t) af_release_array(r_t);
    if (tau) af_release_array(tau);
    error("af_transpose failed for R");
  }

  arrayfire_debug_stage("materialize_q");
  PROTECT(q_r = arrayfire_result_to_r_matrix(q_t));
  arrayfire_debug_stage("materialize_r");
  PROTECT(r_r = arrayfire_result_to_r_matrix(r_t));
  arrayfire_debug_stage("build_result");
  result = arrayfire_named_list2("q", q_r, "r", r_r);
  UNPROTECT(2);

  arrayfire_debug_stage("release");
  af_release_array(ax);
  af_release_array(ax_t);
  af_release_array(q);
  af_release_array(r);
  af_release_array(q_t);
  af_release_array(r_t);
  af_release_array(tau);
  arrayfire_debug_stage("done");
  return result;
#else
  error("arrayfire qr bridge requires arrayfire");
#endif
}
