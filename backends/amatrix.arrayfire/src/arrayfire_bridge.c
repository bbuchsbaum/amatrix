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

/* ── Resident registry ─────────────────────────────────────────── */

typedef struct {
  char*    key;
  af_array array;
  int      in_use;
} amatrix_af_resident_entry;

static amatrix_af_resident_entry* amatrix_af_registry          = NULL;
static size_t                     amatrix_af_registry_capacity = 0;

static void amatrix_af_registry_init(void) {
  if (amatrix_af_registry != NULL) return;
  amatrix_af_registry_capacity = 128;
  amatrix_af_registry = (amatrix_af_resident_entry*)
    calloc(amatrix_af_registry_capacity, sizeof(amatrix_af_resident_entry));
  if (amatrix_af_registry == NULL)
    error("failed to allocate arrayfire residency registry");
}

static amatrix_af_resident_entry* amatrix_af_registry_find(const char* key) {
  amatrix_af_registry_init();
  for (size_t i = 0; i < amatrix_af_registry_capacity; ++i) {
    if (amatrix_af_registry[i].in_use &&
        strcmp(amatrix_af_registry[i].key, key) == 0)
      return &amatrix_af_registry[i];
  }
  return NULL;
}

static amatrix_af_resident_entry* amatrix_af_registry_reserve(const char* key) {
  amatrix_af_resident_entry* existing = amatrix_af_registry_find(key);
  if (existing != NULL) {
    if (existing->array) af_release_array(existing->array);
    existing->array = 0;
    return existing;
  }
  amatrix_af_registry_init();
  for (size_t i = 0; i < amatrix_af_registry_capacity; ++i) {
    if (!amatrix_af_registry[i].in_use) {
      amatrix_af_registry[i].in_use = 1;
      amatrix_af_registry[i].key = strdup(key);
      amatrix_af_registry[i].array = 0;
      if (amatrix_af_registry[i].key == NULL)
        error("failed to allocate arrayfire resident key");
      return &amatrix_af_registry[i];
    }
  }
  error("arrayfire residency registry is full");
  return NULL;
}

static void amatrix_af_registry_drop(const char* key) {
  amatrix_af_resident_entry* e = amatrix_af_registry_find(key);
  if (e == NULL) return;
  if (e->key)   { free(e->key); e->key = NULL; }
  if (e->array) { af_release_array(e->array); e->array = 0; }
  e->in_use = 0;
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

SEXP amatrix_arrayfire_resident_store_bridge(SEXP key, SEXP x) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(key) || LENGTH(key) != 1)
    error("key must be a scalar string");
  if (!isReal(x) || !isMatrix(x))
    error("x must be a numeric matrix");
  const char* k = CHAR(STRING_ELT(key, 0));
  amatrix_af_resident_entry* e = amatrix_af_registry_reserve(k);
  e->array = arrayfire_matrix_from_r(x);
  return ScalarLogical(1);
#else
  error("arrayfire resident store requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_resident_has_bridge(SEXP key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(key) || LENGTH(key) != 1)
    error("key must be a scalar string");
  return ScalarLogical(amatrix_af_registry_find(CHAR(STRING_ELT(key, 0))) != NULL);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_arrayfire_resident_drop_bridge(SEXP key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(key) || LENGTH(key) != 1)
    error("key must be a scalar string");
  amatrix_af_registry_drop(CHAR(STRING_ELT(key, 0)));
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_arrayfire_resident_materialize_bridge(SEXP key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(key) || LENGTH(key) != 1)
    error("key must be a scalar string");
  amatrix_af_resident_entry* e = amatrix_af_registry_find(CHAR(STRING_ELT(key, 0)));
  if (e == NULL)
    error("arrayfire resident key not found: %s", CHAR(STRING_ELT(key, 0)));
  return arrayfire_result_to_r_matrix(e->array);
#else
  error("arrayfire resident materialize requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || !isString(y_key) || !isString(out_key))
    error("keys must be scalar strings");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  amatrix_af_resident_entry* ey = amatrix_af_registry_find(CHAR(STRING_ELT(y_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  if (ey == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(y_key, 0)));
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  af_err err = af_matmul(&eout->array, ex->array, ey->array, AF_MAT_NONE, AF_MAT_NONE);
  if (err != AF_SUCCESS) error("af_matmul (resident) failed");
  return ScalarLogical(1);
#else
  error("arrayfire matmul_resident requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || !isString(out_key))
    error("keys must be scalar strings");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  /* y_key NULL => crossprod(x, x) */
  af_array ay = 0;
  if (isNull(y_key)) {
    ay = ex->array;
  } else {
    amatrix_af_resident_entry* ey = amatrix_af_registry_find(CHAR(STRING_ELT(y_key, 0)));
    if (ey == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(y_key, 0)));
    ay = ey->array;
  }
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  /* crossprod(x,y) = t(x) %*% y.  arrayfire_matrix_from_r stores the R-transpose,
     so the AF array ax = t(x_orig).  To compute t(x_orig) %*% y_orig we need
     ax %*% t(ay) = AF_MAT_NONE on lhs, AF_MAT_TRANS on rhs. */
  af_err err = af_matmul(&eout->array, ex->array, ay, AF_MAT_NONE, AF_MAT_TRANS);
  if (err != AF_SUCCESS) error("af_matmul (crossprod resident) failed");
  return ScalarLogical(1);
#else
  error("arrayfire crossprod_resident requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || !isString(out_key))
    error("keys must be scalar strings");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  af_array ay = 0;
  if (isNull(y_key)) {
    ay = ex->array;
  } else {
    amatrix_af_resident_entry* ey = amatrix_af_registry_find(CHAR(STRING_ELT(y_key, 0)));
    if (ey == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(y_key, 0)));
    ay = ey->array;
  }
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  /* tcrossprod(x,y) = x %*% t(y).  arrayfire_matrix_from_r stores the R-transpose,
     so ax = t(x_orig), ay = t(y_orig).  To compute x_orig %*% t(y_orig) we need
     t(ax) %*% ay = AF_MAT_TRANS on lhs, AF_MAT_NONE on rhs. */
  af_err err = af_matmul(&eout->array, ex->array, ay, AF_MAT_TRANS, AF_MAT_NONE);
  if (err != AF_SUCCESS) error("af_matmul (tcrossprod resident) failed");
  return ScalarLogical(1);
#else
  error("arrayfire tcrossprod_resident requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(lhs_key) || !isString(op) || !isString(out_key))
    error("keys and op must be scalar strings");
  amatrix_af_resident_entry* elhs = amatrix_af_registry_find(CHAR(STRING_ELT(lhs_key, 0)));
  if (elhs == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(lhs_key, 0)));
  const char* op_name = CHAR(STRING_ELT(op, 0));

  af_array arhs = 0;
  int own_arhs = 0;
  if (isString(rhs) && LENGTH(rhs) == 1) {
    /* rhs is a resident key */
    amatrix_af_resident_entry* erhs = amatrix_af_registry_find(CHAR(STRING_ELT(rhs, 0)));
    if (erhs == NULL) error("arrayfire resident rhs key not found");
    arhs = erhs->array;
  } else if (isReal(rhs) && XLENGTH(rhs) == 1) {
    af_err err = arrayfire_scalar_like(&arhs, REAL(rhs)[0], elhs->array);
    if (err != AF_SUCCESS) error("af_constant failed for scalar rhs");
    own_arhs = 1;
  } else {
    error("rhs must be a resident key (string) or a numeric scalar");
  }

  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  af_err err = AF_SUCCESS;
  if      (strcmp(op_name, "+") == 0) err = af_add(&eout->array, elhs->array, arhs, false);
  else if (strcmp(op_name, "-") == 0) err = af_sub(&eout->array, elhs->array, arhs, false);
  else if (strcmp(op_name, "*") == 0) err = af_mul(&eout->array, elhs->array, arhs, false);
  else { if (own_arhs) af_release_array(arhs); error("unsupported ewise op: %s", op_name); }

  if (own_arhs) af_release_array(arhs);
  if (err != AF_SUCCESS) error("arrayfire ewise (resident) failed");
  return ScalarLogical(1);
#else
  error("arrayfire ewise_resident requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_sum_axis_resident_bridge(SEXP x_key, SEXP axis) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || LENGTH(x_key) != 1)
    error("x_key must be a scalar string");
  if (!isInteger(axis) || LENGTH(axis) != 1)
    error("axis must be a scalar integer");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  af_array out = 0;
  af_err err = af_sum(&out, ex->array, INTEGER(axis)[0]);
  if (err != AF_SUCCESS) error("af_sum (resident) failed");
  SEXP result = PROTECT(arrayfire_result_to_r_matrix(out));
  af_release_array(out);
  SEXP vec = PROTECT(allocVector(REALSXP, XLENGTH(result)));
  for (R_xlen_t i = 0; i < XLENGTH(result); ++i) REAL(vec)[i] = REAL(result)[i];
  UNPROTECT(2);
  return vec;
#else
  error("arrayfire sum_axis_resident requires arrayfire");
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
   Native GPU Lanczos Bidiagonalization
   ───────────────────────────────────────────────────────────────────────────
   Uses a *separate* convention from the row-major resident bridge: data is
   passed to af_create_array directly in column-major order (same as R), so
   AF sees the correct matrices without any transpose.  Arrays created here
   must NOT be mixed with the resident registry (different storage convention).
   ═══════════════════════════════════════════════════════════════════════════ */

#ifdef HAVE_ARRAYFIRE

/* ── LBZ globals (visible to all helpers below) ────────────────── */
static af_array g_lbz_A    = 0;
static int      g_lbz_A_m  = 0;
static int      g_lbz_A_n  = 0;
static af_dtype g_lbz_dtype = f32;  /* set to f64 if device supports it */

/* Upload n-length R double vector as an n×1 AF column vector.
 * Uses g_lbz_dtype (f64 if device supports it, else f32). */
static af_array lbz_upload_vec(const double *src, int n) {
  dim_t dims[2] = {(dim_t)n, 1};
  af_array out = 0;
  if (g_lbz_dtype == f64) {
    af_create_array(&out, src, 2, dims, f64);
  } else {
    size_t sz = (size_t)n;
    float *buf = (float *)arrayfire_xmalloc(sz * sizeof(float));
    for (size_t i = 0; i < sz; i++) buf[i] = (float)src[i];
    af_create_array(&out, buf, 2, dims, f32);
    free(buf);
  }
  return out;
}

/* Upload m×n R double matrix (column-major) as AF matrix. */
static af_array lbz_upload_mat(const double *src, int nrow, int ncol) {
  dim_t dims[2] = {(dim_t)nrow, (dim_t)ncol};
  af_array out = 0;
  if (g_lbz_dtype == f64) {
    af_create_array(&out, src, 2, dims, f64);
  } else {
    size_t sz = (size_t)nrow * (size_t)ncol;
    float *buf = (float *)arrayfire_xmalloc(sz * sizeof(float));
    for (size_t i = 0; i < sz; i++) buf[i] = (float)src[i];
    af_create_array(&out, buf, 2, dims, f32);
    free(buf);
  }
  return out;
}

/* Download AF matrix → new R double matrix (column-major).
 * Handles both f32 and f64 device arrays. */
static SEXP lbz_download_mat(const af_array arr) {
  dim_t d0, d1, d2, d3;
  af_get_dims(&d0, &d1, &d2, &d3, arr);
  int nrow = (int)d0, ncol = (int)d1;
  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  if (g_lbz_dtype == f64) {
    af_get_data_ptr(REAL(out), arr);
  } else {
    size_t sz = (size_t)nrow * (size_t)ncol;
    float *buf = (float *)arrayfire_xmalloc(sz * sizeof(float));
    af_get_data_ptr(buf, arr);
    double *dst = REAL(out);
    for (size_t i = 0; i < sz; i++) dst[i] = (double)buf[i];
    free(buf);
  }
  UNPROTECT(1);
  return out;
}

/* out = arr * scalar */
static af_err lbz_scale(af_array *out, const af_array arr, double s) {
  dim_t d0, d1, d2, d3;
  af_get_dims(&d0, &d1, &d2, &d3, arr);
  dim_t dims[2] = {d0, d1};
  af_array sc = 0;
  af_err err = af_constant(&sc, s, 2, dims, g_lbz_dtype);
  if (err != AF_SUCCESS) return err;
  err = af_mul(out, arr, sc, false);
  af_release_array(sc);
  return err;
}

/* Normalize to unit L2 norm; store norm in *nrm. */
static af_err lbz_normalize(af_array *out, double *nrm, const af_array arr) {
  af_err err = af_norm(nrm, arr, AF_NORM_VECTOR_2, 0, 0);
  if (err != AF_SUCCESS) return err;
  if (*nrm < 1e-14) return af_retain_array(out, arr);
  return lbz_scale(out, arr, 1.0 / *nrm);
}

/*
 * CGS2 reorthogonalization: out = vec - basis*(basis^T*vec), two passes.
 *
 * With column-major convention: basis is m×j, vec is m×1.
 *   pass1: proj = basis^T * vec  (AF_MAT_TRANS on lhs = t(basis))  → j×1
 *          corr = basis * proj                                       → m×1
 *          tmp  = vec - corr
 *   pass2: same with tmp
 */
static af_err lbz_cgs2(af_array *out, const af_array vec, const af_array basis) {
  af_array p1 = 0, c1 = 0, t1 = 0, p2 = 0, c2 = 0;
  af_err err;

  err = af_matmul(&p1, basis, vec, AF_MAT_TRANS, AF_MAT_NONE); if (err) goto done;
  err = af_matmul(&c1, basis, p1,  AF_MAT_NONE,  AF_MAT_NONE); if (err) goto done;
  err = af_sub(&t1, vec, c1, false);                            if (err) goto done;
  err = af_matmul(&p2, basis, t1,  AF_MAT_TRANS, AF_MAT_NONE); if (err) goto done;
  err = af_matmul(&c2, basis, p2,  AF_MAT_NONE,  AF_MAT_NONE); if (err) goto done;
  err = af_sub(out, t1, c2, false);

done:
  if (p1) af_release_array(p1);
  if (c1) af_release_array(c1);
  if (t1) af_release_array(t1);
  if (p2) af_release_array(p2);
  if (c2) af_release_array(c2);
  return err;
}

/* Append column vec to matrix *pbasis along dim-1 (grow right). */
static af_err lbz_col_append(af_array *pbasis, const af_array vec) {
  if (*pbasis == 0) return af_retain_array(pbasis, vec);
  af_array nb = 0;
  af_err err = af_join(&nb, 1, *pbasis, vec);
  if (err != AF_SUCCESS) return err;
  af_release_array(*pbasis);
  *pbasis = nb;
  return AF_SUCCESS;
}

#endif /* HAVE_ARRAYFIRE */

/*
 * am_af_lanczos_bidiag_bridge(A, v0, k)
 *
 * Run k steps of Golub-Kahan Lanczos bidiagonalization fully on GPU.
 * All matvecs and CGS2 reorthogonalization stay on device; only 2k
 * scalars (alpha, beta) and the two basis matrices cross PCIe at the end.
 *
 * A   : m×n double R matrix (column-major)
 * v0  : n-length double R vector (starting vector; normalised internally)
 * k   : scalar integer, Lanczos steps (1 ≤ k ≤ n)
 *
 * Returns list(U = m×k, V = n×(k+1), alpha = k, beta = k) where
 *   B = diag(alpha) + superdiag(beta[1..k-1])  is the k×k upper bidiagonal,
 *   beta[k] is the final residual norm for convergence testing.
 */
SEXP am_af_lanczos_bidiag_bridge(SEXP A_r, SEXP v0_r, SEXP k_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isMatrix(A_r) || !isReal(A_r))
    error("A must be a real matrix");
  if (!isReal(v0_r))
    error("v0 must be a real numeric vector");
  if (!isInteger(k_r) || LENGTH(k_r) != 1)
    error("k must be a scalar integer");

  SEXP dim_r = getAttrib(A_r, R_DimSymbol);
  int m = INTEGER(dim_r)[0];
  int n = INTEGER(dim_r)[1];
  int k = INTEGER(k_r)[0];

  if (LENGTH(v0_r) != n)
    error("v0 must have length ncol(A) = %d", n);
  if (k < 1 || k > n)
    error("k must satisfy 1 <= k <= ncol(A)");

  /* Scalars (R-managed, freed automatically) */
  double *alpha_h = (double *)R_alloc(k, sizeof(double));
  double *beta_h  = (double *)R_alloc(k, sizeof(double));

  /* ── Upload A ────────────────────────────────────────────────── */
  af_array d_A = lbz_upload_mat(REAL(A_r), m, n);

  /* ── Normalise v0 → d_v_cur = v_0 ───────────────────────────── */
  af_array d_tmp_v0 = lbz_upload_vec(REAL(v0_r), n);
  double init_nrm = 0.0;
  af_array d_v_cur = 0;
  if (lbz_normalize(&d_v_cur, &init_nrm, d_tmp_v0) != AF_SUCCESS
      || init_nrm < 1e-14) {
    af_release_array(d_tmp_v0);
    af_release_array(d_A);
    if (d_v_cur) af_release_array(d_v_cur);
    error("v0 is essentially zero or normalization failed");
  }
  af_release_array(d_tmp_v0);

  /* Growing basis matrices and recurrence state */
  af_array d_Ubasis = 0;  /* m × j   (grows to m × k)     */
  af_array d_Vbasis = 0;  /* n × j+1 (grows to n × (k+1)) */
  af_array d_u_prev = 0;  /* u_{j-1} for the recurrence   */

  /* ── Lanczos bidiagonalization loop ─────────────────────────── */
  for (int j = 0; j < k; j++) {

    /* Append v_j to V basis */
    lbz_col_append(&d_Vbasis, d_v_cur);

    /* Forward matvec: u_raw = A * v_j */
    af_array d_u_raw = 0;
    af_matmul(&d_u_raw, d_A, d_v_cur, AF_MAT_NONE, AF_MAT_NONE);

    /* Subtract prior contribution: u_raw -= beta_{j-1} * u_{j-1} */
    if (j > 0 && d_u_prev != 0) {
      af_array d_sc = 0, d_sub = 0;
      lbz_scale(&d_sc, d_u_prev, (float)beta_h[j - 1]);
      af_sub(&d_sub, d_u_raw, d_sc, false);
      af_release_array(d_u_raw);
      af_release_array(d_sc);
      d_u_raw = d_sub;
    }

    /* Reorthogonalise u against U basis (GEMM CGS2) */
    if (d_Ubasis != 0) {
      af_array d_u_orth = 0;
      lbz_cgs2(&d_u_orth, d_u_raw, d_Ubasis);
      af_release_array(d_u_raw);
      d_u_raw = d_u_orth;
    }

    /* alpha_j = ||u_raw||;  d_uj = u_raw / alpha_j */
    af_array d_uj = 0;
    lbz_normalize(&d_uj, &alpha_h[j], d_u_raw);
    af_release_array(d_u_raw);

    /* Append u_j to U basis; advance d_u_prev */
    lbz_col_append(&d_Ubasis, d_uj);
    if (d_u_prev) af_release_array(d_u_prev);
    d_u_prev = d_uj;  /* transfer ownership */

    /* Backward matvec: v_raw = A^T * u_j - alpha_j * v_j */
    af_array d_AtU = 0, d_alphav = 0, d_v_raw = 0;
    af_matmul(&d_AtU, d_A, d_uj, AF_MAT_TRANS, AF_MAT_NONE);
    lbz_scale(&d_alphav, d_v_cur, (float)alpha_h[j]);
    af_sub(&d_v_raw, d_AtU, d_alphav, false);
    af_release_array(d_AtU);
    af_release_array(d_alphav);

    /* Reorthogonalise v_raw against V basis (GEMM CGS2) */
    {
      af_array d_v_orth = 0;
      lbz_cgs2(&d_v_orth, d_v_raw, d_Vbasis);
      af_release_array(d_v_raw);
      d_v_raw = d_v_orth;
    }

    /* beta_j = ||v_raw||;  d_v_next = v_raw / beta_j */
    af_array d_v_next = 0;
    lbz_normalize(&d_v_next, &beta_h[j], d_v_raw);
    af_release_array(d_v_raw);

    /* Advance current v */
    af_release_array(d_v_cur);
    d_v_cur = d_v_next;
  }

  /* Append the trailing v_k (residual direction) to V basis */
  lbz_col_append(&d_Vbasis, d_v_cur);
  af_release_array(d_v_cur);  d_v_cur = 0;

  /* ── Materialise results to R ────────────────────────────────── */
  SEXP U_r     = PROTECT(lbz_download_mat(d_Ubasis));
  SEXP V_r     = PROTECT(lbz_download_mat(d_Vbasis));
  SEXP alpha_r = PROTECT(allocVector(REALSXP, k));
  SEXP beta_r  = PROTECT(allocVector(REALSXP, k));
  for (int i = 0; i < k; i++) {
    REAL(alpha_r)[i] = alpha_h[i];
    REAL(beta_r)[i]  = beta_h[i];
  }
  SEXP result = PROTECT(allocVector(VECSXP, 4));
  SEXP names  = PROTECT(allocVector(STRSXP, 4));
  SET_VECTOR_ELT(result, 0, U_r);      SET_STRING_ELT(names, 0, mkChar("U"));
  SET_VECTOR_ELT(result, 1, V_r);      SET_STRING_ELT(names, 1, mkChar("V"));
  SET_VECTOR_ELT(result, 2, alpha_r);  SET_STRING_ELT(names, 2, mkChar("alpha"));
  SET_VECTOR_ELT(result, 3, beta_r);   SET_STRING_ELT(names, 3, mkChar("beta"));
  setAttrib(result, R_NamesSymbol, names);
  UNPROTECT(6);

  af_release_array(d_A);
  if (d_Ubasis) af_release_array(d_Ubasis);
  if (d_Vbasis) af_release_array(d_Vbasis);
  if (d_u_prev) af_release_array(d_u_prev);

  return result;
#else
  error("am_af_lanczos_bidiag requires ArrayFire");
  return R_NilValue;
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
   LBZ A cache management + warm-start Lanczos for IRLBA thick restart
   ═══════════════════════════════════════════════════════════════════════════ */

/*
 * am_af_lbz_upload_A_bridge(A)
 *
 * Upload A to the LBZ GPU cache (column-major float32).  Must be called once
 * before any am_af_lanczos_warm_bridge calls.  Replaces any previously cached A.
 */
SEXP am_af_lbz_upload_A_bridge(SEXP A_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isMatrix(A_r) || !isReal(A_r))
    error("A must be a real matrix");
  SEXP dim_r = getAttrib(A_r, R_DimSymbol);
  g_lbz_A_m  = INTEGER(dim_r)[0];
  g_lbz_A_n  = INTEGER(dim_r)[1];
  /* Detect float64 support by probing: try creating a tiny f64 array.
   * Metal (Apple Silicon) = f32 only; CUDA = f64 available. */
  double probe_val = 1.0;
  dim_t  probe_dims[1] = {1};
  af_array probe = 0;
  bool dbl_ok = (af_create_array(&probe, &probe_val, 1, probe_dims, f64) == AF_SUCCESS
                 && probe != 0);
  if (probe) af_release_array(probe);
  g_lbz_dtype = dbl_ok ? f64 : f32;
  if (g_lbz_A) { af_release_array(g_lbz_A); g_lbz_A = 0; }
  g_lbz_A = lbz_upload_mat(REAL(A_r), g_lbz_A_m, g_lbz_A_n);
  if (g_lbz_A == 0) error("lbz_upload_mat failed (af_create_array returned 0)");
  return ScalarLogical(dbl_ok ? 2 : 1);  /* 2 = f64, 1 = f32 */
#else
  error("requires ArrayFire"); return R_NilValue;
#endif
}

/*
 * am_af_lbz_drop_A_bridge()
 *
 * Release the cached A from GPU memory.
 */
SEXP am_af_lbz_drop_A_bridge(void) {
#ifdef HAVE_ARRAYFIRE
  if (g_lbz_A) { af_release_array(g_lbz_A); g_lbz_A = 0; }
  g_lbz_A_m = g_lbz_A_n = 0;
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}

/*
 * am_af_lanczos_warm_bridge(V_warm, U_warm, p0, k)
 *
 * Warm-start Golub-Kahan Lanczos bidiagonalization using the cached A
 * (g_lbz_A).  A stays on GPU across all restarts; only the small warm
 * basis and starting vector cross PCIe each restart.
 *
 * V_warm : n × nv_warm real matrix of right Ritz vectors (or NULL for cold)
 * U_warm : m × nv_warm real matrix of left  Ritz vectors (or NULL for cold)
 * p0     : n-vector starting direction; normalized + reorthogonalized vs V_warm
 * k      : scalar integer — number of NEW Lanczos steps to run
 *
 * Returns list(U = m×(nv_warm+k),  V = n×(nv_warm+k+1),
 *              alpha = k,           beta = k)
 * where alpha/beta describe the k NEW steps only.  The full U and V
 * (warm + new) are returned so that the R side can assemble B_full
 * and compute the thick-restart rotation without re-uploading.
 * beta[k] (the last element) is the residual norm used for convergence.
 */
SEXP am_af_lanczos_warm_bridge(SEXP V_warm_r, SEXP U_warm_r, SEXP p0_r, SEXP k_r) {
#ifdef HAVE_ARRAYFIRE
  if (g_lbz_A == 0)
    error("no A cached; call am_af_lbz_upload_A_bridge first");

  int m = g_lbz_A_m, n = g_lbz_A_n;

  if (!isReal(p0_r) || LENGTH(p0_r) != n)
    error("p0 must be a real vector of length %d", n);
  if (!isInteger(k_r) || LENGTH(k_r) != 1)
    error("k must be a scalar integer");

  int k       = INTEGER(k_r)[0];
  int nv_warm = 0;

  if (!isNull(V_warm_r)) {
    if (!isMatrix(V_warm_r) || !isReal(V_warm_r))
      error("V_warm must be a real matrix or NULL");
    SEXP vdim = getAttrib(V_warm_r, R_DimSymbol);
    if (INTEGER(vdim)[0] != n)
      error("V_warm nrow must equal ncol(A)=%d", n);
    nv_warm = INTEGER(vdim)[1];
    if (isNull(U_warm_r) || !isMatrix(U_warm_r) || !isReal(U_warm_r))
      error("U_warm must be a real matrix when V_warm is given");
    SEXP udim = getAttrib(U_warm_r, R_DimSymbol);
    if (INTEGER(udim)[0] != m || INTEGER(udim)[1] != nv_warm)
      error("U_warm must be %d x %d", m, nv_warm);
  }
  if (k < 1) error("k must be >= 1");

  double *alpha_h = (double *)R_alloc(k, sizeof(double));
  double *beta_h  = (double *)R_alloc(k, sizeof(double));

  /* ── Initialise warm bases ─────────────────────────────────── */
  af_array d_Ubasis = 0;  /* m × (nv_warm + j)     */
  af_array d_Vbasis = 0;  /* n × (nv_warm + j + 1) */

  if (nv_warm > 0) {
    d_Ubasis = lbz_upload_mat(REAL(U_warm_r), m, nv_warm);
    d_Vbasis = lbz_upload_mat(REAL(V_warm_r), n, nv_warm);
  }

  /* ── Normalise p0 and reorthogonalise against V_warm ────────── */
  af_array d_tmp = lbz_upload_vec(REAL(p0_r), n);
  double nrm = 0.0;
  af_array d_v_cur = 0;
  lbz_normalize(&d_v_cur, &nrm, d_tmp);
  af_release_array(d_tmp);
  if (nrm < 1e-14) {
    if (d_Ubasis) af_release_array(d_Ubasis);
    if (d_Vbasis) af_release_array(d_Vbasis);
    if (d_v_cur)  af_release_array(d_v_cur);
    error("p0 is essentially zero");
  }

  if (d_Vbasis != 0) {
    /* CGS2 of p0 against V_warm ensures numerical orthogonality */
    af_array d_orth = 0;
    lbz_cgs2(&d_orth, d_v_cur, d_Vbasis);
    af_release_array(d_v_cur);
    double nrm2 = 0.0;
    af_array d_renorm = 0;
    lbz_normalize(&d_renorm, &nrm2, d_orth);
    af_release_array(d_orth);
    if (nrm2 < 1e-14) {
      if (d_Ubasis) af_release_array(d_Ubasis);
      if (d_Vbasis) af_release_array(d_Vbasis);
      if (d_renorm) af_release_array(d_renorm);
      error("p0 lies in span(V_warm) after reorthogonalisation");
    }
    d_v_cur = d_renorm;
  }

  /* ── k new Lanczos steps ────────────────────────────────────────
   * The explicit three-term subtractions (u -= beta*u_prev,
   * v -= alpha*v_cur) are omitted: CGS2 against the full accumulated
   * basis subsumes them in exact arithmetic and is more numerically
   * stable in finite precision.  Saves ~4 dispatch calls per step.
   * ──────────────────────────────────────────────────────────────── */
  for (int j = 0; j < k; j++) {

    /* Append v_j to V basis (before backward step so CGS2 sees v_j) */
    lbz_col_append(&d_Vbasis, d_v_cur);

    /* Forward: u_raw = A * v_j */
    af_array d_u_raw = 0;
    af_matmul(&d_u_raw, g_lbz_A, d_v_cur, AF_MAT_NONE, AF_MAT_NONE);

    /* CGS2 against U_basis — implicitly removes beta_{j-1}*u_{j-1}
     * and all prior U components (foot correction on warm restart). */
    if (d_Ubasis != 0) {
      af_array d_u_orth = 0;
      lbz_cgs2(&d_u_orth, d_u_raw, d_Ubasis);
      af_release_array(d_u_raw);
      d_u_raw = d_u_orth;
    }

    /* alpha_j = ||u||;  u_j = normalise */
    af_array d_uj = 0;
    lbz_normalize(&d_uj, &alpha_h[j], d_u_raw);
    af_release_array(d_u_raw);
    lbz_col_append(&d_Ubasis, d_uj);

    /* Backward: v_raw = A^T * u_j */
    af_array d_v_raw = 0;
    af_matmul(&d_v_raw, g_lbz_A, d_uj, AF_MAT_TRANS, AF_MAT_NONE);
    af_release_array(d_uj);

    /* CGS2 against V_basis (includes v_j) — implicitly removes
     * alpha_j*v_j and all prior V components. */
    {
      af_array d_v_orth = 0;
      lbz_cgs2(&d_v_orth, d_v_raw, d_Vbasis);
      af_release_array(d_v_raw);
      d_v_raw = d_v_orth;
    }

    /* beta_j = ||v||;  v_{j+1} = normalise */
    af_array d_v_next = 0;
    lbz_normalize(&d_v_next, &beta_h[j], d_v_raw);
    af_release_array(d_v_raw);

    af_release_array(d_v_cur);
    d_v_cur = d_v_next;
  }

  /* Append residual direction to V */
  lbz_col_append(&d_Vbasis, d_v_cur);
  af_release_array(d_v_cur); d_v_cur = 0;

  /* ── Materialise ─────────────────────────────────────────────── */
  SEXP U_r     = PROTECT(lbz_download_mat(d_Ubasis));  /* m × (nv_warm+k)   */
  SEXP V_r     = PROTECT(lbz_download_mat(d_Vbasis));  /* n × (nv_warm+k+1) */
  SEXP alpha_r = PROTECT(allocVector(REALSXP, k));
  SEXP beta_r  = PROTECT(allocVector(REALSXP, k));
  for (int i = 0; i < k; i++) {
    REAL(alpha_r)[i] = alpha_h[i];
    REAL(beta_r)[i]  = beta_h[i];
  }
  SEXP result = PROTECT(allocVector(VECSXP, 4));
  SEXP names  = PROTECT(allocVector(STRSXP, 4));
  SET_VECTOR_ELT(result, 0, U_r);     SET_STRING_ELT(names, 0, mkChar("U"));
  SET_VECTOR_ELT(result, 1, V_r);     SET_STRING_ELT(names, 1, mkChar("V"));
  SET_VECTOR_ELT(result, 2, alpha_r); SET_STRING_ELT(names, 2, mkChar("alpha"));
  SET_VECTOR_ELT(result, 3, beta_r);  SET_STRING_ELT(names, 3, mkChar("beta"));
  setAttrib(result, R_NamesSymbol, names);
  UNPROTECT(6);

  af_release_array(d_Ubasis);
  af_release_array(d_Vbasis);

  return result;
#else
  error("requires ArrayFire"); return R_NilValue;
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
   Pairwise Distance and Kernel Matrix Bridges (ArrayFire)
   ───────────────────────────────────────────────────────────────────────────
   These bridges use COLUMN-MAJOR upload/download (not the row-major trick
   used by the existing resident bridges), so af_matmul(X, Y, NONE, TRANS)
   = X %*% Y^T works for arbitrary m×p and n×p inputs without the
   nrow(X)==nrow(Y) constraint of the non-resident bridges.
   ═══════════════════════════════════════════════════════════════════════════ */

#ifdef HAVE_ARRAYFIRE

/* Upload R double matrix [nrow×ncol] to AF in column-major order. */
static af_array af_dist_upload(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  int  nrow = INTEGER(dim)[0], ncol = INTEGER(dim)[1];
  dim_t dims[2] = {(dim_t)nrow, (dim_t)ncol};
  size_t n = (size_t)nrow * (size_t)ncol;
  float* buf = (float*) arrayfire_xmalloc(n * sizeof(float));
  const double* src = REAL(x);
  for (size_t i = 0; i < n; i++) buf[i] = (float)src[i];
  af_array out = 0;
  af_create_array(&out, buf, 2, dims, f32);
  free(buf);
  return out;
}

/* Download AF [nrow×ncol] result to R double matrix (column-major). */
static SEXP af_dist_download(af_array arr) {
  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], arr);
  int nrow = (int)dims[0], ncol = (int)dims[1];
  size_t n = (size_t)nrow * (size_t)ncol;
  float* buf = (float*) arrayfire_xmalloc(n * sizeof(float));
  af_get_data_ptr(buf, arr);
  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  double* dst = REAL(out);
  for (size_t i = 0; i < n; i++) dst[i] = (double)buf[i];
  free(buf);
  UNPROTECT(1);
  return out;
}

/* D²[i,j] = ||x_i-y_j||² = nx[i]+ny[j]-2*(X@Y^T)[i,j]
 * ax [m×p], ay [n×p] in column-major AF.
 * Returns a new [m×n] AF array; caller must release. */
static af_err af_dist_sq_colmaj(af_array *D_out, af_array ax, af_array ay) {
  af_array G = 0, X2 = 0, nx = 0, Y2 = 0, ny = 0, ny_t = 0;
  af_array D = 0, sc2 = 0, G2 = 0, D_raw = 0, sc0 = 0;
  dim_t sc_dims[1] = {1};
  af_err err;

  err = af_matmul(&G, ax, ay, AF_MAT_NONE, AF_MAT_TRANS); /* X %*% Y^T */
  if (err != AF_SUCCESS) goto af_dsq_fail;

  err = af_mul(&X2, ax, ax, false);      if (err != AF_SUCCESS) goto af_dsq_fail;
  err = af_sum(&nx, X2, 1);             if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(X2); X2 = 0;

  err = af_mul(&Y2, ay, ay, false);      if (err != AF_SUCCESS) goto af_dsq_fail;
  err = af_sum(&ny, Y2, 1);             if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(Y2); Y2 = 0;

  err = af_transpose(&ny_t, ny, false);  if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(ny); ny = 0;

  err = af_add(&D, nx, ny_t, true);     if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(nx);   nx   = 0;
  af_release_array(ny_t); ny_t = 0;

  err = af_constant(&sc2, 2.0, 1, sc_dims, f32);  if (err != AF_SUCCESS) goto af_dsq_fail;
  err = af_mul(&G2, G, sc2, true);                 if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(sc2); sc2 = 0;
  af_release_array(G);   G   = 0;

  err = af_sub(&D_raw, D, G2, false);   if (err != AF_SUCCESS) goto af_dsq_fail;
  af_release_array(D);  D  = 0;
  af_release_array(G2); G2 = 0;

  err = af_constant(&sc0, 0.0, 1, sc_dims, f32);  if (err != AF_SUCCESS) goto af_dsq_fail;
  err = af_maxof(D_out, D_raw, sc0, true);
  af_release_array(D_raw); D_raw = 0;
  af_release_array(sc0);   sc0   = 0;
  return err;

af_dsq_fail:
  if (G)     af_release_array(G);
  if (X2)    af_release_array(X2);
  if (nx)    af_release_array(nx);
  if (Y2)    af_release_array(Y2);
  if (ny)    af_release_array(ny);
  if (ny_t)  af_release_array(ny_t);
  if (D)     af_release_array(D);
  if (sc2)   af_release_array(sc2);
  if (G2)    af_release_array(G2);
  if (D_raw) af_release_array(D_raw);
  if (sc0)   af_release_array(sc0);
  return err;
}

#endif /* HAVE_ARRAYFIRE */

SEXP am_af_dist_sq_bridge(SEXP X_r, SEXP Y_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isReal(X_r) || !isMatrix(X_r)) error("X must be a numeric matrix");
  SEXP dimX = getAttrib(X_r, R_DimSymbol);
  int pX = INTEGER(dimX)[1];
  int own_Y = !isNull(Y_r);
  if (own_Y) {
    if (!isReal(Y_r) || !isMatrix(Y_r)) error("Y must be a numeric matrix or NULL");
    if (INTEGER(getAttrib(Y_r, R_DimSymbol))[1] != pX)
      error("X and Y must have the same number of columns");
  }
  af_array ax = af_dist_upload(X_r);
  af_array ay = own_Y ? af_dist_upload(Y_r) : ax;
  af_array D_sq = 0;
  af_err err = af_dist_sq_colmaj(&D_sq, ax, ay);
  if (own_Y) af_release_array(ay);
  af_release_array(ax);
  if (err != AF_SUCCESS)
    error("am_af_dist_sq: failed (AF error %d)", (int)err);
  SEXP result = PROTECT(af_dist_download(D_sq));
  af_release_array(D_sq);
  UNPROTECT(1);
  return result;
#else
  error("am_af_dist_sq requires ArrayFire"); return R_NilValue;
#endif
}

SEXP am_af_kernel_bridge(SEXP X_r, SEXP Y_r, SEXP kernel_r,
                          SEXP sigma_r, SEXP degree_r, SEXP coef_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(kernel_r) || LENGTH(kernel_r) != 1)
    error("kernel must be a scalar string");
  const char* kern  = CHAR(STRING_ELT(kernel_r, 0));
  double sigma      = asReal(sigma_r);
  int    degree     = asInteger(degree_r);
  double coef_val   = asReal(coef_r);
  dim_t  sc_dims[1] = {1};

  if (!isReal(X_r) || !isMatrix(X_r)) error("X must be a numeric matrix");
  int pX = INTEGER(getAttrib(X_r, R_DimSymbol))[1];
  int own_Y = !isNull(Y_r);
  if (own_Y) {
    if (!isReal(Y_r) || !isMatrix(Y_r)) error("Y must be a numeric matrix or NULL");
    if (INTEGER(getAttrib(Y_r, R_DimSymbol))[1] != pX)
      error("X and Y must have the same number of columns");
  }

  af_array ax = af_dist_upload(X_r);
  af_array ay = own_Y ? af_dist_upload(Y_r) : ax;
  af_array K  = 0;
  af_err err  = AF_SUCCESS;

  if (strcmp(kern, "linear") == 0) {
    err = af_matmul(&K, ax, ay, AF_MAT_NONE, AF_MAT_TRANS);

  } else if (strcmp(kern, "rbf") == 0) {
    af_array D_sq = 0, sc = 0, neg_D = 0;
    err = af_dist_sq_colmaj(&D_sq, ax, ay);
    if (err == AF_SUCCESS) err = af_constant(&sc, -1.0/(2.0*sigma*sigma), 1, sc_dims, f32);
    if (err == AF_SUCCESS) err = af_mul(&neg_D, D_sq, sc, true);
    if (D_sq) af_release_array(D_sq); if (sc) af_release_array(sc);
    if (err == AF_SUCCESS) err = af_exp(&K, neg_D);
    if (neg_D) af_release_array(neg_D);

  } else if (strcmp(kern, "laplacian") == 0) {
    af_array D_sq = 0, D = 0, sc = 0, neg_D = 0;
    err = af_dist_sq_colmaj(&D_sq, ax, ay);
    if (err == AF_SUCCESS) err = af_sqrt(&D, D_sq);
    if (D_sq) af_release_array(D_sq);
    if (err == AF_SUCCESS) err = af_constant(&sc, -1.0/sigma, 1, sc_dims, f32);
    if (err == AF_SUCCESS) err = af_mul(&neg_D, D, sc, true);
    if (D) af_release_array(D); if (sc) af_release_array(sc);
    if (err == AF_SUCCESS) err = af_exp(&K, neg_D);
    if (neg_D) af_release_array(neg_D);

  } else if (strcmp(kern, "polynomial") == 0) {
    af_array G = 0, Gs = 0, ca = 0, da = 0;
    err = af_matmul(&G, ax, ay, AF_MAT_NONE, AF_MAT_TRANS);
    if (err == AF_SUCCESS) err = af_constant(&ca, coef_val, 1, sc_dims, f32);
    if (err == AF_SUCCESS) err = af_add(&Gs, G, ca, true);
    if (G) af_release_array(G); if (ca) af_release_array(ca);
    if (err == AF_SUCCESS) err = af_constant(&da, (double)degree, 1, sc_dims, f32);
    if (err == AF_SUCCESS) err = af_pow(&K, Gs, da, true);
    if (Gs) af_release_array(Gs); if (da) af_release_array(da);

  } else if (strcmp(kern, "cosine") == 0) {
    af_array G = 0, X2 = 0, nx = 0, Y2 = 0, ny = 0, ny_t = 0;
    af_array nxny = 0, nrm = 0, ea = 0, nrm2 = 0;
    err = af_matmul(&G, ax, ay, AF_MAT_NONE, AF_MAT_TRANS);
    if (err == AF_SUCCESS) err = af_mul(&X2, ax, ax, false);
    if (err == AF_SUCCESS) err = af_sum(&nx, X2, 1);
    if (X2) af_release_array(X2);
    if (err == AF_SUCCESS) err = af_mul(&Y2, ay, ay, false);
    if (err == AF_SUCCESS) err = af_sum(&ny, Y2, 1);
    if (Y2) af_release_array(Y2);
    if (err == AF_SUCCESS) err = af_transpose(&ny_t, ny, false);
    if (ny) af_release_array(ny);
    if (err == AF_SUCCESS) err = af_mul(&nxny, nx, ny_t, true);
    if (nx) af_release_array(nx); if (ny_t) af_release_array(ny_t);
    if (err == AF_SUCCESS) err = af_sqrt(&nrm, nxny);
    if (nxny) af_release_array(nxny);
    if (err == AF_SUCCESS) err = af_constant(&ea, 1e-12, 1, sc_dims, f32);
    if (err == AF_SUCCESS) err = af_add(&nrm2, nrm, ea, true);
    if (nrm) af_release_array(nrm); if (ea) af_release_array(ea);
    if (err == AF_SUCCESS) err = af_div(&K, G, nrm2, false);
    if (G) af_release_array(G); if (nrm2) af_release_array(nrm2);

  } else {
    if (own_Y) af_release_array(ay);
    af_release_array(ax);
    error("unknown kernel '%s': use linear/rbf/laplacian/polynomial/cosine", kern);
    return R_NilValue;
  }

  if (own_Y) af_release_array(ay);
  af_release_array(ax);
  if (err != AF_SUCCESS) error("am_af_kernel: computation failed (AF error %d)", (int)err);
  SEXP result = PROTECT(af_dist_download(K));
  af_release_array(K);
  UNPROTECT(1);
  return result;
#else
  error("am_af_kernel requires ArrayFire"); return R_NilValue;
#endif
}
