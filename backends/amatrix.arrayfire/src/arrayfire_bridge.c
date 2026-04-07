#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
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

/* ── Correct column-major helpers (used by resident path and rsvd) ────── *
 * R and ArrayFire both use column-major storage, so we can pass R data    *
 * directly without the row-major staging buffer in arrayfire_matrix_from_r. *
 * Use these everywhere you need correct non-square matrix support.         */

/* Cast R double column-major → float32, preserving element order */
static float* amatrix_af_r_to_f32(const double *src, int n) {
  float *buf = (float *) arrayfire_xmalloc((size_t)n * sizeof(float));
  for (int k = 0; k < n; k++) buf[k] = (float)src[k];
  return buf;
}

/* Build an AF array from an R matrix (column-major, double→float32) */
static af_array amatrix_af_from_r(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  int m = INTEGER(dim)[0], n = INTEGER(dim)[1];
  dim_t dims[2] = {(dim_t)m, (dim_t)n};
  float *buf = amatrix_af_r_to_f32(REAL(x), m * n);
  af_array out = 0;
  af_create_array(&out, buf, 2, dims, f32);
  free(buf);
  return out;
}

/* Convert AF array → R matrix (column-major, float32→double) */
static SEXP amatrix_af_to_r(af_array arr) {
  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], arr);
  int m = (int)dims[0], n = (int)dims[1];
  float *buf = (float *) arrayfire_xmalloc((size_t)m * n * sizeof(float));
  af_get_data_ptr(buf, arr);
  SEXP out = PROTECT(allocMatrix(REALSXP, m, n));
  double *dst = REAL(out);
  for (int k = 0; k < m * n; k++) dst[k] = (double)buf[k];
  free(buf);
  UNPROTECT(1);
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
  } else if (strcmp(op_name, "/") == 0) {
    err = af_div(&out, alhs, arhs, false);
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
  e->array = amatrix_af_from_r(x);   /* column-major, correct for all shapes */
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
  return amatrix_af_to_r(e->array);   /* column-major, correct for all shapes */
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
  /* crossprod(x,y) = t(x) %*% y.  Resident arrays use column-major (amatrix_af_from_r),
     so ax = x_orig.  AF_MAT_TRANS on lhs gives t(x_orig) %*% y_orig — correct for
     all shapes including cross-covariance (p != q). */
  af_err err = af_matmul(&eout->array, ex->array, ay, AF_MAT_TRANS, AF_MAT_NONE);
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
  /* tcrossprod(x,y) = x %*% t(y).  Resident arrays use column-major (amatrix_af_from_r),
     so ax = x_orig, ay = y_orig.  AF_MAT_TRANS on rhs gives x_orig %*% t(y_orig). */
  af_err err = af_matmul(&eout->array, ex->array, ay, AF_MAT_NONE, AF_MAT_TRANS);
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
  else if (strcmp(op_name, "/") == 0) err = af_div(&eout->array, elhs->array, arhs, false);
  else { if (own_arhs) af_release_array(arhs); error("unsupported ewise op: %s", op_name); }

  if (own_arhs) af_release_array(arhs);
  if (err != AF_SUCCESS) error("arrayfire ewise (resident) failed");
  return ScalarLogical(1);
#else
  error("arrayfire ewise_resident requires arrayfire");
#endif
}

SEXP amatrix_arrayfire_scatter_mean_bridge(SEXP lhs_key, SEXP labels_r, SEXP K_r) {
#ifdef HAVE_ARRAYFIRE
  /* Returns K×p group-sum matrix (R divides by counts to get means) */
  if (!isString(lhs_key) || LENGTH(lhs_key) != 1)
    error("lhs_key must be a scalar string");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(lhs_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(lhs_key, 0)));

  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], ex->array);
  int n = (int)dims[0], p = (int)dims[1];
  int K = INTEGER(K_r)[0];
  const int* labels = INTEGER(labels_r); /* 1-indexed */

  /* Build W (n×K) in column-major: W[i, k] = 1.0 if labels[i]-1 == k */
  float* w_buf = (float*)arrayfire_xmalloc((size_t)n * (size_t)K * sizeof(float));
  memset(w_buf, 0, (size_t)n * (size_t)K * sizeof(float));
  for (int i = 0; i < n; i++) {
    int k = labels[i] - 1; /* 0-indexed */
    if (k >= 0 && k < K) w_buf[(size_t)k * n + i] = 1.0f; /* col-major: col k, row i */
  }

  dim_t wdims[2] = {(dim_t)n, (dim_t)K};
  af_array W = 0;
  af_err err = af_create_array(&W, w_buf, 2, wdims, f32);
  free(w_buf);
  if (err != AF_SUCCESS) error("af_create_array failed for W in scatter_mean");

  /* t(W) %*% X = [K,n] × [n,p] → [K,p] */
  af_array result = 0;
  err = af_matmul(&result, W, ex->array, AF_MAT_TRANS, AF_MAT_NONE);
  af_release_array(W);
  if (err != AF_SUCCESS) error("af_matmul failed in scatter_mean");

  SEXP out = PROTECT(amatrix_af_to_r(result));
  af_release_array(result);
  UNPROTECT(1);
  return out;
#else
  error("arrayfire scatter_mean requires arrayfire");
  return R_NilValue;
#endif
}

/* ── segment_sum / segment_mean (amatrix-ylo) ────────────────────────────── */

SEXP amatrix_arrayfire_segment_sum_bridge(SEXP x_key, SEXP labels_r, SEXP K_r, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || LENGTH(x_key) != 1)
    error("x_key must be a scalar string");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire segment_sum: resident key not found");

  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], ex->array);
  int n = (int)dims[0], p = (int)dims[1];
  int K = INTEGER(K_r)[0];
  const int* labels = INTEGER(labels_r);

  /* W (n×K) col-major: W[i, k] = 1.0 if labels[i]-1 == k */
  float* w_buf = (float*)arrayfire_xmalloc((size_t)n * (size_t)K * sizeof(float));
  memset(w_buf, 0, (size_t)n * (size_t)K * sizeof(float));
  for (int i = 0; i < n; i++) {
    int k = labels[i] - 1;
    if (k >= 0 && k < K) w_buf[(size_t)k * n + i] = 1.0f;
  }
  dim_t wdims[2] = {(dim_t)n, (dim_t)K};
  af_array W = 0;
  af_err err = af_create_array(&W, w_buf, 2, wdims, f32);
  free(w_buf);
  if (err != AF_SUCCESS) error("arrayfire segment_sum: af_create_array failed");

  /* sums = W^T %*% X : [K,n] × [n,p] → [K,p] */
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  err = af_matmul(&eout->array, W, ex->array, AF_MAT_TRANS, AF_MAT_NONE);
  af_release_array(W);
  if (err != AF_SUCCESS) error("arrayfire segment_sum: af_matmul failed");
  return ScalarLogical(1);
#else
  error("arrayfire segment_sum requires arrayfire");
  return R_NilValue;
#endif
}

SEXP amatrix_arrayfire_segment_mean_bridge(SEXP x_key, SEXP labels_r, SEXP K_r, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || LENGTH(x_key) != 1)
    error("x_key must be a scalar string");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire segment_mean: resident key not found");

  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], ex->array);
  int n = (int)dims[0], p = (int)dims[1];
  int K = INTEGER(K_r)[0];
  const int* labels = INTEGER(labels_r);

  /* W (n×K) + counts (K) */
  float* w_buf   = (float*)arrayfire_xmalloc((size_t)n * (size_t)K * sizeof(float));
  float* cnt_buf = (float*)arrayfire_xmalloc((size_t)K * sizeof(float));
  memset(w_buf,   0, (size_t)n * (size_t)K * sizeof(float));
  memset(cnt_buf, 0, (size_t)K * sizeof(float));
  for (int i = 0; i < n; i++) {
    int k = labels[i] - 1;
    if (k >= 0 && k < K) {
      w_buf[(size_t)k * n + i] = 1.0f;
      cnt_buf[k] += 1.0f;
    }
  }
  dim_t wdims[2] = {(dim_t)n, (dim_t)K};
  af_array W = 0;
  af_err err = af_create_array(&W, w_buf, 2, wdims, f32);
  free(w_buf);
  if (err != AF_SUCCESS) { free(cnt_buf); error("arrayfire segment_mean: W create failed"); }

  /* sums = W^T %*% X : [K,p] */
  af_array sums = 0;
  err = af_matmul(&sums, W, ex->array, AF_MAT_TRANS, AF_MAT_NONE);
  af_release_array(W);
  if (err != AF_SUCCESS) { free(cnt_buf); error("arrayfire segment_mean: af_matmul failed"); }

  /* counts as K×1 AF array */
  dim_t cdims[2] = {(dim_t)K, 1};
  af_array counts = 0;
  err = af_create_array(&counts, cnt_buf, 2, cdims, f32);
  free(cnt_buf);
  if (err != AF_SUCCESS) { af_release_array(sums); error("arrayfire segment_mean: counts create failed"); }

  /* means = sums / counts (K×p / K×1 broadcast; 0/0 → NaN for empty clusters) */
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  err = af_div(&eout->array, sums, counts, true);
  af_release_array(sums);
  af_release_array(counts);
  if (err != AF_SUCCESS) error("arrayfire segment_mean: af_div failed");
  return ScalarLogical(1);
#else
  error("arrayfire segment_mean requires arrayfire");
  return R_NilValue;
#endif
}

SEXP amatrix_arrayfire_argreduce_bridge(SEXP lhs_key, SEXP axis_r, SEXP is_max_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(lhs_key) || LENGTH(lhs_key) != 1)
    error("lhs_key must be a scalar string");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(lhs_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(lhs_key, 0)));

  int axis   = INTEGER(axis_r)[0];
  int is_max = LOGICAL(is_max_r)[0];

  /* AF column-major: axis=1 reduces cols → n row results; axis=0 reduces rows → K col results */
  dim_t dims[4] = {0, 0, 0, 0};
  af_get_dims(&dims[0], &dims[1], &dims[2], &dims[3], ex->array);
  int len = (axis == 0) ? (int)dims[1] : (int)dims[0];

  af_array out_val = 0, out_idx = 0;
  af_err err = is_max
    ? af_imax(&out_val, &out_idx, ex->array, axis)
    : af_imin(&out_val, &out_idx, ex->array, axis);
  if (err != AF_SUCCESS) error("af_imax/imin failed");
  af_release_array(out_val);

  uint32_t* buf = (uint32_t*)arrayfire_xmalloc((size_t)len * sizeof(uint32_t));
  af_get_data_ptr(buf, out_idx);
  af_release_array(out_idx);

  SEXP result = PROTECT(allocVector(INTSXP, len));
  int* ires = INTEGER(result);
  for (int i = 0; i < len; i++) ires[i] = (int)buf[i] + 1; /* 0-indexed → 1-indexed */
  free(buf);
  UNPROTECT(1);
  return result;
#else
  error("arrayfire argreduce requires arrayfire");
  return R_NilValue;
#endif
}

SEXP amatrix_arrayfire_broadcast_ewise_resident_bridge(SEXP lhs_key, SEXP v, SEXP margin_r, SEXP op, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(lhs_key) || !isString(op) || !isString(out_key))
    error("keys and op must be scalar strings");
  if (!isReal(v) && !isInteger(v))
    error("v must be a numeric vector");
  if (!isInteger(margin_r) || LENGTH(margin_r) != 1)
    error("margin must be a scalar integer");

  amatrix_af_resident_entry* elhs = amatrix_af_registry_find(CHAR(STRING_ELT(lhs_key, 0)));
  if (elhs == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(lhs_key, 0)));

  const char* op_name = CHAR(STRING_ELT(op, 0));
  int margin = INTEGER(margin_r)[0];
  int len_v = (int)XLENGTH(v);

  /* Build float32 buffer from v (1-D, column-major = identity for 1-D) */
  float* fbuf = (float*)arrayfire_xmalloc((size_t)len_v * sizeof(float));
  if (isReal(v)) {
    const double* src = REAL(v);
    for (int i = 0; i < len_v; i++) fbuf[i] = (float)src[i];
  } else {
    const int* src = INTEGER(v);
    for (int i = 0; i < len_v; i++) fbuf[i] = (float)src[i];
  }

  /* margin=1 → [n,1]; margin=2 → [1,K] — AF broadcasts with batch=true */
  dim_t vdims[2];
  if (margin == 1) { vdims[0] = (dim_t)len_v; vdims[1] = 1; }
  else             { vdims[0] = 1;             vdims[1] = (dim_t)len_v; }

  af_array av = 0;
  af_err err = af_create_array(&av, fbuf, 2, vdims, f32);
  free(fbuf);
  if (err != AF_SUCCESS) error("af_create_array failed for broadcast vector");

  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  if      (strcmp(op_name, "+") == 0) err = af_add(&eout->array, elhs->array, av, true);
  else if (strcmp(op_name, "-") == 0) err = af_sub(&eout->array, elhs->array, av, true);
  else if (strcmp(op_name, "*") == 0) err = af_mul(&eout->array, elhs->array, av, true);
  else if (strcmp(op_name, "/") == 0) err = af_div(&eout->array, elhs->array, av, true);
  else { af_release_array(av); error("unsupported broadcast ewise op: %s", op_name); }

  af_release_array(av);
  if (err != AF_SUCCESS) error("arrayfire broadcast ewise failed");
  return ScalarLogical(1);
#else
  error("arrayfire broadcast_ewise_resident requires arrayfire");
  return R_NilValue;
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

/* ── Cholesky factorization ────────────────────────────────────────────── */

SEXP amatrix_arrayfire_chol_bridge(SEXP x) {
  if (!isReal(x) || !isMatrix(x))
    error("x must be a numeric matrix");
#ifdef HAVE_ARRAYFIRE
  bool lapack = false;
  af_is_lapack_available(&lapack);
  if (!lapack)
    error("amatrix_arrayfire_chol: LAPACK not available in this ArrayFire build");
  af_array ax = 0, out = 0;
  int info = 0;
  ax = arrayfire_matrix_from_r(x);
  af_err err = af_cholesky(&out, &info, ax, true /* upper */);
  af_release_array(ax);
  if (err != AF_SUCCESS || info != 0) {
    if (out) af_release_array(out);
    error("af_cholesky failed (info=%d)", info);
  }
  SEXP result = arrayfire_result_to_r_matrix(out);
  af_release_array(out);
  return result;
#else
  error("amatrix_arrayfire_chol requires arrayfire");
  return R_NilValue;
#endif
}

SEXP amatrix_arrayfire_chol_resident_bridge(SEXP x_key, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || LENGTH(x_key) != 1 || !isString(out_key) || LENGTH(out_key) != 1)
    error("keys must be scalar strings");
  bool lapack = false;
  af_is_lapack_available(&lapack);
  if (!lapack)
    error("amatrix_arrayfire_chol_resident: LAPACK not available");
  amatrix_af_resident_entry* ex = amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  af_array out = 0;
  int info = 0;
  af_err err = af_cholesky(&out, &info, ex->array, true /* upper */);
  if (err != AF_SUCCESS || info != 0) {
    if (out) af_release_array(out);
    error("af_cholesky (resident) failed (info=%d)", info);
  }
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  eout->array = out;
  return ScalarLogical(1);
#else
  error("arrayfire chol_resident requires arrayfire");
  return R_NilValue;
#endif
}

/* ── Linear solve ──────────────────────────────────────────────────────── */

SEXP amatrix_arrayfire_solve_bridge(SEXP a, SEXP b) {
  if (!isReal(a) || !isMatrix(a))
    error("a must be a numeric matrix");
#ifdef HAVE_ARRAYFIRE
  bool lapack = false;
  af_is_lapack_available(&lapack);
  if (!lapack)
    error("amatrix_arrayfire_solve: LAPACK not available in this ArrayFire build");
  af_array aa = 0, ab = 0, out = 0;
  aa = arrayfire_matrix_from_r(a);
  if (isNull(b)) {
    SEXP dim = getAttrib(a, R_DimSymbol);
    int n = INTEGER(dim)[0];
    dim_t dims[2] = {(dim_t)n, (dim_t)n};
    af_identity(&ab, 2, dims, f32);
  } else {
    if (!isReal(b) || !isMatrix(b))
      error("b must be a numeric matrix or NULL");
    ab = arrayfire_matrix_from_r(b);
  }
  af_err err = af_solve(&out, aa, ab, AF_MAT_NONE);
  af_release_array(aa);
  af_release_array(ab);
  if (err != AF_SUCCESS) {
    if (out) af_release_array(out);
    error("af_solve failed");
  }
  SEXP result = arrayfire_result_to_r_matrix(out);
  af_release_array(out);
  return result;
#else
  error("amatrix_arrayfire_solve requires arrayfire");
  return R_NilValue;
#endif
}

SEXP amatrix_arrayfire_solve_resident_bridge(SEXP a_key, SEXP b_key, SEXP out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(a_key) || LENGTH(a_key) != 1 || !isString(out_key) || LENGTH(out_key) != 1)
    error("a_key and out_key must be scalar strings");
  bool lapack = false;
  af_is_lapack_available(&lapack);
  if (!lapack)
    error("amatrix_arrayfire_solve_resident: LAPACK not available");
  amatrix_af_resident_entry* ea = amatrix_af_registry_find(CHAR(STRING_ELT(a_key, 0)));
  if (ea == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(a_key, 0)));
  af_array ab = 0;
  int own_ab = 0;
  if (isNull(b_key)) {
    dim_t d[4] = {0, 0, 0, 0};
    af_get_dims(&d[0], &d[1], &d[2], &d[3], ea->array);
    dim_t dims[2] = {d[0], d[0]};
    af_identity(&ab, 2, dims, f32);
    own_ab = 1;
  } else {
    if (!isString(b_key) || LENGTH(b_key) != 1)
      error("b_key must be a scalar string or NULL");
    amatrix_af_resident_entry* eb = amatrix_af_registry_find(CHAR(STRING_ELT(b_key, 0)));
    if (eb == NULL) error("arrayfire resident key not found: %s", CHAR(STRING_ELT(b_key, 0)));
    ab = eb->array;
  }
  af_array out = 0;
  af_err err = af_solve(&out, ea->array, ab, AF_MAT_NONE);
  if (own_ab) af_release_array(ab);
  if (err != AF_SUCCESS) {
    if (out) af_release_array(out);
    error("af_solve (resident) failed");
  }
  amatrix_af_resident_entry* eout = amatrix_af_registry_reserve(CHAR(STRING_ELT(out_key, 0)));
  eout->array = out;
  return ScalarLogical(1);
#else
  error("arrayfire solve_resident requires arrayfire");
  return R_NilValue;
#endif
}

/* ── QR-Q resident: input key → thin Q stored as new resident key ────── */
SEXP amatrix_arrayfire_qr_Q_resident_bridge(SEXP x_key, SEXP q_out_key) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(x_key) || LENGTH(x_key) != 1 ||
      !isString(q_out_key) || LENGTH(q_out_key) != 1)
    error("keys must be scalar strings");
  amatrix_af_resident_entry* ex =
    amatrix_af_registry_find(CHAR(STRING_ELT(x_key, 0)));
  if (ex == NULL)
    error("arrayfire resident key not found: %s", CHAR(STRING_ELT(x_key, 0)));
  af_array Q = 0, R = 0, tau = 0;
  af_err err = af_qr(&Q, &R, &tau, ex->array);
  af_release_array(R);
  af_release_array(tau);
  if (err != AF_SUCCESS) {
    if (Q) af_release_array(Q);
    error("af_qr failed in qr_Q_resident (%d)", (int)err);
  }
  amatrix_af_resident_entry* eout =
    amatrix_af_registry_reserve(CHAR(STRING_ELT(q_out_key, 0)));
  eout->array = Q;
  return ScalarLogical(1);
#else
  error("amatrix_arrayfire_qr_Q_resident requires ArrayFire");
  return R_NilValue;
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

/* ── am_af_kernel_resident_bridge: compute kernel, store in resident registry ─
 * Same computation as am_af_kernel_bridge but stores K as a resident af_array
 * instead of downloading to R.  Eliminates the CPU round-trip when the caller
 * will immediately use the result as a GPU operand.
 *
 * zero_diag_r: logical; if TRUE and Y is NULL (self-kernel), zeros the diagonal
 * on GPU before storing (avoids a separate diag<-0 + re-upload cycle).
 */
SEXP am_af_kernel_resident_bridge(SEXP out_key_r, SEXP X_r, SEXP Y_r,
                                   SEXP kernel_r, SEXP sigma_r,
                                   SEXP degree_r, SEXP coef_r,
                                   SEXP zero_diag_r) {
#ifdef HAVE_ARRAYFIRE
  if (!isString(out_key_r) || LENGTH(out_key_r) != 1)
    error("out_key must be a scalar string");
  if (!isString(kernel_r) || LENGTH(kernel_r) != 1)
    error("kernel must be a scalar string");
  const char* out_key = CHAR(STRING_ELT(out_key_r, 0));
  const char* kern    = CHAR(STRING_ELT(kernel_r, 0));
  double sigma        = asReal(sigma_r);
  int    degree       = asInteger(degree_r);
  double coef_val     = asReal(coef_r);
  int    do_zero_diag = asLogical(zero_diag_r);
  dim_t  sc_dims[1]   = {1};

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
  if (err != AF_SUCCESS) {
    if (K) af_release_array(K);
    error("am_af_kernel_resident: computation failed (AF error %d)", (int)err);
  }

  /* Zero diagonal on GPU — avoids the diag(W)<-0 + re-upload round-trip */
  if (do_zero_diag && !own_Y && K != 0) {
    dim_t k_dims[4] = {0, 0, 0, 0};
    af_get_dims(&k_dims[0], &k_dims[1], &k_dims[2], &k_dims[3], K);
    dim_t eye_d[2] = {k_dims[0], k_dims[1]};
    af_array eye_m = 0, one_sc = 0, mask = 0, K_out = 0;
    dim_t sc1[1] = {1};
    af_err ze = af_identity(&eye_m, 2, eye_d, f32);
    if (ze == AF_SUCCESS) ze = af_constant(&one_sc, 1.0, 1, sc1, f32);
    if (ze == AF_SUCCESS) ze = af_sub(&mask, one_sc, eye_m, true);
    if (one_sc) af_release_array(one_sc);
    if (eye_m)  af_release_array(eye_m);
    if (ze == AF_SUCCESS) {
      ze = af_mul(&K_out, K, mask, false);
      af_release_array(mask);
      af_release_array(K);
      K = K_out;
    } else {
      if (mask) af_release_array(mask);
      /* zero_diag failed — continue without zeroing rather than error */
    }
  }

  /* Store in resident registry — caller owns the key and must drop it */
  amatrix_af_resident_entry* e = amatrix_af_registry_reserve(out_key);
  e->array = K;
  return ScalarLogical(1);
#else
  error("am_af_kernel_resident requires ArrayFire"); return R_NilValue;
#endif
}

/* ── amatrix_arrayfire_matmul_correct_bridge(A, B) → A %*% B ─────────── */
SEXP amatrix_arrayfire_matmul_correct_bridge(SEXP A_r, SEXP B_r) {
  if (!isReal(A_r) || !isMatrix(A_r) || !isReal(B_r) || !isMatrix(B_r))
    error("inputs must be real matrices");
#ifdef HAVE_ARRAYFIRE
  af_array A = 0, B = 0, C = 0;
  A = amatrix_af_from_r(A_r);
  B = amatrix_af_from_r(B_r);
  af_err err = af_matmul(&C, A, B, AF_MAT_NONE, AF_MAT_NONE);
  af_release_array(A); af_release_array(B);
  if (err != AF_SUCCESS) error("amatrix_arrayfire_matmul_correct: af_matmul failed (%d)", (int)err);
  SEXP result = PROTECT(amatrix_af_to_r(C));
  af_release_array(C);
  UNPROTECT(1);
  return result;
#else
  error("amatrix_arrayfire_matmul_correct requires ArrayFire");
  return R_NilValue;
#endif
}

/* ── amatrix_arrayfire_crossprod_correct_bridge(A, B) → t(A) %*% B ───── */
SEXP amatrix_arrayfire_crossprod_correct_bridge(SEXP A_r, SEXP B_r) {
  if (!isReal(A_r) || !isMatrix(A_r) || !isReal(B_r) || !isMatrix(B_r))
    error("inputs must be real matrices");
#ifdef HAVE_ARRAYFIRE
  af_array A = 0, B = 0, C = 0;
  A = amatrix_af_from_r(A_r);
  B = amatrix_af_from_r(B_r);
  af_err err = af_matmul(&C, A, B, AF_MAT_TRANS, AF_MAT_NONE);
  af_release_array(A); af_release_array(B);
  if (err != AF_SUCCESS) error("amatrix_arrayfire_crossprod_correct: af_matmul failed (%d)", (int)err);
  SEXP result = PROTECT(amatrix_af_to_r(C));
  af_release_array(C);
  UNPROTECT(1);
  return result;
#else
  error("amatrix_arrayfire_crossprod_correct requires ArrayFire");
  return R_NilValue;
#endif
}

/* ── amatrix_arrayfire_tcrossprod_correct_bridge(A, B) → A %*% t(B) ─── */
SEXP amatrix_arrayfire_tcrossprod_correct_bridge(SEXP A_r, SEXP B_r) {
  if (!isReal(A_r) || !isMatrix(A_r) || !isReal(B_r) || !isMatrix(B_r))
    error("inputs must be real matrices");
#ifdef HAVE_ARRAYFIRE
  af_array A = 0, B = 0, C = 0;
  A = amatrix_af_from_r(A_r);
  B = amatrix_af_from_r(B_r);
  af_err err = af_matmul(&C, A, B, AF_MAT_NONE, AF_MAT_TRANS);
  af_release_array(A); af_release_array(B);
  if (err != AF_SUCCESS) error("amatrix_arrayfire_tcrossprod_correct: af_matmul failed (%d)", (int)err);
  SEXP result = PROTECT(amatrix_af_to_r(C));
  af_release_array(C);
  UNPROTECT(1);
  return result;
#else
  error("amatrix_arrayfire_tcrossprod_correct requires ArrayFire");
  return R_NilValue;
#endif
}

/* ── Full thin SVD  (af_svd) → list(d, u, v) matching base::svd ────────
 *
 * af_svd(U, s, Vt, A) for A [m×n] returns:
 *   U  [m × k]   left  singular vectors  (thin, k = min(m,n))
 *   s  [k]        singular values in descending order
 *   Vt [k × n]   right singular vectors, TRANSPOSED
 *
 * R's base::svd returns V = Vt^T [n × k], so we call af_transpose on Vt.
 * nu / nv control how many columns of U / V to return (columns 1..nu, 1..nv).
 * d always has length k regardless of nu/nv — identical to base::svd behaviour.
 * ─────────────────────────────────────────────────────────────────────────── */
SEXP amatrix_arrayfire_svd_bridge(SEXP x, SEXP nu_r, SEXP nv_r) {
  if (!isReal(x) || !isMatrix(x))
    error("x must be a real numeric matrix");
  int nu = asInteger(nu_r);
  int nv = asInteger(nv_r);
  if (nu < 0 || nv < 0)
    error("nu and nv must be non-negative integers");
#ifdef HAVE_ARRAYFIRE
  bool lapack = false;
  af_is_lapack_available(&lapack);
  if (!lapack)
    error("amatrix_arrayfire_svd: LAPACK not available in this ArrayFire build");

  SEXP dim = getAttrib(x, R_DimSymbol);
  int m = INTEGER(dim)[0], n = INTEGER(dim)[1];
  int k = m < n ? m : n;   /* thin rank */

  af_array ax = 0, U_a = 0, s_a = 0, Vt_a = 0, V_a = 0;
  ax = amatrix_af_from_r(x);
  if (!ax) error("amatrix_arrayfire_svd: af_from_r returned NULL");
  af_err err = af_svd(&U_a, &s_a, &Vt_a, ax);
  af_release_array(ax);
  if (err != AF_SUCCESS) {
    if (U_a)  af_release_array(U_a);
    if (s_a)  af_release_array(s_a);
    if (Vt_a) af_release_array(Vt_a);
    error("amatrix_arrayfire_svd: af_svd failed (%d)", (int)err);
  }
  if (!U_a || !s_a || !Vt_a) {
    if (U_a)  af_release_array(U_a);
    if (s_a)  af_release_array(s_a);
    if (Vt_a) af_release_array(Vt_a);
    error("amatrix_arrayfire_svd: af_svd returned NULL output(s) U=%p s=%p Vt=%p",
          (void*)U_a, (void*)s_a, (void*)Vt_a);
  }

  /* Transpose Vt → V to match base::svd's $v convention */
  err = af_transpose(&V_a, Vt_a, false);
  af_release_array(Vt_a);
  if (err != AF_SUCCESS || !V_a) {
    if (U_a) af_release_array(U_a);
    if (s_a) af_release_array(s_a);
    if (V_a) af_release_array(V_a);
    error("amatrix_arrayfire_svd: af_transpose failed (%d)", (int)err);
  }

  /* Read actual output shapes from AF (af_svd may return full or thin matrices
   * depending on the backend/version — never assume thin without checking). */
  dim_t u_dims[4] = {0,0,0,0}, v_dims[4] = {0,0,0,0}, s_dims[4] = {0,0,0,0};
  af_get_dims(&u_dims[0], &u_dims[1], &u_dims[2], &u_dims[3], U_a);
  af_get_dims(&v_dims[0], &v_dims[1], &v_dims[2], &v_dims[3], V_a);
  af_get_dims(&s_dims[0], &s_dims[1], &s_dims[2], &s_dims[3], s_a);
  /* s_a is 1-D of length k; u_dims[0]=m, u_dims[1]=cols_u;
   * V_a = transpose(Vt) has dims [n, cols_vt] or possibly [cols_vt, n].
   * Actual k from s: */
  int k_actual = (int)s_dims[0];
  int u_rows = (int)u_dims[0], u_cols = (int)u_dims[1];
  int v_rows = (int)v_dims[0], v_cols = (int)v_dims[1];

  /* --- Singular values --- */
  size_t s_total = (size_t)k_actual;
  float *s_buf = (float *) arrayfire_xmalloc(s_total * sizeof(float));
  af_get_data_ptr(s_buf, s_a);
  af_release_array(s_a);
  SEXP d_r = PROTECT(allocVector(REALSXP, k_actual));
  for (int i = 0; i < k_actual; i++) REAL(d_r)[i] = (double)s_buf[i];
  free(s_buf);

  /* --- Left singular vectors U: copy first nu_eff columns --- */
  int nu_eff = (nu > k_actual) ? k_actual : nu;
  SEXP u_r = PROTECT(allocMatrix(REALSXP, m, nu_eff));  /* (1) */
  if (nu_eff > 0) {
    size_t u_total = (size_t)u_rows * u_cols;
    float *u_buf = (float *) arrayfire_xmalloc(u_total * sizeof(float));
    af_get_data_ptr(u_buf, U_a);
    double *up = REAL(u_r);
    /* AF column-major: U[i,j] = u_buf[i + u_rows * j]. */
    for (int j = 0; j < nu_eff; j++)
      for (int i = 0; i < m; i++)
        up[i + m * j] = (double)u_buf[i + u_rows * j];
    free(u_buf);
  }
  af_release_array(U_a);

  /* --- Right singular vectors V: copy first nv_eff columns ---
   * V_a = af_transpose(Vt_a) has shape [n × k_actual] (or [n × n] if full).
   * We want the first nv_eff columns of V.                                  */
  int nv_eff = (nv > k_actual) ? k_actual : nv;
  SEXP v_r = PROTECT(allocMatrix(REALSXP, n, nv_eff));  /* (2) */
  if (nv_eff > 0) {
    size_t v_total = (size_t)v_rows * v_cols;
    float *v_buf = (float *) arrayfire_xmalloc(v_total * sizeof(float));
    af_get_data_ptr(v_buf, V_a);
    double *vp = REAL(v_r);
    /* V_a is [n × v_cols] column-major: V[i,j] = v_buf[i + v_rows * j]. */
    for (int j = 0; j < nv_eff; j++)
      for (int i = 0; i < n; i++)
        vp[i + n * j] = (double)v_buf[i + v_rows * j];
    free(v_buf);
  }
  af_release_array(V_a);

  /* --- Assemble list(d, u, v) --- */
  SEXP names = PROTECT(allocVector(STRSXP, 3));   /* (3) */
  SET_STRING_ELT(names, 0, mkChar("d"));
  SET_STRING_ELT(names, 1, mkChar("u"));
  SET_STRING_ELT(names, 2, mkChar("v"));
  SEXP result = PROTECT(allocVector(VECSXP, 3));  /* (4) */
  SET_VECTOR_ELT(result, 0, d_r);
  SET_VECTOR_ELT(result, 1, u_r);
  SET_VECTOR_ELT(result, 2, v_r);
  setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(5);  /* d_r u_r v_r names result */
  return result;
#else
  error("amatrix_arrayfire_svd requires arrayfire");
  return R_NilValue;
#endif
}

/* ── amatrix_arrayfire_svd_safe_bridge() → logical(1) ───────────────────
 *
 * Returns TRUE if af_svd is known safe on the current active backend:
 *   CUDA  (AF_BACKEND_CUDA)    — cuBLAS/NVBLAS path, stable
 *   oneAPI (AF_BACKEND_ONEAPI) — MKL path, stable
 * Metal/OpenCL/CPU backends are NOT guaranteed stable here; use the R-level
 * subprocess probe or fall back to QR→SVD(R) for those.
 *
 * Compile-time override: define AMATRIX_AF_NATIVE_SVD_SAFE (e.g. for a known-
 * good CUDA build) to skip the runtime check and always return TRUE.
 * ─────────────────────────────────────────────────────────────────────────── */
SEXP amatrix_arrayfire_svd_safe_bridge(void) {
#ifdef AMATRIX_AF_NATIVE_SVD_SAFE
  return ScalarLogical(1);
#elif defined(HAVE_ARRAYFIRE)
  af_backend active = AF_BACKEND_DEFAULT;
  af_get_active_backend(&active);
  int safe = (active == AF_BACKEND_CUDA);
#  ifdef AF_BACKEND_ONEAPI
  if (active == AF_BACKEND_ONEAPI) safe = 1;
#  endif
  return ScalarLogical(safe);
#else
  return ScalarLogical(0);
#endif
}

/* ── amatrix_arrayfire_qr_q_correct_bridge(A) → thin Q [m×min(m,n)] ──── */
SEXP amatrix_arrayfire_qr_q_correct_bridge(SEXP A_r) {
  if (!isReal(A_r) || !isMatrix(A_r))
    error("input must be a real matrix");
#ifdef HAVE_ARRAYFIRE
  af_array A = 0, Q = 0, R = 0, tau = 0;
  A = amatrix_af_from_r(A_r);
  af_err err = af_qr(&Q, &R, &tau, A);
  af_release_array(A); af_release_array(R); af_release_array(tau);
  if (err != AF_SUCCESS) error("amatrix_arrayfire_qr_q_correct: af_qr failed (%d)", (int)err);
  SEXP result = PROTECT(amatrix_af_to_r(Q));
  af_release_array(Q);
  UNPROTECT(1);
  return result;
#else
  error("amatrix_arrayfire_qr_q_correct requires ArrayFire");
  return R_NilValue;
#endif
}

/* ── LAPACK dgebrd: bidiagonal reduction ────────────────────────────────────
 *
 * Reduces an m×n real matrix A to bidiagonal form B via orthogonal
 * transformations:  Q^T * A * P = B
 *   m >= n → B is n×n upper bidiagonal
 *   m <  n → B is m×m lower bidiagonal
 *
 * Returns a named list: list(a, d, e, tauq, taup)
 *   a    — m×n packed reflectors (overwrites input)
 *   d    — k = min(m,n) diagonal elements of B
 *   e    — k-1 off-diagonal elements of B
 *   tauq — k scalars for the left Householder reflectors (Q)
 *   taup — k scalars for the right Householder reflectors (P)
 * ─────────────────────────────────────────────────────────────────────────── */
SEXP amatrix_arrayfire_bdc_bidiag_bridge(SEXP A_r) {
  if (!isReal(A_r) || !isMatrix(A_r))
    error("bdc_bidiag: A must be a real matrix");
  SEXP dim = getAttrib(A_r, R_DimSymbol);
  int m = INTEGER(dim)[0], n = INTEGER(dim)[1];
  int k = (m < n) ? m : n;

  /* Working copy of A — dgebrd overwrites it in place */
  SEXP a_sexp = PROTECT(allocMatrix(REALSXP, m, n));
  memcpy(REAL(a_sexp), REAL(A_r), (size_t)m * n * sizeof(double));

  /* Output vectors */
  int e_len = (k > 1) ? k - 1 : 0;
  SEXP d_sexp    = PROTECT(allocVector(REALSXP, k));
  SEXP e_sexp    = PROTECT(allocVector(REALSXP, e_len));
  SEXP tauq_sexp = PROTECT(allocVector(REALSXP, k));
  SEXP taup_sexp = PROTECT(allocVector(REALSXP, k));

  /* Workspace query */
  int lwork = -1, info = 0;
  double work_query = 0.0;
  F77_CALL(dgebrd)(&m, &n, REAL(a_sexp), &m,
                   REAL(d_sexp), REAL(e_sexp),
                   REAL(tauq_sexp), REAL(taup_sexp),
                   &work_query, &lwork, &info);
  lwork = (work_query > 1.0) ? (int)work_query : (m + n) * 64;
  double *work = (double *) R_alloc((size_t)lwork, sizeof(double));

  /* dgebrd call */
  F77_CALL(dgebrd)(&m, &n, REAL(a_sexp), &m,
                   REAL(d_sexp), REAL(e_sexp),
                   REAL(tauq_sexp), REAL(taup_sexp),
                   work, &lwork, &info);
  if (info != 0)
    error("dgebrd failed with INFO = %d", info);

  /* Build named return list */
  SEXP result = PROTECT(allocVector(VECSXP, 5));
  SEXP names  = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names, 0, mkChar("a"));
  SET_STRING_ELT(names, 1, mkChar("d"));
  SET_STRING_ELT(names, 2, mkChar("e"));
  SET_STRING_ELT(names, 3, mkChar("tauq"));
  SET_STRING_ELT(names, 4, mkChar("taup"));
  setAttrib(result, R_NamesSymbol, names);
  SET_VECTOR_ELT(result, 0, a_sexp);
  SET_VECTOR_ELT(result, 1, d_sexp);
  SET_VECTOR_ELT(result, 2, e_sexp);
  SET_VECTOR_ELT(result, 3, tauq_sexp);
  SET_VECTOR_ELT(result, 4, taup_sexp);
  UNPROTECT(7); /* a, d, e, tauq, taup, result, names */
  return result;
}

/* ── LAPACK dorgbr: form Q or P^T from dgebrd output ───────────────────────
 *
 * Arguments:
 *   vect_r — "Q" to generate the left orthogonal factor Q,
 *             "P" to generate the right orthogonal factor P^T
 *   A_r    — the packed-reflector matrix returned by bdc_bidiag (m×n)
 *   tau_r  — tauq (for "Q") or taup (for "P") from bdc_bidiag
 *   M_r    — number of rows of the matrix to generate
 *   N_r    — number of columns of the matrix to generate
 *   K_r    — min(orig_m, orig_n) from the dgebrd call
 *
 * Returns: M×N real matrix (Q or truncated P^T).
 *
 * Standard parameter choices for an orig_m × orig_n input matrix with
 * k = min(orig_m, orig_n):
 *   dorgbr("Q", M=orig_m, N=k, K=k, ...)  → orig_m × k  Q
 *   dorgbr("P", M=k,      N=orig_n, K=k, ...) → k × orig_n P^T
 * ─────────────────────────────────────────────────────────────────────────── */
SEXP amatrix_arrayfire_bdc_orgbr_bridge(SEXP vect_r, SEXP A_r, SEXP tau_r,
                                         SEXP M_r, SEXP N_r, SEXP K_r) {
  if (TYPEOF(vect_r) != STRSXP || length(vect_r) < 1)
    error("bdc_orgbr: vect must be a length-1 character string");
  const char *vect = CHAR(STRING_ELT(vect_r, 0));
  if (vect[0] != 'Q' && vect[0] != 'P')
    error("bdc_orgbr: vect must be \"Q\" or \"P\"");
  if (!isReal(A_r) || !isMatrix(A_r))
    error("bdc_orgbr: A must be a real matrix");
  if (!isReal(tau_r))
    error("bdc_orgbr: tau must be a real vector");

  SEXP dim = getAttrib(A_r, R_DimSymbol);
  int lda   = INTEGER(dim)[0];   /* leading dimension = original m from dgebrd */
  int n_col = INTEGER(dim)[1];   /* original n from dgebrd */
  int M = asInteger(M_r);
  int N = asInteger(N_r);
  int K = asInteger(K_r);

  if (M < 0 || N < 0 || K < 0)
    error("bdc_orgbr: M, N, K must be non-negative");

  /* Working copy of the full m×n packed-reflector matrix.
   * dorgbr reads the embedded reflectors from A and overwrites the first
   * M×N block; we keep the full array so the leading dimension stays lda. */
  size_t full_sz = (size_t)lda * n_col;
  double *a_work = (double *) R_alloc(full_sz, sizeof(double));
  memcpy(a_work, REAL(A_r), full_sz * sizeof(double));

  double *tau = REAL(tau_r);

  /* Workspace query */
  int lwork = -1, info = 0;
  double work_query = 0.0;
  F77_CALL(dorgbr)(vect, &M, &N, &K, a_work, &lda, tau,
                   &work_query, &lwork, &info, (size_t)1);
  lwork = (work_query > 1.0) ? (int)work_query : (M + N + K) * 64;
  double *work = (double *) R_alloc((size_t)lwork, sizeof(double));

  /* dorgbr call — result overwrites the leading M×N block of a_work */
  F77_CALL(dorgbr)(vect, &M, &N, &K, a_work, &lda, tau,
                   work, &lwork, &info, (size_t)1);
  if (info != 0)
    error("dorgbr failed with INFO = %d", info);

  /* Extract the M×N result from the leading portion (a_work has stride lda) */
  SEXP out = PROTECT(allocMatrix(REALSXP, M, N));
  double *out_data = REAL(out);
  for (int j = 0; j < N; j++)
    for (int i = 0; i < M; i++)
      out_data[i + (size_t)j * M] = a_work[i + (size_t)j * lda];

  UNPROTECT(1);
  return out;
}

/* ── LAPACK dbdsdc: bidiagonal divide-and-conquer SVD ───────────────────────
 *
 * Computes the SVD of a real bidiagonal matrix B directly from its diagonal
 * and off-diagonal vectors, without requiring B to be materialised as a full
 * matrix.  This is 5-10× faster than calling dgesdd on the equivalent dense
 * bidiagonal matrix.
 *
 * Arguments:
 *   d_r    — real vector, length N: main diagonal of B (overwritten with
 *             singular values in decreasing order on output)
 *   e_r    — real vector, length N-1: off-diagonal of B
 *             (upper bidiagonal → superdiagonal; lower bidiagonal → subdiag)
 *   uplo_r — "U" (upper bidiagonal) or "L" (lower bidiagonal)
 *
 * Returns a named list: list(d, u, vt)
 *   d  — N singular values (decreasing)
 *   u  — N×N left  singular vector matrix (column per vector)
 *   vt — N×N right singular vector matrix as V^T (row per vector)
 * ─────────────────────────────────────────────────────────────────────────── */
SEXP amatrix_arrayfire_bdc_dbdsdc_bridge(SEXP d_r, SEXP e_r, SEXP uplo_r) {
  if (!isReal(d_r))
    error("bdc_dbdsdc: d must be a real vector");
  if (TYPEOF(uplo_r) != STRSXP || length(uplo_r) < 1)
    error("bdc_dbdsdc: uplo must be a length-1 character string");

  const char *uplo = CHAR(STRING_ELT(uplo_r, 0));
  int N = (int) length(d_r);

  /* Working copy of d (dbdsdc overwrites it with singular values) */
  SEXP d_out = PROTECT(allocVector(REALSXP, N));
  memcpy(REAL(d_out), REAL(d_r), (size_t)N * sizeof(double));

  /* Working copy of e (dbdsdc overwrites it) */
  int e_len = (N > 1) ? N - 1 : 0;
  double *e_work = (double *) R_alloc(e_len > 0 ? e_len : 1, sizeof(double));
  if (e_len > 0) memcpy(e_work, REAL(e_r), (size_t)e_len * sizeof(double));

  /* Output matrices: U (N×N) and VT (N×N) */
  SEXP U_out  = PROTECT(allocMatrix(REALSXP, N, N));
  SEXP VT_out = PROTECT(allocMatrix(REALSXP, N, N));

  /* Workspace: COMPQ='I' requires LWORK >= max(1, 3*N^2 + 4*N),
   * LIWORK >= 8*N (from LAPACK dbdsdc documentation). */
  int lwork  = (N > 0) ? (3 * N * N + 4 * N) : 1;
  int liwork = (N > 0) ? (8 * N) : 1;
  double *work  = (double *) R_alloc((size_t)lwork,  sizeof(double));
  int    *iwork = (int *)    R_alloc((size_t)liwork, sizeof(int));

  /* Dummy Q / IQ (only used when COMPQ='P') */
  double dum_q  = 0.0;
  int    dum_iq = 0;
  int    info   = 0;

  F77_CALL(dbdsdc)(uplo, "I", &N,
                   REAL(d_out), e_work,
                   REAL(U_out),  &N,
                   REAL(VT_out), &N,
                   &dum_q, &dum_iq, work, iwork, &info,
                   (size_t)1, (size_t)1);
  if (info != 0)
    error("dbdsdc failed with INFO = %d", info);

  /* Build return list */
  SEXP result = PROTECT(allocVector(VECSXP, 3));
  SEXP names  = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("d"));
  SET_STRING_ELT(names, 1, mkChar("u"));
  SET_STRING_ELT(names, 2, mkChar("vt"));
  setAttrib(result, R_NamesSymbol, names);
  SET_VECTOR_ELT(result, 0, d_out);
  SET_VECTOR_ELT(result, 1, U_out);
  SET_VECTOR_ELT(result, 2, VT_out);
  UNPROTECT(5); /* d_out, U_out, VT_out, result, names */
  return result;
}

/* ── Sparse×Dense matrix multiply (SpMM) ────────────────────────────────────
 *
 * Accepts a dgCMatrix (CSC: values/@x, col-ptrs/@p, row-idx/@i) plus a dense
 * RHS matrix B, then computes:
 *   trans_lhs=FALSE :  X %*% B
 *   trans_lhs=TRUE  :  t(X) %*% B
 * via ArrayFire af_sparse_matmul (sparse×dense, GPU-accelerated).
 *
 * Arguments
 *   values_r    REALSXP  — NNZ values      (dgCMatrix @x)
 *   p_r         INTSXP   — col pointers    (dgCMatrix @p, length ncol+1)
 *   i_r         INTSXP   — row indices     (dgCMatrix @i, length NNZ)
 *   dim_r       INTSXP   — c(nrow, ncol)   (dgCMatrix @Dim)
 *   B_r         REALSXP  — dense RHS matrix
 *   trans_lhs_r LGLSXP   — TRUE → compute t(X) %*% B
 */
SEXP amatrix_arrayfire_spmm_bridge(SEXP values_r, SEXP p_r, SEXP i_r,
                                    SEXP dim_r,    SEXP B_r, SEXP trans_lhs_r) {
  if (!isReal(values_r))
    error("spmm: values must be a real vector");
  if (TYPEOF(i_r) != INTSXP)
    error("spmm: row indices must be integer");
  if (TYPEOF(p_r) != INTSXP)
    error("spmm: col pointers must be integer");
  if (TYPEOF(dim_r) != INTSXP || length(dim_r) != 2)
    error("spmm: dim must be integer[2]");
  if (!isReal(B_r) || !isMatrix(B_r))
    error("spmm: B must be a real matrix");

#ifdef HAVE_ARRAYFIRE
  int    nrow      = INTEGER(dim_r)[0];
  int    ncol      = INTEGER(dim_r)[1];
  int    nnz       = (int)length(values_r);
  int    trans_lhs = asLogical(trans_lhs_r);

  /* float32 copy of values */
  float *fval = (float *) arrayfire_xmalloc((size_t)nnz * sizeof(float));
  const double *dval = REAL(values_r);
  for (int k = 0; k < nnz; k++) fval[k] = (float)dval[k];

  /* Build AF sparse array (CSC layout matches dgCMatrix exactly):
   *   rowIdx = row indices (NNZ elements)
   *   colIdx = col pointers (ncol+1 elements)                        */
  af_array sp_arr = 0;
  af_err err = af_create_sparse_array_from_ptr(
      &sp_arr,
      (dim_t)nrow, (dim_t)ncol, (dim_t)nnz,
      fval,
      INTEGER(i_r),   /* row indices */
      INTEGER(p_r),   /* col pointers */
      AF_STORAGE_CSC,
      f32,
      afHost);
  free(fval);
  if (err != AF_SUCCESS) {
    error("af_create_sparse_array_from_ptr failed (err=%d)", (int)err);
  }

  /* Dense RHS — column-major path (amatrix_af_from_r) */
  af_array B_af = amatrix_af_from_r(B_r);

  /* SpMM */
  af_array out = 0;
  af_mat_prop opt_lhs = trans_lhs ? AF_MAT_TRANS : AF_MAT_NONE;
  err = af_sparse_matmul(&out, sp_arr, B_af, opt_lhs, AF_MAT_NONE);
  af_release_array(sp_arr);
  af_release_array(B_af);
  if (err != AF_SUCCESS) {
    if (out) af_release_array(out);
    error("af_sparse_matmul failed (err=%d)", (int)err);
  }

  SEXP result = amatrix_af_to_r(out);
  af_release_array(out);
  return result;

#else
  /* ── No-ArrayFire fallback: plain sparse CSC × dense ── */
  SEXP B_dim  = getAttrib(B_r, R_DimSymbol);
  int  B_nrow = INTEGER(B_dim)[0];
  int  B_ncol = INTEGER(B_dim)[1];
  int  X_nrow = INTEGER(dim_r)[0];
  int  X_ncol = INTEGER(dim_r)[1];
  int  trans  = asLogical(trans_lhs_r);
  int  out_nrow = trans ? X_ncol : X_nrow;

  (void)B_nrow;  /* suppress unused warning */

  SEXP out_r = PROTECT(allocMatrix(REALSXP, out_nrow, B_ncol));
  double *res   = REAL(out_r);
  for (int k = 0; k < out_nrow * B_ncol; k++) res[k] = 0.0;

  const double *xdata = REAL(values_r);
  const double *bdata = REAL(B_r);
  const int    *xi    = INTEGER(i_r);
  const int    *xp    = INTEGER(p_r);

  if (!trans) {
    /* X %*% B : iterate CSC columns */
    for (int j = 0; j < X_ncol; j++) {
      for (int sp = xp[j]; sp < xp[j + 1]; sp++) {
        int    ri = xi[sp];
        double v  = xdata[sp];
        for (int cb = 0; cb < B_ncol; cb++)
          res[ri + (size_t)out_nrow * cb] += v * bdata[j + (size_t)X_ncol * cb];
      }
    }
  } else {
    /* t(X) %*% B : col j of X → row j of t(X) */
    for (int j = 0; j < X_ncol; j++) {
      for (int sp = xp[j]; sp < xp[j + 1]; sp++) {
        int    ri = xi[sp];
        double v  = xdata[sp];
        for (int cb = 0; cb < B_ncol; cb++)
          res[j + (size_t)out_nrow * cb] += v * bdata[ri + (size_t)X_nrow * cb];
      }
    }
  }

  UNPROTECT(1);
  return out_r;
#endif
}
