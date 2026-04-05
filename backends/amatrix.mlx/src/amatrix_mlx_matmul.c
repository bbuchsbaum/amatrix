#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MLXC
#include <mlx/c/mlx.h>
#include <mlx/c/linalg.h>

typedef struct {
  char* key;
  mlx_array array;
  bool in_use;
} amatrix_mlx_resident_entry;

static amatrix_mlx_resident_entry* amatrix_mlx_registry = NULL;
static size_t amatrix_mlx_registry_capacity = 0;

static void amatrix_mlx_free_array_if_needed(mlx_array arr);

static void amatrix_mlx_registry_init(void) {
  if (amatrix_mlx_registry != NULL) {
    return;
  }

  amatrix_mlx_registry_capacity = 128;
  amatrix_mlx_registry = (amatrix_mlx_resident_entry*) calloc(amatrix_mlx_registry_capacity, sizeof(amatrix_mlx_resident_entry));
  if (amatrix_mlx_registry == NULL) {
    error("failed to allocate mlx residency registry");
  }
}

static amatrix_mlx_resident_entry* amatrix_mlx_registry_find(const char* key) {
  amatrix_mlx_registry_init();

  for (size_t idx = 0; idx < amatrix_mlx_registry_capacity; ++idx) {
    if (amatrix_mlx_registry[idx].in_use && strcmp(amatrix_mlx_registry[idx].key, key) == 0) {
      return &amatrix_mlx_registry[idx];
    }
  }

  return NULL;
}

static amatrix_mlx_resident_entry* amatrix_mlx_registry_reserve(const char* key) {
  amatrix_mlx_resident_entry* existing = amatrix_mlx_registry_find(key);
  if (existing != NULL) {
    amatrix_mlx_free_array_if_needed(existing->array);
    existing->array = mlx_array_new();
    return existing;
  }

  amatrix_mlx_registry_init();
  for (size_t idx = 0; idx < amatrix_mlx_registry_capacity; ++idx) {
    if (!amatrix_mlx_registry[idx].in_use) {
      amatrix_mlx_registry[idx].in_use = true;
      amatrix_mlx_registry[idx].key = strdup(key);
      amatrix_mlx_registry[idx].array = mlx_array_new();
      if (amatrix_mlx_registry[idx].key == NULL) {
        error("failed to allocate mlx resident key");
      }
      return &amatrix_mlx_registry[idx];
    }
  }

  error("mlx residency registry is full");
  return NULL;
}

static void amatrix_mlx_registry_drop(const char* key) {
  amatrix_mlx_resident_entry* entry = amatrix_mlx_registry_find(key);
  if (entry == NULL) {
    return;
  }

  if (entry->key != NULL) {
    free(entry->key);
    entry->key = NULL;
  }
  amatrix_mlx_free_array_if_needed(entry->array);
  entry->array = mlx_array_new();
  entry->in_use = false;
}
#endif

static void copy_r_to_row_major_float(float* out, const double* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[i * ncol + j] = (float)in[i + nrow * j];
    }
  }
}

static void copy_row_major_float_to_r(double* out, const float* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[i + nrow * j] = (double)in[i * ncol + j];
    }
  }
}

static void copy_r_block_to_row_major_float(float* out, const double* in, int nrow, int ncol, int row_start, int block_nrow) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < block_nrow; ++i) {
      out[i * ncol + j] = (float)in[(row_start + i) + nrow * j];
    }
  }
}

static SEXP make_r_numeric_vector_from_float(const float* data, int n) {
  SEXP out = PROTECT(allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) {
    REAL(out)[i] = (double)data[i];
  }
  UNPROTECT(1);
  return out;
}

static SEXP amatrix_named_list2(const char* name1, SEXP value1, const char* name2, SEXP value2) {
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

static SEXP amatrix_named_list3(
    const char* name1,
    SEXP value1,
    const char* name2,
    SEXP value2,
    const char* name3,
    SEXP value3) {
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));

  SET_VECTOR_ELT(out, 0, value1);
  SET_VECTOR_ELT(out, 1, value2);
  SET_VECTOR_ELT(out, 2, value3);
  SET_STRING_ELT(names, 0, mkChar(name1));
  SET_STRING_ELT(names, 1, mkChar(name2));
  SET_STRING_ELT(names, 2, mkChar(name3));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

static SEXP amatrix_named_list5(
    const char* name1,
    SEXP value1,
    const char* name2,
    SEXP value2,
    const char* name3,
    SEXP value3,
    const char* name4,
    SEXP value4,
    const char* name5,
    SEXP value5) {
  SEXP out = PROTECT(allocVector(VECSXP, 5));
  SEXP names = PROTECT(allocVector(STRSXP, 5));

  SET_VECTOR_ELT(out, 0, value1);
  SET_VECTOR_ELT(out, 1, value2);
  SET_VECTOR_ELT(out, 2, value3);
  SET_VECTOR_ELT(out, 3, value4);
  SET_VECTOR_ELT(out, 4, value5);
  SET_STRING_ELT(names, 0, mkChar(name1));
  SET_STRING_ELT(names, 1, mkChar(name2));
  SET_STRING_ELT(names, 2, mkChar(name3));
  SET_STRING_ELT(names, 3, mkChar(name4));
  SET_STRING_ELT(names, 4, mkChar(name5));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

static SEXP amatrix_named_list6(
    const char* name1,
    SEXP value1,
    const char* name2,
    SEXP value2,
    const char* name3,
    SEXP value3,
    const char* name4,
    SEXP value4,
    const char* name5,
    SEXP value5,
    const char* name6,
    SEXP value6) {
  SEXP out = PROTECT(allocVector(VECSXP, 6));
  SEXP names = PROTECT(allocVector(STRSXP, 6));

  SET_VECTOR_ELT(out, 0, value1);
  SET_VECTOR_ELT(out, 1, value2);
  SET_VECTOR_ELT(out, 2, value3);
  SET_VECTOR_ELT(out, 3, value4);
  SET_VECTOR_ELT(out, 4, value5);
  SET_VECTOR_ELT(out, 5, value6);
  SET_STRING_ELT(names, 0, mkChar(name1));
  SET_STRING_ELT(names, 1, mkChar(name2));
  SET_STRING_ELT(names, 2, mkChar(name3));
  SET_STRING_ELT(names, 3, mkChar(name4));
  SET_STRING_ELT(names, 4, mkChar(name5));
  SET_STRING_ELT(names, 5, mkChar(name6));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

#ifdef HAVE_MLXC

static void amatrix_mlx_noop_error_handler(const char* msg, void* data) {
  (void)msg;
  (void)data;
}

static void amatrix_mlx_install_error_handler(void) {
  mlx_set_error_handler(amatrix_mlx_noop_error_handler, NULL, NULL);
}

static bool amatrix_mlx_gpu_stream_ok(mlx_stream* stream) {
  *stream = mlx_default_gpu_stream_new();
  return stream->ctx != NULL;
}

static void amatrix_mlx_free_array_if_needed(mlx_array arr) {
  if (arr.ctx != NULL) {
    mlx_array_free(arr);
  }
}

static SEXP amatrix_mlx_result_to_r_matrix(const mlx_array arr) {
  if (mlx_array_dtype(arr) != MLX_FLOAT32) {
    error("mlx result must be float32");
  }

  if (mlx_array_eval(arr) != 0) {
    error("mlx_array_eval failed");
  }

  const int* shape = mlx_array_shape(arr);
  int nrow = shape[0];
  int ncol = shape[1];
  const float* data = mlx_array_data_float32(arr);

  if (data == NULL) {
    error("mlx result data is unavailable");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  copy_row_major_float_to_r(REAL(out), data, nrow, ncol);
  UNPROTECT(1);
  return out;
}

static mlx_array amatrix_mlx_matrix_from_r(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  int shape[2] = {nrow, ncol};
  size_t size = (size_t)nrow * (size_t)ncol;
  float* buffer = (float*) R_alloc(size, sizeof(float));

  copy_r_to_row_major_float(buffer, REAL(x), nrow, ncol);
  return mlx_array_new_data(buffer, shape, 2, MLX_FLOAT32);
}

static mlx_array amatrix_mlx_array_from_r_value(SEXP x) {
  if (isMatrix(x)) {
    return amatrix_mlx_matrix_from_r(x);
  }

  if (isReal(x) && XLENGTH(x) == 1) {
    return mlx_array_new_float32((float)REAL(x)[0]);
  }

  if (isInteger(x) && XLENGTH(x) == 1) {
    return mlx_array_new_float32((float)INTEGER(x)[0]);
  }

  error("unsupported ewise operand");
}

static mlx_array amatrix_mlx_array_from_resident_key(SEXP key) {
  const char* key_str = CHAR(asChar(key));
  amatrix_mlx_resident_entry* entry = amatrix_mlx_registry_find(key_str);
  if (entry == NULL || !entry->in_use || entry->array.ctx == NULL) {
    error("unknown resident mlx key");
  }
  return entry->array;
}

static SEXP amatrix_mlx_matmul_real(SEXP x, SEXP y) {
  mlx_stream stream;
  mlx_array ax = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();

  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ax = amatrix_mlx_matrix_from_r(x);
  ay = amatrix_mlx_matrix_from_r(y);

  if (mlx_matmul(&out, ax, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ax);
  amatrix_mlx_free_array_if_needed(ay);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_crossprod_real(SEXP x, SEXP y) {
  mlx_stream stream;
  mlx_array ax = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array at = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();

  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ax = amatrix_mlx_matrix_from_r(x);
  ay = isNull(y) ? ax : amatrix_mlx_matrix_from_r(y);

  if (mlx_transpose(&at, ax, stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_transpose failed");
  }

  if (mlx_matmul(&out, at, ay, stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  if (!isNull(y)) {
    amatrix_mlx_free_array_if_needed(ay);
  }
  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ax);
  amatrix_mlx_free_array_if_needed(at);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_tcrossprod_real(SEXP x, SEXP y) {
  mlx_stream stream;
  mlx_array ax = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array bt = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();

  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ax = amatrix_mlx_matrix_from_r(x);
  ay = isNull(y) ? ax : amatrix_mlx_matrix_from_r(y);

  if (mlx_transpose(&bt, ay, stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_transpose failed");
  }

  if (mlx_matmul(&out, ax, bt, stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    if (!isNull(y)) {
      amatrix_mlx_free_array_if_needed(ay);
    }
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  if (!isNull(y)) {
    amatrix_mlx_free_array_if_needed(ay);
  }
  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ax);
  amatrix_mlx_free_array_if_needed(bt);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_solve_triangular_real(SEXP a, SEXP b, SEXP upper) {
  mlx_stream stream = {0};
  mlx_stream cpu_stream = {0};
  mlx_array aa = mlx_array_new();
  mlx_array bb = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;
  bool use_upper = asLogical(upper) == TRUE;
  bool used_gpu = false;

  amatrix_mlx_install_error_handler();

  cpu_stream = mlx_default_cpu_stream_new();
  if (amatrix_mlx_gpu_stream_ok(&stream)) {
    used_gpu = true;
  } else {
    stream = cpu_stream;
  }

  aa = amatrix_mlx_matrix_from_r(a);
  bb = amatrix_mlx_matrix_from_r(b);

  if (mlx_linalg_solve_triangular(&out, aa, bb, use_upper, stream) != 0) {
    if (used_gpu) {
      mlx_stream_free(stream);
      stream = cpu_stream;
      out = mlx_array_new();
      if (mlx_linalg_solve_triangular(&out, aa, bb, use_upper, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(aa);
        amatrix_mlx_free_array_if_needed(bb);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx_linalg_solve_triangular failed");
      }
      used_gpu = false;
    } else {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(aa);
      amatrix_mlx_free_array_if_needed(bb);
      amatrix_mlx_free_array_if_needed(out);
      error("mlx_linalg_solve_triangular failed");
    }
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(aa);
    amatrix_mlx_free_array_if_needed(bb);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  mlx_stream_free(stream);
  if (used_gpu && cpu_stream.ctx != NULL) {
    mlx_stream_free(cpu_stream);
  }
  amatrix_mlx_free_array_if_needed(aa);
  amatrix_mlx_free_array_if_needed(bb);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_qr_real(SEXP x, SEXP q_key) {
  mlx_stream stream = mlx_default_cpu_stream_new();
  mlx_array ax = mlx_array_new();
  mlx_array q = mlx_array_new();
  mlx_array r = mlx_array_new();
  SEXP r_r = R_NilValue;
  SEXP q_key_r = R_NilValue;
  SEXP result = R_NilValue;
  amatrix_mlx_resident_entry* entry = NULL;

  amatrix_mlx_install_error_handler();

  if (stream.ctx == NULL) {
    error("mlx CPU stream is unavailable");
  }

  ax = amatrix_mlx_matrix_from_r(x);

  if (mlx_linalg_qr(&q, &r, ax, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(q);
    amatrix_mlx_free_array_if_needed(r);
    error("mlx_linalg_qr failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(q);
    amatrix_mlx_free_array_if_needed(r);
    error("mlx_synchronize failed");
  }

  entry = amatrix_mlx_registry_reserve(CHAR(asChar(q_key)));
  entry->array = q;

  PROTECT(q_key_r = ScalarString(asChar(q_key)));
  PROTECT(r_r = amatrix_mlx_result_to_r_matrix(r));
  result = amatrix_named_list3("q", R_NilValue, "q_key", q_key_r, "r", r_r);
  UNPROTECT(2);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ax);
  amatrix_mlx_free_array_if_needed(r);
  return result;
}

static SEXP amatrix_mlx_qr_qty_key_real(SEXP q_key, SEXP y) {
  mlx_stream stream;
  mlx_array q = amatrix_mlx_array_from_resident_key(q_key);
  mlx_array qt = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ay = amatrix_mlx_matrix_from_r(y);

  if (mlx_transpose(&qt, q, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident qr transpose failed");
  }

  if (mlx_matmul(&out, qt, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident qr qty failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident qr qty synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(qt);
  amatrix_mlx_free_array_if_needed(ay);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_qr_qy_key_real(SEXP q_key, SEXP y) {
  mlx_stream stream;
  mlx_array q = amatrix_mlx_array_from_resident_key(q_key);
  mlx_array ay = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ay = amatrix_mlx_matrix_from_r(y);

  if (mlx_matmul(&out, q, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident qr qy failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident qr qy synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ay);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_qr_coef_key_real(SEXP q_key, SEXP r, SEXP y) {
  mlx_stream stream = {0};
  mlx_stream cpu_stream = {0};
  mlx_array q = amatrix_mlx_array_from_resident_key(q_key);
  mlx_array qt = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array rr = mlx_array_new();
  mlx_array qty = mlx_array_new();
  mlx_array coef = mlx_array_new();
  SEXP result = R_NilValue;
  bool used_gpu = false;

  amatrix_mlx_install_error_handler();
  cpu_stream = mlx_default_cpu_stream_new();
  if (amatrix_mlx_gpu_stream_ok(&stream)) {
    used_gpu = true;
  } else {
    stream = cpu_stream;
  }

  ay = amatrix_mlx_matrix_from_r(y);
  rr = amatrix_mlx_matrix_from_r(r);

  if (mlx_transpose(&qt, q, stream) != 0) {
    mlx_stream_free(stream);
    if (used_gpu && cpu_stream.ctx != NULL) {
      mlx_stream_free(cpu_stream);
    }
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx resident qr coef transpose failed");
  }

  if (mlx_matmul(&qty, qt, ay, stream) != 0) {
    mlx_stream_free(stream);
    if (used_gpu && cpu_stream.ctx != NULL) {
      mlx_stream_free(cpu_stream);
    }
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx resident qr coef qty failed");
  }

  if (mlx_linalg_solve_triangular(&coef, rr, qty, true, stream) != 0) {
    if (used_gpu) {
      mlx_stream_free(stream);
      stream = cpu_stream;
      coef = mlx_array_new();
      if (mlx_linalg_solve_triangular(&coef, rr, qty, true, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(rr);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx resident qr coef solve failed");
      }
      used_gpu = false;
    } else {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx resident qr coef solve failed");
    }
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx resident qr coef synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(coef);

  mlx_stream_free(stream);
  if (used_gpu && cpu_stream.ctx != NULL) {
    mlx_stream_free(cpu_stream);
  }
  amatrix_mlx_free_array_if_needed(qt);
  amatrix_mlx_free_array_if_needed(ay);
  amatrix_mlx_free_array_if_needed(rr);
  amatrix_mlx_free_array_if_needed(qty);
  amatrix_mlx_free_array_if_needed(coef);
  return result;
}

static SEXP amatrix_mlx_tsqr_coef_key_real(SEXP q_keys, SEXP block_rows, SEXP top_q_key, SEXP r, SEXP y) {
  const int nblocks = (int) XLENGTH(q_keys);
  SEXP y_dim = getAttrib(y, R_DimSymbol);
  const double* y_data = REAL(y);
  const int nrow = INTEGER(y_dim)[0];
  const int nrhs = INTEGER(y_dim)[1];
  SEXP r_dim = getAttrib(r, R_DimSymbol);
  const int p = INTEGER(r_dim)[1];
  const int head_nrow = nblocks * p;
  mlx_stream stream = {0};
  mlx_stream cpu_stream = {0};
  bool used_gpu = false;
  float* block_buf = NULL;
  float* head_buf = NULL;
  mlx_array q_top = mlx_array_new();
  mlx_array qt = mlx_array_new();
  mlx_array ay = mlx_array_new();
  mlx_array rr = mlx_array_new();
  mlx_array qty = mlx_array_new();
  mlx_array coef = mlx_array_new();
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();
  cpu_stream = mlx_default_cpu_stream_new();
  if (amatrix_mlx_gpu_stream_ok(&stream)) {
    used_gpu = true;
  } else {
    stream = cpu_stream;
  }

  block_buf = (float*) R_alloc((size_t)nrow * (size_t)nrhs, sizeof(float));
  head_buf = (float*) R_alloc((size_t)head_nrow * (size_t)nrhs, sizeof(float));

  for (int block_idx = 0, row_start = 0; block_idx < nblocks; ++block_idx) {
    const int block_nrow = INTEGER(block_rows)[block_idx];
    const int block_shape[2] = {block_nrow, nrhs};
    mlx_array q_block = mlx_array_new();
    mlx_array q_block_t = mlx_array_new();
    mlx_array y_block = mlx_array_new();
    mlx_array qty_block = mlx_array_new();
    const float* qty_data = NULL;

    q_block = amatrix_mlx_array_from_resident_key(STRING_ELT(q_keys, block_idx));
    copy_r_block_to_row_major_float(block_buf, y_data, nrow, nrhs, row_start, block_nrow);
    y_block = mlx_array_new_data(block_buf, block_shape, 2, MLX_FLOAT32);

    if (mlx_transpose(&q_block_t, q_block, stream) != 0) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(q_block_t);
      amatrix_mlx_free_array_if_needed(y_block);
      amatrix_mlx_free_array_if_needed(qty_block);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef local transpose failed");
    }

    if (mlx_matmul(&qty_block, q_block_t, y_block, stream) != 0) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(q_block_t);
      amatrix_mlx_free_array_if_needed(y_block);
      amatrix_mlx_free_array_if_needed(qty_block);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef local qty failed");
    }

    if (mlx_array_eval(qty_block) != 0) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(q_block_t);
      amatrix_mlx_free_array_if_needed(y_block);
      amatrix_mlx_free_array_if_needed(qty_block);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef local eval failed");
    }

    qty_data = mlx_array_data_float32(qty_block);
    if (qty_data == NULL) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(q_block_t);
      amatrix_mlx_free_array_if_needed(y_block);
      amatrix_mlx_free_array_if_needed(qty_block);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef local data is unavailable");
    }

    memcpy(
      head_buf + ((size_t)block_idx * (size_t)p * (size_t)nrhs),
      qty_data,
      (size_t)p * (size_t)nrhs * sizeof(float)
    );

    amatrix_mlx_free_array_if_needed(q_block_t);
    amatrix_mlx_free_array_if_needed(y_block);
    amatrix_mlx_free_array_if_needed(qty_block);
    row_start += block_nrow;
  }

  q_top = amatrix_mlx_array_from_resident_key(top_q_key);
  {
    const int head_shape[2] = {head_nrow, nrhs};
    ay = mlx_array_new_data(head_buf, head_shape, 2, MLX_FLOAT32);
  }
  rr = amatrix_mlx_matrix_from_r(r);

  if (mlx_transpose(&qt, q_top, stream) != 0) {
    mlx_stream_free(stream);
    if (used_gpu && cpu_stream.ctx != NULL) {
      mlx_stream_free(cpu_stream);
    }
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx tsqr coef top transpose failed");
  }

  if (mlx_matmul(&qty, qt, ay, stream) != 0) {
    mlx_stream_free(stream);
    if (used_gpu && cpu_stream.ctx != NULL) {
      mlx_stream_free(cpu_stream);
    }
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx tsqr coef top qty failed");
  }

  if (mlx_linalg_solve_triangular(&coef, rr, qty, true, stream) != 0) {
    if (used_gpu) {
      mlx_stream_free(stream);
      stream = cpu_stream;
      coef = mlx_array_new();
      if (mlx_linalg_solve_triangular(&coef, rr, qty, true, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(rr);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef solve failed");
      }
      used_gpu = false;
    } else {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(rr);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef solve failed");
    }
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(rr);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    error("mlx tsqr coef synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(coef);

  mlx_stream_free(stream);
  if (used_gpu && cpu_stream.ctx != NULL) {
    mlx_stream_free(cpu_stream);
  }
  amatrix_mlx_free_array_if_needed(qt);
  amatrix_mlx_free_array_if_needed(ay);
  amatrix_mlx_free_array_if_needed(rr);
  amatrix_mlx_free_array_if_needed(qty);
  amatrix_mlx_free_array_if_needed(coef);
  return result;
}

static int amatrix_mlx_rank_from_upper_triangular_array(mlx_array r) {
  const int* shape = mlx_array_shape(r);
  const int nrow = shape[0];
  const int ncol = shape[1];
  const int diag_len = nrow < ncol ? nrow : ncol;
  const float* data = mlx_array_data_float32(r);
  float max_diag = 0.0f;
  int rank = 0;

  if (data == NULL || diag_len <= 0) {
    return 0;
  }

  for (int i = 0; i < diag_len; ++i) {
    const float v = data[i * ncol + i];
    const float a = v < 0.0f ? -v : v;
    if (a > max_diag) {
      max_diag = a;
    }
  }

  {
    const double tol = ((double) (nrow > ncol ? nrow : ncol)) * 2.2204460492503131e-16 * (double) max_diag;
    for (int i = 0; i < diag_len; ++i) {
      const float v = data[i * ncol + i];
      const double a = (double) (v < 0.0f ? -v : v);
      if (a > tol) {
        rank += 1;
      }
    }
  }

  return rank;
}

static SEXP amatrix_mlx_tsqr_build_real(SEXP x, SEXP block_rows, SEXP q_keys, SEXP top_q_key, SEXP top_r_key, SEXP r_stack_key) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  const double* x_data = REAL(x);
  const int nrow = INTEGER(dim)[0];
  const int ncol = INTEGER(dim)[1];
  const int block_rows_val = INTEGER(block_rows)[0];
  const int nblocks = (nrow + block_rows_val - 1) / block_rows_val;
  const int top_nrow = nblocks * ncol;
  mlx_stream stream = mlx_default_cpu_stream_new();
  float* block_buf = NULL;
  float* r_stack_buf = NULL;
  SEXP result = R_NilValue;
  SEXP block_rows_r = R_NilValue;
  SEXP q_keys_r = R_NilValue;
  SEXP top_q_key_r = R_NilValue;
  SEXP top_r_key_r = R_NilValue;
  SEXP r_stack_key_r = R_NilValue;
  SEXP rank_r = R_NilValue;

  amatrix_mlx_install_error_handler();

  if (stream.ctx == NULL) {
    error("mlx CPU stream is unavailable");
  }

  block_buf = (float*) R_alloc((size_t)block_rows_val * (size_t)ncol, sizeof(float));
  r_stack_buf = (float*) R_alloc((size_t)top_nrow * (size_t)ncol, sizeof(float));

  PROTECT(block_rows_r = allocVector(INTSXP, nblocks));
  PROTECT(q_keys_r = duplicate(q_keys));

  for (int block_idx = 0; block_idx < nblocks; ++block_idx) {
    const int row_start = block_idx * block_rows_val;
    const int block_nrow = (row_start + block_rows_val <= nrow) ? block_rows_val : (nrow - row_start);
    const int block_shape[2] = {block_nrow, ncol};
    mlx_array ax = mlx_array_new();
    mlx_array q = mlx_array_new();
    mlx_array r = mlx_array_new();
    const float* r_data = NULL;
    amatrix_mlx_resident_entry* entry = NULL;

    INTEGER(block_rows_r)[block_idx] = block_nrow;
    copy_r_block_to_row_major_float(block_buf, x_data, nrow, ncol, row_start, block_nrow);

    ax = mlx_array_new_data(block_buf, block_shape, 2, MLX_FLOAT32);

    if (mlx_linalg_qr(&q, &r, ax, stream) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_linalg_qr failed in tsqr leaf");
    }

    if (mlx_synchronize(stream) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_synchronize failed in tsqr leaf");
    }

    if (mlx_array_eval(r) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_array_eval failed in tsqr leaf");
    }

    r_data = mlx_array_data_float32(r);
    if (r_data == NULL) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx tsqr leaf R data is unavailable");
    }

    memcpy(
      r_stack_buf + ((size_t)block_idx * (size_t)ncol * (size_t)ncol),
      r_data,
      (size_t)ncol * (size_t)ncol * sizeof(float)
    );

    entry = amatrix_mlx_registry_reserve(CHAR(STRING_ELT(q_keys_r, block_idx)));
    entry->array = q;

    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(r);
  }

  {
    const int top_shape[2] = {top_nrow, ncol};
    mlx_array ax = mlx_array_new();
    mlx_array q = mlx_array_new();
    mlx_array r = mlx_array_new();
    amatrix_mlx_resident_entry* entry = NULL;

    ax = mlx_array_new_data(r_stack_buf, top_shape, 2, MLX_FLOAT32);

    if (mlx_linalg_qr(&q, &r, ax, stream) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_linalg_qr failed in tsqr top reduction");
    }

    if (mlx_synchronize(stream) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_synchronize failed in tsqr top reduction");
    }

    if (mlx_array_eval(r) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(ax);
      amatrix_mlx_free_array_if_needed(q);
      amatrix_mlx_free_array_if_needed(r);
      UNPROTECT(2);
      error("mlx_array_eval failed in tsqr top reduction");
    }

    entry = amatrix_mlx_registry_reserve(CHAR(asChar(top_q_key)));
    entry->array = q;

    entry = amatrix_mlx_registry_reserve(CHAR(asChar(r_stack_key)));
    entry->array = ax;
    entry = amatrix_mlx_registry_reserve(CHAR(asChar(top_r_key)));
    entry->array = r;

    PROTECT(top_q_key_r = ScalarString(asChar(top_q_key)));
    PROTECT(top_r_key_r = ScalarString(asChar(top_r_key)));
    PROTECT(r_stack_key_r = ScalarString(asChar(r_stack_key)));
    PROTECT(rank_r = ScalarInteger(amatrix_mlx_rank_from_upper_triangular_array(r)));
    result = amatrix_named_list6(
      "block_rows", block_rows_r,
      "block_q_keys", q_keys_r,
      "r_stack_key", r_stack_key_r,
      "top_q_key", top_q_key_r,
      "top_r_key", top_r_key_r,
      "rank", rank_r
    );
    UNPROTECT(4);

    mlx_stream_free(stream);
  }

  UNPROTECT(2);
  return result;
}

static SEXP amatrix_mlx_ewise_real(SEXP lhs, SEXP rhs, SEXP op) {
  mlx_stream stream;
  mlx_array a = mlx_array_new();
  mlx_array b = mlx_array_new();
  mlx_array out = mlx_array_new();
  const char* op_name = CHAR(asChar(op));
  SEXP result = R_NilValue;
  bool rhs_was_null = isNull(rhs);

  amatrix_mlx_install_error_handler();

  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  a = amatrix_mlx_array_from_r_value(lhs);

  if (rhs_was_null) {
    if (strcmp(op_name, "-") == 0) {
      if (mlx_negative(&out, a, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(a);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx_negative failed");
      }
    } else {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(a);
      error("unsupported unary op for mlx bridge");
    }
  } else {
    b = amatrix_mlx_array_from_r_value(rhs);
    if (strcmp(op_name, "+") == 0) {
      if (mlx_add(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(a);
        amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx_add failed");
      }
    } else if (strcmp(op_name, "-") == 0) {
      if (mlx_subtract(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(a);
        amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx_subtract failed");
      }
    } else if (strcmp(op_name, "*") == 0) {
      if (mlx_multiply(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(a);
        amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx_multiply failed");
      }
    } else {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(a);
      amatrix_mlx_free_array_if_needed(b);
      error("unsupported binary op for mlx bridge");
    }
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(a);
    amatrix_mlx_free_array_if_needed(b);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  result = amatrix_mlx_result_to_r_matrix(out);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(a);
  amatrix_mlx_free_array_if_needed(b);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_sum_axis_real(SEXP x, SEXP axis) {
  mlx_stream stream;
  mlx_array ax = mlx_array_new();
  mlx_array out = mlx_array_new();
  int axis_val = INTEGER(axis)[0];
  SEXP result = R_NilValue;

  amatrix_mlx_install_error_handler();

  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  ax = amatrix_mlx_matrix_from_r(x);

  if (mlx_sum_axis(&out, ax, axis_val, false, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_sum_axis failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_synchronize failed");
  }

  if (mlx_array_dtype(out) != MLX_FLOAT32) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx sum result must be float32");
  }

  if (mlx_array_eval(out) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(ax);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx_array_eval failed");
  }

  result = make_r_numeric_vector_from_float(mlx_array_data_float32(out), (int) mlx_array_size(out));

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(ax);
  amatrix_mlx_free_array_if_needed(out);
  return result;
}

static SEXP amatrix_mlx_resident_store_real(SEXP key, SEXP x) {
  amatrix_mlx_resident_entry* entry = amatrix_mlx_registry_reserve(CHAR(asChar(key)));
  entry->array = amatrix_mlx_matrix_from_r(x);
  return ScalarLogical(1);
}

static SEXP amatrix_mlx_resident_materialize_real(SEXP key) {
  mlx_array arr = amatrix_mlx_array_from_resident_key(key);
  return amatrix_mlx_result_to_r_matrix(arr);
}

static SEXP amatrix_mlx_matmul_resident_real(SEXP x_key, SEXP y_key, SEXP out_key) {
  mlx_stream stream;
  mlx_array ax = amatrix_mlx_array_from_resident_key(x_key);
  mlx_array ay = amatrix_mlx_array_from_resident_key(y_key);
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;
  amatrix_mlx_resident_entry* entry = NULL;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  if (mlx_matmul(&out, ax, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident matmul synchronize failed");
  }

  entry = amatrix_mlx_registry_reserve(CHAR(asChar(out_key)));
  entry->array = out;
  result = amatrix_mlx_result_to_r_matrix(entry->array);
  mlx_stream_free(stream);
  return result;
}

static SEXP amatrix_mlx_crossprod_resident_real(SEXP x_key, SEXP y_key, SEXP out_key) {
  mlx_stream stream;
  mlx_array ax = amatrix_mlx_array_from_resident_key(x_key);
  mlx_array ay = isNull(y_key) ? ax : amatrix_mlx_array_from_resident_key(y_key);
  mlx_array at = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;
  amatrix_mlx_resident_entry* entry = NULL;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  if (mlx_transpose(&at, ax, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident crossprod transpose failed");
  }

  if (mlx_matmul(&out, at, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident crossprod matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(at);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident crossprod synchronize failed");
  }

  entry = amatrix_mlx_registry_reserve(CHAR(asChar(out_key)));
  entry->array = out;
  result = amatrix_mlx_result_to_r_matrix(entry->array);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(at);
  return result;
}

static SEXP amatrix_mlx_tcrossprod_resident_real(SEXP x_key, SEXP y_key, SEXP out_key) {
  mlx_stream stream;
  mlx_array ax = amatrix_mlx_array_from_resident_key(x_key);
  mlx_array ay = isNull(y_key) ? ax : amatrix_mlx_array_from_resident_key(y_key);
  mlx_array bt = mlx_array_new();
  mlx_array out = mlx_array_new();
  SEXP result = R_NilValue;
  amatrix_mlx_resident_entry* entry = NULL;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  if (mlx_transpose(&bt, ay, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident tcrossprod transpose failed");
  }

  if (mlx_matmul(&out, ax, bt, stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident tcrossprod matmul failed");
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    amatrix_mlx_free_array_if_needed(bt);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident tcrossprod synchronize failed");
  }

  entry = amatrix_mlx_registry_reserve(CHAR(asChar(out_key)));
  entry->array = out;
  result = amatrix_mlx_result_to_r_matrix(entry->array);

  mlx_stream_free(stream);
  amatrix_mlx_free_array_if_needed(bt);
  return result;
}

static SEXP amatrix_mlx_ewise_resident_real(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key) {
  mlx_stream stream;
  mlx_array a = amatrix_mlx_array_from_resident_key(lhs_key);
  mlx_array b = mlx_array_new();
  mlx_array out = mlx_array_new();
  const char* op_name = CHAR(asChar(op));
  SEXP result = R_NilValue;
  bool rhs_was_null = isNull(rhs);
  bool owns_rhs = false;
  amatrix_mlx_resident_entry* entry = NULL;

  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    error("mlx GPU stream is unavailable");
  }

  if (rhs_was_null) {
    if (strcmp(op_name, "-") == 0) {
      if (mlx_negative(&out, a, stream) != 0) {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx resident negative failed");
      }
    } else {
      mlx_stream_free(stream);
      error("unsupported unary op for resident mlx bridge");
    }
  } else {
    if (isString(rhs) && XLENGTH(rhs) == 1) {
      b = amatrix_mlx_array_from_resident_key(rhs);
    } else {
      b = amatrix_mlx_array_from_r_value(rhs);
      owns_rhs = true;
    }

    if (strcmp(op_name, "+") == 0) {
      if (mlx_add(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        if (owns_rhs) amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx resident add failed");
      }
    } else if (strcmp(op_name, "-") == 0) {
      if (mlx_subtract(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        if (owns_rhs) amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx resident subtract failed");
      }
    } else if (strcmp(op_name, "*") == 0) {
      if (mlx_multiply(&out, a, b, stream) != 0) {
        mlx_stream_free(stream);
        if (owns_rhs) amatrix_mlx_free_array_if_needed(b);
        amatrix_mlx_free_array_if_needed(out);
        error("mlx resident multiply failed");
      }
    } else {
      mlx_stream_free(stream);
      if (owns_rhs) amatrix_mlx_free_array_if_needed(b);
      error("unsupported binary op for resident mlx bridge");
    }
  }

  if (mlx_synchronize(stream) != 0) {
    mlx_stream_free(stream);
    if (owns_rhs) amatrix_mlx_free_array_if_needed(b);
    amatrix_mlx_free_array_if_needed(out);
    error("mlx resident synchronize failed");
  }

  entry = amatrix_mlx_registry_reserve(CHAR(asChar(out_key)));
  entry->array = out;
  result = amatrix_mlx_result_to_r_matrix(entry->array);

  if (owns_rhs) {
    amatrix_mlx_free_array_if_needed(b);
  }
  mlx_stream_free(stream);
  return result;
}

#endif

SEXP amatrix_mlx_native_available_bridge(void) {
#ifdef HAVE_MLXC
  mlx_stream stream;
  amatrix_mlx_install_error_handler();
  if (!amatrix_mlx_gpu_stream_ok(&stream)) {
    return ScalarLogical(0);
  }
  mlx_stream_free(stream);
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_mlx_bridge_info_bridge(void) {
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, mkChar("compiled"));
  SET_STRING_ELT(names, 1, mkChar("native"));
  SET_STRING_ELT(names, 2, mkChar("engine"));

  SET_VECTOR_ELT(out, 0, ScalarLogical(1));
#ifdef HAVE_MLXC
  SET_VECTOR_ELT(out, 1, ScalarLogical(1));
  SET_VECTOR_ELT(out, 2, mkString("mlx-c"));
#else
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
  SET_VECTOR_ELT(out, 2, mkString("mock-c-bridge"));
#endif

  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP amatrix_mlx_resident_has_bridge(SEXP key) {
#ifdef HAVE_MLXC
  return ScalarLogical(amatrix_mlx_registry_find(CHAR(asChar(key))) != NULL);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_mlx_resident_store_bridge(SEXP key, SEXP x) {
#ifdef HAVE_MLXC
  return amatrix_mlx_resident_store_real(key, x);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_resident_drop_bridge(SEXP key) {
#ifdef HAVE_MLXC
  amatrix_mlx_registry_drop(CHAR(asChar(key)));
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_mlx_resident_materialize_bridge(SEXP key) {
#ifdef HAVE_MLXC
  return amatrix_mlx_resident_materialize_real(key);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_MLXC
  return amatrix_mlx_matmul_resident_real(x_key, y_key, out_key);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_MLXC
  return amatrix_mlx_crossprod_resident_real(x_key, y_key, out_key);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_MLXC
  return amatrix_mlx_tcrossprod_resident_real(x_key, y_key, out_key);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key) {
#ifdef HAVE_MLXC
  return amatrix_mlx_ewise_resident_real(lhs_key, rhs, op, out_key);
#else
  error("mlx residency requires mlx-c");
#endif
}

SEXP amatrix_mlx_matmul_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_matmul_real(x, y);
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

SEXP amatrix_mlx_crossprod_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

  int has_y = !isNull(y);
  if (has_y && (!isReal(y) || !isMatrix(y))) {
    error("y must be NULL or a numeric matrix");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_crossprod_real(x, y);
#else

  SEXP x_dim = getAttrib(x, R_DimSymbol);
  int x_nrow = INTEGER(x_dim)[0];
  int x_ncol = INTEGER(x_dim)[1];

  int y_nrow = x_nrow;
  int y_ncol = x_ncol;
  SEXP y_dim = R_NilValue;

  if (has_y) {
    y_dim = getAttrib(y, R_DimSymbol);
    y_nrow = INTEGER(y_dim)[0];
    y_ncol = INTEGER(y_dim)[1];

    if (x_nrow != y_nrow) {
      error("non-conformable arguments");
    }
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, x_ncol, y_ncol));
  double *x_ptr = REAL(x);
  double *y_ptr = has_y ? REAL(y) : REAL(x);
  double *out_ptr = REAL(out);

  for (int j = 0; j < y_ncol; ++j) {
    for (int i = 0; i < x_ncol; ++i) {
      double acc = 0.0;
      for (int k = 0; k < x_nrow; ++k) {
        acc += x_ptr[k + x_nrow * i] * y_ptr[k + y_nrow * j];
      }
      out_ptr[i + x_ncol * j] = acc;
    }
  }

  UNPROTECT(1);
  return out;
#endif
}

SEXP amatrix_mlx_tcrossprod_bridge(SEXP x, SEXP y) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

  int has_y = !isNull(y);
  if (has_y && (!isReal(y) || !isMatrix(y))) {
    error("y must be NULL or a numeric matrix");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_tcrossprod_real(x, y);
#else

  SEXP x_dim = getAttrib(x, R_DimSymbol);
  int x_nrow = INTEGER(x_dim)[0];
  int x_ncol = INTEGER(x_dim)[1];

  int y_nrow = x_nrow;
  int y_ncol = x_ncol;
  SEXP y_dim = R_NilValue;

  if (has_y) {
    y_dim = getAttrib(y, R_DimSymbol);
    y_nrow = INTEGER(y_dim)[0];
    y_ncol = INTEGER(y_dim)[1];

    if (x_ncol != y_ncol) {
      error("non-conformable arguments");
    }
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, x_nrow, y_nrow));
  double *x_ptr = REAL(x);
  double *y_ptr = has_y ? REAL(y) : REAL(x);
  double *out_ptr = REAL(out);

  for (int j = 0; j < y_nrow; ++j) {
    for (int i = 0; i < x_nrow; ++i) {
      double acc = 0.0;
      for (int k = 0; k < x_ncol; ++k) {
        acc += x_ptr[i + x_nrow * k] * y_ptr[j + y_nrow * k];
      }
      out_ptr[i + x_nrow * j] = acc;
    }
  }

  UNPROTECT(1);
  return out;
#endif
}

SEXP amatrix_mlx_solve_triangular_bridge(SEXP a, SEXP b, SEXP upper) {
  if (!isReal(a) || !isMatrix(a)) {
    error("a must be a numeric matrix");
  }

  if (!isReal(b) || !isMatrix(b)) {
    error("b must be a numeric matrix");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_solve_triangular_real(a, b, upper);
#else
  error("mlx solve_triangular bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op) {
  if ((!isReal(lhs) && !isInteger(lhs)) || (!isMatrix(lhs) && XLENGTH(lhs) != 1)) {
    error("lhs must be a numeric matrix or scalar");
  }

  if (!isNull(rhs) && ((!isReal(rhs) && !isInteger(rhs)) || (!isMatrix(rhs) && XLENGTH(rhs) != 1))) {
    error("rhs must be NULL, a numeric matrix, or a scalar");
  }

  if (!isString(op) || XLENGTH(op) != 1) {
    error("op must be a single string");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_ewise_real(lhs, rhs, op);
#else
  error("mlx ewise bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_sum_axis_bridge(SEXP x, SEXP axis) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

  if (!isInteger(axis) || XLENGTH(axis) != 1) {
    error("axis must be a single integer");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_sum_axis_real(x, axis);
#else
  error("mlx sum bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_qr_bridge(SEXP x, SEXP q_key) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }

  if (!isString(q_key) || XLENGTH(q_key) != 1) {
    error("q_key must be a single string");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_qr_real(x, q_key);
#else
  error("mlx qr bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_qr_qty_key_bridge(SEXP q_key, SEXP y) {
  if (!isString(q_key) || XLENGTH(q_key) != 1) {
    error("q_key must be a single string");
  }
  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }
#ifdef HAVE_MLXC
  return amatrix_mlx_qr_qty_key_real(q_key, y);
#else
  error("mlx qr qty key bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_qr_qy_key_bridge(SEXP q_key, SEXP y) {
  if (!isString(q_key) || XLENGTH(q_key) != 1) {
    error("q_key must be a single string");
  }
  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }
#ifdef HAVE_MLXC
  return amatrix_mlx_qr_qy_key_real(q_key, y);
#else
  error("mlx qr qy key bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_qr_coef_key_bridge(SEXP q_key, SEXP r, SEXP y) {
  if (!isString(q_key) || XLENGTH(q_key) != 1) {
    error("q_key must be a single string");
  }
  if (!isReal(r) || !isMatrix(r)) {
    error("r must be a numeric matrix");
  }
  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }
#ifdef HAVE_MLXC
  return amatrix_mlx_qr_coef_key_real(q_key, r, y);
#else
  error("mlx qr coef key bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_tsqr_coef_key_bridge(SEXP q_keys, SEXP block_rows, SEXP top_q_key, SEXP r, SEXP y) {
  if (!isString(q_keys)) {
    error("q_keys must be a character vector");
  }
  if (!isInteger(block_rows) || XLENGTH(block_rows) != XLENGTH(q_keys)) {
    error("block_rows must be an integer vector matching q_keys");
  }
  if (!isString(top_q_key) || XLENGTH(top_q_key) != 1) {
    error("top_q_key must be a single string");
  }
  if (!isReal(r) || !isMatrix(r)) {
    error("r must be a numeric matrix");
  }
  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_tsqr_coef_key_real(q_keys, block_rows, top_q_key, r, y);
#else
  error("mlx tsqr coef key bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_tsqr_coef_resident_bridge(SEXP q_keys, SEXP block_rows, SEXP top_q_key, SEXP top_r_key, SEXP y) {
  if (!isString(q_keys)) {
    error("q_keys must be a character vector");
  }
  if (!isInteger(block_rows) || XLENGTH(block_rows) != XLENGTH(q_keys)) {
    error("block_rows must be an integer vector matching q_keys");
  }
  if (!isString(top_q_key) || XLENGTH(top_q_key) != 1) {
    error("top_q_key must be a single string");
  }
  if (!isString(top_r_key) || XLENGTH(top_r_key) != 1) {
    error("top_r_key must be a single string");
  }
  if (!isReal(y) || !isMatrix(y)) {
    error("y must be a numeric matrix");
  }

#ifdef HAVE_MLXC
  {
    const int nblocks = (int) XLENGTH(q_keys);
    SEXP y_dim = getAttrib(y, R_DimSymbol);
    const double* y_data = REAL(y);
    const int nrow = INTEGER(y_dim)[0];
    const int nrhs = INTEGER(y_dim)[1];
    mlx_array rr_ref = amatrix_mlx_array_from_resident_key(top_r_key);
    const int* r_shape = mlx_array_shape(rr_ref);
    const int p = r_shape[1];
    const int head_nrow = nblocks * p;
    mlx_stream stream = {0};
    mlx_stream cpu_stream = {0};
    bool used_gpu = false;
    float* block_buf = NULL;
    float* head_buf = NULL;
    mlx_array q_top = mlx_array_new();
    mlx_array qt = mlx_array_new();
    mlx_array ay = mlx_array_new();
    mlx_array qty = mlx_array_new();
    mlx_array coef = mlx_array_new();
    SEXP result = R_NilValue;

    amatrix_mlx_install_error_handler();
    cpu_stream = mlx_default_cpu_stream_new();
    if (amatrix_mlx_gpu_stream_ok(&stream)) {
      used_gpu = true;
    } else {
      stream = cpu_stream;
    }

    block_buf = (float*) R_alloc((size_t)nrow * (size_t)nrhs, sizeof(float));
    head_buf = (float*) R_alloc((size_t)head_nrow * (size_t)nrhs, sizeof(float));

    for (int block_idx = 0, row_start = 0; block_idx < nblocks; ++block_idx) {
      const int block_nrow = INTEGER(block_rows)[block_idx];
      const int block_shape[2] = {block_nrow, nrhs};
      mlx_array q_block = mlx_array_new();
      mlx_array q_block_t = mlx_array_new();
      mlx_array y_block = mlx_array_new();
      mlx_array qty_block = mlx_array_new();
      const float* qty_data = NULL;

      q_block = amatrix_mlx_array_from_resident_key(STRING_ELT(q_keys, block_idx));
      copy_r_block_to_row_major_float(block_buf, y_data, nrow, nrhs, row_start, block_nrow);
      y_block = mlx_array_new_data(block_buf, block_shape, 2, MLX_FLOAT32);

      if (mlx_transpose(&q_block_t, q_block, stream) != 0) {
        mlx_stream_free(stream);
        if (used_gpu && cpu_stream.ctx != NULL) {
          mlx_stream_free(cpu_stream);
        }
        amatrix_mlx_free_array_if_needed(q_block_t);
        amatrix_mlx_free_array_if_needed(y_block);
        amatrix_mlx_free_array_if_needed(qty_block);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef local transpose failed");
      }

      if (mlx_matmul(&qty_block, q_block_t, y_block, stream) != 0) {
        mlx_stream_free(stream);
        if (used_gpu && cpu_stream.ctx != NULL) {
          mlx_stream_free(cpu_stream);
        }
        amatrix_mlx_free_array_if_needed(q_block_t);
        amatrix_mlx_free_array_if_needed(y_block);
        amatrix_mlx_free_array_if_needed(qty_block);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef local qty failed");
      }

      if (mlx_array_eval(qty_block) != 0) {
        mlx_stream_free(stream);
        if (used_gpu && cpu_stream.ctx != NULL) {
          mlx_stream_free(cpu_stream);
        }
        amatrix_mlx_free_array_if_needed(q_block_t);
        amatrix_mlx_free_array_if_needed(y_block);
        amatrix_mlx_free_array_if_needed(qty_block);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef local eval failed");
      }

      qty_data = mlx_array_data_float32(qty_block);
      if (qty_data == NULL) {
        mlx_stream_free(stream);
        if (used_gpu && cpu_stream.ctx != NULL) {
          mlx_stream_free(cpu_stream);
        }
        amatrix_mlx_free_array_if_needed(q_block_t);
        amatrix_mlx_free_array_if_needed(y_block);
        amatrix_mlx_free_array_if_needed(qty_block);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef local data is unavailable");
      }

      memcpy(
        head_buf + ((size_t)block_idx * (size_t)p * (size_t)nrhs),
        qty_data,
        (size_t)p * (size_t)nrhs * sizeof(float)
      );

      amatrix_mlx_free_array_if_needed(q_block_t);
      amatrix_mlx_free_array_if_needed(y_block);
      amatrix_mlx_free_array_if_needed(qty_block);
      row_start += block_nrow;
    }

    q_top = amatrix_mlx_array_from_resident_key(top_q_key);
    {
      const int head_shape[2] = {head_nrow, nrhs};
      ay = mlx_array_new_data(head_buf, head_shape, 2, MLX_FLOAT32);
    }

    if (mlx_transpose(&qt, q_top, stream) != 0) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef top transpose failed");
    }

    if (mlx_matmul(&qty, qt, ay, stream) != 0) {
      mlx_stream_free(stream);
      if (used_gpu && cpu_stream.ctx != NULL) {
        mlx_stream_free(cpu_stream);
      }
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef top qty failed");
    }

    if (mlx_linalg_solve_triangular(&coef, rr_ref, qty, true, stream) != 0) {
      if (used_gpu) {
        mlx_stream_free(stream);
        stream = cpu_stream;
        coef = mlx_array_new();
        if (mlx_linalg_solve_triangular(&coef, rr_ref, qty, true, stream) != 0) {
          mlx_stream_free(stream);
          amatrix_mlx_free_array_if_needed(qt);
          amatrix_mlx_free_array_if_needed(ay);
          amatrix_mlx_free_array_if_needed(qty);
          amatrix_mlx_free_array_if_needed(coef);
          error("mlx tsqr coef solve failed");
        }
        used_gpu = false;
      } else {
        mlx_stream_free(stream);
        amatrix_mlx_free_array_if_needed(qt);
        amatrix_mlx_free_array_if_needed(ay);
        amatrix_mlx_free_array_if_needed(qty);
        amatrix_mlx_free_array_if_needed(coef);
        error("mlx tsqr coef solve failed");
      }
    }

    if (mlx_synchronize(stream) != 0) {
      mlx_stream_free(stream);
      amatrix_mlx_free_array_if_needed(qt);
      amatrix_mlx_free_array_if_needed(ay);
      amatrix_mlx_free_array_if_needed(qty);
      amatrix_mlx_free_array_if_needed(coef);
      error("mlx tsqr coef synchronize failed");
    }

    result = amatrix_mlx_result_to_r_matrix(coef);

    mlx_stream_free(stream);
    if (used_gpu && cpu_stream.ctx != NULL) {
      mlx_stream_free(cpu_stream);
    }
    amatrix_mlx_free_array_if_needed(qt);
    amatrix_mlx_free_array_if_needed(ay);
    amatrix_mlx_free_array_if_needed(qty);
    amatrix_mlx_free_array_if_needed(coef);
    return result;
  }
#else
  error("mlx tsqr coef resident bridge requires mlx-c");
#endif
}

SEXP amatrix_mlx_tsqr_build_bridge(SEXP x, SEXP block_rows, SEXP q_keys, SEXP top_q_key, SEXP top_r_key, SEXP r_stack_key) {
  if (!isReal(x) || !isMatrix(x)) {
    error("x must be a numeric matrix");
  }
  if (!isInteger(block_rows) || XLENGTH(block_rows) != 1) {
    error("block_rows must be a single integer");
  }
  if (INTEGER(block_rows)[0] < 1) {
    error("block_rows must be positive");
  }
  if (!isString(q_keys)) {
    error("q_keys must be a character vector");
  }
  if (!isString(top_q_key) || XLENGTH(top_q_key) != 1) {
    error("top_q_key must be a single string");
  }
  if (!isString(top_r_key) || XLENGTH(top_r_key) != 1) {
    error("top_r_key must be a single string");
  }
  if (!isString(r_stack_key) || XLENGTH(r_stack_key) != 1) {
    error("r_stack_key must be a single string");
  }

  {
    SEXP dim = getAttrib(x, R_DimSymbol);
    const int nrow = INTEGER(dim)[0];
    const int ncol = INTEGER(dim)[1];
    const int block_rows_val = INTEGER(block_rows)[0];
    const int nblocks = (nrow + block_rows_val - 1) / block_rows_val;
    if ((int) XLENGTH(q_keys) != nblocks) {
      error("q_keys length must match the TSQR block count");
    }
    if (ncol < 1) {
      error("x must have at least one column");
    }
  }

#ifdef HAVE_MLXC
  return amatrix_mlx_tsqr_build_real(x, block_rows, q_keys, top_q_key, top_r_key, r_stack_key);
#else
  error("mlx tsqr build bridge requires mlx-c");
#endif
}
