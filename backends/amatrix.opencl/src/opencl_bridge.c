#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_OPENCL
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#endif

#ifdef HAVE_CLBLAST
#include <clblast_c.h>
#endif

#define AMATRIX_OPENCL_MAX_RESIDENT 256
#define AMATRIX_OPENCL_MAX_SPARSE_RESIDENT 256
#define AMATRIX_OPENCL_REDUCE_WG 64
#define AMATRIX_OPENCL_CHOL_BLOCK 32
#define AMATRIX_OPENCL_CHOL_PANEL_MAX 64

typedef struct {
  char key[64];
  int in_use;
  int nrow;
  int ncol;
  SEXP host_value;
#ifdef HAVE_OPENCL
  cl_mem buffer;
  int on_device;
#endif
} amatrix_opencl_entry;

typedef struct {
  char key[64];
  int in_use;
  int nrow;
  int ncol;
  int nnz;
  double *values;
  double *csr_values;
  int *p;
  int *i;
  int *csr_row_ptr;
  int *csr_col_idx;
#ifdef HAVE_OPENCL
  cl_mem csr_row_ptr_buffer;
  cl_mem csr_col_idx_buffer;
  cl_mem csr_values_buffer;
  cl_mem csc_col_ptr_buffer;
  cl_mem csc_row_idx_buffer;
  cl_mem csc_values_buffer;
  int on_device;
#endif
} amatrix_opencl_sparse_entry;

#ifdef HAVE_OPENCL
typedef struct {
  cl_mem buffer;
  size_t elements;
} amatrix_opencl_workspace;
#endif

static amatrix_opencl_entry g_entries[AMATRIX_OPENCL_MAX_RESIDENT];
static amatrix_opencl_sparse_entry g_sparse_entries[AMATRIX_OPENCL_MAX_SPARSE_RESIDENT];

#ifdef HAVE_OPENCL
static cl_platform_id g_platform = NULL;
static cl_device_id g_device = NULL;
static cl_context g_context = NULL;
static cl_command_queue g_queue = NULL;
static int g_runtime_attempted = 0;
static int g_runtime_available = 0;
static char g_device_name[256] = "";
static cl_program g_custom_program = NULL;
static cl_kernel g_ewise_add_kernel = NULL;
static cl_kernel g_ewise_sub_kernel = NULL;
static cl_kernel g_ewise_div_kernel = NULL;
static cl_kernel g_scalar_mul_kernel = NULL;
static cl_kernel g_broadcast_sweep_kernel = NULL;
static cl_kernel g_row_sum_kernel = NULL;
static cl_kernel g_col_sum_kernel = NULL;
static cl_kernel g_chol_panel_kernel = NULL;
static cl_kernel g_spmm_csr_kernel = NULL;
static cl_kernel g_spmm_csc_trans_kernel = NULL;
static cl_program g_sym_fill_program = NULL;
static cl_kernel g_sym_fill_kernel = NULL;
static cl_kernel g_zero_strict_lower_kernel = NULL;
static amatrix_opencl_workspace g_factor_workspace = {NULL, 0};
static amatrix_opencl_workspace g_status_workspace = {NULL, 0};
static float *g_chol_panel_workspace = NULL;
static size_t g_chol_panel_workspace_len = 0;
#endif

static int amatrix_opencl_nrow(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) {
    Rf_error("expected a matrix with 2 dimensions");
  }
  return INTEGER(dim)[0];
}

static int amatrix_opencl_ncol(SEXP x) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) {
    Rf_error("expected a matrix with 2 dimensions");
  }
  return INTEGER(dim)[1];
}

static void amatrix_opencl_require_matrix(SEXP x, const char *name) {
  if (TYPEOF(x) != REALSXP) {
    Rf_error("%s must be a double matrix", name);
  }
  if (TYPEOF(getAttrib(x, R_DimSymbol)) != INTSXP || XLENGTH(getAttrib(x, R_DimSymbol)) != 2) {
    Rf_error("%s must be a matrix", name);
  }
}

static SEXP amatrix_opencl_named_list(int n, const char **names) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; ++i) {
    SET_STRING_ELT(out_names, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  UNPROTECT(2);
  return out;
}

static void amatrix_opencl_copy_r_to_f32(float *out, const double *in, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    out[i] = (float)in[i];
  }
}

static void amatrix_opencl_copy_f32_to_r(double *out, const float *in, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    out[i] = (double)in[i];
  }
}

static int amatrix_opencl_chol_block_size(void) {
  const char *raw = getenv("AMATRIX_OPENCL_CHOL_BLOCK");
  long value = AMATRIX_OPENCL_CHOL_BLOCK;
  char *end = NULL;

  if (raw != NULL && raw[0] != '\0') {
    value = strtol(raw, &end, 10);
    if (end == raw || *end != '\0') {
      value = AMATRIX_OPENCL_CHOL_BLOCK;
    }
  }

  if (value < 8L) {
    value = 8L;
  }
  if (value > AMATRIX_OPENCL_CHOL_PANEL_MAX) {
    value = AMATRIX_OPENCL_CHOL_PANEL_MAX;
  }

  return (int)value;
}

static int amatrix_opencl_find_entry(const char *key) {
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_RESIDENT; ++idx) {
    if (g_entries[idx].in_use && strcmp(g_entries[idx].key, key) == 0) {
      return idx;
    }
  }
  return -1;
}

static int amatrix_opencl_find_free_entry(void) {
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_RESIDENT; ++idx) {
    if (!g_entries[idx].in_use) {
      return idx;
    }
  }
  return -1;
}

static amatrix_opencl_entry *amatrix_opencl_lookup_entry(const char *key) {
  int idx = amatrix_opencl_find_entry(key);
  if (idx < 0) {
    Rf_error("resident key '%s' not found", key);
  }
  return &g_entries[idx];
}

static int amatrix_opencl_find_sparse_entry(const char *key) {
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_SPARSE_RESIDENT; ++idx) {
    if (g_sparse_entries[idx].in_use && strcmp(g_sparse_entries[idx].key, key) == 0) {
      return idx;
    }
  }
  return -1;
}

static int amatrix_opencl_find_free_sparse_entry(void) {
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_SPARSE_RESIDENT; ++idx) {
    if (!g_sparse_entries[idx].in_use) {
      return idx;
    }
  }
  return -1;
}

static amatrix_opencl_sparse_entry *amatrix_opencl_lookup_sparse_entry(const char *key) {
  int idx = amatrix_opencl_find_sparse_entry(key);
  if (idx < 0) {
    Rf_error("sparse resident key '%s' not found", key);
  }
  return &g_sparse_entries[idx];
}

static void amatrix_opencl_release_entry(amatrix_opencl_entry *entry) {
  if (entry == NULL) {
    return;
  }

  if (entry->host_value != NULL && entry->host_value != R_NilValue) {
    R_ReleaseObject(entry->host_value);
  }
  entry->host_value = NULL;

#ifdef HAVE_OPENCL
  if (entry->buffer != NULL) {
    clReleaseMemObject(entry->buffer);
  }
  entry->buffer = NULL;
  entry->on_device = 0;
#endif

  entry->nrow = 0;
  entry->ncol = 0;
}

static void amatrix_opencl_release_sparse_entry(amatrix_opencl_sparse_entry *entry) {
  if (entry == NULL) {
    return;
  }

  if (entry->values != NULL) {
    free(entry->values);
    entry->values = NULL;
  }
  if (entry->csr_values != NULL) {
    free(entry->csr_values);
    entry->csr_values = NULL;
  }
  if (entry->p != NULL) {
    free(entry->p);
    entry->p = NULL;
  }
  if (entry->i != NULL) {
    free(entry->i);
    entry->i = NULL;
  }
  if (entry->csr_row_ptr != NULL) {
    free(entry->csr_row_ptr);
    entry->csr_row_ptr = NULL;
  }
  if (entry->csr_col_idx != NULL) {
    free(entry->csr_col_idx);
    entry->csr_col_idx = NULL;
  }

#ifdef HAVE_OPENCL
  if (entry->csr_row_ptr_buffer != NULL) {
    clReleaseMemObject(entry->csr_row_ptr_buffer);
    entry->csr_row_ptr_buffer = NULL;
  }
  if (entry->csr_col_idx_buffer != NULL) {
    clReleaseMemObject(entry->csr_col_idx_buffer);
    entry->csr_col_idx_buffer = NULL;
  }
  if (entry->csr_values_buffer != NULL) {
    clReleaseMemObject(entry->csr_values_buffer);
    entry->csr_values_buffer = NULL;
  }
  if (entry->csc_col_ptr_buffer != NULL) {
    clReleaseMemObject(entry->csc_col_ptr_buffer);
    entry->csc_col_ptr_buffer = NULL;
  }
  if (entry->csc_row_idx_buffer != NULL) {
    clReleaseMemObject(entry->csc_row_idx_buffer);
    entry->csc_row_idx_buffer = NULL;
  }
  if (entry->csc_values_buffer != NULL) {
    clReleaseMemObject(entry->csc_values_buffer);
    entry->csc_values_buffer = NULL;
  }
  entry->on_device = 0;
#endif

  entry->nrow = 0;
  entry->ncol = 0;
  entry->nnz = 0;
}

static void amatrix_opencl_commit_entry(
  int idx,
  const char *key,
  int nrow,
  int ncol,
  SEXP host_value
#ifdef HAVE_OPENCL
  , cl_mem buffer,
  int on_device
#endif
) {
  amatrix_opencl_entry *entry = &g_entries[idx];

  if (entry->in_use) {
    amatrix_opencl_release_entry(entry);
  }

  entry->in_use = 1;
  strncpy(entry->key, key, sizeof(entry->key) - 1);
  entry->key[sizeof(entry->key) - 1] = '\0';
  entry->nrow = nrow;
  entry->ncol = ncol;
  entry->host_value = host_value;

#ifdef HAVE_OPENCL
  entry->buffer = buffer;
  entry->on_device = on_device;
#endif
}

static void amatrix_opencl_commit_sparse_entry(int idx, const char *key, amatrix_opencl_sparse_entry *value) {
  amatrix_opencl_sparse_entry *entry = &g_sparse_entries[idx];

  if (entry->in_use) {
    amatrix_opencl_release_sparse_entry(entry);
  }

  *entry = *value;
  entry->in_use = 1;
  strncpy(entry->key, key, sizeof(entry->key) - 1);
  entry->key[sizeof(entry->key) - 1] = '\0';
}

static int amatrix_opencl_resident_count(void) {
  int total = 0;
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_RESIDENT; ++idx) {
    if (g_entries[idx].in_use) {
      total += 1;
    }
  }
  return total;
}

static int amatrix_opencl_host_resident_count(void) {
  int total = 0;
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_RESIDENT; ++idx) {
    if (!g_entries[idx].in_use) {
      continue;
    }
#ifdef HAVE_OPENCL
    if (!g_entries[idx].on_device && g_entries[idx].host_value != NULL && g_entries[idx].host_value != R_NilValue) {
      total += 1;
    }
#else
    if (g_entries[idx].host_value != NULL && g_entries[idx].host_value != R_NilValue) {
      total += 1;
    }
#endif
  }
  return total;
}

#ifdef HAVE_OPENCL
static int amatrix_opencl_device_resident_count(void) {
  int total = 0;
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_RESIDENT; ++idx) {
    if (g_entries[idx].in_use && g_entries[idx].on_device && g_entries[idx].buffer != NULL) {
      total += 1;
    }
  }
  return total;
}

static int amatrix_opencl_probe_enabled(void) {
  const char *probe = getenv("AMATRIX_OPENCL_PROBE_GPU");
  return probe != NULL && strcmp(probe, "1") == 0;
}

static void amatrix_opencl_runtime_clear(void) {
  if (g_ewise_add_kernel != NULL) {
    clReleaseKernel(g_ewise_add_kernel);
    g_ewise_add_kernel = NULL;
  }
  if (g_ewise_sub_kernel != NULL) {
    clReleaseKernel(g_ewise_sub_kernel);
    g_ewise_sub_kernel = NULL;
  }
  if (g_ewise_div_kernel != NULL) {
    clReleaseKernel(g_ewise_div_kernel);
    g_ewise_div_kernel = NULL;
  }
  if (g_scalar_mul_kernel != NULL) {
    clReleaseKernel(g_scalar_mul_kernel);
    g_scalar_mul_kernel = NULL;
  }
  if (g_broadcast_sweep_kernel != NULL) {
    clReleaseKernel(g_broadcast_sweep_kernel);
    g_broadcast_sweep_kernel = NULL;
  }
  if (g_row_sum_kernel != NULL) {
    clReleaseKernel(g_row_sum_kernel);
    g_row_sum_kernel = NULL;
  }
  if (g_col_sum_kernel != NULL) {
    clReleaseKernel(g_col_sum_kernel);
    g_col_sum_kernel = NULL;
  }
  if (g_chol_panel_kernel != NULL) {
    clReleaseKernel(g_chol_panel_kernel);
    g_chol_panel_kernel = NULL;
  }
  if (g_spmm_csr_kernel != NULL) {
    clReleaseKernel(g_spmm_csr_kernel);
    g_spmm_csr_kernel = NULL;
  }
  if (g_spmm_csc_trans_kernel != NULL) {
    clReleaseKernel(g_spmm_csc_trans_kernel);
    g_spmm_csc_trans_kernel = NULL;
  }
  if (g_custom_program != NULL) {
    clReleaseProgram(g_custom_program);
    g_custom_program = NULL;
  }
  if (g_sym_fill_kernel != NULL) {
    clReleaseKernel(g_sym_fill_kernel);
    g_sym_fill_kernel = NULL;
  }
  if (g_zero_strict_lower_kernel != NULL) {
    clReleaseKernel(g_zero_strict_lower_kernel);
    g_zero_strict_lower_kernel = NULL;
  }
  if (g_sym_fill_program != NULL) {
    clReleaseProgram(g_sym_fill_program);
    g_sym_fill_program = NULL;
  }
  if (g_factor_workspace.buffer != NULL) {
    clReleaseMemObject(g_factor_workspace.buffer);
    g_factor_workspace.buffer = NULL;
    g_factor_workspace.elements = 0;
  }
  if (g_status_workspace.buffer != NULL) {
    clReleaseMemObject(g_status_workspace.buffer);
    g_status_workspace.buffer = NULL;
    g_status_workspace.elements = 0;
  }
  if (g_chol_panel_workspace != NULL) {
    free(g_chol_panel_workspace);
    g_chol_panel_workspace = NULL;
    g_chol_panel_workspace_len = 0;
  }
  for (int idx = 0; idx < AMATRIX_OPENCL_MAX_SPARSE_RESIDENT; ++idx) {
    if (!g_sparse_entries[idx].in_use) {
      continue;
    }
    if (g_sparse_entries[idx].csr_row_ptr_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csr_row_ptr_buffer);
      g_sparse_entries[idx].csr_row_ptr_buffer = NULL;
    }
    if (g_sparse_entries[idx].csr_col_idx_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csr_col_idx_buffer);
      g_sparse_entries[idx].csr_col_idx_buffer = NULL;
    }
    if (g_sparse_entries[idx].csr_values_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csr_values_buffer);
      g_sparse_entries[idx].csr_values_buffer = NULL;
    }
    if (g_sparse_entries[idx].csc_col_ptr_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csc_col_ptr_buffer);
      g_sparse_entries[idx].csc_col_ptr_buffer = NULL;
    }
    if (g_sparse_entries[idx].csc_row_idx_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csc_row_idx_buffer);
      g_sparse_entries[idx].csc_row_idx_buffer = NULL;
    }
    if (g_sparse_entries[idx].csc_values_buffer != NULL) {
      clReleaseMemObject(g_sparse_entries[idx].csc_values_buffer);
      g_sparse_entries[idx].csc_values_buffer = NULL;
    }
    g_sparse_entries[idx].on_device = 0;
  }
  if (g_queue != NULL) {
    clReleaseCommandQueue(g_queue);
    g_queue = NULL;
  }
  if (g_context != NULL) {
    clReleaseContext(g_context);
    g_context = NULL;
  }
  g_platform = NULL;
  g_device = NULL;
  g_runtime_available = 0;
  g_device_name[0] = '\0';
}

static int amatrix_opencl_try_init(void) {
  cl_int err = CL_SUCCESS;
  cl_uint num_platforms = 0;
  cl_platform_id platforms[16];

  if (g_runtime_available) {
    return 1;
  }
  if (!amatrix_opencl_probe_enabled()) {
    return 0;
  }
  if (g_runtime_attempted) {
    return 0;
  }

  g_runtime_attempted = 1;

  err = clGetPlatformIDs((cl_uint)(sizeof(platforms) / sizeof(platforms[0])), platforms, &num_platforms);
  if (err != CL_SUCCESS || num_platforms == 0) {
    amatrix_opencl_runtime_clear();
    return 0;
  }

  for (cl_uint p = 0; p < num_platforms; ++p) {
    cl_uint num_devices = 0;
    err = clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
    if (err != CL_SUCCESS || num_devices == 0) {
      continue;
    }

    cl_device_id devices[16];
    if (num_devices > (cl_uint)(sizeof(devices) / sizeof(devices[0]))) {
      num_devices = (cl_uint)(sizeof(devices) / sizeof(devices[0]));
    }

    err = clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
    if (err != CL_SUCCESS || num_devices == 0) {
      continue;
    }

    g_platform = platforms[p];
    g_device = devices[0];
    break;
  }

  if (g_device == NULL) {
    amatrix_opencl_runtime_clear();
    return 0;
  }

  g_context = clCreateContext(NULL, 1, &g_device, NULL, NULL, &err);
  if (err != CL_SUCCESS || g_context == NULL) {
    amatrix_opencl_runtime_clear();
    return 0;
  }

  g_queue = clCreateCommandQueue(g_context, g_device, 0, &err);
  if (err != CL_SUCCESS || g_queue == NULL) {
    amatrix_opencl_runtime_clear();
    return 0;
  }

  {
    size_t name_size = 0;
    err = clGetDeviceInfo(g_device, CL_DEVICE_NAME, sizeof(g_device_name), g_device_name, &name_size);
    if (err != CL_SUCCESS || name_size == 0) {
      strncpy(g_device_name, "unknown", sizeof(g_device_name) - 1);
      g_device_name[sizeof(g_device_name) - 1] = '\0';
    }
  }

  g_runtime_available = 1;
  return 1;
}

static cl_mem amatrix_opencl_buffer_from_r(SEXP x) {
  size_t n = (size_t)amatrix_opencl_nrow(x) * (size_t)amatrix_opencl_ncol(x);
  size_t bytes = n * sizeof(float);
  float *host = (float *)R_alloc(n, sizeof(float));
  cl_mem buffer = NULL;
  cl_int err = CL_SUCCESS;

  amatrix_opencl_copy_r_to_f32(host, REAL(x), n);

  buffer = clCreateBuffer(g_context, CL_MEM_READ_WRITE, bytes, NULL, &err);
  if (err != CL_SUCCESS || buffer == NULL) {
    return NULL;
  }

  err = clEnqueueWriteBuffer(g_queue, buffer, CL_TRUE, 0, bytes, host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    clReleaseMemObject(buffer);
    return NULL;
  }

  return buffer;
}

static cl_mem amatrix_opencl_buffer_from_host(const void *host, size_t bytes) {
  cl_mem buffer = NULL;
  cl_int err = CL_SUCCESS;

  if (bytes == 0) {
    return NULL;
  }

  buffer = clCreateBuffer(g_context, CL_MEM_READ_WRITE, bytes, NULL, &err);
  if (err != CL_SUCCESS || buffer == NULL) {
    return NULL;
  }

  err = clEnqueueWriteBuffer(g_queue, buffer, CL_TRUE, 0, bytes, host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    clReleaseMemObject(buffer);
    return NULL;
  }

  return buffer;
}

static SEXP amatrix_opencl_matrix_from_buffer(cl_mem buffer, int nrow, int ncol) {
  size_t n = (size_t)nrow * (size_t)ncol;
  size_t bytes = n * sizeof(float);
  float *host = (float *)R_alloc(n, sizeof(float));
  cl_int err = clEnqueueReadBuffer(g_queue, buffer, CL_TRUE, 0, bytes, host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    Rf_error("failed to read OpenCL resident buffer");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
  amatrix_opencl_copy_f32_to_r(REAL(out), host, n);
  UNPROTECT(1);
  return out;
}

static cl_mem amatrix_opencl_buffer_from_vector(SEXP x, int *length_out) {
  size_t n = (size_t)XLENGTH(x);
  size_t bytes = n * sizeof(float);
  float *host = NULL;
  cl_mem buffer = NULL;
  cl_int err = CL_SUCCESS;

  if (TYPEOF(x) != REALSXP) {
    return NULL;
  }

  host = (float *)R_alloc(n, sizeof(float));
  amatrix_opencl_copy_r_to_f32(host, REAL(x), n);

  buffer = clCreateBuffer(g_context, CL_MEM_READ_WRITE, bytes, NULL, &err);
  if (err != CL_SUCCESS || buffer == NULL) {
    return NULL;
  }

  err = clEnqueueWriteBuffer(g_queue, buffer, CL_TRUE, 0, bytes, host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    clReleaseMemObject(buffer);
    return NULL;
  }

  if (length_out != NULL) {
    *length_out = (int)n;
  }
  return buffer;
}

static SEXP amatrix_opencl_vector_from_buffer(cl_mem buffer, int length) {
  size_t n = (size_t)length;
  size_t bytes = n * sizeof(float);
  float *host = (float *)R_alloc(n, sizeof(float));
  cl_int err = clEnqueueReadBuffer(g_queue, buffer, CL_TRUE, 0, bytes, host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    Rf_error("failed to read OpenCL vector buffer");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, length));
  amatrix_opencl_copy_f32_to_r(REAL(out), host, n);
  UNPROTECT(1);
  return out;
}

static int amatrix_opencl_entry_vector_length(const amatrix_opencl_entry *entry) {
  if (entry == NULL) {
    return 0;
  }
  if (entry->ncol == 1) {
    return entry->nrow;
  }
  if (entry->nrow == 1) {
    return entry->ncol;
  }
  return 0;
}

static SEXP amatrix_opencl_get_entry_vector(const char *key) {
  int idx = amatrix_opencl_find_entry(key);
  amatrix_opencl_entry *entry = NULL;
  SEXP materialized = R_NilValue;
  int length = 0;

  if (idx < 0) {
    Rf_error("resident key '%s' not found", key);
  }

  entry = &g_entries[idx];
  length = amatrix_opencl_entry_vector_length(entry);
  if (length <= 0) {
    Rf_error("resident key '%s' is not vector-shaped", key);
  }

#ifdef HAVE_OPENCL
  if (entry->on_device && entry->buffer != NULL) {
    if (!amatrix_opencl_try_init()) {
      Rf_error("OpenCL runtime is not available for resident key '%s'", key);
    }
    return amatrix_opencl_vector_from_buffer(entry->buffer, length);
  }
#endif

  if (entry->host_value == NULL || entry->host_value == R_NilValue) {
    Rf_error("resident key '%s' has no materializable value", key);
  }

  materialized = entry->host_value;
  if (TYPEOF(materialized) == REALSXP && TYPEOF(getAttrib(materialized, R_DimSymbol)) != INTSXP) {
    return Rf_duplicate(materialized);
  }

  if (TYPEOF(materialized) == REALSXP &&
      TYPEOF(getAttrib(materialized, R_DimSymbol)) == INTSXP &&
      XLENGTH(getAttrib(materialized, R_DimSymbol)) == 2) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, length));
    memcpy(REAL(out), REAL(materialized), (size_t)length * sizeof(double));
    UNPROTECT(1);
    return out;
  }

  Rf_error("resident key '%s' has unsupported vector materialization type", key);
  return R_NilValue;
}

static cl_mem amatrix_opencl_vector_buffer_from_arg(SEXP arg, int *length_out, int *release_buffer) {
#ifdef HAVE_OPENCL
  if (TYPEOF(arg) == STRSXP && XLENGTH(arg) == 1) {
    amatrix_opencl_entry *entry = amatrix_opencl_lookup_entry(CHAR(STRING_ELT(arg, 0)));
    int length = amatrix_opencl_entry_vector_length(entry);
    if (length > 0 && entry->on_device && entry->buffer != NULL) {
      if (length_out != NULL) {
        *length_out = length;
      }
      if (release_buffer != NULL) {
        *release_buffer = 0;
      }
      return entry->buffer;
    }
    return NULL;
  }

  if (TYPEOF(arg) == REALSXP && TYPEOF(getAttrib(arg, R_DimSymbol)) != INTSXP) {
    cl_mem buffer = amatrix_opencl_buffer_from_vector(arg, length_out);
    if (release_buffer != NULL) {
      *release_buffer = (buffer != NULL);
    }
    return buffer;
  }
#endif
  return NULL;
}

#ifdef HAVE_CLBLAST
static int amatrix_opencl_make_product_dims(
  int ar, int ac, int br, int bc, int trans_a, int trans_b, int *m, int *n, int *k
) {
  int out_m = trans_a ? ac : ar;
  int inner_a = trans_a ? ar : ac;
  int inner_b = trans_b ? bc : br;
  int out_n = trans_b ? br : bc;

  if (inner_a != inner_b) {
    return 0;
  }

  *m = out_m;
  *n = out_n;
  *k = inner_a;
  return 1;
}

static cl_mem amatrix_opencl_alloc_buffer(size_t elements) {
  cl_int err = CL_SUCCESS;
  cl_mem buffer = clCreateBuffer(g_context, CL_MEM_READ_WRITE, elements * sizeof(float), NULL, &err);
  if (err != CL_SUCCESS || buffer == NULL) {
    return NULL;
  }
  return buffer;
}

static int amatrix_opencl_copy_buffer_into(cl_mem dest, cl_mem source, size_t elements) {
  cl_int err = clEnqueueCopyBuffer(g_queue, source, dest, 0, 0, elements * sizeof(float), 0, NULL, NULL);
  return err == CL_SUCCESS;
}

static int amatrix_opencl_workspace_ensure(amatrix_opencl_workspace *workspace, size_t elements) {
  cl_mem buffer = NULL;

  if (workspace == NULL) {
    return 0;
  }
  if (workspace->buffer != NULL && workspace->elements >= elements) {
    return 1;
  }

  if (workspace->buffer != NULL) {
    clReleaseMemObject(workspace->buffer);
    workspace->buffer = NULL;
    workspace->elements = 0;
  }

  buffer = amatrix_opencl_alloc_buffer(elements);
  if (buffer == NULL) {
    return 0;
  }

  workspace->buffer = buffer;
  workspace->elements = elements;
  return 1;
}

static int amatrix_opencl_copy_buffer(cl_mem source, size_t elements, cl_mem *out_buffer) {
  cl_int err = CL_SUCCESS;
  cl_mem out = amatrix_opencl_alloc_buffer(elements);
  if (out == NULL) {
    return 0;
  }

  err = clEnqueueCopyBuffer(g_queue, source, out, 0, 0, elements * sizeof(float), 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  return 1;
}

static float *amatrix_opencl_panel_workspace(size_t elements) {
  if (elements == 0) {
    return NULL;
  }

  if (elements > g_chol_panel_workspace_len) {
    float *next = (float *)realloc(g_chol_panel_workspace, elements * sizeof(float));
    if (next == NULL) {
      return NULL;
    }
    g_chol_panel_workspace = next;
    g_chol_panel_workspace_len = elements;
  }

  return g_chol_panel_workspace;
}

static int amatrix_opencl_read_block(
  cl_mem buffer,
  int ld,
  int row_offset,
  int col_offset,
  int rows,
  int cols,
  float *host
) {
  const size_t buffer_origin[3] = {
    (size_t)row_offset * sizeof(float),
    (size_t)col_offset,
    0
  };
  const size_t host_origin[3] = {0, 0, 0};
  const size_t region[3] = {
    (size_t)rows * sizeof(float),
    (size_t)cols,
    1
  };
  const size_t buffer_row_pitch = (size_t)ld * sizeof(float);
  const size_t host_row_pitch = (size_t)rows * sizeof(float);
  const size_t buffer_slice_pitch = buffer_row_pitch * (size_t)cols;
  const size_t host_slice_pitch = host_row_pitch * (size_t)cols;
  cl_int err = clEnqueueReadBufferRect(
    g_queue,
    buffer,
    CL_TRUE,
    buffer_origin,
    host_origin,
    region,
    buffer_row_pitch,
    buffer_slice_pitch,
    host_row_pitch,
    host_slice_pitch,
    host,
    0,
    NULL,
    NULL
  );
  return err == CL_SUCCESS;
}

static int amatrix_opencl_write_block(
  cl_mem buffer,
  int ld,
  int row_offset,
  int col_offset,
  int rows,
  int cols,
  const float *host
) {
  const size_t buffer_origin[3] = {
    (size_t)row_offset * sizeof(float),
    (size_t)col_offset,
    0
  };
  const size_t host_origin[3] = {0, 0, 0};
  const size_t region[3] = {
    (size_t)rows * sizeof(float),
    (size_t)cols,
    1
  };
  const size_t buffer_row_pitch = (size_t)ld * sizeof(float);
  const size_t host_row_pitch = (size_t)rows * sizeof(float);
  const size_t buffer_slice_pitch = buffer_row_pitch * (size_t)cols;
  const size_t host_slice_pitch = host_row_pitch * (size_t)cols;
  cl_int err = clEnqueueWriteBufferRect(
    g_queue,
    buffer,
    CL_TRUE,
    buffer_origin,
    host_origin,
    region,
    buffer_row_pitch,
    buffer_slice_pitch,
    host_row_pitch,
    host_slice_pitch,
    host,
    0,
    NULL,
    NULL
  );
  return err == CL_SUCCESS;
}

static int amatrix_opencl_zero_strict_lower_buffer(cl_mem buffer, int n);
static int amatrix_opencl_ensure_custom_kernels(void);
static int amatrix_opencl_run_chol_panel_upper_inplace(cl_mem buffer, cl_mem status_buffer, int ld, int offset, int block);

static int amatrix_opencl_host_chol_upper(float *block, int n) {
  for (int j = 0; j < n; ++j) {
    double diag = (double)block[j + j * n];

    for (int k = 0; k < j; ++k) {
      double u = (double)block[k + j * n];
      diag -= u * u;
    }

    if (!(diag > 0.0) || !R_FINITE(diag)) {
      return 0;
    }

    block[j + j * n] = (float)sqrt(diag);

    for (int col = j + 1; col < n; ++col) {
      double value = (double)block[j + col * n];
      for (int k = 0; k < j; ++k) {
        value -= (double)block[k + j * n] * (double)block[k + col * n];
      }
      block[j + col * n] = (float)(value / (double)block[j + j * n]);
    }

    for (int row = j + 1; row < n; ++row) {
      block[row + j * n] = 0.0f;
    }
  }

  return 1;
}

static int amatrix_opencl_run_trsm_left(
  cl_mem a_buffer,
  size_t a_offset,
  size_t a_ld,
  int lower,
  int transpose_a,
  cl_mem b_buffer,
  size_t b_offset,
  size_t b_ld,
  int m,
  int n
) {
  CLBlastStatusCode status = CLBlastStrsm(
    CLBlastLayoutColMajor,
    CLBlastSideLeft,
    lower ? CLBlastTriangleLower : CLBlastTriangleUpper,
    transpose_a ? CLBlastTransposeYes : CLBlastTransposeNo,
    CLBlastDiagonalNonUnit,
    (size_t)m,
    (size_t)n,
    1.0f,
    a_buffer,
    a_offset,
    a_ld,
    b_buffer,
    b_offset,
    b_ld,
    &g_queue,
    NULL
  );

  return status == CLBlastSuccess;
}

static int amatrix_opencl_run_triangular_solve(
  cl_mem factor_buffer,
  int n,
  int lower,
  int transpose,
  cl_mem rhs_buffer,
  int nrhs,
  cl_mem *out_buffer
) {
  cl_mem result_buffer = NULL;

  if (!amatrix_opencl_copy_buffer(rhs_buffer, (size_t)n * (size_t)nrhs, &result_buffer)) {
    return 0;
  }

  if (!amatrix_opencl_run_trsm_left(
        factor_buffer, 0, (size_t)n, lower, transpose,
        result_buffer, 0, (size_t)n,
        n, nrhs
      )) {
    clReleaseMemObject(result_buffer);
    return 0;
  }

  *out_buffer = result_buffer;
  return 1;
}

static int amatrix_opencl_run_syrk_update_upper(
  cl_mem a_buffer,
  size_t a_offset,
  size_t a_ld,
  cl_mem c_buffer,
  size_t c_offset,
  size_t c_ld,
  int n,
  int k
) {
  CLBlastStatusCode status = CLBlastSsyrk(
    CLBlastLayoutColMajor,
    CLBlastTriangleUpper,
    CLBlastTransposeYes,
    (size_t)n,
    (size_t)k,
    -1.0f,
    a_buffer,
    a_offset,
    a_ld,
    1.0f,
    c_buffer,
    c_offset,
    c_ld,
    &g_queue,
    NULL
  );

  return status == CLBlastSuccess;
}

static cl_mem amatrix_opencl_identity_buffer(int n) {
  cl_mem buffer = NULL;
  float *host = NULL;
  size_t total = (size_t)n * (size_t)n;
  cl_int err = CL_SUCCESS;

  host = (float *)R_alloc(total, sizeof(float));
  memset(host, 0, total * sizeof(float));
  for (int j = 0; j < n; ++j) {
    host[j + j * n] = 1.0f;
  }

  buffer = amatrix_opencl_alloc_buffer(total);
  if (buffer == NULL) {
    return NULL;
  }

  err = clEnqueueWriteBuffer(g_queue, buffer, CL_TRUE, 0, total * sizeof(float), host, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    clReleaseMemObject(buffer);
    return NULL;
  }

  return buffer;
}

static int amatrix_opencl_run_chol_upper_inplace(cl_mem buffer, int n) {
  int block_size = amatrix_opencl_chol_block_size();
  cl_mem status_buffer = NULL;
  cl_int status_value = 1;

  if (n <= 0) {
    return 1;
  }

  if (amatrix_opencl_ensure_custom_kernels() && g_chol_panel_kernel != NULL) {
    if (!amatrix_opencl_workspace_ensure(&g_status_workspace, 1)) {
      return 0;
    }
    status_buffer = g_status_workspace.buffer;
    if (status_buffer == NULL ||
        clEnqueueWriteBuffer(g_queue, status_buffer, CL_FALSE, 0, sizeof(cl_int), &status_value, 0, NULL, NULL) != CL_SUCCESS) {
      return 0;
    }
  }

  for (int offset = 0; offset < n; offset += block_size) {
    int block = n - offset;
    if (block > block_size) {
      block = block_size;
    }

    if (!amatrix_opencl_run_chol_panel_upper_inplace(buffer, status_buffer, n, offset, block)) {
      return 0;
    }

    if (offset + block < n) {
      int trailing = n - offset - block;
      size_t diag_offset = (size_t)offset + (size_t)offset * (size_t)n;
      size_t row_offset = (size_t)offset + (size_t)(offset + block) * (size_t)n;
      size_t trail_offset = (size_t)(offset + block) + (size_t)(offset + block) * (size_t)n;

      if (!amatrix_opencl_run_trsm_left(
            buffer, diag_offset, (size_t)n, 0, 1,
            buffer, row_offset, (size_t)n,
            block, trailing
          )) {
        return 0;
      }

      if (!amatrix_opencl_run_syrk_update_upper(
            buffer, row_offset, (size_t)n,
            buffer, trail_offset, (size_t)n,
            trailing, block
          )) {
        return 0;
      }
    }
  }

  if (!amatrix_opencl_zero_strict_lower_buffer(buffer, n)) {
    return 0;
  }

  if (status_buffer != NULL &&
      (clEnqueueReadBuffer(g_queue, status_buffer, CL_TRUE, 0, sizeof(cl_int), &status_value, 0, NULL, NULL) != CL_SUCCESS ||
       status_value != 1)) {
    return 0;
  }

  return 1;
}

static int amatrix_opencl_run_chol_solve(
  cl_mem a_buffer,
  int n,
  cl_mem b_buffer,
  int nrhs,
  cl_mem *out_buffer
) {
  cl_mem factor_buffer = NULL;
  cl_mem result_buffer = NULL;
  size_t factor_elements = (size_t)n * (size_t)n;

  if (!amatrix_opencl_workspace_ensure(&g_factor_workspace, factor_elements)) {
    return 0;
  }
  factor_buffer = g_factor_workspace.buffer;
  if (!amatrix_opencl_copy_buffer_into(factor_buffer, a_buffer, factor_elements)) {
    return 0;
  }
  if (!amatrix_opencl_run_chol_upper_inplace(factor_buffer, n)) {
    return 0;
  }

  if (b_buffer == NULL) {
    result_buffer = amatrix_opencl_identity_buffer(n);
    nrhs = n;
  } else if (!amatrix_opencl_copy_buffer(b_buffer, (size_t)n * (size_t)nrhs, &result_buffer)) {
    return 0;
  }

  if (result_buffer == NULL) {
    return 0;
  }

  if (!amatrix_opencl_run_trsm_left(
        factor_buffer, 0, (size_t)n, 0, 1,
        result_buffer, 0, (size_t)n,
        n, nrhs
      ) ||
      !amatrix_opencl_run_trsm_left(
        factor_buffer, 0, (size_t)n, 0, 0,
        result_buffer, 0, (size_t)n,
        n, nrhs
      )) {
    clReleaseMemObject(result_buffer);
    return 0;
  }

  *out_buffer = result_buffer;
  return 1;
}

static int amatrix_opencl_op_code(const char *op) {
  if (strcmp(op, "+") == 0) return 1;
  if (strcmp(op, "-") == 0) return 2;
  if (strcmp(op, "*") == 0) return 3;
  if (strcmp(op, "/") == 0) return 4;
  return 0;
}

static int amatrix_opencl_ensure_custom_kernels(void) {
  static const char *source =
    "__kernel void ewise_add(__global const float* a, __global const float* b, const float scalar, const int use_scalar, __global float* out, const int n) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  if (i < n) { const float rhs = use_scalar ? scalar : b[i]; out[i] = a[i] + rhs; }\n"
    "}\n"
    "__kernel void ewise_sub(__global const float* a, __global const float* b, const float scalar, const int use_scalar, __global float* out, const int n) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  if (i < n) { const float rhs = use_scalar ? scalar : b[i]; out[i] = a[i] - rhs; }\n"
    "}\n"
    "__kernel void ewise_div(__global const float* a, __global const float* b, const float scalar, const int use_scalar, __global float* out, const int n) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  if (i < n) { const float rhs = use_scalar ? scalar : b[i]; out[i] = a[i] / rhs; }\n"
    "}\n"
    "__kernel void scalar_mul(__global const float* a, const float scalar, __global float* out, const int n) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  if (i < n) { out[i] = a[i] * scalar; }\n"
    "}\n"
    "__kernel void broadcast_sweep(__global const float* a, __global const float* v, const int nrow, const int ncol, const int margin, const int op_code, __global float* out) {\n"
    "  const int idx = (int)get_global_id(0);\n"
    "  const int n = nrow * ncol;\n"
    "  if (idx < n) {\n"
    "    const int i = idx % nrow;\n"
    "    const int j = idx / nrow;\n"
    "    const float lhs = a[idx];\n"
    "    const float rhs = (margin == 1) ? v[i] : v[j];\n"
    "    if (op_code == 1) out[idx] = lhs + rhs;\n"
    "    else if (op_code == 2) out[idx] = lhs - rhs;\n"
    "    else if (op_code == 3) out[idx] = lhs * rhs;\n"
    "    else if (op_code == 4) out[idx] = lhs / rhs;\n"
    "  }\n"
    "}\n"
    "#define AMATRIX_REDUCE_WG 64\n"
    "__kernel void row_sum(__global const float* a, const int nrow, const int ncol, __global float* out) {\n"
    "  const int row = (int)get_group_id(0);\n"
    "  const int lid = (int)get_local_id(0);\n"
    "  const int lsize = (int)get_local_size(0);\n"
    "  __local float scratch[AMATRIX_REDUCE_WG];\n"
    "  float acc = 0.0f;\n"
    "  if (row < nrow) {\n"
    "    for (int j = lid; j < ncol; j += lsize) acc += a[row + j * nrow];\n"
    "  }\n"
    "  scratch[lid] = acc;\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  for (int stride = lsize / 2; stride > 0; stride >>= 1) {\n"
    "    if (lid < stride) scratch[lid] += scratch[lid + stride];\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  }\n"
    "  if (lid == 0 && row < nrow) out[row] = scratch[0];\n"
    "}\n"
    "__kernel void col_sum(__global const float* a, const int nrow, const int ncol, __global float* out) {\n"
    "  const int col = (int)get_group_id(0);\n"
    "  const int lid = (int)get_local_id(0);\n"
    "  const int lsize = (int)get_local_size(0);\n"
    "  __local float scratch[AMATRIX_REDUCE_WG];\n"
    "  float acc = 0.0f;\n"
    "  if (col < ncol) {\n"
    "    const int base = col * nrow;\n"
    "    for (int i = lid; i < nrow; i += lsize) acc += a[base + i];\n"
    "  }\n"
    "  scratch[lid] = acc;\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  for (int stride = lsize / 2; stride > 0; stride >>= 1) {\n"
    "    if (lid < stride) scratch[lid] += scratch[lid + stride];\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  }\n"
    "  if (lid == 0 && col < ncol) out[col] = scratch[0];\n"
    "}\n"
    "#define AMATRIX_CHOL_PANEL_MAX 64\n"
    "__kernel void chol_panel_upper(__global float* x, const int ld, const int offset, const int block, __global int* status) {\n"
    "  const int lid = (int)get_local_id(0);\n"
    "  const int lsize = (int)get_local_size(0);\n"
    "  __local float tile[AMATRIX_CHOL_PANEL_MAX * AMATRIX_CHOL_PANEL_MAX];\n"
    "  __local int ok;\n"
    "  if (block > AMATRIX_CHOL_PANEL_MAX) { if (lid == 0) status[0] = 0; return; }\n"
    "  for (int idx = lid; idx < block * block; idx += lsize) {\n"
    "    const int row = idx % block;\n"
    "    const int col = idx / block;\n"
    "    tile[idx] = x[(offset + row) + (offset + col) * ld];\n"
    "  }\n"
    "  if (lid == 0) ok = 1;\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  for (int j = 0; j < block; ++j) {\n"
    "    if (lid == 0) {\n"
    "      float diag = tile[j + j * block];\n"
    "      for (int k = 0; k < j; ++k) diag -= tile[k + j * block] * tile[k + j * block];\n"
    "      if (!(diag > 0.0f)) { ok = 0; status[0] = 0; }\n"
    "      else tile[j + j * block] = sqrt(diag);\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "    if (!ok) return;\n"
    "    if (lid > j && lid < block) {\n"
    "      const int col = lid;\n"
    "      float value = tile[j + col * block];\n"
    "      for (int k = 0; k < j; ++k) value -= tile[k + j * block] * tile[k + col * block];\n"
    "      tile[j + col * block] = value / tile[j + j * block];\n"
    "      tile[col + j * block] = 0.0f;\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  }\n"
    "  for (int idx = lid; idx < block * block; idx += lsize) {\n"
    "    const int row = idx % block;\n"
    "    const int col = idx / block;\n"
    "    x[(offset + row) + (offset + col) * ld] = tile[idx];\n"
    "  }\n"
    "}\n"
    "__kernel void sparse_spmm_csr(__global const int* row_ptr, __global const int* col_idx, __global const float* values, __global const float* b, const int out_nrow, const int rhs_ncol, const int rhs_ld, __global float* out) {\n"
    "  const int row = (int)get_global_id(0);\n"
    "  const int col = (int)get_global_id(1);\n"
    "  if (row >= out_nrow || col >= rhs_ncol) return;\n"
    "  float acc = 0.0f;\n"
    "  for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx) {\n"
    "    acc += values[idx] * b[col_idx[idx] + col * rhs_ld];\n"
    "  }\n"
    "  out[row + col * out_nrow] = acc;\n"
    "}\n"
    "__kernel void sparse_spmm_csc_trans(__global const int* col_ptr, __global const int* row_idx, __global const float* values, __global const float* b, const int out_nrow, const int rhs_ncol, const int rhs_ld, __global float* out) {\n"
    "  const int xcol = (int)get_global_id(0);\n"
    "  const int col = (int)get_global_id(1);\n"
    "  if (xcol >= out_nrow || col >= rhs_ncol) return;\n"
    "  float acc = 0.0f;\n"
    "  for (int idx = col_ptr[xcol]; idx < col_ptr[xcol + 1]; ++idx) {\n"
    "    acc += values[idx] * b[row_idx[idx] + col * rhs_ld];\n"
    "  }\n"
    "  out[xcol + col * out_nrow] = acc;\n"
    "}\n";

  if (g_custom_program != NULL) {
    if (g_ewise_add_kernel != NULL &&
        g_ewise_sub_kernel != NULL &&
        g_ewise_div_kernel != NULL &&
        g_scalar_mul_kernel != NULL &&
        g_broadcast_sweep_kernel != NULL &&
        g_row_sum_kernel != NULL &&
        g_col_sum_kernel != NULL &&
        g_chol_panel_kernel != NULL &&
        g_spmm_csr_kernel != NULL &&
        g_spmm_csc_trans_kernel != NULL) {
      return 1;
    }

    if (g_chol_panel_kernel == NULL) {
      cl_int err = CL_SUCCESS;
      g_chol_panel_kernel = clCreateKernel(g_custom_program, "chol_panel_upper", &err);
      if (err != CL_SUCCESS || g_chol_panel_kernel == NULL) {
        if (g_chol_panel_kernel != NULL) {
          clReleaseKernel(g_chol_panel_kernel);
          g_chol_panel_kernel = NULL;
        }
        return 0;
      }
    }

    if (g_spmm_csr_kernel == NULL) {
      cl_int err = CL_SUCCESS;
      g_spmm_csr_kernel = clCreateKernel(g_custom_program, "sparse_spmm_csr", &err);
      if (err != CL_SUCCESS || g_spmm_csr_kernel == NULL) {
        if (g_spmm_csr_kernel != NULL) {
          clReleaseKernel(g_spmm_csr_kernel);
          g_spmm_csr_kernel = NULL;
        }
        return 0;
      }
    }

    if (g_spmm_csc_trans_kernel == NULL) {
      cl_int err = CL_SUCCESS;
      g_spmm_csc_trans_kernel = clCreateKernel(g_custom_program, "sparse_spmm_csc_trans", &err);
      if (err != CL_SUCCESS || g_spmm_csc_trans_kernel == NULL) {
        if (g_spmm_csc_trans_kernel != NULL) {
          clReleaseKernel(g_spmm_csc_trans_kernel);
          g_spmm_csc_trans_kernel = NULL;
        }
        return 0;
      }
    }

    if (g_ewise_add_kernel != NULL &&
        g_ewise_sub_kernel != NULL &&
        g_ewise_div_kernel != NULL &&
        g_scalar_mul_kernel != NULL &&
        g_broadcast_sweep_kernel != NULL &&
        g_row_sum_kernel != NULL &&
        g_col_sum_kernel != NULL &&
        g_chol_panel_kernel != NULL &&
        g_spmm_csr_kernel != NULL &&
        g_spmm_csc_trans_kernel != NULL) {
      return 1;
    }

    return 0;
  }

  {
    cl_int err = CL_SUCCESS;
    size_t length = strlen(source);

    g_custom_program = clCreateProgramWithSource(g_context, 1, &source, &length, &err);
    if (err != CL_SUCCESS || g_custom_program == NULL) {
      g_custom_program = NULL;
      return 0;
    }

    err = clBuildProgram(g_custom_program, 1, &g_device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
      clReleaseProgram(g_custom_program);
      g_custom_program = NULL;
      return 0;
    }

    g_ewise_add_kernel = clCreateKernel(g_custom_program, "ewise_add", &err);
    if (err != CL_SUCCESS || g_ewise_add_kernel == NULL) return 0;
    g_ewise_sub_kernel = clCreateKernel(g_custom_program, "ewise_sub", &err);
    if (err != CL_SUCCESS || g_ewise_sub_kernel == NULL) return 0;
    g_ewise_div_kernel = clCreateKernel(g_custom_program, "ewise_div", &err);
    if (err != CL_SUCCESS || g_ewise_div_kernel == NULL) return 0;
    g_scalar_mul_kernel = clCreateKernel(g_custom_program, "scalar_mul", &err);
    if (err != CL_SUCCESS || g_scalar_mul_kernel == NULL) return 0;
    g_broadcast_sweep_kernel = clCreateKernel(g_custom_program, "broadcast_sweep", &err);
    if (err != CL_SUCCESS || g_broadcast_sweep_kernel == NULL) return 0;
    g_row_sum_kernel = clCreateKernel(g_custom_program, "row_sum", &err);
    if (err != CL_SUCCESS || g_row_sum_kernel == NULL) return 0;
    g_col_sum_kernel = clCreateKernel(g_custom_program, "col_sum", &err);
    if (err != CL_SUCCESS || g_col_sum_kernel == NULL) return 0;
    g_chol_panel_kernel = clCreateKernel(g_custom_program, "chol_panel_upper", &err);
    if (err != CL_SUCCESS || g_chol_panel_kernel == NULL) return 0;
    g_spmm_csr_kernel = clCreateKernel(g_custom_program, "sparse_spmm_csr", &err);
    if (err != CL_SUCCESS || g_spmm_csr_kernel == NULL) return 0;
    g_spmm_csc_trans_kernel = clCreateKernel(g_custom_program, "sparse_spmm_csc_trans", &err);
    if (err != CL_SUCCESS || g_spmm_csc_trans_kernel == NULL) return 0;
  }

  return 1;
}

static int amatrix_opencl_run_opencl_kernel_1d(cl_kernel kernel, size_t global_size) {
  cl_int err = clEnqueueNDRangeKernel(g_queue, kernel, 1, NULL, &global_size, NULL, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    return 0;
  }
  return 1;
}

static int amatrix_opencl_run_opencl_kernel_2d(cl_kernel kernel, size_t global_x, size_t global_y) {
  cl_int err = CL_SUCCESS;
  size_t global[2];

  global[0] = global_x;
  global[1] = global_y;
  err = clEnqueueNDRangeKernel(g_queue, kernel, 2, NULL, global, NULL, 0, NULL, NULL);
  return err == CL_SUCCESS;
}

static int amatrix_opencl_run_sparse_spmm_device(
  const amatrix_opencl_sparse_entry *sparse,
  cl_mem rhs_buffer,
  int rhs_nrow,
  int rhs_ncol,
  int trans_lhs,
  cl_mem *out_buffer,
  int *out_nrow,
  int *out_ncol
) {
  cl_mem out = NULL;
  cl_kernel kernel = NULL;
  cl_int err = CL_SUCCESS;
  int result_nrow = trans_lhs ? sparse->ncol : sparse->nrow;
  int expected_rhs_nrow = trans_lhs ? sparse->nrow : sparse->ncol;
  int rhs_ld = rhs_nrow;

  if (!amatrix_opencl_ensure_custom_kernels() || rhs_buffer == NULL || sparse == NULL) {
    return 0;
  }
  if (!sparse->on_device) {
    return 0;
  }
  if (rhs_nrow != expected_rhs_nrow) {
    return -1;
  }

  out = amatrix_opencl_alloc_buffer((size_t)result_nrow * (size_t)rhs_ncol);
  if (out == NULL) {
    return 0;
  }

  kernel = trans_lhs ? g_spmm_csc_trans_kernel : g_spmm_csr_kernel;
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), trans_lhs ? (const void *)&sparse->csc_col_ptr_buffer : (const void *)&sparse->csr_row_ptr_buffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), trans_lhs ? (const void *)&sparse->csc_row_idx_buffer : (const void *)&sparse->csr_col_idx_buffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), trans_lhs ? (const void *)&sparse->csc_values_buffer : (const void *)&sparse->csr_values_buffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &rhs_buffer);
  err |= clSetKernelArg(kernel, 4, sizeof(int), &result_nrow);
  err |= clSetKernelArg(kernel, 5, sizeof(int), &rhs_ncol);
  err |= clSetKernelArg(kernel, 6, sizeof(int), &rhs_ld);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &out);
  if (err != CL_SUCCESS ||
      !amatrix_opencl_run_opencl_kernel_2d(kernel, (size_t)result_nrow, (size_t)rhs_ncol)) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  *out_nrow = result_nrow;
  *out_ncol = rhs_ncol;
  return 1;
}

static int amatrix_opencl_run_chol_panel_upper_inplace(cl_mem buffer, cl_mem status_buffer, int ld, int offset, int block) {
  cl_int err = CL_SUCCESS;
  size_t global_size = (size_t)block;
  size_t local_size = (size_t)block;
  float *panel = NULL;

  if (status_buffer == NULL || !amatrix_opencl_ensure_custom_kernels() || g_chol_panel_kernel == NULL) {
    panel = amatrix_opencl_panel_workspace((size_t)block * (size_t)block);
    if (panel == NULL) {
      return 0;
    }
    if (!amatrix_opencl_read_block(buffer, ld, offset, offset, block, block, panel)) {
      return 0;
    }
    if (!amatrix_opencl_host_chol_upper(panel, block)) {
      return 0;
    }
    return amatrix_opencl_write_block(buffer, ld, offset, offset, block, block, panel);
  }

  err = clSetKernelArg(g_chol_panel_kernel, 0, sizeof(cl_mem), &buffer);
  err |= clSetKernelArg(g_chol_panel_kernel, 1, sizeof(int), &ld);
  err |= clSetKernelArg(g_chol_panel_kernel, 2, sizeof(int), &offset);
  err |= clSetKernelArg(g_chol_panel_kernel, 3, sizeof(int), &block);
  err |= clSetKernelArg(g_chol_panel_kernel, 4, sizeof(cl_mem), &status_buffer);
  if (err != CL_SUCCESS ||
      clEnqueueNDRangeKernel(g_queue, g_chol_panel_kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL) != CL_SUCCESS) {
    return 0;
  }

  return 1;
}

static int amatrix_opencl_run_hadamard(cl_mem lhs, cl_mem rhs, size_t n, cl_mem *out_buffer) {
  cl_mem out = amatrix_opencl_alloc_buffer(n);
  CLBlastStatusCode status;

  if (out == NULL) {
    return 0;
  }

  status = CLBlastShad(
    n, 1.0f,
    lhs, 0, 1,
    rhs, 0, 1,
    0.0f,
    out, 0, 1,
    &g_queue, NULL
  );

  if (status != CLBlastSuccess) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  return 1;
}

static int amatrix_opencl_run_ewise(cl_mem lhs, cl_mem rhs, float scalar, int use_scalar, const char *op, size_t n, cl_mem *out_buffer) {
  cl_mem out = NULL;
  cl_kernel kernel = NULL;
  cl_int err = CL_SUCCESS;
  int n_int = (int)n;

  if (!amatrix_opencl_ensure_custom_kernels()) {
    return 0;
  }

  if (strcmp(op, "*") == 0) {
    if (use_scalar) {
      kernel = g_scalar_mul_kernel;
    } else {
      return amatrix_opencl_run_hadamard(lhs, rhs, n, out_buffer);
    }
  } else if (strcmp(op, "+") == 0) {
    kernel = g_ewise_add_kernel;
  } else if (strcmp(op, "-") == 0) {
    kernel = g_ewise_sub_kernel;
  } else if (strcmp(op, "/") == 0) {
    kernel = g_ewise_div_kernel;
  } else {
    return 0;
  }

  out = amatrix_opencl_alloc_buffer(n);
  if (out == NULL) {
    return 0;
  }

  if (kernel == g_scalar_mul_kernel) {
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &lhs);
    err |= clSetKernelArg(kernel, 1, sizeof(float), &scalar);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &out);
    err |= clSetKernelArg(kernel, 3, sizeof(int), &n_int);
  } else {
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &lhs);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &rhs);
    err |= clSetKernelArg(kernel, 2, sizeof(float), &scalar);
    err |= clSetKernelArg(kernel, 3, sizeof(int), &use_scalar);
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &out);
    err |= clSetKernelArg(kernel, 5, sizeof(int), &n_int);
  }
  if (err != CL_SUCCESS || !amatrix_opencl_run_opencl_kernel_1d(kernel, n)) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  return 1;
}

static int amatrix_opencl_run_broadcast(cl_mem lhs, int nrow, int ncol, cl_mem v, int margin, const char *op, cl_mem *out_buffer) {
  cl_mem out = NULL;
  cl_int err = CL_SUCCESS;
  int op_code = amatrix_opencl_op_code(op);
  int n = nrow * ncol;

  if (!amatrix_opencl_ensure_custom_kernels() || op_code == 0) {
    return 0;
  }

  out = amatrix_opencl_alloc_buffer((size_t)n);
  if (out == NULL) {
    return 0;
  }

  err = clSetKernelArg(g_broadcast_sweep_kernel, 0, sizeof(cl_mem), &lhs);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 1, sizeof(cl_mem), &v);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 2, sizeof(int), &nrow);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 3, sizeof(int), &ncol);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 4, sizeof(int), &margin);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 5, sizeof(int), &op_code);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 6, sizeof(cl_mem), &out);
  if (err != CL_SUCCESS || !amatrix_opencl_run_opencl_kernel_1d(g_broadcast_sweep_kernel, (size_t)n)) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  return 1;
}

static int amatrix_opencl_run_broadcast_into(cl_mem lhs, int nrow, int ncol, cl_mem v, int margin, const char *op, cl_mem out_buffer) {
  cl_int err = CL_SUCCESS;
  int op_code = amatrix_opencl_op_code(op);
  int n = nrow * ncol;

  if (!amatrix_opencl_ensure_custom_kernels() || op_code == 0 || out_buffer == NULL) {
    return 0;
  }

  err = clSetKernelArg(g_broadcast_sweep_kernel, 0, sizeof(cl_mem), &lhs);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 1, sizeof(cl_mem), &v);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 2, sizeof(int), &nrow);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 3, sizeof(int), &ncol);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 4, sizeof(int), &margin);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 5, sizeof(int), &op_code);
  err |= clSetKernelArg(g_broadcast_sweep_kernel, 6, sizeof(cl_mem), &out_buffer);
  if (err != CL_SUCCESS || !amatrix_opencl_run_opencl_kernel_1d(g_broadcast_sweep_kernel, (size_t)n)) {
    return 0;
  }

  return 1;
}

static int amatrix_opencl_run_axis_sum(cl_mem x, int nrow, int ncol, int axis, cl_mem *out_buffer, int *out_length) {
  cl_mem out = NULL;
  cl_kernel kernel = NULL;
  cl_int err = CL_SUCCESS;
  int length = 0;
  size_t global_size = 0;
  size_t local_size = (size_t)AMATRIX_OPENCL_REDUCE_WG;

  if (!amatrix_opencl_ensure_custom_kernels()) {
    return 0;
  }

  if (axis == 0) {
    kernel = g_row_sum_kernel;
    length = nrow;
  } else if (axis == 1) {
    kernel = g_col_sum_kernel;
    length = ncol;
  } else {
    return 0;
  }

  out = amatrix_opencl_alloc_buffer((size_t)length);
  if (out == NULL) {
    return 0;
  }

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &x);
  err |= clSetKernelArg(kernel, 1, sizeof(int), &nrow);
  err |= clSetKernelArg(kernel, 2, sizeof(int), &ncol);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &out);
  global_size = (size_t)length * local_size;
  if (err != CL_SUCCESS ||
      clEnqueueNDRangeKernel(g_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL) != CL_SUCCESS) {
    clReleaseMemObject(out);
    return 0;
  }

  *out_buffer = out;
  *out_length = length;
  return 1;
}

static int amatrix_opencl_ensure_sym_fill_kernel(void) {
  static const char *source =
    "__kernel void sym_fill(__global float* x, const int n, const int upper_to_lower) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  const int j = (int)get_global_id(1);\n"
    "  if (i >= n || j >= n || i == j) return;\n"
    "  if (upper_to_lower) {\n"
    "    if (i > j) x[i + j * n] = x[j + i * n];\n"
    "  } else {\n"
    "    if (i < j) x[i + j * n] = x[j + i * n];\n"
    "  }\n"
    "}\n"
    "__kernel void zero_strict_lower(__global float* x, const int n) {\n"
    "  const int i = (int)get_global_id(0);\n"
    "  const int j = (int)get_global_id(1);\n"
    "  if (i >= n || j >= n) return;\n"
    "  if (i > j) x[i + j * n] = 0.0f;\n"
    "}\n";

  if (g_sym_fill_kernel != NULL && g_zero_strict_lower_kernel != NULL) {
    return 1;
  }

  {
    cl_int err = CL_SUCCESS;
    size_t length = strlen(source);
    g_sym_fill_program = clCreateProgramWithSource(g_context, 1, &source, &length, &err);
    if (err != CL_SUCCESS || g_sym_fill_program == NULL) {
      g_sym_fill_program = NULL;
      return 0;
    }

    err = clBuildProgram(g_sym_fill_program, 1, &g_device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
      clReleaseProgram(g_sym_fill_program);
      g_sym_fill_program = NULL;
      return 0;
    }

    g_sym_fill_kernel = clCreateKernel(g_sym_fill_program, "sym_fill", &err);
    if (err != CL_SUCCESS || g_sym_fill_kernel == NULL) {
      if (g_sym_fill_kernel != NULL) {
        clReleaseKernel(g_sym_fill_kernel);
        g_sym_fill_kernel = NULL;
      }
      clReleaseProgram(g_sym_fill_program);
      g_sym_fill_program = NULL;
      return 0;
    }

    g_zero_strict_lower_kernel = clCreateKernel(g_sym_fill_program, "zero_strict_lower", &err);
    if (err != CL_SUCCESS || g_zero_strict_lower_kernel == NULL) {
      if (g_zero_strict_lower_kernel != NULL) {
        clReleaseKernel(g_zero_strict_lower_kernel);
        g_zero_strict_lower_kernel = NULL;
      }
      clReleaseKernel(g_sym_fill_kernel);
      g_sym_fill_kernel = NULL;
      clReleaseProgram(g_sym_fill_program);
      g_sym_fill_program = NULL;
      return 0;
    }
  }

  return 1;
}

static int amatrix_opencl_sym_fill_buffer(cl_mem buffer, int n, int upper_to_lower) {
  cl_int err = CL_SUCCESS;
  size_t global[2];

  if (!amatrix_opencl_ensure_sym_fill_kernel()) {
    return 0;
  }

  err = clSetKernelArg(g_sym_fill_kernel, 0, sizeof(cl_mem), &buffer);
  err |= clSetKernelArg(g_sym_fill_kernel, 1, sizeof(int), &n);
  err |= clSetKernelArg(g_sym_fill_kernel, 2, sizeof(int), &upper_to_lower);
  if (err != CL_SUCCESS) {
    return 0;
  }

  global[0] = (size_t)n;
  global[1] = (size_t)n;
  err = clEnqueueNDRangeKernel(g_queue, g_sym_fill_kernel, 2, NULL, global, NULL, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    return 0;
  }

  return 1;
}

static int amatrix_opencl_zero_strict_lower_buffer(cl_mem buffer, int n) {
  cl_int err = CL_SUCCESS;
  size_t global[2];

  if (!amatrix_opencl_ensure_sym_fill_kernel()) {
    return 0;
  }

  err = clSetKernelArg(g_zero_strict_lower_kernel, 0, sizeof(cl_mem), &buffer);
  err |= clSetKernelArg(g_zero_strict_lower_kernel, 1, sizeof(int), &n);
  if (err != CL_SUCCESS) {
    return 0;
  }

  global[0] = (size_t)n;
  global[1] = (size_t)n;
  err = clEnqueueNDRangeKernel(g_queue, g_zero_strict_lower_kernel, 2, NULL, global, NULL, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    return 0;
  }

  return 1;
}

static int amatrix_opencl_run_gemm(
  cl_mem a_buffer, int ar, int ac, int trans_a,
  cl_mem b_buffer, int br, int bc, int trans_b,
  cl_mem *out_buffer, int *out_nrow, int *out_ncol
) {
  int m = 0;
  int n = 0;
  int k = 0;
  cl_mem c_buffer = NULL;
  CLBlastStatusCode status;

  if (!amatrix_opencl_make_product_dims(ar, ac, br, bc, trans_a, trans_b, &m, &n, &k)) {
    return -1;
  }

  c_buffer = amatrix_opencl_alloc_buffer((size_t)m * (size_t)n);
  if (c_buffer == NULL) {
    return 0;
  }

  status = CLBlastSgemm(
    CLBlastLayoutColMajor,
    trans_a ? CLBlastTransposeYes : CLBlastTransposeNo,
    trans_b ? CLBlastTransposeYes : CLBlastTransposeNo,
    (size_t)m,
    (size_t)n,
    (size_t)k,
    1.0f,
    a_buffer,
    0,
    (size_t)ar,
    b_buffer,
    0,
    (size_t)br,
    0.0f,
    c_buffer,
    0,
    (size_t)m,
    &g_queue,
    NULL
  );

  if (status != CLBlastSuccess) {
    clReleaseMemObject(c_buffer);
    return 0;
  }

  *out_buffer = c_buffer;
  *out_nrow = m;
  *out_ncol = n;
  return 1;
}

static int amatrix_opencl_run_syrk(
  cl_mem a_buffer, int ar, int ac, int use_crossprod,
  cl_mem *out_buffer, int *out_nrow, int *out_ncol
) {
  int n = use_crossprod ? ac : ar;
  int k = use_crossprod ? ar : ac;
  cl_mem c_buffer = NULL;
  CLBlastTriangle triangle = use_crossprod ? CLBlastTriangleUpper : CLBlastTriangleLower;
  CLBlastTranspose transpose = use_crossprod ? CLBlastTransposeYes : CLBlastTransposeNo;
  int upper_to_lower = use_crossprod ? 1 : 0;
  CLBlastStatusCode status;

  c_buffer = amatrix_opencl_alloc_buffer((size_t)n * (size_t)n);
  if (c_buffer == NULL) {
    return 0;
  }

  status = CLBlastSsyrk(
    CLBlastLayoutColMajor,
    triangle,
    transpose,
    (size_t)n,
    (size_t)k,
    1.0f,
    a_buffer,
    0,
    (size_t)ar,
    0.0f,
    c_buffer,
    0,
    (size_t)n,
    &g_queue,
    NULL
  );

  if (status != CLBlastSuccess || !amatrix_opencl_sym_fill_buffer(c_buffer, n, upper_to_lower)) {
    clReleaseMemObject(c_buffer);
    return 0;
  }

  *out_buffer = c_buffer;
  *out_nrow = n;
  *out_ncol = n;
  return 1;
}
#endif

static SEXP amatrix_opencl_sparse_spmm_host_impl(const amatrix_opencl_sparse_entry *entry, SEXP b, int trans_lhs) {
  int b_nrow = amatrix_opencl_nrow(b);
  int b_ncol = amatrix_opencl_ncol(b);
  int out_nrow = trans_lhs ? entry->ncol : entry->nrow;
  int expected_rhs_nrow = trans_lhs ? entry->nrow : entry->ncol;
  const double *b_ptr = REAL(b);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, out_nrow, b_ncol));
  double *out_ptr = REAL(out);

  if (b_nrow != expected_rhs_nrow) {
    UNPROTECT(1);
    Rf_error("spmm: dimension mismatch");
  }

  memset(out_ptr, 0, (size_t)out_nrow * (size_t)b_ncol * sizeof(double));

  if (!trans_lhs) {
    for (int col = 0; col < b_ncol; ++col) {
      const double *b_col = b_ptr + (size_t)b_nrow * (size_t)col;
      double *out_col = out_ptr + (size_t)out_nrow * (size_t)col;
      for (int xcol = 0; xcol < entry->ncol; ++xcol) {
        double rhs = b_col[xcol];
        if (rhs == 0.0) {
          continue;
        }
        for (int sp = entry->p[xcol]; sp < entry->p[xcol + 1]; ++sp) {
          out_col[entry->i[sp]] += entry->values[sp] * rhs;
        }
      }
    }
  } else {
    for (int col = 0; col < b_ncol; ++col) {
      const double *b_col = b_ptr + (size_t)b_nrow * (size_t)col;
      double *out_col = out_ptr + (size_t)out_nrow * (size_t)col;
      for (int xcol = 0; xcol < entry->ncol; ++xcol) {
        double acc = 0.0;
        for (int sp = entry->p[xcol]; sp < entry->p[xcol + 1]; ++sp) {
          acc += entry->values[sp] * b_col[entry->i[sp]];
        }
        out_col[xcol] = acc;
      }
    }
  }

  UNPROTECT(1);
  return out;
}

static int amatrix_opencl_store_sparse_entry(const char *key, amatrix_opencl_sparse_entry *value) {
  int idx = amatrix_opencl_find_sparse_entry(key);
  if (idx < 0) {
    idx = amatrix_opencl_find_free_sparse_entry();
  }
  if (idx < 0) {
    Rf_error("sparse resident registry is full");
  }

  amatrix_opencl_commit_sparse_entry(idx, key, value);
  return 1;
}

static amatrix_opencl_sparse_entry *amatrix_opencl_sparse_entry_from_slots(
  const char *key,
  const double *values,
  int nnz,
  const int *p,
  int np,
  const int *i,
  int nrow,
  int ncol
) {
  amatrix_opencl_sparse_entry *entry = NULL;
  int *next = NULL;

  if (np != (ncol + 1)) {
    Rf_error("sparse slots have incompatible column pointer length");
  }

  entry = (amatrix_opencl_sparse_entry *)calloc(1, sizeof(amatrix_opencl_sparse_entry));
  if (entry == NULL) {
    Rf_error("failed to allocate sparse entry");
  }

  strncpy(entry->key, key, sizeof(entry->key) - 1);
  entry->key[sizeof(entry->key) - 1] = '\0';
  entry->nrow = nrow;
  entry->ncol = ncol;
  entry->nnz = nnz;
  entry->values = (double *)malloc((size_t)nnz * sizeof(double));
  entry->csr_values = (double *)malloc((size_t)nnz * sizeof(double));
  entry->p = (int *)malloc((size_t)np * sizeof(int));
  entry->i = (int *)malloc((size_t)nnz * sizeof(int));
  entry->csr_row_ptr = (int *)calloc((size_t)nrow + 1U, sizeof(int));
  entry->csr_col_idx = (int *)malloc((size_t)nnz * sizeof(int));

  if (entry->values == NULL || entry->csr_values == NULL || entry->p == NULL ||
      entry->i == NULL || entry->csr_row_ptr == NULL || entry->csr_col_idx == NULL) {
    amatrix_opencl_release_sparse_entry(entry);
    free(entry);
    Rf_error("failed to allocate sparse entry slots");
  }

  memcpy(entry->values, values, (size_t)nnz * sizeof(double));
  memcpy(entry->p, p, (size_t)np * sizeof(int));
  memcpy(entry->i, i, (size_t)nnz * sizeof(int));

  for (int idx = 0; idx < nnz; ++idx) {
    entry->csr_row_ptr[(size_t)i[idx] + 1U] += 1;
  }
  for (int row = 0; row < nrow; ++row) {
    entry->csr_row_ptr[(size_t)row + 1U] += entry->csr_row_ptr[(size_t)row];
  }

  next = (int *)malloc((size_t)nrow * sizeof(int));
  if (next == NULL) {
    amatrix_opencl_release_sparse_entry(entry);
    free(entry);
    Rf_error("failed to allocate sparse row cursor");
  }
  memcpy(next, entry->csr_row_ptr, (size_t)nrow * sizeof(int));

  for (int col = 0; col < ncol; ++col) {
    for (int sp = p[col]; sp < p[col + 1]; ++sp) {
      int row = i[sp];
      int dest = next[row]++;
      entry->csr_col_idx[dest] = col;
      entry->csr_values[dest] = values[sp];
    }
  }

  free(next);
  return entry;
}

#ifdef HAVE_OPENCL
static int amatrix_opencl_sparse_entry_upload_buffers(amatrix_opencl_sparse_entry *entry) {
  float *csr_values_f32 = NULL;
  float *csc_values_f32 = NULL;

  if (entry == NULL) {
    return 0;
  }
  if (!amatrix_opencl_try_init()) {
    return 0;
  }
  if (entry->on_device &&
      entry->csr_row_ptr_buffer != NULL &&
      entry->csr_col_idx_buffer != NULL &&
      entry->csr_values_buffer != NULL &&
      entry->csc_col_ptr_buffer != NULL &&
      entry->csc_row_idx_buffer != NULL &&
      entry->csc_values_buffer != NULL) {
    return 1;
  }

  if (entry->csr_row_ptr_buffer != NULL) {
    clReleaseMemObject(entry->csr_row_ptr_buffer);
    entry->csr_row_ptr_buffer = NULL;
  }
  if (entry->csr_col_idx_buffer != NULL) {
    clReleaseMemObject(entry->csr_col_idx_buffer);
    entry->csr_col_idx_buffer = NULL;
  }
  if (entry->csr_values_buffer != NULL) {
    clReleaseMemObject(entry->csr_values_buffer);
    entry->csr_values_buffer = NULL;
  }
  if (entry->csc_col_ptr_buffer != NULL) {
    clReleaseMemObject(entry->csc_col_ptr_buffer);
    entry->csc_col_ptr_buffer = NULL;
  }
  if (entry->csc_row_idx_buffer != NULL) {
    clReleaseMemObject(entry->csc_row_idx_buffer);
    entry->csc_row_idx_buffer = NULL;
  }
  if (entry->csc_values_buffer != NULL) {
    clReleaseMemObject(entry->csc_values_buffer);
    entry->csc_values_buffer = NULL;
  }

  csr_values_f32 = (float *)malloc((size_t)entry->nnz * sizeof(float));
  csc_values_f32 = (float *)malloc((size_t)entry->nnz * sizeof(float));
  if (csr_values_f32 == NULL || csc_values_f32 == NULL) {
    free(csr_values_f32);
    free(csc_values_f32);
    return 0;
  }

  amatrix_opencl_copy_r_to_f32(csr_values_f32, entry->csr_values, (size_t)entry->nnz);
  amatrix_opencl_copy_r_to_f32(csc_values_f32, entry->values, (size_t)entry->nnz);

  entry->csr_row_ptr_buffer = amatrix_opencl_buffer_from_host(
    entry->csr_row_ptr,
    ((size_t)entry->nrow + 1U) * sizeof(int)
  );
  entry->csr_col_idx_buffer = amatrix_opencl_buffer_from_host(
    entry->csr_col_idx,
    (size_t)entry->nnz * sizeof(int)
  );
  entry->csr_values_buffer = amatrix_opencl_buffer_from_host(
    csr_values_f32,
    (size_t)entry->nnz * sizeof(float)
  );
  entry->csc_col_ptr_buffer = amatrix_opencl_buffer_from_host(
    entry->p,
    ((size_t)entry->ncol + 1U) * sizeof(int)
  );
  entry->csc_row_idx_buffer = amatrix_opencl_buffer_from_host(
    entry->i,
    (size_t)entry->nnz * sizeof(int)
  );
  entry->csc_values_buffer = amatrix_opencl_buffer_from_host(
    csc_values_f32,
    (size_t)entry->nnz * sizeof(float)
  );

  free(csr_values_f32);
  free(csc_values_f32);

  entry->on_device =
    entry->csr_row_ptr_buffer != NULL &&
    entry->csr_col_idx_buffer != NULL &&
    entry->csr_values_buffer != NULL &&
    entry->csc_col_ptr_buffer != NULL &&
    entry->csc_row_idx_buffer != NULL &&
    entry->csc_values_buffer != NULL;

  return entry->on_device;
}
#endif

static int amatrix_opencl_store_device_buffer(const char *key, cl_mem buffer, int nrow, int ncol) {
  int idx = amatrix_opencl_find_entry(key);
  if (idx < 0) {
    idx = amatrix_opencl_find_free_entry();
  }
  if (idx < 0) {
#ifdef HAVE_OPENCL
    if (buffer != NULL) {
      clReleaseMemObject(buffer);
    }
#endif
    Rf_error("resident registry is full");
  }

  amatrix_opencl_commit_entry(idx, key, nrow, ncol, NULL, buffer, 1);
  return 1;
}
#endif

static void amatrix_opencl_store_host_entry(const char *key, SEXP value) {
  int idx = amatrix_opencl_find_entry(key);
  SEXP copy = PROTECT(Rf_duplicate(value));
  R_PreserveObject(copy);

  if (idx < 0) {
    idx = amatrix_opencl_find_free_entry();
  }
  if (idx < 0) {
    R_ReleaseObject(copy);
    UNPROTECT(1);
    Rf_error("resident registry is full");
  }

  amatrix_opencl_commit_entry(idx, key, amatrix_opencl_nrow(value), amatrix_opencl_ncol(value), copy
#ifdef HAVE_OPENCL
    , NULL, 0
#endif
  );
  UNPROTECT(1);
}

static void amatrix_opencl_store_host_vector_entry(const char *key, SEXP value) {
  int idx = amatrix_opencl_find_entry(key);
  SEXP copy = PROTECT(Rf_duplicate(value));
  R_PreserveObject(copy);

  if (idx < 0) {
    idx = amatrix_opencl_find_free_entry();
  }
  if (idx < 0) {
    R_ReleaseObject(copy);
    UNPROTECT(1);
    Rf_error("resident registry is full");
  }

  amatrix_opencl_commit_entry(idx, key, (int)XLENGTH(value), 1, copy
#ifdef HAVE_OPENCL
    , NULL, 0
#endif
  );
  UNPROTECT(1);
}

static void amatrix_opencl_store_entry(const char *key, SEXP value) {
#ifdef HAVE_OPENCL
  if (amatrix_opencl_try_init()) {
    int idx = amatrix_opencl_find_entry(key);
    cl_mem buffer = amatrix_opencl_buffer_from_r(value);

    if (buffer == NULL) {
      Rf_error("failed to allocate or upload OpenCL resident buffer");
    }

    if (idx < 0) {
      idx = amatrix_opencl_find_free_entry();
    }
    if (idx < 0) {
      clReleaseMemObject(buffer);
      Rf_error("resident registry is full");
    }

    amatrix_opencl_commit_entry(idx, key, amatrix_opencl_nrow(value), amatrix_opencl_ncol(value), NULL, buffer, 1);
    return;
  }
#endif
  amatrix_opencl_store_host_entry(key, value);
}

static SEXP amatrix_opencl_get_entry_materialized(const char *key) {
  int idx = amatrix_opencl_find_entry(key);
  amatrix_opencl_entry *entry = NULL;

  if (idx < 0) {
    Rf_error("resident key '%s' not found", key);
  }

  entry = &g_entries[idx];

#ifdef HAVE_OPENCL
  if (entry->on_device && entry->buffer != NULL) {
    if (!amatrix_opencl_try_init()) {
      Rf_error("OpenCL runtime is not available for resident key '%s'", key);
    }
    return amatrix_opencl_matrix_from_buffer(entry->buffer, entry->nrow, entry->ncol);
  }
#endif

  if (entry->host_value == NULL || entry->host_value == R_NilValue) {
    Rf_error("resident key '%s' has no materializable value", key);
  }

  return Rf_duplicate(entry->host_value);
}

static SEXP amatrix_opencl_matmul_impl(SEXP a, SEXP b, int trans_a, int trans_b) {
  int ar = amatrix_opencl_nrow(a);
  int ac = amatrix_opencl_ncol(a);
  int br = amatrix_opencl_nrow(b);
  int bc = amatrix_opencl_ncol(b);
  int m = trans_a ? ac : ar;
  int k_a = trans_a ? ar : ac;
  int k_b = trans_b ? bc : br;
  int n = trans_b ? br : bc;

  if (k_a != k_b) {
    Rf_error("non-conformable arguments");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, m, n));
  double *a_ptr = REAL(a);
  double *b_ptr = REAL(b);
  double *out_ptr = REAL(out);

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      double acc = 0.0;
      for (int k = 0; k < k_a; ++k) {
        double a_val = trans_a ? a_ptr[k + i * ar] : a_ptr[i + k * ar];
        double b_val = trans_b ? b_ptr[j + k * br] : b_ptr[k + j * br];
        acc += a_val * b_val;
      }
      out_ptr[i + j * m] = acc;
    }
  }

  UNPROTECT(1);
  return out;
}

static double amatrix_opencl_apply_binary(double lhs, double rhs, const char *op) {
  if (strcmp(op, "+") == 0) return lhs + rhs;
  if (strcmp(op, "-") == 0) return lhs - rhs;
  if (strcmp(op, "*") == 0) return lhs * rhs;
  if (strcmp(op, "/") == 0) return lhs / rhs;
  Rf_error("unsupported ewise op '%s'", op);
  return 0.0;
}

static SEXP amatrix_opencl_ewise_impl(SEXP lhs, SEXP rhs, const char *op) {
  int nr = amatrix_opencl_nrow(lhs);
  int nc = amatrix_opencl_ncol(lhs);
  R_xlen_t n = XLENGTH(lhs);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nr, nc));
  double *lhs_ptr = REAL(lhs);
  double *out_ptr = REAL(out);

  if (Rf_isNull(rhs)) {
    Rf_error("rhs cannot be NULL");
  } else if (TYPEOF(rhs) == REALSXP && Rf_isMatrix(rhs)) {
    if (amatrix_opencl_nrow(rhs) != nr || amatrix_opencl_ncol(rhs) != nc) {
      Rf_error("ewise matrix rhs must match lhs dimensions");
    }
    {
      double *rhs_ptr = REAL(rhs);
      for (R_xlen_t i = 0; i < n; ++i) {
        out_ptr[i] = amatrix_opencl_apply_binary(lhs_ptr[i], rhs_ptr[i], op);
      }
    }
  } else if (TYPEOF(rhs) == REALSXP && XLENGTH(rhs) == 1) {
    double scalar = REAL(rhs)[0];
    for (R_xlen_t i = 0; i < n; ++i) {
      out_ptr[i] = amatrix_opencl_apply_binary(lhs_ptr[i], scalar, op);
    }
  } else {
    Rf_error("rhs must be a matrix or scalar");
  }

  UNPROTECT(1);
  return out;
}

static SEXP amatrix_opencl_sum_axis_impl(SEXP x, int axis) {
  int nr = amatrix_opencl_nrow(x);
  int nc = amatrix_opencl_ncol(x);
  double *x_ptr = REAL(x);

  if (axis == 0) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, nr));
    double *out_ptr = REAL(out);
    for (int i = 0; i < nr; ++i) {
      out_ptr[i] = 0.0;
    }
    for (int j = 0; j < nc; ++j) {
      for (int i = 0; i < nr; ++i) {
        out_ptr[i] += x_ptr[i + j * nr];
      }
    }
    UNPROTECT(1);
    return out;
  }

  if (axis == 1) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, nc));
    double *out_ptr = REAL(out);
    for (int j = 0; j < nc; ++j) {
      double acc = 0.0;
      for (int i = 0; i < nr; ++i) {
        acc += x_ptr[i + j * nr];
      }
      out_ptr[j] = acc;
    }
    UNPROTECT(1);
    return out;
  }

  Rf_error("axis must be 0 (rows) or 1 (cols)");
  return R_NilValue;
}

static SEXP amatrix_opencl_broadcast_ewise_impl(SEXP lhs, SEXP v, int margin, const char *op) {
  int nr = amatrix_opencl_nrow(lhs);
  int nc = amatrix_opencl_ncol(lhs);
  double *lhs_ptr = REAL(lhs);

  if (TYPEOF(v) != REALSXP) {
    Rf_error("broadcast vector must be double");
  }
  if (margin == 1 && XLENGTH(v) != nr) {
    Rf_error("row broadcast vector length must match nrow(lhs)");
  }
  if (margin == 2 && XLENGTH(v) != nc) {
    Rf_error("column broadcast vector length must match ncol(lhs)");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nr, nc));
  double *out_ptr = REAL(out);
  double *v_ptr = REAL(v);

  for (int j = 0; j < nc; ++j) {
    for (int i = 0; i < nr; ++i) {
      double rhs = (margin == 1) ? v_ptr[i] : v_ptr[j];
      out_ptr[i + j * nr] = amatrix_opencl_apply_binary(lhs_ptr[i + j * nr], rhs, op);
    }
  }

  UNPROTECT(1);
  return out;
}

static const char *amatrix_opencl_engine_name(void) {
#ifdef HAVE_CLBLAST
  return "opencl-clblast-scaffold";
#elif defined(HAVE_OPENCL)
  return "opencl-runtime";
#else
  return "mock-c-bridge";
#endif
}

SEXP amatrix_opencl_native_available_bridge(void) {
#ifdef HAVE_OPENCL
  return ScalarLogical(amatrix_opencl_try_init());
#else
  return ScalarLogical(0);
#endif
}

SEXP amatrix_opencl_bridge_info_bridge(void) {
  static const char *names[] = {"compiled", "clblast", "native", "engine"};
  SEXP out = PROTECT(amatrix_opencl_named_list(4, names));

#ifdef HAVE_OPENCL
  SET_VECTOR_ELT(out, 0, ScalarLogical(1));
#else
  SET_VECTOR_ELT(out, 0, ScalarLogical(0));
#endif

#ifdef HAVE_CLBLAST
  SET_VECTOR_ELT(out, 1, ScalarLogical(1));
#else
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
#endif

#ifdef HAVE_OPENCL
  SET_VECTOR_ELT(out, 2, ScalarLogical(amatrix_opencl_try_init()));
#else
  SET_VECTOR_ELT(out, 2, ScalarLogical(0));
#endif
  SET_VECTOR_ELT(out, 3, Rf_mkString(amatrix_opencl_engine_name()));
  UNPROTECT(1);
  return out;
}

SEXP amatrix_opencl_diagnostics_bridge(void) {
  static const char *names[] = {
    "compiled", "clblast", "native", "engine", "probe_enabled",
    "resident_entries", "resident_device_entries", "resident_host_entries", "device_name"
  };
  SEXP out = PROTECT(amatrix_opencl_named_list(9, names));

#ifdef HAVE_OPENCL
  SET_VECTOR_ELT(out, 0, ScalarLogical(1));
  SET_VECTOR_ELT(out, 2, ScalarLogical(amatrix_opencl_try_init()));
  SET_VECTOR_ELT(out, 4, ScalarLogical(amatrix_opencl_probe_enabled()));
  SET_VECTOR_ELT(out, 6, ScalarInteger(amatrix_opencl_device_resident_count()));
  SET_VECTOR_ELT(out, 8, Rf_mkString(g_device_name[0] != '\0' ? g_device_name : ""));
#else
  SET_VECTOR_ELT(out, 0, ScalarLogical(0));
  SET_VECTOR_ELT(out, 2, ScalarLogical(0));
  SET_VECTOR_ELT(out, 4, ScalarLogical(0));
  SET_VECTOR_ELT(out, 6, ScalarInteger(0));
  SET_VECTOR_ELT(out, 8, Rf_mkString(""));
#endif

#ifdef HAVE_CLBLAST
  SET_VECTOR_ELT(out, 1, ScalarLogical(1));
#else
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
#endif

  SET_VECTOR_ELT(out, 3, Rf_mkString(amatrix_opencl_engine_name()));
  SET_VECTOR_ELT(out, 5, ScalarInteger(amatrix_opencl_resident_count()));
  SET_VECTOR_ELT(out, 7, ScalarInteger(amatrix_opencl_host_resident_count()));
  UNPROTECT(1);
  return out;
}

SEXP amatrix_opencl_sparse_store_bridge(SEXP key, SEXP values, SEXP p, SEXP i, SEXP dim) {
  amatrix_opencl_sparse_entry *entry = NULL;
  const char *key_c = NULL;

  if (TYPEOF(key) != STRSXP || XLENGTH(key) != 1) {
    Rf_error("sparse_store: key must be a scalar character");
  }
  if (TYPEOF(values) != REALSXP || TYPEOF(p) != INTSXP || TYPEOF(i) != INTSXP ||
      TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) {
    Rf_error("sparse_store: invalid sparse slots");
  }

  key_c = CHAR(asChar(key));
  entry = amatrix_opencl_sparse_entry_from_slots(
    key_c,
    REAL(values),
    (int)XLENGTH(values),
    INTEGER(p),
    (int)XLENGTH(p),
    INTEGER(i),
    INTEGER(dim)[0],
    INTEGER(dim)[1]
  );

#ifdef HAVE_OPENCL
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_sparse_entry_upload_buffers(entry);
  }
#endif

  amatrix_opencl_store_sparse_entry(key_c, entry);
  free(entry);
  return ScalarLogical(1);
}

SEXP amatrix_opencl_sparse_has_bridge(SEXP key) {
  if (TYPEOF(key) != STRSXP || XLENGTH(key) != 1) {
    return ScalarLogical(0);
  }
  return ScalarLogical(amatrix_opencl_find_sparse_entry(CHAR(asChar(key))) >= 0);
}

SEXP amatrix_opencl_sparse_drop_bridge(SEXP key) {
  int idx = -1;

  if (TYPEOF(key) != STRSXP || XLENGTH(key) != 1) {
    return ScalarLogical(1);
  }

  idx = amatrix_opencl_find_sparse_entry(CHAR(asChar(key)));
  if (idx >= 0) {
    amatrix_opencl_release_sparse_entry(&g_sparse_entries[idx]);
    g_sparse_entries[idx].in_use = 0;
    g_sparse_entries[idx].key[0] = '\0';
  }
  return ScalarLogical(1);
}

SEXP amatrix_opencl_spmm_bridge(SEXP values, SEXP p, SEXP i, SEXP dim, SEXP b, SEXP trans_lhs) {
  amatrix_opencl_sparse_entry *entry = NULL;
  SEXP out = R_NilValue;
  int trans = asLogical(trans_lhs);

  amatrix_opencl_require_matrix(b, "b");
  if (TYPEOF(values) != REALSXP || TYPEOF(p) != INTSXP || TYPEOF(i) != INTSXP ||
      TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) {
    Rf_error("spmm: invalid sparse slots");
  }

  entry = amatrix_opencl_sparse_entry_from_slots(
    "bridge",
    REAL(values),
    (int)XLENGTH(values),
    INTEGER(p),
    (int)XLENGTH(p),
    INTEGER(i),
    INTEGER(dim)[0],
    INTEGER(dim)[1]
  );

#ifdef HAVE_OPENCL
  if (amatrix_opencl_try_init() && amatrix_opencl_sparse_entry_upload_buffers(entry)) {
    cl_mem rhs_buffer = amatrix_opencl_buffer_from_r(b);
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;
    int ok = 0;

    if (rhs_buffer != NULL) {
      ok = amatrix_opencl_run_sparse_spmm_device(
        entry,
        rhs_buffer,
        amatrix_opencl_nrow(b),
        amatrix_opencl_ncol(b),
        trans,
        &out_buffer,
        &out_nrow,
        &out_ncol
      );
    }

    if (rhs_buffer != NULL) {
      clReleaseMemObject(rhs_buffer);
    }

    if (ok > 0) {
      out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
      clReleaseMemObject(out_buffer);
      amatrix_opencl_release_sparse_entry(entry);
      free(entry);
      UNPROTECT(1);
      return out;
    }

    if (out_buffer != NULL) {
      clReleaseMemObject(out_buffer);
    }
  }
#endif

  out = PROTECT(amatrix_opencl_sparse_spmm_host_impl(entry, b, trans));
  amatrix_opencl_release_sparse_entry(entry);
  free(entry);
  UNPROTECT(1);
  return out;
}

SEXP amatrix_opencl_spmm_resident_bridge(SEXP key, SEXP b, SEXP trans_lhs) {
  amatrix_opencl_sparse_entry *entry = NULL;
  int trans = asLogical(trans_lhs);

  if (TYPEOF(key) != STRSXP || XLENGTH(key) != 1) {
    Rf_error("spmm_resident: sparse key must be a scalar character");
  }
  amatrix_opencl_require_matrix(b, "b");
  entry = amatrix_opencl_lookup_sparse_entry(CHAR(asChar(key)));

#ifdef HAVE_OPENCL
  if (amatrix_opencl_try_init() && amatrix_opencl_sparse_entry_upload_buffers(entry)) {
    cl_mem rhs_buffer = amatrix_opencl_buffer_from_r(b);
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;
    int ok = 0;

    if (rhs_buffer != NULL) {
      ok = amatrix_opencl_run_sparse_spmm_device(
        entry,
        rhs_buffer,
        amatrix_opencl_nrow(b),
        amatrix_opencl_ncol(b),
        trans,
        &out_buffer,
        &out_nrow,
        &out_ncol
      );
    }

    if (rhs_buffer != NULL) {
      clReleaseMemObject(rhs_buffer);
    }

    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (out_buffer != NULL) {
      clReleaseMemObject(out_buffer);
    }
  }
#endif

  return amatrix_opencl_sparse_spmm_host_impl(entry, b, trans);
}

SEXP amatrix_opencl_spmm_resident_key_bridge(SEXP sp_key, SEXP y_key, SEXP out_key, SEXP trans_lhs, SEXP defer) {
  amatrix_opencl_sparse_entry *sparse = NULL;
  amatrix_opencl_entry *dense = NULL;
  int trans = asLogical(trans_lhs);
  int should_defer = asLogical(defer);

  if (TYPEOF(sp_key) != STRSXP || XLENGTH(sp_key) != 1 ||
      TYPEOF(y_key) != STRSXP || XLENGTH(y_key) != 1 ||
      TYPEOF(out_key) != STRSXP || XLENGTH(out_key) != 1) {
    Rf_error("spmm_resident_key: invalid arguments");
  }

  sparse = amatrix_opencl_lookup_sparse_entry(CHAR(asChar(sp_key)));
  dense = amatrix_opencl_lookup_entry(CHAR(asChar(y_key)));

#ifdef HAVE_OPENCL
  if (amatrix_opencl_try_init() &&
      dense->on_device &&
      dense->buffer != NULL &&
      amatrix_opencl_sparse_entry_upload_buffers(sparse)) {
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;
    int ok = amatrix_opencl_run_sparse_spmm_device(
      sparse,
      dense->buffer,
      dense->nrow,
      dense->ncol,
      trans,
      &out_buffer,
      &out_nrow,
      &out_ncol
    );

    if (ok > 0) {
      amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
      if (should_defer) {
        return R_NilValue;
      }
      return amatrix_opencl_get_entry_materialized(CHAR(asChar(out_key)));
    }

    if (out_buffer != NULL) {
      clReleaseMemObject(out_buffer);
    }
  }
#endif

  {
    SEXP rhs = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(y_key))));
    SEXP host_out = PROTECT(amatrix_opencl_sparse_spmm_host_impl(sparse, rhs, trans));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), host_out);
    if (should_defer) {
      UNPROTECT(2);
      return R_NilValue;
    }
    UNPROTECT(2);
    return host_out;
  }
}

SEXP amatrix_opencl_matmul_bridge(SEXP x, SEXP y) {
  amatrix_opencl_require_matrix(x, "x");
  amatrix_opencl_require_matrix(y, "y");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem a_buffer = amatrix_opencl_buffer_from_r(x);
    cl_mem b_buffer = amatrix_opencl_buffer_from_r(y);
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;

    if (a_buffer != NULL && b_buffer != NULL &&
        amatrix_opencl_run_gemm(
          a_buffer, amatrix_opencl_nrow(x), amatrix_opencl_ncol(x), 0,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        ) > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
      clReleaseMemObject(a_buffer);
      clReleaseMemObject(b_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (a_buffer != NULL) clReleaseMemObject(a_buffer);
    if (b_buffer != NULL) clReleaseMemObject(b_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  return amatrix_opencl_matmul_impl(x, y, 0, 0);
}

SEXP amatrix_opencl_crossprod_bridge(SEXP x, SEXP y) {
  amatrix_opencl_require_matrix(x, "x");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem a_buffer = amatrix_opencl_buffer_from_r(x);
    cl_mem b_buffer = NULL;
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;
    int ok = 0;

    if (!Rf_isNull(y)) {
      amatrix_opencl_require_matrix(y, "y");
      b_buffer = amatrix_opencl_buffer_from_r(y);
      if (a_buffer != NULL && b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          a_buffer, amatrix_opencl_nrow(x), amatrix_opencl_ncol(x), 1,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }
    } else if (a_buffer != NULL) {
      ok = amatrix_opencl_run_syrk(
        a_buffer, amatrix_opencl_nrow(x), amatrix_opencl_ncol(x), 1,
        &out_buffer, &out_nrow, &out_ncol
      );
    }

    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
      clReleaseMemObject(a_buffer);
      if (b_buffer != NULL) clReleaseMemObject(b_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (a_buffer != NULL) clReleaseMemObject(a_buffer);
    if (b_buffer != NULL) clReleaseMemObject(b_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  if (Rf_isNull(y)) {
    return amatrix_opencl_matmul_impl(x, x, 1, 0);
  }
  amatrix_opencl_require_matrix(y, "y");
  return amatrix_opencl_matmul_impl(x, y, 1, 0);
}

SEXP amatrix_opencl_tcrossprod_bridge(SEXP x, SEXP y) {
  amatrix_opencl_require_matrix(x, "x");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem a_buffer = amatrix_opencl_buffer_from_r(x);
    cl_mem b_buffer = NULL;
    cl_mem out_buffer = NULL;
    int out_nrow = 0;
    int out_ncol = 0;
    int ok = 0;

    if (!Rf_isNull(y)) {
      amatrix_opencl_require_matrix(y, "y");
      b_buffer = amatrix_opencl_buffer_from_r(y);
      if (a_buffer != NULL && b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          a_buffer, amatrix_opencl_nrow(x), amatrix_opencl_ncol(x), 0,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 1,
          &out_buffer, &out_nrow, &out_ncol
        );
      }
    } else if (a_buffer != NULL) {
      ok = amatrix_opencl_run_syrk(
        a_buffer, amatrix_opencl_nrow(x), amatrix_opencl_ncol(x), 0,
        &out_buffer, &out_nrow, &out_ncol
      );
    }

    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
      clReleaseMemObject(a_buffer);
      if (b_buffer != NULL) clReleaseMemObject(b_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (a_buffer != NULL) clReleaseMemObject(a_buffer);
    if (b_buffer != NULL) clReleaseMemObject(b_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  if (Rf_isNull(y)) {
    return amatrix_opencl_matmul_impl(x, x, 0, 1);
  }
  amatrix_opencl_require_matrix(y, "y");
  return amatrix_opencl_matmul_impl(x, y, 0, 1);
}

SEXP amatrix_opencl_ewise_bridge(SEXP lhs, SEXP rhs, SEXP op) {
  amatrix_opencl_require_matrix(lhs, "lhs");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem lhs_buffer = amatrix_opencl_buffer_from_r(lhs);
    cl_mem rhs_buffer = NULL;
    cl_mem out_buffer = NULL;
    int ok = 0;
    float scalar = 0.0f;
    int use_scalar = 0;
    size_t n = (size_t)amatrix_opencl_nrow(lhs) * (size_t)amatrix_opencl_ncol(lhs);

    if (TYPEOF(rhs) == REALSXP && Rf_isMatrix(rhs)) {
      if (amatrix_opencl_nrow(rhs) != amatrix_opencl_nrow(lhs) || amatrix_opencl_ncol(rhs) != amatrix_opencl_ncol(lhs)) {
        if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
        Rf_error("ewise matrix rhs must match lhs dimensions");
      }
      rhs_buffer = amatrix_opencl_buffer_from_r(rhs);
    } else if (TYPEOF(rhs) == REALSXP && XLENGTH(rhs) == 1) {
      scalar = (float)REAL(rhs)[0];
      use_scalar = 1;
    } else {
      if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
      Rf_error("rhs must be a matrix or scalar");
    }

    if (lhs_buffer != NULL) {
      ok = amatrix_opencl_run_ewise(lhs_buffer, rhs_buffer, scalar, use_scalar, CHAR(asChar(op)), n, &out_buffer);
    }
    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, amatrix_opencl_nrow(lhs), amatrix_opencl_ncol(lhs)));
      clReleaseMemObject(lhs_buffer);
      if (rhs_buffer != NULL) clReleaseMemObject(rhs_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
    if (rhs_buffer != NULL) clReleaseMemObject(rhs_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  return amatrix_opencl_ewise_impl(lhs, rhs, CHAR(asChar(op)));
}

SEXP amatrix_opencl_broadcast_ewise_bridge(SEXP lhs, SEXP v, SEXP margin, SEXP op) {
  amatrix_opencl_require_matrix(lhs, "lhs");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem lhs_buffer = amatrix_opencl_buffer_from_r(lhs);
    cl_mem v_buffer = NULL;
    cl_mem out_buffer = NULL;
    int vec_length = 0;
    int ok = 0;
    int margin_int = INTEGER(margin)[0];
    int nrow = amatrix_opencl_nrow(lhs);
    int ncol = amatrix_opencl_ncol(lhs);

    if (TYPEOF(v) != REALSXP) {
      if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
      Rf_error("broadcast vector must be double");
    }
    if ((margin_int == 1 && XLENGTH(v) != nrow) || (margin_int == 2 && XLENGTH(v) != ncol)) {
      if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
      Rf_error("broadcast vector length must match the selected margin");
    }

    v_buffer = amatrix_opencl_buffer_from_vector(v, &vec_length);
    if (lhs_buffer != NULL && v_buffer != NULL) {
      ok = amatrix_opencl_run_broadcast(lhs_buffer, nrow, ncol, v_buffer, margin_int, CHAR(asChar(op)), &out_buffer);
    }
    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, nrow, ncol));
      clReleaseMemObject(lhs_buffer);
      clReleaseMemObject(v_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (lhs_buffer != NULL) clReleaseMemObject(lhs_buffer);
    if (v_buffer != NULL) clReleaseMemObject(v_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  return amatrix_opencl_broadcast_ewise_impl(lhs, v, INTEGER(margin)[0], CHAR(asChar(op)));
}

SEXP amatrix_opencl_sum_axis_bridge(SEXP x, SEXP axis) {
  amatrix_opencl_require_matrix(x, "x");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    cl_mem x_buffer = amatrix_opencl_buffer_from_r(x);
    cl_mem out_buffer = NULL;
    int out_length = 0;
    int ok = 0;

    if (x_buffer != NULL) {
      ok = amatrix_opencl_run_axis_sum(
        x_buffer,
        amatrix_opencl_nrow(x),
        amatrix_opencl_ncol(x),
        INTEGER(axis)[0],
        &out_buffer,
        &out_length
      );
    }

    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_vector_from_buffer(out_buffer, out_length));
      clReleaseMemObject(x_buffer);
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }

    if (x_buffer != NULL) clReleaseMemObject(x_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  return amatrix_opencl_sum_axis_impl(x, INTEGER(axis)[0]);
}

SEXP amatrix_opencl_resident_store_bridge(SEXP key, SEXP x) {
  amatrix_opencl_require_matrix(x, "x");
  amatrix_opencl_store_entry(CHAR(asChar(key)), x);
  return R_NilValue;
}

SEXP amatrix_opencl_resident_has_bridge(SEXP key) {
  return ScalarLogical(amatrix_opencl_find_entry(CHAR(asChar(key))) >= 0);
}

SEXP amatrix_opencl_resident_drop_bridge(SEXP key) {
  int idx = amatrix_opencl_find_entry(CHAR(asChar(key)));
  if (idx >= 0) {
    amatrix_opencl_entry *entry = &g_entries[idx];
    amatrix_opencl_release_entry(entry);
    entry->in_use = 0;
    entry->key[0] = '\0';
  }
  return R_NilValue;
}

SEXP amatrix_opencl_resident_materialize_bridge(SEXP key) {
  return amatrix_opencl_get_entry_materialized(CHAR(asChar(key)));
}

SEXP amatrix_opencl_chol_resident_bridge(SEXP x_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    cl_mem out_buffer = NULL;

    if (x_entry == NULL) {
      Rf_error("resident key '%s' was not found", CHAR(asChar(x_key)));
    }
    if (!x_entry->on_device || x_entry->buffer == NULL) {
      Rf_error("resident key '%s' is not device-backed", CHAR(asChar(x_key)));
    }
    if (x_entry->nrow != x_entry->ncol) {
      Rf_error("chol requires a square matrix");
    }

    if (!amatrix_opencl_copy_buffer(x_entry->buffer, (size_t)x_entry->nrow * (size_t)x_entry->ncol, &out_buffer) ||
        !amatrix_opencl_run_chol_upper_inplace(out_buffer, x_entry->nrow)) {
      if (out_buffer != NULL) {
        clReleaseMemObject(out_buffer);
      }
      Rf_error("OpenCL resident chol failed");
    }

    amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, x_entry->nrow, x_entry->ncol);
    return R_NilValue;
  }
#endif
  Rf_error("OpenCL resident chol is unavailable");
  return R_NilValue;
}

SEXP amatrix_opencl_solve_resident_bridge(SEXP a_key, SEXP b_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *a_entry = amatrix_opencl_lookup_entry(CHAR(asChar(a_key)));
    amatrix_opencl_entry *b_entry = Rf_isNull(b_key) ? NULL : amatrix_opencl_lookup_entry(CHAR(asChar(b_key)));
    cl_mem out_buffer = NULL;
    int out_ncol = 0;

    if (a_entry == NULL) {
      Rf_error("resident key '%s' was not found", CHAR(asChar(a_key)));
    }
    if (!a_entry->on_device || a_entry->buffer == NULL) {
      Rf_error("resident key '%s' is not device-backed", CHAR(asChar(a_key)));
    }
    if (a_entry->nrow != a_entry->ncol) {
      Rf_error("solve requires a square matrix");
    }
    out_ncol = a_entry->ncol;

    if (b_entry != NULL) {
      if (!b_entry->on_device || b_entry->buffer == NULL) {
        Rf_error("resident rhs key '%s' is not device-backed", CHAR(asChar(b_key)));
      }
      if (b_entry->nrow != a_entry->nrow) {
        Rf_error("solve rhs has incompatible dimensions");
      }
      out_ncol = b_entry->ncol;
    }

    if (!amatrix_opencl_run_chol_solve(
          a_entry->buffer,
          a_entry->nrow,
          b_entry == NULL ? NULL : b_entry->buffer,
          b_entry == NULL ? a_entry->ncol : b_entry->ncol,
          &out_buffer
        )) {
      if (out_buffer != NULL) {
        clReleaseMemObject(out_buffer);
      }
      Rf_error("OpenCL resident solve failed");
    }

    amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, a_entry->nrow, out_ncol);
    return R_NilValue;
  }
#endif
  Rf_error("OpenCL resident solve is unavailable");
  return R_NilValue;
}

SEXP amatrix_opencl_solve_triangular_resident_bridge(SEXP factor_key, SEXP rhs_key, SEXP out_key, SEXP lower, SEXP transpose) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *factor_entry = amatrix_opencl_lookup_entry(CHAR(asChar(factor_key)));
    amatrix_opencl_entry *rhs_entry = amatrix_opencl_lookup_entry(CHAR(asChar(rhs_key)));
    cl_mem out_buffer = NULL;
    int lower_flag = asLogical(lower);
    int transpose_flag = asLogical(transpose);

    if (factor_entry == NULL) {
      Rf_error("resident factor key '%s' was not found", CHAR(asChar(factor_key)));
    }
    if (rhs_entry == NULL) {
      Rf_error("resident rhs key '%s' was not found", CHAR(asChar(rhs_key)));
    }
    if (!factor_entry->on_device || factor_entry->buffer == NULL) {
      Rf_error("resident factor key '%s' is not device-backed", CHAR(asChar(factor_key)));
    }
    if (!rhs_entry->on_device || rhs_entry->buffer == NULL) {
      Rf_error("resident rhs key '%s' is not device-backed", CHAR(asChar(rhs_key)));
    }
    if (factor_entry->nrow != factor_entry->ncol) {
      Rf_error("triangular solve requires a square factor");
    }
    if (rhs_entry->nrow != factor_entry->nrow) {
      Rf_error("triangular solve rhs has incompatible dimensions");
    }

    if (!amatrix_opencl_run_triangular_solve(
          factor_entry->buffer,
          factor_entry->nrow,
          lower_flag,
          transpose_flag,
          rhs_entry->buffer,
          rhs_entry->ncol,
          &out_buffer
        )) {
      if (out_buffer != NULL) {
        clReleaseMemObject(out_buffer);
      }
      Rf_error("OpenCL resident triangular solve failed");
    }

    amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, rhs_entry->nrow, rhs_entry->ncol);
    return R_NilValue;
  }
#endif
  Rf_error("OpenCL resident triangular solve is unavailable");
  return R_NilValue;
}

SEXP amatrix_opencl_chol_solve_resident_bridge(SEXP factor_key, SEXP rhs_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *factor_entry = amatrix_opencl_lookup_entry(CHAR(asChar(factor_key)));
    amatrix_opencl_entry *rhs_entry = amatrix_opencl_lookup_entry(CHAR(asChar(rhs_key)));
    cl_mem out_buffer = NULL;

    if (factor_entry == NULL) {
      Rf_error("resident factor key '%s' was not found", CHAR(asChar(factor_key)));
    }
    if (rhs_entry == NULL) {
      Rf_error("resident rhs key '%s' was not found", CHAR(asChar(rhs_key)));
    }
    if (!factor_entry->on_device || factor_entry->buffer == NULL) {
      Rf_error("resident factor key '%s' is not device-backed", CHAR(asChar(factor_key)));
    }
    if (!rhs_entry->on_device || rhs_entry->buffer == NULL) {
      Rf_error("resident rhs key '%s' is not device-backed", CHAR(asChar(rhs_key)));
    }
    if (factor_entry->nrow != factor_entry->ncol) {
      Rf_error("chol_solve requires a square factor");
    }
    if (rhs_entry->nrow != factor_entry->nrow) {
      Rf_error("chol_solve rhs has incompatible dimensions");
    }

    if (!amatrix_opencl_copy_buffer(rhs_entry->buffer, (size_t)rhs_entry->nrow * (size_t)rhs_entry->ncol, &out_buffer) ||
        !amatrix_opencl_run_trsm_left(
          factor_entry->buffer, 0, (size_t)factor_entry->nrow, 0, 1,
          out_buffer, 0, (size_t)rhs_entry->nrow,
          factor_entry->nrow, rhs_entry->ncol
        ) ||
        !amatrix_opencl_run_trsm_left(
          factor_entry->buffer, 0, (size_t)factor_entry->nrow, 0, 0,
          out_buffer, 0, (size_t)rhs_entry->nrow,
          factor_entry->nrow, rhs_entry->ncol
        )) {
      if (out_buffer != NULL) {
        clReleaseMemObject(out_buffer);
      }
      Rf_error("OpenCL resident chol_solve failed");
    }

    amatrix_opencl_store_device_buffer(
      CHAR(asChar(out_key)),
      out_buffer,
      rhs_entry->nrow,
      rhs_entry->ncol
    );
    return R_NilValue;
  }
#endif
  Rf_error("OpenCL resident chol_solve is unavailable");
  return R_NilValue;
}

SEXP amatrix_opencl_matmul_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    amatrix_opencl_entry *y_entry = amatrix_opencl_lookup_entry(CHAR(asChar(y_key)));

    if (x_entry->on_device && y_entry->on_device && x_entry->buffer != NULL && y_entry->buffer != NULL) {
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = amatrix_opencl_run_gemm(
        x_entry->buffer, x_entry->nrow, x_entry->ncol, 0,
        y_entry->buffer, y_entry->nrow, y_entry->ncol, 0,
        &out_buffer, &out_nrow, &out_ncol
      );
      if (ok > 0) {
        amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
        return R_NilValue;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP y = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(y_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 0, 0));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(3);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_crossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    amatrix_opencl_entry *y_entry = Rf_isNull(y_key) ? NULL : amatrix_opencl_lookup_entry(CHAR(asChar(y_key)));

    if (x_entry->on_device && x_entry->buffer != NULL &&
        (y_entry == NULL || (y_entry->on_device && y_entry->buffer != NULL))) {
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (y_entry == NULL) {
        ok = amatrix_opencl_run_syrk(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 1,
          &out_buffer, &out_nrow, &out_ncol
        );
      } else {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 1,
          y_entry->buffer, y_entry->nrow, y_entry->ncol, 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (ok > 0) {
        amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
        return R_NilValue;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP y = PROTECT(Rf_isNull(y_key) ? Rf_duplicate(x) : amatrix_opencl_get_entry_materialized(CHAR(asChar(y_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 1, 0));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(3);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_tcrossprod_resident_bridge(SEXP x_key, SEXP y_key, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    amatrix_opencl_entry *y_entry = Rf_isNull(y_key) ? NULL : amatrix_opencl_lookup_entry(CHAR(asChar(y_key)));

    if (x_entry->on_device && x_entry->buffer != NULL &&
        (y_entry == NULL || (y_entry->on_device && y_entry->buffer != NULL))) {
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (y_entry == NULL) {
        ok = amatrix_opencl_run_syrk(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      } else {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 0,
          y_entry->buffer, y_entry->nrow, y_entry->ncol, 1,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (ok > 0) {
        amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
        return R_NilValue;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP y = PROTECT(Rf_isNull(y_key) ? Rf_duplicate(x) : amatrix_opencl_get_entry_materialized(CHAR(asChar(y_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 0, 1));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(3);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_matmul_resident_host_bridge(SEXP x_key, SEXP y) {
  amatrix_opencl_require_matrix(y, "y");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));

    if (x_entry != NULL && x_entry->on_device && x_entry->buffer != NULL) {
      cl_mem b_buffer = amatrix_opencl_buffer_from_r(y);
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 0,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (b_buffer != NULL) clReleaseMemObject(b_buffer);

      if (ok > 0) {
        SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
        clReleaseMemObject(out_buffer);
        UNPROTECT(1);
        return out;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 0, 0));
    UNPROTECT(2);
    return out;
  }
}

SEXP amatrix_opencl_matmul_resident_host_into_bridge(SEXP x_key, SEXP y, SEXP out_key) {
  amatrix_opencl_require_matrix(y, "y");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));

    if (x_entry != NULL && x_entry->on_device && x_entry->buffer != NULL) {
      cl_mem b_buffer = amatrix_opencl_buffer_from_r(y);
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 0,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (b_buffer != NULL) clReleaseMemObject(b_buffer);

      if (ok > 0) {
        amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
        return R_NilValue;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 0, 0));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(2);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_crossprod_resident_host_bridge(SEXP x_key, SEXP y) {
  amatrix_opencl_require_matrix(y, "y");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));

    if (x_entry != NULL && x_entry->on_device && x_entry->buffer != NULL) {
      cl_mem b_buffer = amatrix_opencl_buffer_from_r(y);
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 1,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (b_buffer != NULL) clReleaseMemObject(b_buffer);

      if (ok > 0) {
        SEXP out = PROTECT(amatrix_opencl_matrix_from_buffer(out_buffer, out_nrow, out_ncol));
        clReleaseMemObject(out_buffer);
        UNPROTECT(1);
        return out;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 1, 0));
    UNPROTECT(2);
    return out;
  }
}

SEXP amatrix_opencl_crossprod_resident_host_into_bridge(SEXP x_key, SEXP y, SEXP out_key) {
  amatrix_opencl_require_matrix(y, "y");
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));

    if (x_entry != NULL && x_entry->on_device && x_entry->buffer != NULL) {
      cl_mem b_buffer = amatrix_opencl_buffer_from_r(y);
      cl_mem out_buffer = NULL;
      int out_nrow = 0;
      int out_ncol = 0;
      int ok = 0;

      if (b_buffer != NULL) {
        ok = amatrix_opencl_run_gemm(
          x_entry->buffer, x_entry->nrow, x_entry->ncol, 1,
          b_buffer, amatrix_opencl_nrow(y), amatrix_opencl_ncol(y), 0,
          &out_buffer, &out_nrow, &out_ncol
        );
      }

      if (b_buffer != NULL) clReleaseMemObject(b_buffer);

      if (ok > 0) {
        amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_nrow, out_ncol);
        return R_NilValue;
      }
      if (out_buffer != NULL) clReleaseMemObject(out_buffer);
    }
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = PROTECT(amatrix_opencl_matmul_impl(x, y, 1, 0));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(2);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_ewise_resident_bridge(SEXP lhs_key, SEXP rhs, SEXP op, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *lhs_entry = amatrix_opencl_lookup_entry(CHAR(asChar(lhs_key)));
    cl_mem rhs_buffer = NULL;
    cl_mem out_buffer = NULL;
    int rhs_is_temp = 0;
    int ok = 0;
    float scalar = 0.0f;
    int use_scalar = 0;
    size_t n = (size_t)lhs_entry->nrow * (size_t)lhs_entry->ncol;

    if (lhs_entry->on_device && lhs_entry->buffer != NULL) {
      if (TYPEOF(rhs) == STRSXP && XLENGTH(rhs) == 1) {
        amatrix_opencl_entry *rhs_entry = amatrix_opencl_lookup_entry(CHAR(STRING_ELT(rhs, 0)));
        if (rhs_entry->on_device && rhs_entry->buffer != NULL) {
          rhs_buffer = rhs_entry->buffer;
          ok = amatrix_opencl_run_ewise(lhs_entry->buffer, rhs_buffer, 0.0f, 0, CHAR(asChar(op)), n, &out_buffer);
        }
      } else if (TYPEOF(rhs) == REALSXP && XLENGTH(rhs) == 1) {
        scalar = (float)REAL(rhs)[0];
        use_scalar = 1;
        ok = amatrix_opencl_run_ewise(lhs_entry->buffer, NULL, scalar, use_scalar, CHAR(asChar(op)), n, &out_buffer);
      } else if (TYPEOF(rhs) == REALSXP && Rf_isMatrix(rhs) &&
                 amatrix_opencl_nrow(rhs) == lhs_entry->nrow &&
                 amatrix_opencl_ncol(rhs) == lhs_entry->ncol) {
        rhs_buffer = amatrix_opencl_buffer_from_r(rhs);
        rhs_is_temp = 1;
        if (rhs_buffer != NULL) {
          ok = amatrix_opencl_run_ewise(lhs_entry->buffer, rhs_buffer, 0.0f, 0, CHAR(asChar(op)), n, &out_buffer);
        }
      }
    }

    if (ok > 0) {
      amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, lhs_entry->nrow, lhs_entry->ncol);
      if (rhs_is_temp && rhs_buffer != NULL) clReleaseMemObject(rhs_buffer);
      return R_NilValue;
    }

    if (rhs_is_temp && rhs_buffer != NULL) clReleaseMemObject(rhs_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  {
    SEXP lhs = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(lhs_key))));
    SEXP rhs_host = R_NilValue;
    SEXP out = NULL;

    if (TYPEOF(rhs) == STRSXP && XLENGTH(rhs) == 1) {
      rhs_host = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(STRING_ELT(rhs, 0))));
      out = PROTECT(amatrix_opencl_ewise_impl(lhs, rhs_host, CHAR(asChar(op))));
      amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
      UNPROTECT(3);
      return R_NilValue;
    }

    out = PROTECT(amatrix_opencl_ewise_impl(lhs, rhs, CHAR(asChar(op))));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(2);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_broadcast_ewise_resident_bridge(SEXP lhs_key, SEXP v, SEXP margin, SEXP op, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *lhs_entry = amatrix_opencl_lookup_entry(CHAR(asChar(lhs_key)));
    cl_mem v_buffer = NULL;
    cl_mem out_buffer = NULL;
    int vec_length = 0;
    int ok = 0;
    int release_v_buffer = 0;

    if (lhs_entry->on_device && lhs_entry->buffer != NULL) {
      v_buffer = amatrix_opencl_vector_buffer_from_arg(v, &vec_length, &release_v_buffer);
      if (v_buffer != NULL) {
        ok = amatrix_opencl_run_broadcast(
          lhs_entry->buffer,
          lhs_entry->nrow,
          lhs_entry->ncol,
          v_buffer,
          INTEGER(margin)[0],
          CHAR(asChar(op)),
          &out_buffer
        );
      }
    }

    if (ok > 0) {
      amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, lhs_entry->nrow, lhs_entry->ncol);
      if (release_v_buffer && v_buffer != NULL) clReleaseMemObject(v_buffer);
      return R_NilValue;
    }

    if (release_v_buffer && v_buffer != NULL) clReleaseMemObject(v_buffer);
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  {
    SEXP lhs = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(lhs_key))));
    SEXP v_host = PROTECT((TYPEOF(v) == STRSXP && XLENGTH(v) == 1) ? amatrix_opencl_get_entry_vector(CHAR(STRING_ELT(v, 0))) : Rf_duplicate(v));
    SEXP out = PROTECT(amatrix_opencl_broadcast_ewise_impl(lhs, v_host, INTEGER(margin)[0], CHAR(asChar(op))));
    amatrix_opencl_store_host_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(3);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_broadcast_ewise_resident_inplace_bridge(SEXP lhs_key, SEXP v, SEXP margin, SEXP op) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *lhs_entry = amatrix_opencl_lookup_entry(CHAR(asChar(lhs_key)));
    cl_mem v_buffer = NULL;
    int vec_length = 0;
    int ok = 0;
    int release_v_buffer = 0;

    if (lhs_entry->on_device && lhs_entry->buffer != NULL) {
      v_buffer = amatrix_opencl_vector_buffer_from_arg(v, &vec_length, &release_v_buffer);
      if (v_buffer != NULL) {
        ok = amatrix_opencl_run_broadcast_into(
          lhs_entry->buffer,
          lhs_entry->nrow,
          lhs_entry->ncol,
          v_buffer,
          INTEGER(margin)[0],
          CHAR(asChar(op)),
          lhs_entry->buffer
        );
      }
    }

    if (release_v_buffer && v_buffer != NULL) {
      clReleaseMemObject(v_buffer);
    }
    if (ok > 0) {
      return R_NilValue;
    }
  }
#endif
  {
    SEXP lhs = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(lhs_key))));
    SEXP v_host = PROTECT((TYPEOF(v) == STRSXP && XLENGTH(v) == 1) ? amatrix_opencl_get_entry_vector(CHAR(STRING_ELT(v, 0))) : Rf_duplicate(v));
    SEXP out = PROTECT(amatrix_opencl_broadcast_ewise_impl(lhs, v_host, INTEGER(margin)[0], CHAR(asChar(op))));
    amatrix_opencl_store_entry(CHAR(asChar(lhs_key)), out);
    UNPROTECT(3);
    return R_NilValue;
  }
}

SEXP amatrix_opencl_sum_axis_resident_bridge(SEXP x_key, SEXP axis) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    cl_mem out_buffer = NULL;
    int out_length = 0;
    int ok = 0;

    if (x_entry->on_device && x_entry->buffer != NULL) {
      ok = amatrix_opencl_run_axis_sum(
        x_entry->buffer,
        x_entry->nrow,
        x_entry->ncol,
        INTEGER(axis)[0],
        &out_buffer,
        &out_length
      );
    }
    if (ok > 0) {
      SEXP out = PROTECT(amatrix_opencl_vector_from_buffer(out_buffer, out_length));
      clReleaseMemObject(out_buffer);
      UNPROTECT(1);
      return out;
    }
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = amatrix_opencl_sum_axis_impl(x, INTEGER(axis)[0]);
    UNPROTECT(1);
    return out;
  }
}

SEXP amatrix_opencl_sum_axis_resident_key_bridge(SEXP x_key, SEXP axis, SEXP out_key) {
#ifdef HAVE_CLBLAST
  if (amatrix_opencl_try_init()) {
    amatrix_opencl_entry *x_entry = amatrix_opencl_lookup_entry(CHAR(asChar(x_key)));
    cl_mem out_buffer = NULL;
    int out_length = 0;
    int ok = 0;

    if (x_entry->on_device && x_entry->buffer != NULL) {
      ok = amatrix_opencl_run_axis_sum(
        x_entry->buffer,
        x_entry->nrow,
        x_entry->ncol,
        INTEGER(axis)[0],
        &out_buffer,
        &out_length
      );
    }
    if (ok > 0) {
      amatrix_opencl_store_device_buffer(CHAR(asChar(out_key)), out_buffer, out_length, 1);
      return R_NilValue;
    }
    if (out_buffer != NULL) clReleaseMemObject(out_buffer);
  }
#endif
  {
    SEXP x = PROTECT(amatrix_opencl_get_entry_materialized(CHAR(asChar(x_key))));
    SEXP out = PROTECT(amatrix_opencl_sum_axis_impl(x, INTEGER(axis)[0]));
    amatrix_opencl_store_host_vector_entry(CHAR(asChar(out_key)), out);
    UNPROTECT(2);
    return R_NilValue;
  }
}
