// C++ and Objective-C headers MUST precede the R headers: Rinternals.h
// (without R_NO_REMAP) defines function-like macros such as length() and
// error() that poison libc++/Foundation headers on newer SDKs
// ("too many arguments provided to function-like macro invocation" in
// libc++ __locale, seen with the macOS 26 SDK).
#include <algorithm>
#include <chrono>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef HAVE_METAL
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>

#ifdef error
#undef error
#endif

namespace {

struct SparseEntry {
  std::string key;
  int nrow = 0;
  int ncol = 0;
  int nnz = 0;
  std::vector<double> values;
  std::vector<int> p;
  std::vector<int> i;
  std::vector<int> csr_row_ptr;
  std::vector<int> csr_col_idx;
  std::vector<float> csr_values;
  std::vector<float> csc_values;
#ifdef HAVE_METAL
  id<MTLBuffer> csr_row_ptr_buffer = nil;
  id<MTLBuffer> csr_col_idx_buffer = nil;
  id<MTLBuffer> csr_values_buffer = nil;
  id<MTLBuffer> csc_col_ptr_buffer = nil;
  id<MTLBuffer> csc_row_idx_buffer = nil;
  id<MTLBuffer> csc_values_buffer = nil;
#endif

  ~SparseEntry() {
#ifdef HAVE_METAL
    if (csr_row_ptr_buffer != nil) {
      [csr_row_ptr_buffer release];
      csr_row_ptr_buffer = nil;
    }
    if (csr_col_idx_buffer != nil) {
      [csr_col_idx_buffer release];
      csr_col_idx_buffer = nil;
    }
    if (csr_values_buffer != nil) {
      [csr_values_buffer release];
      csr_values_buffer = nil;
    }
    if (csc_col_ptr_buffer != nil) {
      [csc_col_ptr_buffer release];
      csc_col_ptr_buffer = nil;
    }
    if (csc_row_idx_buffer != nil) {
      [csc_row_idx_buffer release];
      csc_row_idx_buffer = nil;
    }
    if (csc_values_buffer != nil) {
      [csc_values_buffer release];
      csc_values_buffer = nil;
    }
#endif
  }
};

struct DenseEntry {
  std::string key;
  int nrow = 0;
  int ncol = 0;
  bool host_cache_valid = false;
  std::vector<float> row_major_cache;
#ifdef HAVE_METAL
  id<MTLBuffer> buffer = nil;
  id<MTLCommandBuffer> pending_command_buffer = nil;
#endif

  ~DenseEntry() {
#ifdef HAVE_METAL
    if (pending_command_buffer != nil) {
      [pending_command_buffer release];
      pending_command_buffer = nil;
    }
    if (buffer != nil) {
      [buffer release];
      buffer = nil;
    }
#endif
  }
};

static std::unordered_map<std::string, SparseEntry*> sparse_registry;
static std::unordered_map<std::string, DenseEntry*> dense_registry;

struct ProfileCounters {
  double sparse_upload_ms = 0.0;
  long long sparse_upload_count = 0;
  long long sparse_upload_reuse_count = 0;
  double dense_upload_ms = 0.0;
  long long dense_upload_count = 0;
  double spmm_submit_ms = 0.0;
  long long spmm_submit_count = 0;
  double spmm_wait_ms = 0.0;
  long long spmm_wait_count = 0;
  double dense_sparse_submit_ms = 0.0;
  long long dense_sparse_submit_count = 0;
  double dense_sparse_wait_ms = 0.0;
  long long dense_sparse_wait_count = 0;
  double transpose_submit_ms = 0.0;
  long long transpose_submit_count = 0;
  double transpose_wait_ms = 0.0;
  long long transpose_wait_count = 0;
  double pending_wait_ms = 0.0;
  long long pending_wait_count = 0;
  double materialize_ms = 0.0;
  long long materialize_count = 0;
};

static bool g_profile_enabled = false;
static ProfileCounters g_profile;

static double profile_time_ms() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double, std::milli>(clock::now().time_since_epoch()).count();
}

static double profile_start_ms() {
  return g_profile_enabled ? profile_time_ms() : 0.0;
}

static void profile_add(double* total_ms, long long* count, double start_ms) {
  if (!g_profile_enabled) {
    return;
  }
  *total_ms += profile_time_ms() - start_ms;
  *count += 1;
}

static void profile_inc(long long* count) {
  if (!g_profile_enabled) {
    return;
  }
  *count += 1;
}

static void copy_r_to_row_major_float(float* out, const double* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[(size_t)i * (size_t)ncol + (size_t)j] = (float) in[i + nrow * j];
    }
  }
}

static void copy_row_major_float_to_r(double* out, const float* in, int nrow, int ncol) {
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      out[i + nrow * j] = (double) in[(size_t)i * (size_t)ncol + (size_t)j];
    }
  }
}

static void transpose_row_major_float(float* out, const float* in, int nrow, int ncol) {
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      out[(size_t)j * (size_t)nrow + (size_t)i] = in[(size_t)i * (size_t)ncol + (size_t)j];
    }
  }
}

static void spmm_cpu_compute(double* out,
                             const double* xdata,
                             const int* xi,
                             const int* xp,
                             int x_nrow,
                             int x_ncol,
                             const double* bdata,
                             int b_nrow,
                             int b_ncol,
                             bool trans_lhs) {
  int out_nrow = trans_lhs ? x_ncol : x_nrow;
  std::fill(out, out + (size_t)out_nrow * (size_t)b_ncol, 0.0);

  if (!trans_lhs) {
    for (int cb = 0; cb < b_ncol; ++cb) {
      const double* bcol = bdata + (size_t)x_ncol * (size_t)cb;
      double* outcol = out + (size_t)out_nrow * (size_t)cb;
      for (int col = 0; col < x_ncol; ++col) {
        double bj = bcol[col];
        if (bj == 0.0) {
          continue;
        }
        for (int sp = xp[col]; sp < xp[col + 1]; ++sp) {
          outcol[xi[sp]] += xdata[sp] * bj;
        }
      }
    }
    return;
  }

  for (int cb = 0; cb < b_ncol; ++cb) {
    const double* bcol = bdata + (size_t)b_nrow * (size_t)cb;
    double* outcol = out + (size_t)out_nrow * (size_t)cb;
    for (int col = 0; col < x_ncol; ++col) {
      double acc = 0.0;
      for (int sp = xp[col]; sp < xp[col + 1]; ++sp) {
        acc += xdata[sp] * bcol[xi[sp]];
      }
      outcol[col] = acc;
    }
  }
}

static SparseEntry* sparse_registry_find(const char* key) {
  auto it = sparse_registry.find(std::string(key));
  return it == sparse_registry.end() ? nullptr : it->second;
}

static void sparse_registry_drop(const char* key) {
  auto it = sparse_registry.find(std::string(key));
  if (it == sparse_registry.end()) {
    return;
  }
  delete it->second;
  sparse_registry.erase(it);
}

static DenseEntry* dense_registry_find(const char* key) {
  auto it = dense_registry.find(std::string(key));
  return it == dense_registry.end() ? nullptr : it->second;
}

static void dense_registry_drop(const char* key) {
  auto it = dense_registry.find(std::string(key));
  if (it == dense_registry.end()) {
    return;
  }
  delete it->second;
  dense_registry.erase(it);
}

static void dense_registry_store(const char* key, DenseEntry* entry) {
  dense_registry_drop(key);
  dense_registry[std::string(key)] = entry;
}

static SparseEntry* sparse_entry_from_slots(const char* key,
                                            const double* values,
                                            int nnz,
                                            const int* p,
                                            int np,
                                            const int* i,
                                            int nrow,
                                            int ncol) {
  SparseEntry* entry = new SparseEntry();
  entry->key = key;
  entry->nrow = nrow;
  entry->ncol = ncol;
  entry->nnz = nnz;
  entry->values.assign(values, values + nnz);
  entry->p.assign(p, p + np);
  entry->i.assign(i, i + nnz);
  entry->csr_row_ptr.assign((size_t)nrow + 1U, 0);
  entry->csr_col_idx.assign((size_t)nnz, 0);
  entry->csr_values.assign((size_t)nnz, 0.0f);
  entry->csc_values.assign((size_t)nnz, 0.0f);

  for (int idx = 0; idx < nnz; ++idx) {
    entry->csc_values[(size_t)idx] = (float) values[idx];
    entry->csr_row_ptr[(size_t)i[idx] + 1U] += 1;
  }
  for (int row = 0; row < nrow; ++row) {
    entry->csr_row_ptr[(size_t)row + 1U] += entry->csr_row_ptr[(size_t)row];
  }

  std::vector<int> next(entry->csr_row_ptr.begin(), entry->csr_row_ptr.end() - 1);
  for (int col = 0; col < ncol; ++col) {
    for (int sp = p[col]; sp < p[col + 1]; ++sp) {
      int row = i[sp];
      int dest = next[(size_t)row]++;
      entry->csr_col_idx[(size_t)dest] = col;
      entry->csr_values[(size_t)dest] = (float) values[sp];
    }
  }

  return entry;
}

static DenseEntry* dense_entry_alloc(const char* key, int nrow, int ncol, bool with_cache) {
  DenseEntry* entry = new DenseEntry();
  entry->key = key;
  entry->nrow = nrow;
  entry->ncol = ncol;
  entry->host_cache_valid = with_cache;
  if (with_cache) {
    entry->row_major_cache.assign((size_t)nrow * (size_t)ncol, 0.0f);
  }
  return entry;
}

static DenseEntry* dense_entry_from_matrix(const char* key, SEXP x_r) {
  if (!isReal(x_r) || !isMatrix(x_r)) {
    Rf_error("dense_store: expected numeric matrix");
  }

  SEXP dim = getAttrib(x_r, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  DenseEntry* entry = dense_entry_alloc(key, nrow, ncol, true);
  copy_r_to_row_major_float(entry->row_major_cache.data(), REAL(x_r), nrow, ncol);
  return entry;
}

static SEXP dense_entry_materialize_sexp(DenseEntry* entry);

#ifdef HAVE_METAL
static id<MTLDevice> g_device = nil;
static id<MTLCommandQueue> g_queue = nil;
static id<MTLComputePipelineState> g_spmm_pipeline = nil;
static id<MTLComputePipelineState> g_spmm_trans_pipeline = nil;
static id<MTLComputePipelineState> g_dense_sparse_pipeline = nil;
static id<MTLComputePipelineState> g_transpose_pipeline = nil;
static bool g_metal_initialized = false;
static bool g_metal_available = false;

static bool metal_runtime_ready();

static void dense_entry_clear_pending(DenseEntry* entry) {
  if (entry->pending_command_buffer != nil) {
    [entry->pending_command_buffer release];
    entry->pending_command_buffer = nil;
  }
}

static bool dense_entry_wait_pending(DenseEntry* entry, const char* context) {
  if (entry->pending_command_buffer == nil) {
    return true;
  }

  double start_ms = profile_start_ms();
  [entry->pending_command_buffer waitUntilCompleted];
  profile_add(&g_profile.pending_wait_ms, &g_profile.pending_wait_count, start_ms);
  if (entry->pending_command_buffer.error != nil) {
    dense_entry_clear_pending(entry);
    Rf_error("%s: Metal command buffer failed", context);
  }

  dense_entry_clear_pending(entry);
  return true;
}

static const char* metal_kernel_source =
  "#include <metal_stdlib>\n"
  "using namespace metal;\n"
  "kernel void amatrix_sparse_spmm_csr(\n"
  "  device const int* row_ptr [[buffer(0)]],\n"
  "  device const int* col_idx [[buffer(1)]],\n"
  "  device const float* values [[buffer(2)]],\n"
  "  device const float* B [[buffer(3)]],\n"
  "  device float* out [[buffer(4)]],\n"
  "  constant uint& rows [[buffer(5)]],\n"
  "  constant uint& cols [[buffer(6)]],\n"
  "  uint2 gid [[thread_position_in_grid]]) {\n"
  "    if (gid.x >= rows || gid.y >= cols) return;\n"
  "    uint row = gid.x;\n"
  "    uint col = gid.y;\n"
  "    float acc = 0.0f;\n"
  "    for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx) {\n"
  "      acc += values[idx] * B[(size_t)col_idx[idx] * cols + col];\n"
  "    }\n"
  "    out[(size_t)row * cols + col] = acc;\n"
  "  }\n"
  "kernel void amatrix_sparse_spmm_csc_trans(\n"
  "  device const int* col_ptr [[buffer(0)]],\n"
  "  device const int* row_idx [[buffer(1)]],\n"
  "  device const float* values [[buffer(2)]],\n"
  "  device const float* B [[buffer(3)]],\n"
  "  device float* out [[buffer(4)]],\n"
  "  constant uint& rows [[buffer(5)]],\n"
  "  constant uint& cols [[buffer(6)]],\n"
  "  uint2 gid [[thread_position_in_grid]]) {\n"
  "    if (gid.x >= rows || gid.y >= cols) return;\n"
  "    uint xcol = gid.x;\n"
  "    uint col = gid.y;\n"
  "    float acc = 0.0f;\n"
  "    for (int idx = col_ptr[xcol]; idx < col_ptr[xcol + 1]; ++idx) {\n"
  "      acc += values[idx] * B[(size_t)row_idx[idx] * cols + col];\n"
  "    }\n"
  "    out[(size_t)xcol * cols + col] = acc;\n"
  "  }\n"
  "kernel void amatrix_dense_sparse_matmul_csc(\n"
  "  device const int* col_ptr [[buffer(0)]],\n"
  "  device const int* row_idx [[buffer(1)]],\n"
  "  device const float* values [[buffer(2)]],\n"
  "  device const float* A [[buffer(3)]],\n"
  "  device float* out [[buffer(4)]],\n"
  "  constant uint& rows [[buffer(5)]],\n"
  "  constant uint& cols [[buffer(6)]],\n"
  "  constant uint& inner [[buffer(7)]],\n"
  "  uint2 gid [[thread_position_in_grid]]) {\n"
  "    if (gid.x >= rows || gid.y >= cols) return;\n"
  "    uint row = gid.x;\n"
  "    uint col = gid.y;\n"
  "    float acc = 0.0f;\n"
  "    for (int idx = col_ptr[col]; idx < col_ptr[col + 1]; ++idx) {\n"
  "      acc += A[(size_t)row * inner + row_idx[idx]] * values[idx];\n"
  "    }\n"
  "    out[(size_t)row * cols + col] = acc;\n"
  "  }\n"
  "kernel void amatrix_dense_transpose(\n"
  "  device const float* input [[buffer(0)]],\n"
  "  device float* output [[buffer(1)]],\n"
  "  constant uint& in_rows [[buffer(2)]],\n"
  "  constant uint& in_cols [[buffer(3)]],\n"
  "  uint2 gid [[thread_position_in_grid]]) {\n"
  "    if (gid.x >= in_cols || gid.y >= in_rows) return;\n"
  "    output[(size_t)gid.x * in_rows + gid.y] = input[(size_t)gid.y * in_cols + gid.x];\n"
  "  }\n";

static bool dense_entry_sync_cache_from_buffer(DenseEntry* entry) {
  if (entry->buffer == nil) {
    return false;
  }

  size_t len = (size_t) entry->nrow * (size_t) entry->ncol;
  if (entry->row_major_cache.size() != len) {
    entry->row_major_cache.assign(len, 0.0f);
  }
  std::memcpy(entry->row_major_cache.data(), [entry->buffer contents], len * sizeof(float));
  entry->host_cache_valid = true;
  return true;
}

static bool dense_entry_allocate_buffer(DenseEntry* entry) {
  if (!metal_runtime_ready()) {
    return false;
  }

  dense_entry_clear_pending(entry);
  if (entry->buffer != nil) {
    [entry->buffer release];
    entry->buffer = nil;
  }

  entry->buffer = [g_device newBufferWithLength:(size_t) entry->nrow * (size_t) entry->ncol * sizeof(float)
                                        options:MTLResourceStorageModeShared];
  return entry->buffer != nil;
}

static bool dense_entry_upload_buffer(DenseEntry* entry) {
  if (!metal_runtime_ready()) {
    return false;
  }

  double start_ms = profile_start_ms();
  size_t len = (size_t) entry->nrow * (size_t) entry->ncol;
  if (entry->row_major_cache.size() != len) {
    return false;
  }

  dense_entry_clear_pending(entry);
  if (entry->buffer != nil) {
    [entry->buffer release];
    entry->buffer = nil;
  }

  entry->buffer = [g_device newBufferWithBytes:entry->row_major_cache.data()
                                        length:len * sizeof(float)
                                       options:MTLResourceStorageModeShared];
  bool ok = entry->buffer != nil;
  profile_add(&g_profile.dense_upload_ms, &g_profile.dense_upload_count, start_ms);
  return ok;
}

static bool make_pipeline(id<MTLLibrary> library,
                          const char* function_name,
                          id<MTLComputePipelineState>* out_pipeline) {
  NSString* function_ns = [NSString stringWithUTF8String:function_name];
  id<MTLFunction> function = [library newFunctionWithName:function_ns];
  if (function == nil) {
    return false;
  }

  NSError* metal_error = nil;
  *out_pipeline = [g_device newComputePipelineStateWithFunction:function error:&metal_error];
  [function release];
  return *out_pipeline != nil && metal_error == nil;
}

static bool metal_runtime_ready() {
  if (g_metal_initialized) {
    return g_metal_available;
  }

  g_metal_initialized = true;
  @autoreleasepool {
    g_device = MTLCreateSystemDefaultDevice();
    if (g_device == nil) {
      g_metal_available = false;
      return false;
    }
    [g_device retain];

    g_queue = [g_device newCommandQueue];
    if (g_queue == nil) {
      g_metal_available = false;
      return false;
    }

    NSString* source = [NSString stringWithUTF8String:metal_kernel_source];
    NSError* metal_error = nil;
    id<MTLLibrary> library = [g_device newLibraryWithSource:source options:nil error:&metal_error];
    if (library == nil || metal_error != nil) {
      if (library != nil) {
        [library release];
      }
      g_metal_available = false;
      return false;
    }

    bool ok = make_pipeline(library, "amatrix_sparse_spmm_csr", &g_spmm_pipeline) &&
      make_pipeline(library, "amatrix_sparse_spmm_csc_trans", &g_spmm_trans_pipeline) &&
      make_pipeline(library, "amatrix_dense_sparse_matmul_csc", &g_dense_sparse_pipeline) &&
      make_pipeline(library, "amatrix_dense_transpose", &g_transpose_pipeline);
    [library release];
    if (!ok) {
      g_metal_available = false;
      return false;
    }

    g_metal_available = true;
    return true;
  }
}

static bool sparse_entry_upload_buffers(SparseEntry* entry) {
  if (!metal_runtime_ready()) {
    return false;
  }

  if (entry->csr_row_ptr_buffer != nil &&
      entry->csr_col_idx_buffer != nil &&
      entry->csr_values_buffer != nil &&
      entry->csc_col_ptr_buffer != nil &&
      entry->csc_row_idx_buffer != nil &&
      entry->csc_values_buffer != nil) {
    profile_inc(&g_profile.sparse_upload_reuse_count);
    return true;
  }

  double start_ms = profile_start_ms();

  if (entry->csr_row_ptr_buffer != nil) {
    [entry->csr_row_ptr_buffer release];
    entry->csr_row_ptr_buffer = nil;
  }
  if (entry->csr_col_idx_buffer != nil) {
    [entry->csr_col_idx_buffer release];
    entry->csr_col_idx_buffer = nil;
  }
  if (entry->csr_values_buffer != nil) {
    [entry->csr_values_buffer release];
    entry->csr_values_buffer = nil;
  }
  if (entry->csc_col_ptr_buffer != nil) {
    [entry->csc_col_ptr_buffer release];
    entry->csc_col_ptr_buffer = nil;
  }
  if (entry->csc_row_idx_buffer != nil) {
    [entry->csc_row_idx_buffer release];
    entry->csc_row_idx_buffer = nil;
  }
  if (entry->csc_values_buffer != nil) {
    [entry->csc_values_buffer release];
    entry->csc_values_buffer = nil;
  }

  entry->csr_row_ptr_buffer = [g_device newBufferWithBytes:entry->csr_row_ptr.data()
                                                    length:entry->csr_row_ptr.size() * sizeof(int)
                                                   options:MTLResourceStorageModeShared];
  entry->csr_col_idx_buffer = [g_device newBufferWithBytes:entry->csr_col_idx.data()
                                                    length:entry->csr_col_idx.size() * sizeof(int)
                                                   options:MTLResourceStorageModeShared];
  entry->csr_values_buffer = [g_device newBufferWithBytes:entry->csr_values.data()
                                                   length:entry->csr_values.size() * sizeof(float)
                                                  options:MTLResourceStorageModeShared];
  entry->csc_col_ptr_buffer = [g_device newBufferWithBytes:entry->p.data()
                                                    length:entry->p.size() * sizeof(int)
                                                   options:MTLResourceStorageModeShared];
  entry->csc_row_idx_buffer = [g_device newBufferWithBytes:entry->i.data()
                                                    length:entry->i.size() * sizeof(int)
                                                   options:MTLResourceStorageModeShared];
  entry->csc_values_buffer = [g_device newBufferWithBytes:entry->csc_values.data()
                                                   length:entry->csc_values.size() * sizeof(float)
                                                  options:MTLResourceStorageModeShared];

  bool ok = entry->csr_row_ptr_buffer != nil &&
    entry->csr_col_idx_buffer != nil &&
    entry->csr_values_buffer != nil &&
    entry->csc_col_ptr_buffer != nil &&
    entry->csc_row_idx_buffer != nil &&
    entry->csc_values_buffer != nil;
  profile_add(&g_profile.sparse_upload_ms, &g_profile.sparse_upload_count, start_ms);
  return ok;
}

static bool run_spmm_pipeline(const SparseEntry* sparse,
                              const DenseEntry* dense,
                              DenseEntry* out,
                              bool trans_lhs,
                              bool wait_for_completion) {
  if (!metal_runtime_ready() ||
      dense->buffer == nil ||
      out->buffer == nil) {
    return false;
  }

  id<MTLBuffer> buf0 = trans_lhs ? sparse->csc_col_ptr_buffer : sparse->csr_row_ptr_buffer;
  id<MTLBuffer> buf1 = trans_lhs ? sparse->csc_row_idx_buffer : sparse->csr_col_idx_buffer;
  id<MTLBuffer> buf2 = trans_lhs ? sparse->csc_values_buffer : sparse->csr_values_buffer;
  if (buf0 == nil || buf1 == nil || buf2 == nil) {
    return false;
  }

  uint32_t rows = (uint32_t) out->nrow;
  uint32_t cols = (uint32_t) out->ncol;
  id<MTLComputePipelineState> pipeline = trans_lhs ? g_spmm_trans_pipeline : g_spmm_pipeline;
  if (pipeline == nil) {
    return false;
  }

  @autoreleasepool {
    double submit_start_ms = profile_start_ms();
    id<MTLCommandBuffer> command_buffer = [g_queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:buf0 offset:0 atIndex:0];
    [encoder setBuffer:buf1 offset:0 atIndex:1];
    [encoder setBuffer:buf2 offset:0 atIndex:2];
    [encoder setBuffer:dense->buffer offset:0 atIndex:3];
    [encoder setBuffer:out->buffer offset:0 atIndex:4];
    [encoder setBytes:&rows length:sizeof(uint32_t) atIndex:5];
    [encoder setBytes:&cols length:sizeof(uint32_t) atIndex:6];

    MTLSize grid = MTLSizeMake((NSUInteger) out->nrow, (NSUInteger) out->ncol, 1);
    NSUInteger tg_width = 8;
    NSUInteger tg_height = 8;
    if (pipeline.maxTotalThreadsPerThreadgroup < tg_width * tg_height) {
      tg_height = std::max<NSUInteger>(1, pipeline.maxTotalThreadsPerThreadgroup / tg_width);
    }
    MTLSize threadgroup = MTLSizeMake(tg_width, tg_height, 1);
    [encoder dispatchThreads:grid threadsPerThreadgroup:threadgroup];
    [encoder endEncoding];

    [command_buffer commit];
    profile_add(&g_profile.spmm_submit_ms, &g_profile.spmm_submit_count, submit_start_ms);
    if (wait_for_completion) {
      double wait_start_ms = profile_start_ms();
      [command_buffer waitUntilCompleted];
      profile_add(&g_profile.spmm_wait_ms, &g_profile.spmm_wait_count, wait_start_ms);
      if (command_buffer.error != nil) {
        Rf_error("spmm: Metal command buffer failed");
      }
      dense_entry_clear_pending(out);
    } else {
      dense_entry_clear_pending(out);
      out->pending_command_buffer = [command_buffer retain];
    }
  }

  out->host_cache_valid = false;
  return true;
}

static bool run_dense_sparse_pipeline(const DenseEntry* dense,
                                      const SparseEntry* sparse,
                                      DenseEntry* out,
                                      bool wait_for_completion) {
  if (!metal_runtime_ready() ||
      dense->buffer == nil ||
      out->buffer == nil ||
      g_dense_sparse_pipeline == nil ||
      sparse->csc_col_ptr_buffer == nil ||
      sparse->csc_row_idx_buffer == nil ||
      sparse->csc_values_buffer == nil) {
    return false;
  }

  uint32_t rows = (uint32_t) out->nrow;
  uint32_t cols = (uint32_t) out->ncol;
  uint32_t inner = (uint32_t) dense->ncol;

  @autoreleasepool {
    double submit_start_ms = profile_start_ms();
    id<MTLCommandBuffer> command_buffer = [g_queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
    [encoder setComputePipelineState:g_dense_sparse_pipeline];
    [encoder setBuffer:sparse->csc_col_ptr_buffer offset:0 atIndex:0];
    [encoder setBuffer:sparse->csc_row_idx_buffer offset:0 atIndex:1];
    [encoder setBuffer:sparse->csc_values_buffer offset:0 atIndex:2];
    [encoder setBuffer:dense->buffer offset:0 atIndex:3];
    [encoder setBuffer:out->buffer offset:0 atIndex:4];
    [encoder setBytes:&rows length:sizeof(uint32_t) atIndex:5];
    [encoder setBytes:&cols length:sizeof(uint32_t) atIndex:6];
    [encoder setBytes:&inner length:sizeof(uint32_t) atIndex:7];

    MTLSize grid = MTLSizeMake((NSUInteger) out->nrow, (NSUInteger) out->ncol, 1);
    NSUInteger tg_width = 8;
    NSUInteger tg_height = 8;
    if (g_dense_sparse_pipeline.maxTotalThreadsPerThreadgroup < tg_width * tg_height) {
      tg_height = std::max<NSUInteger>(1, g_dense_sparse_pipeline.maxTotalThreadsPerThreadgroup / tg_width);
    }
    MTLSize threadgroup = MTLSizeMake(tg_width, tg_height, 1);
    [encoder dispatchThreads:grid threadsPerThreadgroup:threadgroup];
    [encoder endEncoding];

    [command_buffer commit];
    profile_add(&g_profile.dense_sparse_submit_ms, &g_profile.dense_sparse_submit_count, submit_start_ms);
    if (wait_for_completion) {
      double wait_start_ms = profile_start_ms();
      [command_buffer waitUntilCompleted];
      profile_add(&g_profile.dense_sparse_wait_ms, &g_profile.dense_sparse_wait_count, wait_start_ms);
      if (command_buffer.error != nil) {
        Rf_error("dense_sparse_matmul: Metal command buffer failed");
      }
      dense_entry_clear_pending(out);
    } else {
      dense_entry_clear_pending(out);
      out->pending_command_buffer = [command_buffer retain];
    }
  }

  out->host_cache_valid = false;
  return true;
}

static bool run_transpose_pipeline(const DenseEntry* input,
                                   DenseEntry* out,
                                   bool wait_for_completion) {
  if (!metal_runtime_ready() || input->buffer == nil || out->buffer == nil || g_transpose_pipeline == nil) {
    return false;
  }

  uint32_t in_rows = (uint32_t) input->nrow;
  uint32_t in_cols = (uint32_t) input->ncol;

  @autoreleasepool {
    double submit_start_ms = profile_start_ms();
    id<MTLCommandBuffer> command_buffer = [g_queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
    [encoder setComputePipelineState:g_transpose_pipeline];
    [encoder setBuffer:input->buffer offset:0 atIndex:0];
    [encoder setBuffer:out->buffer offset:0 atIndex:1];
    [encoder setBytes:&in_rows length:sizeof(uint32_t) atIndex:2];
    [encoder setBytes:&in_cols length:sizeof(uint32_t) atIndex:3];

    MTLSize grid = MTLSizeMake((NSUInteger) input->ncol, (NSUInteger) input->nrow, 1);
    NSUInteger tg_width = 8;
    NSUInteger tg_height = 8;
    if (g_transpose_pipeline.maxTotalThreadsPerThreadgroup < tg_width * tg_height) {
      tg_height = std::max<NSUInteger>(1, g_transpose_pipeline.maxTotalThreadsPerThreadgroup / tg_width);
    }
    MTLSize threadgroup = MTLSizeMake(tg_width, tg_height, 1);
    [encoder dispatchThreads:grid threadsPerThreadgroup:threadgroup];
    [encoder endEncoding];

    [command_buffer commit];
    profile_add(&g_profile.transpose_submit_ms, &g_profile.transpose_submit_count, submit_start_ms);
    if (wait_for_completion) {
      double wait_start_ms = profile_start_ms();
      [command_buffer waitUntilCompleted];
      profile_add(&g_profile.transpose_wait_ms, &g_profile.transpose_wait_count, wait_start_ms);
      if (command_buffer.error != nil) {
        Rf_error("transpose_resident: Metal command buffer failed");
      }
      dense_entry_clear_pending(out);
    } else {
      dense_entry_clear_pending(out);
      out->pending_command_buffer = [command_buffer retain];
    }
  }

  out->host_cache_valid = false;
  return true;
}
#else
static bool dense_entry_sync_cache_from_buffer(DenseEntry* entry) {
  (void) entry;
  return false;
}
#endif

static bool dense_entry_ensure_cache(DenseEntry* entry) {
  size_t len = (size_t) entry->nrow * (size_t) entry->ncol;
  if (entry->host_cache_valid && entry->row_major_cache.size() == len) {
    return true;
  }

#ifdef HAVE_METAL
  dense_entry_wait_pending(entry, "resident_materialize");
  if (dense_entry_sync_cache_from_buffer(entry)) {
    return true;
  }
#endif

  return false;
}

static SEXP dense_entry_materialize_sexp(DenseEntry* entry) {
  double start_ms = profile_start_ms();
  if (!dense_entry_ensure_cache(entry)) {
    Rf_error("resident_materialize: dense key not materializable");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, entry->nrow, entry->ncol));
  copy_row_major_float_to_r(REAL(out), entry->row_major_cache.data(), entry->nrow, entry->ncol);
  profile_add(&g_profile.materialize_ms, &g_profile.materialize_count, start_ms);
  UNPROTECT(1);
  return out;
}

static DenseEntry* dense_entry_transpose_host(const char* out_key, DenseEntry* input) {
  if (!dense_entry_ensure_cache(input)) {
    return nullptr;
  }

  DenseEntry* out = dense_entry_alloc(out_key, input->ncol, input->nrow, true);
  transpose_row_major_float(out->row_major_cache.data(),
                            input->row_major_cache.data(),
                            input->nrow,
                            input->ncol);
  return out;
}

static DenseEntry* sparse_spmm_host_entry(const char* out_key,
                                          const SparseEntry* sparse,
                                          DenseEntry* dense,
                                          bool trans_lhs) {
  if (!dense_entry_ensure_cache(dense)) {
    return nullptr;
  }

  int expected_rows = trans_lhs ? sparse->nrow : sparse->ncol;
  if (dense->nrow != expected_rows) {
    Rf_error("spmm_resident_key: dimension mismatch");
  }

  std::vector<double> b_col_major((size_t) dense->nrow * (size_t) dense->ncol);
  copy_row_major_float_to_r(b_col_major.data(), dense->row_major_cache.data(), dense->nrow, dense->ncol);

  int out_nrow = trans_lhs ? sparse->ncol : sparse->nrow;
  DenseEntry* out = dense_entry_alloc(out_key, out_nrow, dense->ncol, true);
  std::vector<double> out_col_major((size_t) out_nrow * (size_t) dense->ncol);
  spmm_cpu_compute(out_col_major.data(),
                   sparse->values.data(),
                   sparse->i.data(),
                   sparse->p.data(),
                   sparse->nrow,
                   sparse->ncol,
                   b_col_major.data(),
                   dense->nrow,
                   dense->ncol,
                   trans_lhs);
  copy_r_to_row_major_float(out->row_major_cache.data(), out_col_major.data(), out_nrow, dense->ncol);

#ifdef HAVE_METAL
  if (metal_runtime_ready()) {
    dense_entry_upload_buffer(out);
  }
#endif

  return out;
}

static DenseEntry* dense_sparse_host_entry(const char* out_key,
                                           DenseEntry* dense,
                                           const SparseEntry* sparse) {
  if (!dense_entry_ensure_cache(dense)) {
    return nullptr;
  }
  if (dense->ncol != sparse->nrow) {
    Rf_error("dense_sparse_matmul: dimension mismatch");
  }

  DenseEntry* out = dense_entry_alloc(out_key, dense->nrow, sparse->ncol, true);
  std::fill(out->row_major_cache.begin(), out->row_major_cache.end(), 0.0f);

  for (int col = 0; col < sparse->ncol; ++col) {
    for (int sp = sparse->p[(size_t) col]; sp < sparse->p[(size_t) col + 1U]; ++sp) {
      int inner = sparse->i[(size_t) sp];
      float value = sparse->csc_values[(size_t) sp];
      for (int row = 0; row < dense->nrow; ++row) {
        out->row_major_cache[(size_t) row * (size_t) sparse->ncol + (size_t) col] +=
          dense->row_major_cache[(size_t) row * (size_t) dense->ncol + (size_t) inner] * value;
      }
    }
  }

#ifdef HAVE_METAL
  if (metal_runtime_ready()) {
    dense_entry_upload_buffer(out);
  }
#endif

  return out;
}

static SEXP sparse_spmm_host(const SparseEntry* entry, SEXP b_r, bool trans_lhs) {
  DenseEntry* dense = dense_entry_from_matrix("host", b_r);
  DenseEntry* out = sparse_spmm_host_entry("host_out", entry, dense, trans_lhs);
  delete dense;
  if (out == nullptr) {
    Rf_error("spmm: failed host fallback");
  }
  SEXP out_r = PROTECT(dense_entry_materialize_sexp(out));
  delete out;
  UNPROTECT(1);
  return out_r;
}

static SEXP sparse_spmm_direct(const SparseEntry* entry, SEXP b_r, bool trans_lhs) {
  DenseEntry* dense = dense_entry_from_matrix("bridge_rhs", b_r);
#ifdef HAVE_METAL
  bool gpu_ready = metal_runtime_ready() &&
    sparse_entry_upload_buffers(const_cast<SparseEntry*>(entry)) &&
    dense_entry_upload_buffer(dense);
  if (gpu_ready) {
    DenseEntry* out = dense_entry_alloc("bridge_out",
                                        trans_lhs ? entry->ncol : entry->nrow,
                                        dense->ncol,
                                        false);
    if (dense_entry_allocate_buffer(out) &&
        run_spmm_pipeline(entry, dense, out, trans_lhs, true)) {
      SEXP out_r = PROTECT(dense_entry_materialize_sexp(out));
      delete dense;
      delete out;
      UNPROTECT(1);
      return out_r;
    }
    delete out;
  }
#endif

  DenseEntry* out = sparse_spmm_host_entry("bridge_out", entry, dense, trans_lhs);
  delete dense;
  if (out == nullptr) {
    Rf_error("spmm: failed host fallback");
  }
  SEXP out_r = PROTECT(dense_entry_materialize_sexp(out));
  delete out;
  UNPROTECT(1);
  return out_r;
}

static SEXP make_bridge_info() {
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("compiled"));
  SET_STRING_ELT(names, 1, mkChar("native"));
  SET_STRING_ELT(names, 2, mkChar("engine"));
  setAttrib(out, R_NamesSymbol, names);
#ifdef HAVE_METAL
  SET_VECTOR_ELT(out, 0, ScalarLogical(1));
  SET_VECTOR_ELT(out, 1, ScalarLogical(metal_runtime_ready() ? 1 : 0));
  SET_VECTOR_ELT(out, 2, mkString(metal_runtime_ready() ? "metal-runtime" : "metal-unavailable"));
#else
  SET_VECTOR_ELT(out, 0, ScalarLogical(0));
  SET_VECTOR_ELT(out, 1, ScalarLogical(0));
  SET_VECTOR_ELT(out, 2, mkString("mock-metal-bridge"));
#endif
  UNPROTECT(2);
  return out;
}

static void set_profile_metric(SEXP values, SEXP names, int idx, const char* name, double value) {
  REAL(values)[idx] = value;
  SET_STRING_ELT(names, idx, mkChar(name));
}

static SEXP make_profile_info() {
  const int n = 25;
  SEXP out = PROTECT(allocVector(REALSXP, n));
  SEXP names = PROTECT(allocVector(STRSXP, n));
  int idx = 0;

  set_profile_metric(out, names, idx++, "enabled", g_profile_enabled ? 1.0 : 0.0);
  set_profile_metric(out, names, idx++, "sparse_upload_ms", g_profile.sparse_upload_ms);
  set_profile_metric(out, names, idx++, "sparse_upload_count", (double) g_profile.sparse_upload_count);
  set_profile_metric(out, names, idx++, "sparse_upload_reuse_count", (double) g_profile.sparse_upload_reuse_count);
  set_profile_metric(out, names, idx++, "dense_upload_ms", g_profile.dense_upload_ms);
  set_profile_metric(out, names, idx++, "dense_upload_count", (double) g_profile.dense_upload_count);
  set_profile_metric(out, names, idx++, "spmm_submit_ms", g_profile.spmm_submit_ms);
  set_profile_metric(out, names, idx++, "spmm_submit_count", (double) g_profile.spmm_submit_count);
  set_profile_metric(out, names, idx++, "spmm_wait_ms", g_profile.spmm_wait_ms);
  set_profile_metric(out, names, idx++, "spmm_wait_count", (double) g_profile.spmm_wait_count);
  set_profile_metric(out, names, idx++, "dense_sparse_submit_ms", g_profile.dense_sparse_submit_ms);
  set_profile_metric(out, names, idx++, "dense_sparse_submit_count", (double) g_profile.dense_sparse_submit_count);
  set_profile_metric(out, names, idx++, "dense_sparse_wait_ms", g_profile.dense_sparse_wait_ms);
  set_profile_metric(out, names, idx++, "dense_sparse_wait_count", (double) g_profile.dense_sparse_wait_count);
  set_profile_metric(out, names, idx++, "transpose_submit_ms", g_profile.transpose_submit_ms);
  set_profile_metric(out, names, idx++, "transpose_submit_count", (double) g_profile.transpose_submit_count);
  set_profile_metric(out, names, idx++, "transpose_wait_ms", g_profile.transpose_wait_ms);
  set_profile_metric(out, names, idx++, "transpose_wait_count", (double) g_profile.transpose_wait_count);
  set_profile_metric(out, names, idx++, "pending_wait_ms", g_profile.pending_wait_ms);
  set_profile_metric(out, names, idx++, "pending_wait_count", (double) g_profile.pending_wait_count);
  set_profile_metric(out, names, idx++, "materialize_ms", g_profile.materialize_ms);
  set_profile_metric(out, names, idx++, "materialize_count", (double) g_profile.materialize_count);
  set_profile_metric(out, names, idx++, "sparse_resident_count", (double) sparse_registry.size());
  set_profile_metric(out, names, idx++, "dense_resident_count", (double) dense_registry.size());
  set_profile_metric(out, names, idx++, "profile_schema_version", 1.0);

  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

}  // namespace

extern "C" SEXP amatrix_metal_native_available_bridge(void) {
#ifdef HAVE_METAL
  return ScalarLogical(metal_runtime_ready() ? 1 : 0);
#else
  return ScalarLogical(0);
#endif
}

extern "C" SEXP amatrix_metal_bridge_info_bridge(void) {
  return make_bridge_info();
}

extern "C" SEXP amatrix_metal_profile_set_enabled_bridge(SEXP enabled_r) {
  g_profile_enabled = Rf_asLogical(enabled_r) == TRUE;
  return ScalarLogical(g_profile_enabled ? 1 : 0);
}

extern "C" SEXP amatrix_metal_profile_reset_bridge(void) {
  g_profile = ProfileCounters();
  return ScalarLogical(1);
}

extern "C" SEXP amatrix_metal_profile_bridge(void) {
  return make_profile_info();
}

extern "C" SEXP amatrix_metal_sparse_store_bridge(SEXP key_r, SEXP values_r, SEXP p_r, SEXP i_r, SEXP dim_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1) {
    Rf_error("sparse_store: key must be a scalar character");
  }
  if (!isReal(values_r) || TYPEOF(p_r) != INTSXP || TYPEOF(i_r) != INTSXP ||
      TYPEOF(dim_r) != INTSXP || LENGTH(dim_r) != 2) {
    Rf_error("sparse_store: invalid sparse slots");
  }

  const char* key = CHAR(asChar(key_r));
  SparseEntry* entry = sparse_entry_from_slots(
    key,
    REAL(values_r),
    LENGTH(values_r),
    INTEGER(p_r),
    LENGTH(p_r),
    INTEGER(i_r),
    INTEGER(dim_r)[0],
    INTEGER(dim_r)[1]
  );

#ifdef HAVE_METAL
  if (metal_runtime_ready()) {
    sparse_entry_upload_buffers(entry);
  }
#endif

  sparse_registry_drop(key);
  sparse_registry[std::string(key)] = entry;
  return ScalarLogical(1);
}

extern "C" SEXP amatrix_metal_sparse_has_bridge(SEXP key_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1) {
    return ScalarLogical(0);
  }
  return ScalarLogical(sparse_registry_find(CHAR(asChar(key_r))) != nullptr);
}

extern "C" SEXP amatrix_metal_sparse_drop_bridge(SEXP key_r) {
  if (isString(key_r) && LENGTH(key_r) == 1) {
    sparse_registry_drop(CHAR(asChar(key_r)));
  }
  return ScalarLogical(1);
}

extern "C" SEXP amatrix_metal_dense_store_bridge(SEXP key_r, SEXP x_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1) {
    Rf_error("dense_store: key must be a scalar character");
  }

  const char* key = CHAR(asChar(key_r));
  DenseEntry* entry = dense_entry_from_matrix(key, x_r);
#ifdef HAVE_METAL
  if (metal_runtime_ready()) {
    dense_entry_upload_buffer(entry);
  }
#endif
  dense_registry_store(key, entry);
  return ScalarLogical(1);
}

extern "C" SEXP amatrix_metal_dense_has_bridge(SEXP key_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1) {
    return ScalarLogical(0);
  }
  return ScalarLogical(dense_registry_find(CHAR(asChar(key_r))) != nullptr);
}

extern "C" SEXP amatrix_metal_dense_drop_bridge(SEXP key_r) {
  if (isString(key_r) && LENGTH(key_r) == 1) {
    dense_registry_drop(CHAR(asChar(key_r)));
  }
  return ScalarLogical(1);
}

extern "C" SEXP amatrix_metal_dense_materialize_bridge(SEXP key_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1) {
    Rf_error("resident_materialize: invalid dense key");
  }

  DenseEntry* entry = dense_registry_find(CHAR(asChar(key_r)));
  if (entry == nullptr) {
    Rf_error("resident_materialize: dense key not found");
  }
  return dense_entry_materialize_sexp(entry);
}

extern "C" SEXP amatrix_metal_spmm_bridge(SEXP values_r, SEXP p_r, SEXP i_r, SEXP dim_r, SEXP b_r, SEXP trans_lhs_r) {
  if (!isReal(values_r) || TYPEOF(p_r) != INTSXP || TYPEOF(i_r) != INTSXP ||
      TYPEOF(dim_r) != INTSXP || LENGTH(dim_r) != 2 || !isReal(b_r) || !isMatrix(b_r)) {
    Rf_error("spmm: invalid arguments");
  }

  const bool trans_lhs = asLogical(trans_lhs_r);
  SparseEntry* entry = sparse_entry_from_slots(
    "bridge",
    REAL(values_r),
    LENGTH(values_r),
    INTEGER(p_r),
    LENGTH(p_r),
    INTEGER(i_r),
    INTEGER(dim_r)[0],
    INTEGER(dim_r)[1]
  );

  SEXP out = PROTECT(sparse_spmm_direct(entry, b_r, trans_lhs));
  delete entry;
  UNPROTECT(1);
  return out;
}

extern "C" SEXP amatrix_metal_spmm_resident_bridge(SEXP key_r, SEXP b_r, SEXP trans_lhs_r) {
  if (!isString(key_r) || LENGTH(key_r) != 1 || !isReal(b_r) || !isMatrix(b_r)) {
    Rf_error("spmm_resident: invalid arguments");
  }

  const bool trans_lhs = asLogical(trans_lhs_r);
  SparseEntry* entry = sparse_registry_find(CHAR(asChar(key_r)));
  if (entry == nullptr) {
    Rf_error("spmm_resident: sparse key not found");
  }

  return sparse_spmm_direct(entry, b_r, trans_lhs);
}

extern "C" SEXP amatrix_metal_spmm_resident_key_bridge(SEXP sp_key_r,
                                                       SEXP y_key_r,
                                                       SEXP out_key_r,
                                                       SEXP trans_lhs_r,
                                                       SEXP defer_r) {
  if (!isString(sp_key_r) || LENGTH(sp_key_r) != 1 ||
      !isString(y_key_r) || LENGTH(y_key_r) != 1 ||
      !isString(out_key_r) || LENGTH(out_key_r) != 1) {
    Rf_error("spmm_resident_key: invalid arguments");
  }

  const bool trans_lhs = asLogical(trans_lhs_r);
  const bool defer = asLogical(defer_r);
  const char* sp_key = CHAR(asChar(sp_key_r));
  const char* y_key = CHAR(asChar(y_key_r));
  const char* out_key = CHAR(asChar(out_key_r));

  SparseEntry* sparse = sparse_registry_find(sp_key);
  if (sparse == nullptr) {
    Rf_error("spmm_resident_key: sparse key not found");
  }
  DenseEntry* dense = dense_registry_find(y_key);
  if (dense == nullptr) {
    Rf_error("spmm_resident_key: dense rhs key not found");
  }

  int expected_rows = trans_lhs ? sparse->nrow : sparse->ncol;
  if (dense->nrow != expected_rows) {
    Rf_error("spmm_resident_key: dimension mismatch");
  }

  DenseEntry* out = dense_entry_alloc(out_key,
                                      trans_lhs ? sparse->ncol : sparse->nrow,
                                      dense->ncol,
                                      false);

#ifdef HAVE_METAL
  bool gpu_ready = metal_runtime_ready() &&
    sparse_entry_upload_buffers(sparse) &&
    dense->buffer != nil;
  if (gpu_ready && dense_entry_allocate_buffer(out) &&
      run_spmm_pipeline(sparse, dense, out, trans_lhs, !defer)) {
    dense_registry_store(out_key, out);
    if (defer) {
      return R_NilValue;
    }
    return dense_entry_materialize_sexp(out);
  }
#endif

  delete out;
  out = sparse_spmm_host_entry(out_key, sparse, dense, trans_lhs);
  if (out == nullptr) {
    Rf_error("spmm_resident_key: failed host fallback");
  }
  dense_registry_store(out_key, out);
  if (defer) {
    return R_NilValue;
  }
  return dense_entry_materialize_sexp(out);
}

extern "C" SEXP amatrix_metal_dense_sparse_matmul_resident_key_bridge(SEXP x_key_r,
                                                                      SEXP sp_key_r,
                                                                      SEXP out_key_r,
                                                                      SEXP defer_r) {
  if (!isString(x_key_r) || LENGTH(x_key_r) != 1 ||
      !isString(sp_key_r) || LENGTH(sp_key_r) != 1 ||
      !isString(out_key_r) || LENGTH(out_key_r) != 1) {
    Rf_error("dense_sparse_matmul_resident_key: invalid arguments");
  }

  const bool defer = asLogical(defer_r);
  const char* x_key = CHAR(asChar(x_key_r));
  const char* sp_key = CHAR(asChar(sp_key_r));
  const char* out_key = CHAR(asChar(out_key_r));

  DenseEntry* dense = dense_registry_find(x_key);
  if (dense == nullptr) {
    Rf_error("dense_sparse_matmul_resident_key: dense lhs key not found");
  }
  SparseEntry* sparse = sparse_registry_find(sp_key);
  if (sparse == nullptr) {
    Rf_error("dense_sparse_matmul_resident_key: sparse key not found");
  }
  if (dense->ncol != sparse->nrow) {
    Rf_error("dense_sparse_matmul_resident_key: dimension mismatch");
  }

  DenseEntry* out = dense_entry_alloc(out_key, dense->nrow, sparse->ncol, false);

#ifdef HAVE_METAL
  bool gpu_ready = metal_runtime_ready() &&
    dense->buffer != nil &&
    sparse_entry_upload_buffers(sparse);
  if (gpu_ready && dense_entry_allocate_buffer(out) &&
      run_dense_sparse_pipeline(dense, sparse, out, !defer)) {
    dense_registry_store(out_key, out);
    if (defer) {
      return R_NilValue;
    }
    return dense_entry_materialize_sexp(out);
  }
#endif

  delete out;
  out = dense_sparse_host_entry(out_key, dense, sparse);
  if (out == nullptr) {
    Rf_error("dense_sparse_matmul_resident_key: failed host fallback");
  }
  dense_registry_store(out_key, out);
  if (defer) {
    return R_NilValue;
  }
  return dense_entry_materialize_sexp(out);
}

extern "C" SEXP amatrix_metal_transpose_resident_bridge(SEXP x_key_r, SEXP out_key_r) {
  if (!isString(x_key_r) || LENGTH(x_key_r) != 1 ||
      !isString(out_key_r) || LENGTH(out_key_r) != 1) {
    Rf_error("transpose_resident: invalid arguments");
  }

  const char* x_key = CHAR(asChar(x_key_r));
  const char* out_key = CHAR(asChar(out_key_r));

  DenseEntry* input = dense_registry_find(x_key);
  if (input == nullptr) {
    Rf_error("transpose_resident: dense key not found");
  }

  DenseEntry* out = dense_entry_alloc(out_key, input->ncol, input->nrow, false);

#ifdef HAVE_METAL
  bool gpu_ready = metal_runtime_ready() && input->buffer != nil;
  if (gpu_ready && dense_entry_allocate_buffer(out) && run_transpose_pipeline(input, out, false)) {
    dense_registry_store(out_key, out);
    return ScalarLogical(1);
  }
#endif

  delete out;
  out = dense_entry_transpose_host(out_key, input);
  if (out == nullptr) {
    Rf_error("transpose_resident: failed host fallback");
  }
#ifdef HAVE_METAL
  if (metal_runtime_ready()) {
    dense_entry_upload_buffer(out);
  }
#endif
  dense_registry_store(out_key, out);
  return ScalarLogical(1);
}
