#include <arrayfire.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void fail_af(const char* stage, af_err err) {
  fprintf(stderr, "[af_qr_repro] %s failed with error code %d\n", stage, (int)err);
  exit(2);
}

static af_backend parse_backend(const char* name) {
  if (name == NULL || strcmp(name, "default") == 0) return AF_BACKEND_DEFAULT;
  if (strcmp(name, "cpu") == 0) return AF_BACKEND_CPU;
  if (strcmp(name, "opencl") == 0) return AF_BACKEND_OPENCL;
  if (strcmp(name, "cuda") == 0) return AF_BACKEND_CUDA;
  if (strcmp(name, "oneapi") == 0) return AF_BACKEND_ONEAPI;
  fprintf(stderr, "unknown backend: %s\n", name);
  exit(1);
}

int main(int argc, char** argv) {
  const char* backend_name = argc > 1 ? argv[1] : "default";
  int n = argc > 2 ? atoi(argv[2]) : 96;
  af_backend backend = parse_backend(backend_name);
  af_backend active = AF_BACKEND_DEFAULT;
  bool lapack = false;
  dim_t dims[2] = {n, n};
  size_t size = (size_t)n * (size_t)n;
  float* data = (float*)malloc(size * sizeof(float));
  float* q_host = NULL;
  float* r_host = NULL;
  af_array ax = 0, ax_t = 0, q = 0, r = 0, tau = 0, q_t = 0, r_t = 0;
  af_err err = AF_SUCCESS;

  if (data == NULL) {
    fprintf(stderr, "failed to allocate input buffer\n");
    return 1;
  }

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      data[i * n + j] = (float)((i + 1) * 0.01 + (j + 1) * 0.001);
    }
  }

  err = af_init();
  if (err != AF_SUCCESS) fail_af("af_init", err);
  if (backend != AF_BACKEND_DEFAULT) {
    err = af_set_backend(backend);
    if (err != AF_SUCCESS) fail_af("af_set_backend", err);
  }
  err = af_get_active_backend(&active);
  if (err != AF_SUCCESS) fail_af("af_get_active_backend", err);
  err = af_is_lapack_available(&lapack);
  if (err != AF_SUCCESS) fail_af("af_is_lapack_available", err);

  printf("[af_qr_repro] backend=%s active=%d n=%d lapack=%d\n", backend_name, (int)active, n, (int)lapack);
  fflush(stdout);

  err = af_create_array(&ax, data, 2, dims, f32);
  if (err != AF_SUCCESS) fail_af("af_create_array", err);
  printf("[af_qr_repro] created input\n");
  fflush(stdout);

  err = af_transpose(&ax_t, ax, false);
  if (err != AF_SUCCESS) fail_af("af_transpose", err);
  printf("[af_qr_repro] transposed input\n");
  fflush(stdout);

  err = af_qr(&q, &r, &tau, ax_t);
  if (err != AF_SUCCESS) fail_af("af_qr", err);
  printf("[af_qr_repro] qr succeeded\n");
  fflush(stdout);

  err = af_transpose(&q_t, q, false);
  if (err != AF_SUCCESS) fail_af("af_transpose(q)", err);
  err = af_transpose(&r_t, r, false);
  if (err != AF_SUCCESS) fail_af("af_transpose(r)", err);

  q_host = (float*)malloc(size * sizeof(float));
  r_host = (float*)malloc(size * sizeof(float));
  if (q_host == NULL || r_host == NULL) {
    fprintf(stderr, "failed to allocate host output buffers\n");
    return 1;
  }

  err = af_get_data_ptr(q_host, q_t);
  if (err != AF_SUCCESS) fail_af("af_get_data_ptr(q)", err);
  err = af_get_data_ptr(r_host, r_t);
  if (err != AF_SUCCESS) fail_af("af_get_data_ptr(r)", err);

  printf("[af_qr_repro] q0=%g r0=%g\n", q_host[0], r_host[0]);
  fflush(stdout);

  if (tau) af_release_array(tau);
  if (r_t) af_release_array(r_t);
  if (q_t) af_release_array(q_t);
  if (r) af_release_array(r);
  if (q) af_release_array(q);
  if (ax_t) af_release_array(ax_t);
  if (ax) af_release_array(ax);
  free(r_host);
  free(q_host);
  free(data);
  return 0;
}
