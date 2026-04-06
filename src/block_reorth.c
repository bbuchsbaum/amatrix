#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <string.h>
#include <math.h>

static SEXP amatrix_named_list2_local(const char* name1, SEXP value1, const char* name2, SEXP value2) {
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

static SEXP amatrix_block_reorth_impl(SEXP z, SEXP basis, int k_used, int do_return_projection) {
  const double alpha = 1.0;
  const double beta0 = 0.0;
  const double beta1 = 1.0;
  const double neg_alpha = -1.0;
  const double reorth_ratio = 0.5;
  const int inc = 1;
  SEXP z_out = R_NilValue;
  SEXP coeff = R_NilValue;
  SEXP coeff2 = R_NilValue;
  SEXP result = R_NilValue;
  int m, b, k;
  int len_z;
  double z_norm;
  double z_reorth_norm;
  int basis_cols;

  m = INTEGER(getAttrib(z, R_DimSymbol))[0];
  b = INTEGER(getAttrib(z, R_DimSymbol))[1];
  basis_cols = INTEGER(getAttrib(basis, R_DimSymbol))[1];
  k = k_used;
  if (INTEGER(getAttrib(basis, R_DimSymbol))[0] != m) {
    error("basis and z must have the same number of rows");
  }
  if (k < 0 || k > basis_cols) {
    error("basis_cols must be between 0 and ncol(basis)");
  }
  z_out = PROTECT(duplicate(z));

  if (k <= 0) {
    if (do_return_projection) {
      result = amatrix_named_list2_local("z", z_out, "coeff", R_NilValue);
      UNPROTECT(1);
      return result;
    }
    UNPROTECT(1);
    return z_out;
  }

  coeff = PROTECT(allocMatrix(REALSXP, k, b));
  memset(REAL(coeff), 0, (size_t) k * (size_t) b * sizeof(double));

  F77_CALL(dgemm)(
    "T", "N",
    &k, &b, &m,
    &alpha,
    REAL(basis), &m,
    REAL(z), &m,
    &beta0,
    REAL(coeff), &k
    FCONE FCONE
  );

  F77_CALL(dgemm)(
    "N", "N",
    &m, &b, &k,
    &neg_alpha,
    REAL(basis), &m,
    REAL(coeff), &k,
    &beta1,
    REAL(z_out), &m
    FCONE FCONE
  );

  len_z = m * b;
  z_norm = F77_CALL(dnrm2)(&len_z, REAL(z), &inc);
  z_reorth_norm = F77_CALL(dnrm2)(&len_z, REAL(z_out), &inc);

  if (R_FINITE(z_norm) && z_norm > 0.0 && z_reorth_norm <= reorth_ratio * z_norm) {
    R_xlen_t coeff_len;

    coeff2 = PROTECT(allocMatrix(REALSXP, k, b));
    memset(REAL(coeff2), 0, (size_t) k * (size_t) b * sizeof(double));

    F77_CALL(dgemm)(
      "T", "N",
      &k, &b, &m,
      &alpha,
      REAL(basis), &m,
      REAL(z_out), &m,
      &beta0,
      REAL(coeff2), &k
      FCONE FCONE
    );

    F77_CALL(dgemm)(
      "N", "N",
      &m, &b, &k,
      &neg_alpha,
      REAL(basis), &m,
      REAL(coeff2), &k,
      &beta1,
      REAL(z_out), &m
      FCONE FCONE
    );

    coeff_len = (R_xlen_t) k * (R_xlen_t) b;
    for (R_xlen_t idx = 0; idx < coeff_len; ++idx) {
      REAL(coeff)[idx] += REAL(coeff2)[idx];
    }
  }

  if (do_return_projection) {
    result = amatrix_named_list2_local("z", z_out, "coeff", coeff);
  } else {
    result = z_out;
  }

  UNPROTECT(coeff2 == R_NilValue ? 2 : 3);
  return result;
}

SEXP amatrix_block_reorth_bridge(SEXP z, SEXP basis, SEXP return_projection) {
  if (!isReal(z) || !isMatrix(z)) {
    error("z must be a numeric matrix");
  }
  if (!isReal(basis) || !isMatrix(basis)) {
    error("basis must be a numeric matrix");
  }
  if (!isLogical(return_projection) || XLENGTH(return_projection) != 1) {
    error("return_projection must be a single logical");
  }

  return amatrix_block_reorth_impl(
    z,
    basis,
    INTEGER(getAttrib(basis, R_DimSymbol))[1],
    asLogical(return_projection)
  );
}

SEXP amatrix_block_reorth_prefix_bridge(SEXP z, SEXP basis, SEXP basis_cols, SEXP return_projection) {
  if (!isReal(z) || !isMatrix(z)) {
    error("z must be a numeric matrix");
  }
  if (!isReal(basis) || !isMatrix(basis)) {
    error("basis must be a numeric matrix");
  }
  if (!isInteger(basis_cols) || XLENGTH(basis_cols) != 1) {
    error("basis_cols must be a single integer");
  }
  if (!isLogical(return_projection) || XLENGTH(return_projection) != 1) {
    error("return_projection must be a single logical");
  }

  return amatrix_block_reorth_impl(
    z,
    basis,
    INTEGER(basis_cols)[0],
    asLogical(return_projection)
  );
}

SEXP amatrix_block_thin_qr_bridge(SEXP z) {
  SEXP q = R_NilValue;
  SEXP r = R_NilValue;
  SEXP tau = R_NilValue;
  SEXP work = R_NilValue;
  SEXP work2 = R_NilValue;
  SEXP result = R_NilValue;
  int* dims = NULL;
  int m, n, k, info, lwork;
  double work_query;

  if (!isReal(z) || !isMatrix(z)) {
    error("z must be a numeric matrix");
  }

  dims = INTEGER(getAttrib(z, R_DimSymbol));
  m = dims[0];
  n = dims[1];
  k = (m < n) ? m : n;

  q = PROTECT(duplicate(z));
  tau = PROTECT(allocVector(REALSXP, k));

  lwork = -1;
  F77_CALL(dgeqrf)(&m, &n, REAL(q), &m, REAL(tau), &work_query, &lwork, &info);
  if (info != 0) {
    UNPROTECT(2);
    error("dgeqrf workspace query failed");
  }

  lwork = (int) fmax(1.0, work_query);
  work = PROTECT(allocVector(REALSXP, lwork));
  F77_CALL(dgeqrf)(&m, &n, REAL(q), &m, REAL(tau), REAL(work), &lwork, &info);
  if (info != 0) {
    UNPROTECT(3);
    error("dgeqrf failed");
  }

  r = PROTECT(allocMatrix(REALSXP, k, n));
  memset(REAL(r), 0, (size_t) k * (size_t) n * sizeof(double));
  for (int col = 0; col < n; ++col) {
    int row_max = (col < k) ? col : (k - 1);
    for (int row = 0; row <= row_max; ++row) {
      REAL(r)[row + k * col] = REAL(q)[row + m * col];
    }
  }

  lwork = -1;
  F77_CALL(dorgqr)(&m, &k, &k, REAL(q), &m, REAL(tau), &work_query, &lwork, &info);
  if (info != 0) {
    UNPROTECT(4);
    error("dorgqr workspace query failed");
  }

  lwork = (int) fmax(1.0, work_query);
  work2 = PROTECT(allocVector(REALSXP, lwork));
  F77_CALL(dorgqr)(&m, &k, &k, REAL(q), &m, REAL(tau), REAL(work2), &lwork, &info);
  if (info != 0) {
    UNPROTECT(5);
    error("dorgqr failed");
  }

  if (n > k) {
    SEXP q_thin = PROTECT(allocMatrix(REALSXP, m, k));
    for (int col = 0; col < k; ++col) {
      memcpy(
        REAL(q_thin) + ((size_t) m * (size_t) col),
        REAL(q) + ((size_t) m * (size_t) col),
        (size_t) m * sizeof(double)
      );
    }
    q = q_thin;
    result = amatrix_named_list2_local("q", q, "r", r);
    UNPROTECT(6);
    return result;
  }

  result = amatrix_named_list2_local("q", q, "r", r);
  UNPROTECT(5);
  return result;
}
