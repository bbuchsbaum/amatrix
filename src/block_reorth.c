#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
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

SEXP amatrix_block_reorth_bridge(SEXP z, SEXP basis, SEXP return_projection) {
  const double alpha = 1.0;
  const double beta0 = 0.0;
  const double beta1 = 1.0;
  const double neg_alpha = -1.0;
  const double reorth_ratio = 0.717;
  const int inc = 1;
  SEXP z_out = R_NilValue;
  SEXP coeff = R_NilValue;
  SEXP coeff2 = R_NilValue;
  SEXP result = R_NilValue;
  int m, b, k;
  int len_z;
  int do_return_projection;
  double z_norm;
  double z_reorth_norm;

  if (!isReal(z) || !isMatrix(z)) {
    error("z must be a numeric matrix");
  }
  if (!isReal(basis) || !isMatrix(basis)) {
    error("basis must be a numeric matrix");
  }
  if (!isLogical(return_projection) || XLENGTH(return_projection) != 1) {
    error("return_projection must be a single logical");
  }

  m = INTEGER(getAttrib(z, R_DimSymbol))[0];
  b = INTEGER(getAttrib(z, R_DimSymbol))[1];
  k = INTEGER(getAttrib(basis, R_DimSymbol))[1];
  if (INTEGER(getAttrib(basis, R_DimSymbol))[0] != m) {
    error("basis and z must have the same number of rows");
  }

  do_return_projection = asLogical(return_projection);
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
