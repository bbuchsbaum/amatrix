#include <R.h>
#include <Rinternals.h>
#include <string.h>

/* O(NNZ) sparse segment sum for dgCMatrix CSC data.
 *
 * Arguments:
 *   values_r  REALSXP  — NNZ values      (dgCMatrix @x)
 *   p_r       INTSXP   — col pointers    (dgCMatrix @p, length ncol+1)
 *   i_r       INTSXP   — row indices     (dgCMatrix @i, 0-based, length NNZ)
 *   dim_r     INTSXP   — c(nrow, ncol)   (dgCMatrix @Dim)
 *   labels_r  INTSXP   — length nrow, 1-based group labels
 *   K_r       INTSXP   — number of groups (scalar)
 *
 * Returns: K × ncol double matrix (group sums).
 */
SEXP am_sparse_segment_sum_c(SEXP values_r, SEXP p_r, SEXP i_r,
                              SEXP dim_r, SEXP labels_r, SEXP K_r) {
  if (!isReal(values_r))
    error("sparse_segment_sum: values must be real");
  if (TYPEOF(i_r) != INTSXP)
    error("sparse_segment_sum: row indices must be integer");
  if (TYPEOF(p_r) != INTSXP)
    error("sparse_segment_sum: col pointers must be integer");
  if (TYPEOF(dim_r) != INTSXP || length(dim_r) != 2)
    error("sparse_segment_sum: dim must be integer[2]");
  if (TYPEOF(labels_r) != INTSXP)
    error("sparse_segment_sum: labels must be integer");

  int nrow = INTEGER(dim_r)[0];
  int ncol = INTEGER(dim_r)[1];
  int K    = asInteger(K_r);

  if (length(labels_r) != nrow)
    error("sparse_segment_sum: labels length (%d) != nrow (%d)",
          (int)length(labels_r), nrow);

  const double *xdata  = REAL(values_r);
  const int    *xi     = INTEGER(i_r);
  const int    *xp     = INTEGER(p_r);
  const int    *labels = INTEGER(labels_r);

  SEXP out_r = PROTECT(allocMatrix(REALSXP, K, ncol));
  double *res = REAL(out_r);
  memset(res, 0, (size_t)K * (size_t)ncol * sizeof(double));

  for (int j = 0; j < ncol; j++) {
    for (int sp = xp[j]; sp < xp[j + 1]; sp++) {
      int ri = xi[sp];           /* 0-based row index */
      int g  = labels[ri] - 1;   /* convert 1-based label to 0-based */
      if (g >= 0 && g < K) {
        res[g + (size_t)K * j] += xdata[sp];
      }
    }
  }

  UNPROTECT(1);
  return out_r;
}
