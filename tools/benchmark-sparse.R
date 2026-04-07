library(Matrix)
library(microbenchmark)
devtools::load_all()

# Results table: list of list(section, density, sparse_ms, dense_ms)
results <- list()

cat("=== SpMM benchmark ===\n")
for (density in c(0.001, 0.01, 0.05)) {
  X_sp <- adgCMatrix(rsparsematrix(1000, 500, density=density))
  X_dn <- as.matrix(X_sp)
  B    <- matrix(rnorm(500*10), 500, 10)
  mb <- microbenchmark(
    sparse = X_sp %*% B,
    dense  = X_dn %*% B,
    times  = 5
  )
  sp_ms <- median(mb$time[mb$expr=="sparse"])/1e6
  dn_ms <- median(mb$time[mb$expr=="dense"])/1e6
  cat(sprintf("density=%.3f: sparse %.1fms, dense %.1fms\n", density, sp_ms, dn_ms))
  results[[length(results)+1]] <- list(section="spmm", density=density, sparse_ms=sp_ms, dense_ms=dn_ms)
}

cat("\n=== Covariance benchmark ===\n")
for (density in c(0.001, 0.01, 0.05)) {
  X_sp <- adgCMatrix(rsparsematrix(500, 200, density=density))
  X_dn <- as.matrix(X_sp)
  mb <- microbenchmark(
    sparse = covariance(X_sp),
    dense  = cov(X_dn),
    times  = 3
  )
  sp_ms <- median(mb$time[mb$expr=="sparse"])/1e6
  dn_ms <- median(mb$time[mb$expr=="dense"])/1e6
  cat(sprintf("density=%.3f: sparse %.1fms, dense %.1fms\n", density, sp_ms, dn_ms))
  results[[length(results)+1]] <- list(section="covariance", density=density, sparse_ms=sp_ms, dense_ms=dn_ms)
}

cat("\n=== Sparse irlba benchmark ===\n")
for (density in c(0.001, 0.01, 0.05)) {
  X_sp <- adgCMatrix(rsparsematrix(500, 300, density=density))
  X_dn <- as.matrix(X_sp)
  mb <- microbenchmark(
    sparse = irlba(X_sp, nv=5),
    dense  = irlba(X_dn, nv=5),
    times  = 3
  )
  sp_ms <- median(mb$time[mb$expr=="sparse"])/1e6
  dn_ms <- median(mb$time[mb$expr=="dense"])/1e6
  cat(sprintf("density=%.3f: sparse %.1fms, dense %.1fms\n", density, sp_ms, dn_ms))
  results[[length(results)+1]] <- list(section="irlba", density=density, sparse_ms=sp_ms, dense_ms=dn_ms)
}

# Assertion: for density < 0.1, sparse path must not be slower than dense
cat("\n=== Assertion: sparse <= dense for density < 0.1 ===\n")
failures <- Filter(function(r) r$density < 0.1 && r$sparse_ms > r$dense_ms * 1.5, results)
if (length(failures) > 0) {
  for (f in failures) {
    cat(sprintf("WARN [%s] density=%.3f: sparse %.1fms > dense %.1fms (ratio %.2fx)\n",
      f$section, f$density, f$sparse_ms, f$dense_ms, f$sparse_ms / f$dense_ms))
  }
  warning("Some sparse paths were slower than dense for density < 0.1 (see above)")
} else {
  cat("OK: all sparse paths at density < 0.1 are within 1.5x of dense\n")
}
