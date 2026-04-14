.track3_seed_count <- function(default = 8L) {
  raw <- suppressWarnings(as.integer(Sys.getenv("AMATRIX_STRESS_SEEDS", unset = as.character(default))))
  if (is.na(raw) || raw < 1L) {
    default
  } else {
    raw
  }
}

.track3_spd_matrix <- function(seed, n = 6L) {
  set.seed(seed)
  z <- matrix(rnorm(n * n), nrow = n, ncol = n)
  crossprod(z) + diag(seq_len(n)) * 0.25
}

.track3_dense_case <- function(seed, n = 15L, p = 4L, k = 3L) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- matrix(rnorm(n * k), nrow = n, ncol = k)
  shift <- rnorm(p)
  perm <- sample.int(n)

  list(
    x = x,
    y = y,
    shift = shift,
    perm = perm
  )
}

.track3_qr_signs <- function(q_ref, q_cmp) {
  signs <- sign(diag(crossprod(q_ref, q_cmp)))
  signs[!is.finite(signs) | signs == 0] <- 1
  signs
}

test_that("dense workflow invariants hold across stress seeds", {
  seed_count <- .track3_seed_count()

  for (seed in seq_len(seed_count)) {
    case <- .track3_dense_case(2026041300L + seed)
    x <- case$x
    y <- case$y
    shift <- case$shift
    perm <- case$perm

    cov_ref <- covariance(x)
    cov_shifted <- covariance(sweep(x, 2L, shift, "+"))
    expect_equal(
      as.matrix(cov_shifted),
      as.matrix(cov_ref),
      tolerance = 1e-10,
      info = paste("covariance translation invariance seed", seed)
    )

    dist_ref <- dist_matrix(x, method = "euclidean")
    dist_shifted <- dist_matrix(sweep(x, 2L, shift, "+"), method = "euclidean")
    dist_perm <- dist_matrix(x[perm, , drop = FALSE], method = "euclidean")
    expect_equal(
      dist_shifted,
      dist_ref,
      tolerance = 1e-6,
      info = paste("dist translation invariance seed", seed)
    )
    expect_equal(
      dist_perm,
      dist_ref[perm, perm, drop = FALSE],
      tolerance = 1e-6,
      info = paste("dist permutation equivariance seed", seed)
    )

    kernel_ref <- kernel_matrix(x, kernel = "rbf", sigma = 0.7)
    kernel_shifted <- kernel_matrix(sweep(x, 2L, shift, "+"), kernel = "rbf", sigma = 0.7)
    expect_equal(
      kernel_shifted,
      kernel_ref,
      tolerance = 1e-5,
      info = paste("rbf translation invariance seed", seed)
    )

    lm_ref <- many_lm(x, y, method = "qr", cache = FALSE)$coefficients
    lm_perm <- many_lm(x[perm, , drop = FALSE], y[perm, , drop = FALSE], method = "qr", cache = FALSE)$coefficients
    expect_equal(
      as.matrix(lm_perm),
      as.matrix(lm_ref),
      tolerance = 1e-10,
      info = paste("many_lm row permutation invariance seed", seed)
    )

    ridge_ref <- ridge_path(x, y[, 1L], lambdas = c(0.1, 1, 10))$coef
    ridge_neg <- ridge_path(x, -y[, 1L], lambdas = c(0.1, 1, 10))$coef
    expect_equal(
      ridge_neg,
      -ridge_ref,
      tolerance = 1e-10,
      info = paste("ridge response sign symmetry seed", seed)
    )
  }
})

test_that("core algebra invariants hold across stress seeds", {
  seed_count <- .track3_seed_count()

  for (seed in seq_len(seed_count)) {
    set.seed(2026041700L + seed)
    x_host <- matrix(rnorm(24L), nrow = 6L, ncol = 4L)
    y_host <- matrix(rnorm(12L), nrow = 4L, ncol = 3L)
    z_host <- matrix(rnorm(24L), nrow = 6L, ncol = 4L)
    scale <- runif(1L, 0.25, 3)
    shift <- rnorm(1L)

    x <- adgeMatrix(x_host)
    z <- adgeMatrix(z_host)

    expect_equal(
      as.matrix(matmul(scale * x, y_host)),
      scale * (x_host %*% y_host),
      tolerance = 1e-10,
      info = paste("matmul homogeneity seed", seed)
    )

    expect_equal(
      as.matrix(crossprod(scale * x)),
      (scale^2) * crossprod(x_host),
      tolerance = 1e-10,
      info = paste("crossprod homogeneity seed", seed)
    )

    expect_equal(
      as.matrix(tcrossprod(scale * x)),
      (scale^2) * tcrossprod(x_host),
      tolerance = 1e-10,
      info = paste("tcrossprod homogeneity seed", seed)
    )

    expect_equal(
      rowSums(x + shift),
      rowSums(x_host) + ncol(x_host) * shift,
      tolerance = 1e-10,
      info = paste("rowSums shift law seed", seed)
    )

    expect_equal(
      colSums(x + shift),
      colSums(x_host) + nrow(x_host) * shift,
      tolerance = 1e-10,
      info = paste("colSums shift law seed", seed)
    )

    expect_equal(
      as.matrix((x + z) - z),
      x_host,
      tolerance = 1e-10,
      info = paste("ewise cancellation seed", seed)
    )

    expect_equal(
      as.matrix(t(t(x))),
      x_host,
      tolerance = 1e-10,
      info = paste("transpose involution seed", seed)
    )

    spd <- crossprod(x_host) + diag(ncol(x_host))
    rhs <- rnorm(ncol(x_host))
    expect_equal(
      solve(scale * as_adgeMatrix(spd), scale * rhs),
      solve(spd, rhs),
      tolerance = 1e-10,
      info = paste("solve joint scaling invariance seed", seed)
    )
  }
})

test_that("factorization and spectral invariants hold across stress seeds", {
  seed_count <- min(.track3_seed_count(default = 6L), 12L)

  for (seed in seq_len(seed_count)) {
    spd <- .track3_spd_matrix(2026041400L + seed, n = 6L)
    chol_ref <- chol_factor(as_adgeMatrix(spd))
    chol_scaled <- chol_factor(as_adgeMatrix(9 * spd))

    expect_equal(
      as.matrix(chol_scaled),
      3 * as.matrix(chol_ref),
      tolerance = 1e-5,
      info = paste("chol scaling seed", seed)
    )

    set.seed(2026041500L + seed)
    x <- matrix(rnorm(12L * 4L), nrow = 12L, ncol = 4L)
    fit <- block_lanczos(x, nv = 3L, nu = 3L, block_size = 3L, n_steps = 4L)
    fit_scaled <- block_lanczos(2 * x, nv = 3L, nu = 3L, block_size = 3L, n_steps = 4L)

    expect_equal(
      fit_scaled$d,
      2 * fit$d,
      tolerance = 1e-6,
      info = paste("block_lanczos scaling seed", seed)
    )
    expect_true(
      all(diff(fit$d) <= 0),
      info = paste("block_lanczos singular values must be ordered seed", seed)
    )
  }
})

test_that("QR, LU, triangular, and batch helper invariants hold across stress seeds", {
  seed_count <- min(.track3_seed_count(default = 6L), 12L)

  for (seed in seq_len(seed_count)) {
    set.seed(2026041800L + seed)
    x <- matrix(rnorm(10L * 4L), nrow = 10L, ncol = 4L)
    scale <- runif(1L, 0.5, 3)

    qr_ref <- am_qr(x)
    qr_scaled <- am_qr(scale * x)
    info_ref <- qr_info(qr_ref)
    info_scaled <- qr_info(qr_scaled)
    q_ref <- as.matrix(qr.Q(qr_ref))
    q_scaled <- as.matrix(qr.Q(qr_scaled))
    r_ref <- as.matrix(qr.R(qr_ref))
    r_scaled <- as.matrix(qr.R(qr_scaled))
    signs <- .track3_qr_signs(q_ref, q_scaled)
    sign_mat <- diag(signs, nrow = length(signs), ncol = length(signs))

    expect_identical(info_ref$dim, c(10L, 4L))
    expect_identical(info_scaled$dim, c(10L, 4L))
    expect_equal(info_scaled$rank, info_ref$rank, info = paste("qr rank seed", seed))
    expect_equal(
      q_scaled %*% sign_mat,
      q_ref,
      tolerance = 1e-8,
      info = paste("qr.Q positive scaling seed", seed)
    )
    expect_equal(
      sign_mat %*% r_scaled,
      scale * r_ref,
      tolerance = 1e-8,
      info = paste("qr.R positive scaling seed", seed)
    )

    a <- matrix(rnorm(25L), nrow = 5L, ncol = 5L) + diag(5L) * 4
    b <- rnorm(5L)
    expect_equal(
      lu_solve(lu_factor(scale * a), scale * b),
      lu_solve(lu_factor(a), b),
      tolerance = 1e-10,
      info = paste("lu joint scaling seed", seed)
    )

    r_tri <- chol(crossprod(matrix(rnorm(25L), nrow = 5L, ncol = 5L)) + diag(5L))
    rhs <- rnorm(5L)
    expect_equal(
      solve_triangular(scale * r_tri, scale * (r_tri %*% rhs)),
      rhs,
      tolerance = 1e-10,
      info = paste("solve_triangular joint scaling seed", seed)
    )

    batch_size <- 3L
    mats <- lapply(seq_len(batch_size), function(i) .track3_spd_matrix(2026041900L + 10L * seed + i, n = 4L))
    rhs_batch <- lapply(seq_len(batch_size), function(i) matrix(rnorm(8L), nrow = 4L, ncol = 2L))
    perm <- sample.int(batch_size)
    ref_factors <- batch_chol(mats)
    perm_factors <- batch_chol(mats[perm])
    ref_solutions <- batch_solve(ref_factors, rhs_batch)
    perm_solutions <- batch_solve(perm_factors, rhs_batch[perm])

    for (idx in seq_len(batch_size)) {
      expect_equal(
        chol_diag(perm_factors[[idx]]),
        chol_diag(ref_factors[[perm[[idx]]]]),
        tolerance = 1e-10,
        info = paste("batch_chol permutation seed", seed, "index", idx)
      )
      expect_equal(
        perm_solutions[[idx]],
        ref_solutions[[perm[[idx]]]],
        tolerance = 1e-10,
        info = paste("batch_solve permutation seed", seed, "index", idx)
      )
    }
  }
})

test_that("adversarial inputs fail or degrade honestly", {
  one <- adgeMatrix(matrix(5, nrow = 1L, ncol = 1L))

  expect_equal(as.matrix(one %*% matrix(2, nrow = 1L, ncol = 1L)), matrix(10, nrow = 1L, ncol = 1L))
  expect_equal(as.matrix(crossprod(one)), matrix(25, nrow = 1L, ncol = 1L))
  expect_equal(as.matrix(tcrossprod(one)), matrix(25, nrow = 1L, ncol = 1L))
  expect_equal(rowSums(matrix(0, nrow = 3L, ncol = 2L)), c(0, 0, 0))
  expect_equal(colSums(matrix(0, nrow = 3L, ncol = 2L)), c(0, 0))

  expect_equal(
    as.matrix(chol_factor(as_adgeMatrix(matrix(4, nrow = 1L, ncol = 1L)))),
    matrix(2, nrow = 1L, ncol = 1L),
    tolerance = 1e-12
  )

  expect_error(solve(as_adgeMatrix(matrix(c(1, 2, 2, 4), nrow = 2L))))
  expect_error(chol_factor(as_adgeMatrix(matrix(c(1, 2, 2, 4), nrow = 2L))))
  expect_error(chol_factor(as_adgeMatrix(matrix(c(1, 0, 0, NaN), nrow = 2L))))
  expect_true(any(!is.finite(as.matrix(chol_factor(as_adgeMatrix(matrix(c(1, 0, 0, Inf), nrow = 2L)))))))

  zero_fit <- block_lanczos(matrix(0, nrow = 8L, ncol = 5L), nv = 3L, nu = 3L, block_size = 3L, n_steps = 3L)
  expect_equal(zero_fit$d, c(0, 0, 0), tolerance = 0)
  expect_true(all(is.finite(zero_fit$d)))

  one_kernel <- kernel_matrix(matrix(1, nrow = 1L, ncol = 1L), kernel = "rbf", sigma = 0.5)
  expect_identical(dim(one_kernel), c(1L, 1L))
  expect_equal(one_kernel[1, 1], 1, tolerance = 1e-12)

  set.seed(2026041601L)
  x_rank_def <- matrix(rnorm(15L * 4L), nrow = 15L, ncol = 4L)
  x_rank_def[, 4L] <- x_rank_def[, 1L] + x_rank_def[, 2L]
  fit <- ridge_path(x_rank_def, rnorm(15L), lambdas = c(0.1, 1, 10))
  expect_true(all(is.finite(fit$coef)))
})

test_that("QR helper regressions and validation failures are pinned down", {
  set.seed(2026042001L)
  x <- matrix(rnorm(32L), nrow = 8L, ncol = 4L)

  qr_fit <- am_qr(x)
  expect_s3_class(qr_fit, "amQR")
  expect_equal(qr_info(qr_fit)$dim, c(8L, 4L))
  expect_equal(as.matrix(qr.R(qr_fit)), qr.R(qr(x)), tolerance = 1e-8)

  expect_error(qr_downdate(qr_fit, 0L, X = x), "row_idx must be between 1 and 8")
  expect_error(qr_downdate(qr_fit, 9L, X = x), "row_idx must be between 1 and 8")
  expect_error(qr_downdate(qr_fit, 1.5, X = x), "row_idx must be a single positive integer")
  expect_error(qr_downdate(qr_fit, c(1L, 2L), X = x), "row_idx must be a single positive integer")
})
