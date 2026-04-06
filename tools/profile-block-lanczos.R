#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else if (requireNamespace("amatrix", quietly = TRUE)) {
    library(amatrix)
  } else {
    stop("Either the installed 'amatrix' package or 'pkgload' is required", call. = FALSE)
  }
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for these profiling helpers", call. = FALSE)
  }
  if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("backends/amatrix.mlx")) {
    pkgload::load_all("backends/amatrix.mlx", quiet = TRUE)
  } else if (requireNamespace("amatrix.mlx", quietly = TRUE)) {
    library(amatrix.mlx)
  }
})

if (exists("amatrix_mlx_is_available", mode = "function")) {
  options(amatrix.mlx.available = TRUE)
}

.profile_elapsed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  list(value = value, elapsed = unname(elapsed))
}

.profile_add <- function(totals, counts, label, elapsed) {
  if (is.na(elapsed) || !nzchar(label)) {
    return(list(totals = totals, counts = counts))
  }
  if (is.null(totals[[label]])) {
    totals[[label]] <- 0
    counts[[label]] <- 0L
  }
  totals[[label]] <- totals[[label]] + elapsed
  counts[[label]] <- counts[[label]] + 1L
  list(totals = totals, counts = counts)
}

.profile_stage <- function(totals, counts, label, expr) {
  timed <- .profile_elapsed(expr)
  updated <- .profile_add(totals, counts, label, timed$elapsed)
  list(
    value = timed$value,
    totals = updated$totals,
    counts = updated$counts
  )
}

.profile_finalize <- function(totals, counts, total_elapsed, metadata = list()) {
  labels <- names(totals)
  rows <- data.frame(
    stage = labels,
    elapsed = unlist(totals, use.names = FALSE),
    calls = unlist(counts[labels], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  rows$share <- if (total_elapsed > 0) rows$elapsed / total_elapsed else 0
  rows <- rows[order(rows$elapsed, decreasing = TRUE), , drop = FALSE]
  rows$elapsed <- round(rows$elapsed, 6L)
  rows$share <- round(rows$share, 4L)
  c(metadata, list(total_elapsed = round(total_elapsed, 6L), stages = rows))
}

profile_block_lanczos_case <- function(n,
                                       p,
                                       k,
                                       block_size,
                                       n_steps,
                                       seed = 20260406L,
                                       backend = "mlx") {
  set.seed(seed)
  host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x <- adgeMatrix(host, preferred_backend = backend, precision = "fast")
  reorth <- get(".amatrix_block_reorth", envir = asNamespace("amatrix"), inherits = FALSE)
  needs_final_qr <- get(".amatrix_block_basis_needs_final_qr", envir = asNamespace("amatrix"), inherits = FALSE)
  left_product <- get(".amatrix_block_lanczos_left_product", envir = asNamespace("amatrix"), inherits = FALSE)
  right_operator <- get(".amatrix_block_lanczos_right_operator", envir = asNamespace("amatrix"), inherits = FALSE)(x)
  right_product <- get(".amatrix_block_lanczos_right_product", envir = asNamespace("amatrix"), inherits = FALSE)
  right_drop <- get(".amatrix_block_lanczos_drop_right_operator", envir = asNamespace("amatrix"), inherits = FALSE)
  on.exit(right_drop(right_operator), add = TRUE)

  b <- min(as.integer(block_size), n, p)
  J <- as.integer(n_steps)
  totals <- list()
  counts <- list()

  overall <- .profile_elapsed({
    prof <- .profile_stage(
      totals,
      counts,
      "init_random_qr",
      qr(matrix(rnorm(p * b), nrow = p, ncol = b))
    )
    Q_cur <- qr.Q(prof$value)
    totals <- prof$totals
    counts <- prof$counts

    Q_left_blocks <- vector("list", J)
    Q_right_blocks <- vector("list", J)
    B_proj <- matrix(0, nrow = J * b, ncol = J * b)

    for (j in seq_len(J)) {
      prof <- .profile_stage(
        totals,
        counts,
        "bind_prev_left_basis",
        if (j > 1L) do.call(cbind, Q_left_blocks[seq_len(j - 1L)]) else NULL
      )
      QL_prev <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(
        totals,
        counts,
        "bind_prev_right_basis",
        if (j > 1L) do.call(cbind, Q_right_blocks[seq_len(j - 1L)]) else NULL
      )
      QR_prev <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "left_matmul_dispatch", left_product(x, Q_cur))
      Z_left_raw <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(
        totals,
        counts,
        "left_reorth",
        reorth(Z_left_raw, QL_prev, return_projection = TRUE)
      )
      left_reorth <- prof$value
      Z_left <- left_reorth$z
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "left_qr", qr(Z_left))
      left_qr <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "left_q_extract", qr.Q(left_qr))
      QL_j <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(
        totals,
        counts,
        "left_r_extract",
        qr.R(left_qr, complete = FALSE)
      )
      R_left_j <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      storage.mode(QL_j) <- "double"
      storage.mode(R_left_j) <- "double"
      Q_left_blocks[[j]] <- QL_j

      if (j > 1L) {
        col_idx <- ((j - 2L) * b + 1L):((j - 1L) * b)
        if (!is.null(left_reorth$coeff) && nrow(left_reorth$coeff) > 0L) {
          B_proj[seq_len(nrow(left_reorth$coeff)), col_idx] <- left_reorth$coeff
        }
        row_idx <- ((j - 1L) * b + 1L):(j * b)
        B_proj[row_idx, col_idx] <- R_left_j[seq_len(b), seq_len(b), drop = FALSE]
      }

      right_label <- if (is.null(right_operator)) "right_crossprod_dispatch" else "right_matmul_dispatch"
      prof <- .profile_stage(
        totals,
        counts,
        right_label,
        right_product(x, QL_j, operator = right_operator)
      )
      right_prod <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "right_reorth", reorth(right_prod, QR_prev))
      Z_right <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "right_qr", qr(Z_right))
      right_qr <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "right_q_extract", qr.Q(right_qr))
      QR_j <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      storage.mode(QR_j) <- "double"
      Q_right_blocks[[j]] <- QR_j
      Q_cur <- QR_j
    }

    prof <- .profile_stage(totals, counts, "bind_full_left_basis", do.call(cbind, Q_left_blocks))
    Q_L <- prof$value
    totals <- prof$totals
    counts <- prof$counts

    prof <- .profile_stage(totals, counts, "bind_full_right_basis", do.call(cbind, Q_right_blocks))
    Q_R <- prof$value
    totals <- prof$totals
    counts <- prof$counts

    needs_final_left_qr <- J > 1L && needs_final_qr(Q_L)
    if (needs_final_left_qr) {
      prof <- .profile_stage(totals, counts, "final_left_basis_qr", qr.Q(qr(Q_L)))
      Q_L <- prof$value
      totals <- prof$totals
      counts <- prof$counts
    }

    needs_final_right_qr <- J > 1L && needs_final_qr(Q_R)
    if (needs_final_right_qr) {
      prof <- .profile_stage(totals, counts, "final_right_basis_qr", qr.Q(qr(Q_R)))
      Q_R <- prof$value
      totals <- prof$totals
      counts <- prof$counts
    }

    if (needs_final_left_qr || needs_final_right_qr) {
      prof <- .profile_stage(totals, counts, "final_project_dispatch", x %*% Q_R)
      AQ_R_prod <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "final_project_materialize", as.matrix(AQ_R_prod))
      AQ_R <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "projected_crossprod", crossprod(Q_L, AQ_R))
      B <- prof$value
      totals <- prof$totals
      counts <- prof$counts
    } else {
      prof <- .profile_stage(totals, counts, "last_block_dispatch", left_product(x, Q_cur))
      AQ_last <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      prof <- .profile_stage(totals, counts, "last_block_crossprod", crossprod(Q_L, AQ_last))
      last_col <- prof$value
      totals <- prof$totals
      counts <- prof$counts

      B <- B_proj
      last_col_idx <- ((J - 1L) * b + 1L):(J * b)
      B[, last_col_idx] <- last_col
    }

    k_out <- min(as.integer(k), ncol(Q_L), ncol(Q_R), nrow(B), ncol(B))

    prof <- .profile_stage(totals, counts, "projected_svd", base::svd(B, nu = k_out, nv = k_out))
    svd_B <- prof$value
    totals <- prof$totals
    counts <- prof$counts

    prof <- .profile_stage(totals, counts, "lift_left_vectors", Q_L %*% svd_B$u[, seq_len(k_out), drop = FALSE])
    U <- prof$value
    totals <- prof$totals
    counts <- prof$counts

    prof <- .profile_stage(totals, counts, "lift_right_vectors", Q_R %*% svd_B$v[, seq_len(k_out), drop = FALSE])
    V <- prof$value
    totals <- prof$totals
    counts <- prof$counts

    list(
      u = U[, seq_len(k_out), drop = FALSE],
      d = svd_B$d[seq_len(k_out)],
      v = V[, seq_len(k_out), drop = FALSE]
    )
  })

  .profile_finalize(
    totals = totals,
    counts = counts,
    total_elapsed = overall$elapsed,
    metadata = list(
      case = sprintf("%dx%d", n, p),
      k = as.integer(k),
      block_size = b,
      n_steps = J,
      backend = backend
    )
  )
}

benchmark_block_lanczos_products <- function(n,
                                             p,
                                             block_cols,
                                             reps = 5L,
                                             seed = 20260408L,
                                             backend = "mlx") {
  set.seed(seed)
  host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  q_left <- qr.Q(qr(matrix(rnorm(p * block_cols), nrow = p, ncol = block_cols)))
  q_right <- qr.Q(qr(matrix(rnorm(n * block_cols), nrow = n, ncol = block_cols)))
  x <- adgeMatrix(host, preferred_backend = backend, precision = "fast")

  bench <- function(fn) {
    timings <- numeric(reps)
    last <- NULL
    for (idx in seq_len(reps)) {
      gc()
      timings[[idx]] <- system.time(last <- fn())[["elapsed"]]
    }
    c(elapsed = median(timings), nrow = nrow(last), ncol = ncol(last))
  }

  rows <- rbind(
    data.frame(
      operation = "matmul_cpu",
      elapsed = bench(function() host %*% q_left)[["elapsed"]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      operation = "matmul_mlx_materialize",
      elapsed = bench(function() as.matrix(x %*% q_left))[["elapsed"]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      operation = "crossprod_cpu",
      elapsed = bench(function() crossprod(host, q_right))[["elapsed"]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      operation = "crossprod_mlx_materialize",
      elapsed = bench(function() as.matrix(crossprod(x, q_right)))[["elapsed"]],
      stringsAsFactors = FALSE
    )
  )

  rows$case <- sprintf("%dx%d_x_%d", n, p, block_cols)
  rows$elapsed <- round(rows$elapsed, 6L)
  rows
}

print_block_lanczos_profile <- function(n = 3000L,
                                        p = 1200L,
                                        k = 20L,
                                        block_size = 24L,
                                        n_steps = 4L,
                                        seed = 20260406L,
                                        backend = "mlx") {
  profile <- profile_block_lanczos_case(
    n = n,
    p = p,
    k = k,
    block_size = block_size,
    n_steps = n_steps,
    seed = seed,
    backend = backend
  )
  cat(sprintf(
    "Block Lanczos profile: case=%s k=%d block_size=%d n_steps=%d backend=%s total=%.3fs\n\n",
    profile$case,
    profile$k,
    profile$block_size,
    profile$n_steps,
    profile$backend,
    profile$total_elapsed
  ))
  print(profile$stages, row.names = FALSE)
  cat("\nRepresentative product timings:\n")
  print(benchmark_block_lanczos_products(
    n = n,
    p = p,
    block_cols = block_size,
    seed = seed + 1L,
    backend = backend
  ), row.names = FALSE)
  invisible(profile)
}

if (isTRUE(getOption("amatrix.block_lanczos.profile.autorun", FALSE))) {
  print_block_lanczos_profile()
}
