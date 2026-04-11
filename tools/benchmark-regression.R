#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x
r_string_literal <- function(x) encodeString(x, quote = "\"")

timestamp_tag <- function(x = Sys.time()) {
  format(as.POSIXct(x, tz = Sys.timezone()), "%Y%m%d-%H%M%S")
}

parse_args <- function(args) {
  out <- list(
    update = "--update" %in% args,
    warm_only = "--warm-only" %in% args,
    worker = "--worker" %in% args,
    safe_main = "--safe-main" %in% args,
    baseline = "tools/baseline.csv",
    output_dir = file.path("tools", "benchmark-results", timestamp_tag()),
    plan = NULL,
    out = NULL,
    group_id = NULL,
    suites = c("dense", "sparse"),
    backends = NULL,
    include_arrayfire = identical(Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE", unset = ""), "1") ||
      identical(Sys.getenv("AMATRIX_ARRAYFIRE_PROBE_GPU", unset = ""), "1")
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    key <- sub("^--", "", arg)
    if (!grepl("=", key, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(key, "=", fixed = TRUE)[[1L]]
    name <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = "=")
    if (identical(name, "baseline")) out$baseline <- value
    if (identical(name, "output-dir")) out$output_dir <- value
    if (identical(name, "plan")) out$plan <- value
    if (identical(name, "out")) out$out <- value
    if (identical(name, "group-id")) out$group_id <- value
    if (identical(name, "suites")) out$suites <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(name, "backends")) out$backends <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(name, "include-arrayfire")) out$include_arrayfire <- identical(value, "1")
  }

  out$suites <- unique(out$suites[nzchar(out$suites)])
  out
}

relaunch_safe_master_if_needed <- function(args) {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_arg <- grep("^--file=", raw_args, value = TRUE)
  direct_file_entry <- length(direct_file_arg) > 0L

  if (!direct_file_entry || isTRUE(args$worker) || isTRUE(args$safe_main)) {
    return(invisible(FALSE))
  }

  script_path <- normalizePath(sub("^--file=", "", direct_file_arg[[1L]]), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  expr <- sprintf(
    "setwd(%s); source(%s, local = globalenv())",
    r_string_literal(repo_root),
    r_string_literal(script_path)
  )
  relaunch_args <- c("-e", expr, "--args", "--safe-main", commandArgs(trailingOnly = TRUE))
  relaunch_cmd <- paste(
    c(
      shQuote(file.path(R.home("bin"), "Rscript")),
      vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE)
    ),
    collapse = " "
  )
  warned_status <- NULL
  relaunch_output <- withCallingHandlers(
    system(paste(relaunch_cmd, "2>&1"), intern = TRUE),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  status <- attr(relaunch_output, "status") %||% warned_status %||% 0L

  if (length(relaunch_output) > 0L) {
    cat(paste(relaunch_output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(relaunch_output, collapse = "\n"))) {
      cat("\n")
    }
  }

  quit(save = "no", status = status)
}

initialize_regression_benchmark_context <- local({
  initialized <- FALSE

  function() {
    if (initialized) {
      return(invisible(TRUE))
    }

    suppressPackageStartupMessages({
      if (file.exists(file.path("tools", "benchmark-helpers.R"))) {
        source(file.path("tools", "benchmark-helpers.R"), local = FALSE)
      }
      load_benchmark_amatrix()
      if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Matrix package required for sparse benchmark suite", call. = FALSE)
      }
    })

    initialized <<- TRUE
    invisible(TRUE)
  }
})

canonical_backend_specs <- function(include_arrayfire = .benchmark_arrayfire_requested()) {
  available_benchmark_backends(
    include_cpu = TRUE,
    include_mlx = TRUE,
    include_metal = TRUE,
    include_opencl = TRUE,
    include_arrayfire = include_arrayfire
  )
}

filter_backend_specs <- function(specs, names = NULL) {
  if (is.null(names)) {
    return(specs)
  }

  keep <- names[names %in% vapply(specs, `[[`, character(1), "name")]
  specs[vapply(specs, function(spec) spec$name %in% keep, logical(1))]
}

prime_requested_backend <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  if (identical(backend_name, "cpu")) {
    return(invisible(TRUE))
  }

  specs <- .benchmark_optional_backend_specs(include_arrayfire = include_arrayfire)
  spec <- specs[[backend_name]]
  if (is.null(spec)) {
    return(invisible(FALSE))
  }

  invisible(.benchmark_enable_backend(spec))
}

dense_sizes <- function() {
  list(
    small = list(n = 256L, p = 32L, sink_n = 128L),
    medium = list(n = 1024L, p = 128L, sink_n = 512L),
    large = list(n = 4096L, p = 128L, sink_n = 1024L)
  )
}

sparse_sizes <- function() {
  list(
    medium = list(nrow = 4000L, ncol = 1000L),
    large = list(nrow = 8000L, ncol = 2000L)
  )
}

sparse_densities <- function() c(0.001, 0.005, 0.01, 0.05)

spmm_rhs_widths <- function() c(8L, 32L)

sparse_iterative_ops <- function() c("block_lanczos", "svd_factor_subspace")

sparse_iterative_densities <- function() c(0.05)

sparse_iterative_params <- function(op, size_label) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  switch(
    op,
    block_lanczos = list(
      k = 8L,
      block_size = min(16L, size$nrow, size$ncol),
      n_steps = 4L
    ),
    svd_factor_subspace = list(
      k = 8L,
      n_oversamples = 8L,
      n_iter = 1L
    ),
    stop(sprintf("unknown sparse iterative op '%s'", op), call. = FALSE)
  )
}

new_result_row <- function(...) {
  defaults <- list(
    suite = NA_character_,
    op = NA_character_,
    size_label = NA_character_,
    variant = NA_character_,
    requested_backend = NA_character_,
    dispatch_backend = NA_character_,
    dispatch_path = NA_character_,
    status = "ok",
    error_message = NA_character_,
    nrow = NA_integer_,
    ncol = NA_integer_,
    rhs_width = NA_integer_,
    nnz = NA_integer_,
    density = NA_real_,
    density_bucket = NA_character_,
    reps = NA_integer_,
    median_ms = NA_real_
  )
  values <- modifyList(defaults, list(...))
  as.data.frame(values, stringsAsFactors = FALSE)
}

benchmark_time_ms <- function(fn, reps = 7L, warmup = 1L) {
  reps <- as.integer(reps)
  warmup <- as.integer(warmup)

  if (warmup > 0L) {
    for (idx in seq_len(warmup)) {
      fn()
    }
  }

  timings <- vapply(seq_len(reps), function(idx) {
    gc()
    start <- proc.time()[["elapsed"]]
    fn()
    (proc.time()[["elapsed"]] - start) * 1000
  }, numeric(1))

  median(timings)
}

make_dense_operand <- function(host, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  adgeMatrix(host, preferred_backend = backend, precision = precision)
}

make_sparse_operand <- function(host, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  as_adgCMatrix(host, preferred_backend = backend, precision = precision)
}

with_sparse_benchmark_backend_options <- function(requested_backend, rhs_width, code) {
  rhs_width <- as.integer(rhs_width)
  old <- options()
  on.exit(options(old), add = TRUE)

  if (identical(requested_backend, "mlx")) {
    options(
      amatrix.mlx.spmv_min_nnz = 1L,
      amatrix.mlx.spmm_min_nnz = 1L
    )
  }

  if (identical(requested_backend, "metal")) {
    options(
      amatrix.metal.spmv_min_nnz = 1L,
      amatrix.metal.spmm_min_nnz = 1L
    )
  }

  if (identical(requested_backend, "opencl")) {
    options(
      amatrix.opencl.spmv_min_nnz = 1L,
      amatrix.opencl.spmm_min_nnz = 1L
    )
  }

  force(code)
}

release_residency <- function(x) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  entry <- amatrix:::.amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(amatrix:::.amatrix_get_backend(entry$backend), error = function(e) NULL)
  if (is.null(backend)) {
    return(invisible(FALSE))
  }

  if (isTRUE(entry$sparse) && is.function(backend$sparse_resident_drop)) {
    try(backend$sparse_resident_drop(entry$resident_key), silent = TRUE)
  } else if (is.function(backend$resident_drop)) {
    try(backend$resident_drop(entry$resident_key), silent = TRUE)
  }
  amatrix:::.amatrix_drop_resident_binding(x)
  invisible(TRUE)
}

drop_svd_factor_cache <- function(x, k, n_oversamples, n_iter) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  plan <- tryCatch(
    amatrix:::.amatrix_svd_factor_plan(
      X = x,
      k = as.integer(k),
      method = "subspace",
      n_oversamples = as.integer(n_oversamples),
      n_iter = as.integer(n_iter)
    ),
    error = function(e) NULL
  )
  if (is.null(plan)) {
    return(invisible(FALSE))
  }

  key <- tryCatch(amatrix:::.amatrix_svd_cache_key(x, as.integer(k), plan), error = function(e) NULL)
  if (is.null(key) || !nzchar(key)) {
    return(invisible(FALSE))
  }

  if (exists(key, envir = amatrix:::.amatrix_state$model_cache, inherits = FALSE)) {
    rm(list = key, envir = amatrix:::.amatrix_state$model_cache)
  }
  if (exists(key, envir = amatrix:::.amatrix_state$cache_atime, inherits = FALSE)) {
    rm(list = key, envir = amatrix:::.amatrix_state$cache_atime)
  }

  invisible(TRUE)
}

compute_dispatch_info <- function(x, op, y = NULL) {
  plan <- tryCatch(amatrix_backend_plan(x, op, y = y), error = function(e) NULL)
  if (is.null(plan)) {
    return(list(dispatch_backend = NA_character_, dispatch_path = NA_character_))
  }
  list(dispatch_backend = plan$chosen %||% NA_character_, dispatch_path = plan$chosen_path %||% NA_character_)
}

run_dense_case <- function(requested_backend, op, size_label, variant) {
  size <- dense_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown dense size '%s'", size_label), call. = FALSE)
  }

  set.seed(1L)
  X <- matrix(rnorm(size$n * size$p), size$n, size$p)
  set.seed(2L)
  Y <- matrix(rnorm(size$n * 3L), size$n, 3L)
  set.seed(3L)
  B_host <- matrix(rnorm(size$p * size$p), size$p, size$p)
  set.seed(4L)
  Z <- matrix(rnorm(size$p * size$p), size$p, size$p)
  SPD <- crossprod(Z) + diag(size$p) * 2
  SPD_rhs <- matrix(rnorm(size$p * 3L), size$p, 3L)
  set.seed(5L)
  sink_host <- exp(matrix(rnorm(size$sink_n * size$sink_n), size$sink_n, size$sink_n))
  n_dist <- min(size$n, 512L)
  Xs <- X[seq_len(n_dist), , drop = FALSE]

  if (identical(op, "sinkhorn")) {
    probe_x <- make_dense_operand(sink_host, requested_backend)
    dispatch <- compute_dispatch_info(probe_x, "matmul")
    release_residency(probe_x)
  } else {
    probe_x <- make_dense_operand(if (op %in% c("chol", "solve_rhs", "eigen_sym")) SPD else X, requested_backend)
    probe_y <- switch(
      op,
      matmul = make_dense_operand(B_host, requested_backend),
      solve_rhs = make_dense_operand(SPD_rhs, requested_backend),
      NULL
    )
    dispatch_op <- switch(
      op,
      matmul = "matmul",
      crossprod = "crossprod",
      covariance = "crossprod",
      dist = "crossprod",
      chol = "chol",
      solve_rhs = "solve",
      eigen_sym = "eigen",
      many_lm = "matmul",
      rsvd = "svd",
      sinkhorn = "matmul"
    )
    dispatch <- compute_dispatch_info(probe_x, dispatch_op, y = probe_y)
    release_residency(probe_x)
    if (inherits(probe_y, "aMatrix")) {
      release_residency(probe_y)
    }
  }

  reps <- 7L
  runner <- switch(
    variant,
    cold = switch(
      op,
      matmul = function() { aX <- make_dense_operand(X, requested_backend); aB <- make_dense_operand(B_host, requested_backend); on.exit({ release_residency(aX); release_residency(aB) }, add = TRUE); invisible(aX %*% aB) },
      crossprod = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(crossprod(aX)) },
      covariance = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(covariance(aX)) },
      dist = function() { aXs <- make_dense_operand(Xs, requested_backend); on.exit(release_residency(aXs), add = TRUE); invisible(dist_matrix(aXs)) },
      chol = function() { aS <- make_dense_operand(SPD, requested_backend); on.exit(release_residency(aS), add = TRUE); invisible(chol(aS)) },
      solve_rhs = function() { aS <- make_dense_operand(SPD, requested_backend); aR <- make_dense_operand(SPD_rhs, requested_backend); on.exit({ release_residency(aS); release_residency(aR) }, add = TRUE); invisible(solve(aS, aR)) },
      eigen_sym = function() { aS <- make_dense_operand(SPD, requested_backend); on.exit(release_residency(aS), add = TRUE); invisible(eigh(aS)) },
      many_lm = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(many_lm(aX, Y, method = "qr", cache = FALSE)) },
      rsvd = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(rsvd(aX, k = 10L)) },
      sinkhorn = function() { aSink <- make_dense_operand(sink_host, requested_backend); on.exit(release_residency(aSink), add = TRUE); invisible(sinkhorn(aSink, max_iter = 25L, tol = 0, return_info = FALSE)) }
    ),
    warm = {
      aX <- if (op %in% c("chol", "solve_rhs", "eigen_sym")) NULL else make_dense_operand(X, requested_backend)
      aB <- if (identical(op, "matmul")) make_dense_operand(B_host, requested_backend) else NULL
      aXs <- if (identical(op, "dist")) make_dense_operand(Xs, requested_backend) else NULL
      aS <- if (op %in% c("chol", "solve_rhs", "eigen_sym")) make_dense_operand(SPD, requested_backend) else NULL
      aR <- if (identical(op, "solve_rhs")) make_dense_operand(SPD_rhs, requested_backend) else NULL
      aSink <- if (identical(op, "sinkhorn")) make_dense_operand(sink_host, requested_backend) else NULL
      on.exit({
        for (obj in list(aX, aB, aXs, aS, aR, aSink)) {
          release_residency(obj)
        }
      }, add = TRUE)

      tryCatch({
        if (!is.null(aB)) invisible(aX %*% aB)
        if (identical(op, "crossprod")) invisible(crossprod(aX))
        if (identical(op, "covariance")) invisible(covariance(aX))
        if (identical(op, "dist")) invisible(dist_matrix(aXs))
        if (identical(op, "chol")) invisible(chol(aS))
        if (identical(op, "solve_rhs")) invisible(solve(aS, aR))
        if (identical(op, "eigen_sym")) invisible(eigh(aS))
        if (identical(op, "many_lm")) invisible(many_lm(aX, Y, method = "qr", cache = TRUE))
        if (identical(op, "rsvd")) invisible(rsvd(aX, k = 10L))
        if (identical(op, "sinkhorn")) invisible(sinkhorn(aSink, max_iter = 5L, tol = 0, return_info = FALSE))
      }, error = function(e) NULL)

      switch(
        op,
        matmul = function() invisible(aX %*% aB),
        crossprod = function() invisible(crossprod(aX)),
        covariance = function() invisible(covariance(aX)),
        dist = function() invisible(dist_matrix(aXs)),
        chol = function() invisible(chol(aS)),
        solve_rhs = function() invisible(solve(aS, aR)),
        eigen_sym = function() invisible(eigh(aS)),
        many_lm = function() invisible(many_lm(aX, Y, method = "qr", cache = TRUE)),
        rsvd = function() invisible(rsvd(aX, k = 10L)),
        sinkhorn = function() invisible(sinkhorn(aSink, max_iter = 25L, tol = 0, return_info = FALSE))
      )
    },
    stop(sprintf("unknown dense variant '%s'", variant), call. = FALSE)
  )

  median_ms <- benchmark_time_ms(runner, reps = reps)

  new_result_row(
    suite = "dense",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_backend = dispatch$dispatch_backend,
    dispatch_path = dispatch$dispatch_path,
    nrow = if (identical(op, "sinkhorn")) size$sink_n else if (identical(op, "eigen_sym")) size$p else size$n,
    ncol = if (identical(op, "sinkhorn")) size$sink_n else size$p,
    rhs_width = switch(op, matmul = ncol(B_host), solve_rhs = ncol(SPD_rhs), many_lm = ncol(Y), sinkhorn = size$sink_n, 0L),
    nnz = 0L,
    density = 0,
    density_bucket = "dense",
    reps = reps,
    median_ms = median_ms
  )
}

make_sparse_case <- function(nrow, ncol, density, rhs_width, seed) {
  set.seed(seed)
  list(
    X_host = Matrix::rsparsematrix(nrow, ncol, density = density),
    rhs_host = matrix(rnorm(ncol * rhs_width), nrow = ncol, ncol = rhs_width)
  )
}

run_sparse_case <- function(requested_backend, op, size_label, density, rhs_width, variant) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  rhs_width <- as.integer(rhs_width)
  case <- make_sparse_case(
    nrow = size$nrow,
    ncol = size$ncol,
    density = density,
    rhs_width = rhs_width,
    seed = 1000L + rhs_width + as.integer(density * 1000)
  )

  sparse_result <- with_sparse_benchmark_backend_options(requested_backend, rhs_width, {
    probe_x <- make_sparse_operand(case$X_host, requested_backend)
    probe_y <- make_dense_operand(case$rhs_host, requested_backend)

    if (identical(variant, "resident") &&
        requested_backend %in% c("mlx", "metal", "arrayfire")) {
      probe_x <- amatrix::amatrix_bind_resident(probe_x, backend = requested_backend)
      probe_y <- amatrix::amatrix_bind_resident(probe_y, backend = requested_backend)
    }

    dispatch <- compute_dispatch_info(probe_x, "matmul", y = probe_y)
    release_residency(probe_x)
    release_residency(probe_y)

    reps <- 5L
    runner <- switch(
      variant,
      cold = function() {
        x <- make_sparse_operand(case$X_host, requested_backend)
        y <- make_dense_operand(case$rhs_host, requested_backend)
        on.exit({
          release_residency(x)
          release_residency(y)
        }, add = TRUE)
        invisible(x %*% y)
      },
      resident = {
        resident_x <- make_sparse_operand(case$X_host, requested_backend)
        resident_y <- make_dense_operand(case$rhs_host, requested_backend)

        if (requested_backend %in% c("mlx", "metal", "arrayfire")) {
          resident_x <- amatrix::amatrix_bind_resident(resident_x, backend = requested_backend)
          resident_y <- amatrix::amatrix_bind_resident(resident_y, backend = requested_backend)
        }

        on.exit({
          release_residency(resident_x)
          release_residency(resident_y)
        }, add = TRUE)
        invisible(resident_x %*% resident_y)
        function() invisible(resident_x %*% resident_y)
      },
      stop(sprintf("unknown sparse variant '%s'", variant), call. = FALSE)
    )

    list(
      dispatch = dispatch,
      reps = reps,
      median_ms = benchmark_time_ms(runner, reps = reps)
    )
  })

  dispatch <- sparse_result$dispatch
  reps <- sparse_result$reps
  median_ms <- sparse_result$median_ms
  density_value <- amatrix:::.amatrix_sparse_density(case$X_host)

  new_result_row(
    suite = "sparse",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_backend = dispatch$dispatch_backend,
    dispatch_path = dispatch$dispatch_path,
    nrow = size$nrow,
    ncol = size$ncol,
    rhs_width = rhs_width,
    nnz = length(case$X_host@x),
    density = density_value,
    density_bucket = amatrix:::.amatrix_sparse_density_bucket(density_value),
    reps = reps,
    median_ms = median_ms
  )
}

run_sparse_iterative_case <- function(requested_backend, op, size_label, density, variant) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  params <- sparse_iterative_params(op, size_label)
  case <- make_sparse_case(
    nrow = size$nrow,
    ncol = size$ncol,
    density = density,
    rhs_width = 1L,
    seed = 2000L + as.integer(density * 1000)
  )

  iterative_result <- with_sparse_benchmark_backend_options(requested_backend, 32L, {
    probe_x <- make_sparse_operand(case$X_host, requested_backend)
    dispatch_backend <- if (identical(op, "svd_factor_subspace")) {
      plan <- tryCatch(amatrix:::.amatrix_svd_factor_plan(
        probe_x,
        k = params$k,
        method = "subspace",
        n_oversamples = params$n_oversamples %||% 0L,
        n_iter = params$n_iter %||% 0L
      ), error = function(e) NULL)
      plan$subspace_backend %||% requested_backend
    } else {
      tryCatch(amatrix::amatrix_resident_backend_for(probe_x, op = "matmul"), error = function(e) NULL) %||% requested_backend
    }
    release_residency(probe_x)

    reps <- 3L
    runner <- switch(
      variant,
      cold = function() {
        x <- make_sparse_operand(case$X_host, requested_backend)
        on.exit(release_residency(x), add = TRUE)
        if (identical(op, "block_lanczos")) {
          invisible(block_lanczos(
            x,
            nv = params$k,
            nu = params$k,
            block_size = params$block_size,
            n_steps = params$n_steps
          ))
        } else {
          invisible(svd_factor(
            x,
            k = params$k,
            method = "subspace",
            n_oversamples = params$n_oversamples,
            n_iter = params$n_iter
          ))
        }
      },
      warm = {
        resident_x <- make_sparse_operand(case$X_host, requested_backend)
        if (requested_backend %in% c("mlx", "metal", "arrayfire")) {
          resident_x <- amatrix::amatrix_bind_resident(resident_x, backend = requested_backend, op = "matmul")
        }

        on.exit(release_residency(resident_x), add = TRUE)
        if (identical(op, "block_lanczos")) {
          invisible(block_lanczos(
            resident_x,
            nv = params$k,
            nu = params$k,
            block_size = params$block_size,
            n_steps = params$n_steps
          ))
          function() invisible(block_lanczos(
            resident_x,
            nv = params$k,
            nu = params$k,
            block_size = params$block_size,
            n_steps = params$n_steps
          ))
        } else {
          drop_svd_factor_cache(
            resident_x,
            k = params$k,
            n_oversamples = params$n_oversamples,
            n_iter = params$n_iter
          )
          invisible(svd_factor(
            resident_x,
            k = params$k,
            method = "subspace",
            n_oversamples = params$n_oversamples,
            n_iter = params$n_iter
          ))
          function() {
            drop_svd_factor_cache(
              resident_x,
              k = params$k,
              n_oversamples = params$n_oversamples,
              n_iter = params$n_iter
            )
            invisible(svd_factor(
              resident_x,
              k = params$k,
              method = "subspace",
              n_oversamples = params$n_oversamples,
              n_iter = params$n_iter
            ))
          }
        }
      },
      stop(sprintf("unknown sparse iterative variant '%s'", variant), call. = FALSE)
    )

    list(
      dispatch_backend = dispatch_backend,
      reps = reps,
      median_ms = benchmark_time_ms(runner, reps = reps)
    )
  })

  density_value <- amatrix:::.amatrix_sparse_density(case$X_host)

  new_result_row(
    suite = "sparse",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_backend = iterative_result$dispatch_backend,
    dispatch_path = "iterative",
    nrow = size$nrow,
    ncol = size$ncol,
    rhs_width = params$k,
    nnz = length(case$X_host@x),
    density = density_value,
    density_bucket = amatrix:::.amatrix_sparse_density_bucket(density_value),
    reps = iterative_result$reps,
    median_ms = iterative_result$median_ms
  )
}

safe_case <- function(fun, cell_info) {
  tryCatch(
    fun(),
    error = function(e) {
      new_result_row(
        suite = cell_info$suite,
        op = cell_info$op,
        size_label = cell_info$size_label,
        variant = cell_info$variant,
        requested_backend = cell_info$requested_backend,
        nrow = cell_info$nrow %||% NA_integer_,
        ncol = cell_info$ncol %||% NA_integer_,
        rhs_width = cell_info$rhs_width %||% NA_integer_,
        density = cell_info$density %||% NA_real_,
        status = "error",
        error_message = conditionMessage(e)
      )
    }
  )
}

format_backend_diagnostics <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  parts <- c(sprintf("requested_backend=%s", backend_name))

  spec <- .benchmark_optional_backend_specs(include_arrayfire = include_arrayfire)[[backend_name]]
  if (!is.null(spec)) {
    if (!is.null(spec$env)) {
      env_parts <- vapply(
        names(spec$env),
        function(name) sprintf("%s=%s", name, Sys.getenv(name, unset = "<unset>")),
        character(1)
      )
      parts <- c(parts, env_parts)
    }
    if (!is.null(spec$options)) {
      option_parts <- vapply(
        names(spec$options),
        function(name) sprintf("%s=%s", name, as.character(getOption(name))),
        character(1)
      )
      parts <- c(parts, option_parts)
    }
  }

  backend <- tryCatch(amatrix:::.amatrix_get_backend(backend_name), error = function(e) NULL)
  parts <- c(parts, sprintf("registered=%s", !is.null(backend)))
  if (!is.null(backend)) {
    available <- tryCatch(isTRUE(backend$available()), error = function(e) NA)
    parts <- c(parts, sprintf("backend_available=%s", as.character(available)))
  }

  if (identical(backend_name, "opencl")) {
    ns <- tryCatch(
      ensure_optional_backend_namespace("amatrix.opencl", repo_dir = "backends/amatrix.opencl"),
      error = function(e) NULL
    )
    if (!is.null(ns)) {
      info <- tryCatch(get("amatrix_opencl_bridge_info", envir = ns, inherits = FALSE)(), error = function(e) NULL)
      if (is.list(info)) {
        bridge_keys <- intersect(c("compiled", "clblast", "native", "probe_enabled", "engine", "available"), names(info))
        bridge_parts <- vapply(
          bridge_keys,
          function(name) sprintf("bridge_%s=%s", name, as.character(info[[name]])),
          character(1)
        )
        parts <- c(parts, bridge_parts)
      }
    }
  }

  paste(parts, collapse = ", ")
}

group_plan <- function(backends, suites = c("dense", "sparse")) {
  out <- list()

  if ("dense" %in% suites) {
    dense_ops <- c("matmul", "crossprod", "covariance", "dist", "chol", "solve_rhs", "eigen_sym", "many_lm", "rsvd", "sinkhorn")
    for (backend in backends) {
      for (op in dense_ops) {
        out[[length(out) + 1L]] <- list(group_id = sprintf("dense-%s-%s", backend, op), suite = "dense", requested_backend = backend, op = op)
      }
    }
  }

  if ("sparse" %in% suites) {
    sparse_ops <- c("spmv", "spmm", sparse_iterative_ops())
    for (backend in backends) {
      for (op in sparse_ops) {
        out[[length(out) + 1L]] <- list(group_id = sprintf("sparse-%s-%s", backend, op), suite = "sparse", requested_backend = backend, op = op)
      }
    }
  }

  out
}

expand_group_rows <- function(group) {
  rows <- list()
  if (identical(group$suite, "dense")) {
    for (size_label in names(dense_sizes())) {
      for (variant in c("cold", "warm")) {
        rows[[length(rows) + 1L]] <- new_result_row(
          suite = "dense",
          op = group$op,
          size_label = size_label,
          variant = variant,
          requested_backend = group$requested_backend,
          nrow = dense_sizes()[[size_label]]$n,
          ncol = dense_sizes()[[size_label]]$p
        )
      }
    }
  }

  if (identical(group$suite, "sparse")) {
    if (group$op %in% sparse_iterative_ops()) {
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_iterative_densities()) {
          for (variant in c("cold", "warm")) {
            params <- sparse_iterative_params(group$op, size_label)
            rows[[length(rows) + 1L]] <- new_result_row(
              suite = "sparse",
              op = group$op,
              size_label = size_label,
              variant = variant,
              requested_backend = group$requested_backend,
              nrow = sparse_sizes()[[size_label]]$nrow,
              ncol = sparse_sizes()[[size_label]]$ncol,
              rhs_width = params$k,
              density = density
            )
          }
        }
      }
    } else {
      rhs_widths <- if (identical(group$op, "spmv")) 1L else spmm_rhs_widths()
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_densities()) {
          for (rhs_width in rhs_widths) {
            for (variant in c("cold", "resident")) {
              rows[[length(rows) + 1L]] <- new_result_row(
                suite = "sparse",
                op = group$op,
                size_label = size_label,
                variant = variant,
                requested_backend = group$requested_backend,
                nrow = sparse_sizes()[[size_label]]$nrow,
                ncol = sparse_sizes()[[size_label]]$ncol,
                rhs_width = rhs_width,
                density = density
              )
            }
          }
        }
      }
    }
  }

  do.call(rbind, rows)
}

run_group <- function(group) {
  backend_name <- group$requested_backend
  backend <- tryCatch(amatrix:::.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !isTRUE(backend$available())) {
    rows <- expand_group_rows(group)
    rows$status <- "unavailable"
    rows$error_message <- sprintf(
      "backend '%s' is not available [%s]",
      backend_name,
      format_backend_diagnostics(backend_name)
    )
    return(rows)
  }

  rows <- list()
  if (identical(group$suite, "dense")) {
    for (size_label in names(dense_sizes())) {
      for (variant in c("cold", "warm")) {
        cell <- list(suite = "dense", op = group$op, size_label = size_label, variant = variant, requested_backend = backend_name)
        rows[[length(rows) + 1L]] <- safe_case(
          function() run_dense_case(backend_name, group$op, size_label, variant),
          cell
        )
      }
    }
  }

  if (identical(group$suite, "sparse")) {
    if (group$op %in% sparse_iterative_ops()) {
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_iterative_densities()) {
          for (variant in c("cold", "warm")) {
            params <- sparse_iterative_params(group$op, size_label)
            cell <- list(
              suite = "sparse",
              op = group$op,
              size_label = size_label,
              variant = variant,
              requested_backend = backend_name,
              rhs_width = params$k,
              density = density,
              nrow = sparse_sizes()[[size_label]]$nrow,
              ncol = sparse_sizes()[[size_label]]$ncol
            )
            rows[[length(rows) + 1L]] <- safe_case(
              function() run_sparse_iterative_case(backend_name, group$op, size_label, density, variant),
              cell
            )
          }
        }
      }
    } else {
      rhs_widths <- if (identical(group$op, "spmv")) 1L else spmm_rhs_widths()
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_densities()) {
          for (rhs_width in rhs_widths) {
            for (variant in c("cold", "resident")) {
              cell <- list(
                suite = "sparse",
                op = group$op,
                size_label = size_label,
                variant = variant,
                requested_backend = backend_name,
                rhs_width = rhs_width,
                density = density,
                nrow = sparse_sizes()[[size_label]]$nrow,
                ncol = sparse_sizes()[[size_label]]$ncol
              )
              rows[[length(rows) + 1L]] <- safe_case(
                function() run_sparse_case(backend_name, group$op, size_label, density, rhs_width, variant),
                cell
              )
            }
          }
        }
      }
    }
  }

  do.call(rbind, rows)
}

key_columns <- c(
  "suite", "op", "size_label", "variant", "requested_backend",
  "dispatch_backend", "dispatch_path", "nrow", "ncol", "rhs_width",
  "density", "nnz", "density_bucket"
)

add_cpu_reference <- function(results) {
  cpu <- results[
    results$status == "ok" &
      results$requested_backend == "cpu",
    c("suite", "op", "size_label", "variant", "nrow", "ncol", "rhs_width", "density", "nnz", "density_bucket", "median_ms")
  ]
  names(cpu)[names(cpu) == "median_ms"] <- "cpu_reference_ms"

  merge(
    results,
    cpu,
    by = c("suite", "op", "size_label", "variant", "nrow", "ncol", "rhs_width", "density", "nnz", "density_bucket"),
    all.x = TRUE
  )
}

add_baseline_compare <- function(results, baseline_path) {
  results$baseline_ms <- NA_real_
  results$ratio_vs_baseline <- NA_real_

  if (!file.exists(baseline_path)) {
    return(results)
  }

  baseline <- tryCatch(read.csv(baseline_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(baseline) || !"median_ms" %in% names(baseline)) {
    return(results)
  }

  needed <- c(key_columns, "status", "median_ms")
  if (!all(needed %in% names(baseline))) {
    return(results)
  }

  baseline <- baseline[baseline$status == "ok", c(key_columns, "median_ms")]
  names(baseline)[names(baseline) == "median_ms"] <- "baseline_ms"

  merged <- merge(results, baseline, by = key_columns, all.x = TRUE)
  merged$ratio_vs_baseline <- ifelse(
    is.na(merged$baseline_ms) | is.na(merged$median_ms),
    NA_real_,
    merged$median_ms / merged$baseline_ms
  )
  merged
}

summarize_results <- function(results) {
  display <- results[order(results$suite, results$op, results$size_label, results$variant, results$requested_backend), ]
  display$speedup_vs_cpu <- ifelse(
    is.na(display$cpu_reference_ms) | is.na(display$median_ms),
    NA_real_,
    display$cpu_reference_ms / display$median_ms
  )
  display$dispatch_state <- ifelse(
    display$status != "ok",
    display$status,
    ifelse(
      display$requested_backend == "cpu",
      "cpu_baseline",
      ifelse(
        is.na(display$dispatch_backend),
        "unknown_dispatch",
        ifelse(
          display$dispatch_backend == display$requested_backend,
          "accelerated",
          ifelse(display$dispatch_backend == "cpu", "cpu_fallback", "rerouted")
        )
      )
    )
  )
  display
}

write_outputs <- function(results, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  raw_path <- file.path(output_dir, "raw-results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  incidents_path <- file.path(output_dir, "incidents.csv")
  fallbacks_path <- file.path(output_dir, "fallbacks.csv")
  meta_path <- file.path(output_dir, "metadata.rds")

  summary <- summarize_results(results)
  incidents <- summary[summary$status != "ok", , drop = FALSE]
  fallbacks <- summary[
    summary$status == "ok" &
      summary$requested_backend != "cpu" &
      summary$dispatch_state == "cpu_fallback",
    ,
    drop = FALSE
  ]

  write.csv(results, raw_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)
  write.csv(incidents, incidents_path, row.names = FALSE)
  write.csv(fallbacks, fallbacks_path, row.names = FALSE)
  saveRDS(
    list(
      created_at = Sys.time(),
      hostname = Sys.info()[["nodename"]],
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      timezone = Sys.timezone()
    ),
    meta_path
  )

  list(raw = raw_path, summary = summary_path, incidents = incidents_path, fallbacks = fallbacks_path, metadata = meta_path)
}

print_run_summary <- function(results, baseline_path, output_paths, warm_only = FALSE) {
  summary <- summarize_results(results)
  display <- if (warm_only) {
    summary[summary$variant %in% c("warm", "resident"), , drop = FALSE]
  } else {
    summary
  }

  regressions <- display[
    display$status == "ok" &
      display$dispatch_state != "cpu_fallback" &
      !is.na(display$ratio_vs_baseline) &
      display$ratio_vs_baseline > 1.2,
    c("suite", "op", "size_label", "variant", "requested_backend", "median_ms", "baseline_ms", "ratio_vs_baseline")
  ]

  incidents <- display[display$status != "ok", c("suite", "op", "size_label", "variant", "requested_backend", "status", "error_message")]
  fallbacks <- display[
    display$status == "ok" &
      display$requested_backend != "cpu" &
      display$dispatch_state == "cpu_fallback",
    c("suite", "op", "size_label", "variant", "requested_backend", "dispatch_backend", "median_ms", "cpu_reference_ms")
  ]

  if (nrow(regressions) > 0L) {
    message("\nRegressions (>20% slower than baseline):")
    print(regressions, row.names = FALSE)
  } else if (file.exists(baseline_path)) {
    message("\nOK - no regressions vs baseline.")
  } else {
    message("\nNo baseline file found; skipped regression comparison.")
  }

  if (nrow(incidents) > 0L) {
    message("\nIncidents:")
    print(incidents, row.names = FALSE)
  } else {
    message("\nNo incidents recorded.")
  }

  ok_rows <- display[
    display$status == "ok" & display$dispatch_state != "cpu_fallback",
    c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "dispatch_backend", "dispatch_state", "rhs_width", "density_bucket", "median_ms",
      "cpu_reference_ms", "speedup_vs_cpu", "ratio_vs_baseline"
    )
  ]
  message("\nSummary:")
  print(ok_rows, row.names = FALSE)

  if (nrow(fallbacks) > 0L) {
    message("\nCPU Fallbacks:")
    print(fallbacks, row.names = FALSE)
  }

  message("\nArtifacts:")
  message("  raw: ", output_paths$raw)
  message("  summary: ", output_paths$summary)
  message("  incidents: ", output_paths$incidents)
  message("  fallbacks: ", output_paths$fallbacks)
  message("  metadata: ", output_paths$metadata)
}

run_worker <- function(args) {
  if (is.null(args$plan) || is.null(args$out) || is.null(args$group_id)) {
    stop("worker mode requires --plan, --group-id, and --out", call. = FALSE)
  }

  groups <- readRDS(args$plan)
  group <- Filter(function(x) identical(x$group_id, args$group_id), groups)[[1L]]
  prime_requested_backend(group$requested_backend, include_arrayfire = args$include_arrayfire)
  result <- run_group(group)
  saveRDS(result, args$out)
  invisible(NULL)
}

run_master <- function(args) {
  specs <- canonical_backend_specs(include_arrayfire = args$include_arrayfire)
  specs <- filter_backend_specs(specs, args$backends)
  detected_backend_names <- vapply(specs, `[[`, character(1), "name")
  backend_names <- if (is.null(args$backends)) {
    detected_backend_names
  } else {
    unique(args$backends[nzchar(args$backends)])
  }

  if (!is.null(args$backends)) {
    missing_backends <- setdiff(backend_names, detected_backend_names)
    if (length(missing_backends) > 0L) {
      message(
        "Requested backends unresolved during master discovery; worker validation will decide availability: ",
        paste(missing_backends, collapse = ", ")
      )
    }
  }

  groups <- group_plan(backend_names, suites = args$suites)
  plan_path <- tempfile("benchmark-plan-", fileext = ".rds")
  saveRDS(groups, plan_path)
  on.exit(unlink(plan_path), add = TRUE)

  message(
    "Running canonical benchmark harness for backends: ",
    paste(backend_names, collapse = ", "),
    " | suites: ", paste(args$suites, collapse = ", ")
  )

  rows <- list()
  script_path <- normalizePath(sys.frame(1)$ofile %||% "tools/benchmark-regression.R", mustWork = FALSE)

  for (group in groups) {
    out_path <- tempfile(sprintf("%s-", group$group_id), fileext = ".rds")
    message(sprintf("  %s", group$group_id))
    launch <- benchmark_system2_capture(
      file.path(R.home("bin"), "Rscript"),
      benchmark_rscript_source_args(
        script_path,
        args = c("--worker", paste0("--plan=", plan_path), paste0("--group-id=", group$group_id), paste0("--out=", out_path))
      )
    )
    log <- launch$output
    status <- launch$status

    if (status == 0L && file.exists(out_path)) {
      rows[[length(rows) + 1L]] <- readRDS(out_path)
    } else {
      failed <- expand_group_rows(group)
      failed$status <- "crash"
      failed$error_message <- sprintf(
        "worker exited with status %s%s",
        status,
        if (length(log) > 0L) paste0(": ", tail(log, 1L)) else ""
      )
      rows[[length(rows) + 1L]] <- failed
    }

    unlink(out_path)
  }

  results <- do.call(rbind, rows)
  results <- add_cpu_reference(results)
  results$speedup_vs_cpu <- ifelse(
    is.na(results$cpu_reference_ms) | is.na(results$median_ms),
    NA_real_,
    results$cpu_reference_ms / results$median_ms
  )
  results <- add_baseline_compare(results, args$baseline)
  results <- results[order(results$suite, results$op, results$size_label, results$variant, results$requested_backend), ]

  output_paths <- write_outputs(results, args$output_dir)

  if (args$update || !file.exists(args$baseline)) {
    baseline_rows <- results[results$status == "ok", c(key_columns, "status", "median_ms")]
    write.csv(baseline_rows, args$baseline, row.names = FALSE)
    message("\nBaseline written to ", args$baseline)
  }

  print_run_summary(results, args$baseline, output_paths, warm_only = args$warm_only)
  invisible(results)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
relaunch_safe_master_if_needed(args)
initialize_regression_benchmark_context()

if (isTRUE(args$worker)) {
  run_worker(args)
} else {
  run_master(args)
}
