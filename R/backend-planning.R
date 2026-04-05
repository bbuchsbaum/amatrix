amatrix_execution_info <- function(
  x,
  ops = c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums"),
  y_map = list()
) {
  stopifnot(inherits(x, "aMatrix"))

  residency <- amatrix_residency_info(x)
  plans <- amatrix_backend_matrix(x, ops = ops, y_map = y_map)

  list(
    object_id = x@object_id,
    preferred_backend = x@preferred_backend,
    pinned_backend = .amatrix_live_resident_backend(x),
    policy = x@policy,
    precision = x@precision,
    residency = residency,
    plans = plans
  )
}

.amatrix_backend_for <- function(x, op, y = NULL) {
  plan <- amatrix_backend_plan(x, op, y = y)
  chosen <- plan$candidates[[match(TRUE, vapply(plan$candidates, function(candidate) isTRUE(candidate$chosen), logical(1)))]]

  list(name = chosen$name, backend = .amatrix_get_backend(chosen$name))
}

amatrix_backend_plan <- function(x, op, y = NULL) {
  pinned_backend <- .amatrix_live_resident_backend(x)
  preferred <- .amatrix_backend_preference(x, op = op)
  candidates <- vector("list", length(preferred))
  found <- FALSE

  for (idx in seq_along(preferred)) {
    candidate_name <- preferred[[idx]]
    entry <- list(
      name = candidate_name,
      registered = candidate_name %in% amatrix_backend_names(),
      capabilities = character(),
      features = character(),
      precision_modes = character(),
      available = FALSE,
      precision_compatible = FALSE,
      resident_active = FALSE,
      supported_cold = FALSE,
      supported_resident = FALSE,
      supported = FALSE,
      chosen_path = NA_character_,
      chosen = FALSE
    )

    if (entry$registered) {
      backend <- .amatrix_get_backend(candidate_name)
      entry$capabilities <- unique(backend$capabilities())
      entry$features <- unique(backend$features())
      entry$precision_modes <- unique(backend$precision_modes())
      entry$available <- isTRUE(backend$available())
      entry$precision_compatible <- x@precision %in% entry$precision_modes
      if (entry$available && entry$precision_compatible) {
        entry$supported_cold <- isTRUE(backend$supports(op = op, x = x, y = y))
        entry$resident_active <- .amatrix_object_is_resident(x, candidate_name)
        entry$supported_resident <- (
          entry$resident_active &&
            .amatrix_backend_residency_capable(backend) &&
            .amatrix_backend_supports_resident_op(backend, op)
        )
        entry$supported <- isTRUE(entry$supported_cold || entry$supported_resident)
        if (entry$supported) {
          entry$chosen_path <- if (isTRUE(entry$supported_cold)) "cold" else "resident"
        }
      }
    }

    if (!found && entry$registered && entry$available && entry$precision_compatible && entry$supported) {
      entry$chosen <- TRUE
      found <- TRUE
    }

    candidates[[idx]] <- entry
  }

  if (!found) {
    cpu_idx <- match("cpu", vapply(candidates, `[[`, character(1), "name"))
    if (is.na(cpu_idx)) {
      candidates[[length(candidates) + 1L]] <- list(
        name = "cpu",
        registered = TRUE,
        capabilities = amatrix_backend_capabilities("cpu"),
        features = amatrix_backend_features("cpu"),
        precision_modes = amatrix_backend_precision_modes("cpu"),
        available = TRUE,
        precision_compatible = TRUE,
        resident_active = FALSE,
        supported_cold = TRUE,
        supported_resident = FALSE,
        supported = TRUE,
        chosen_path = "cold",
        chosen = TRUE
      )
    } else {
      candidates[[cpu_idx]]$chosen <- TRUE
      candidates[[cpu_idx]]$registered <- TRUE
      candidates[[cpu_idx]]$capabilities <- amatrix_backend_capabilities("cpu")
      candidates[[cpu_idx]]$features <- amatrix_backend_features("cpu")
      candidates[[cpu_idx]]$precision_modes <- amatrix_backend_precision_modes("cpu")
      candidates[[cpu_idx]]$available <- TRUE
      candidates[[cpu_idx]]$precision_compatible <- TRUE
      candidates[[cpu_idx]]$resident_active <- FALSE
      candidates[[cpu_idx]]$supported_cold <- TRUE
      candidates[[cpu_idx]]$supported_resident <- FALSE
      candidates[[cpu_idx]]$supported <- TRUE
      candidates[[cpu_idx]]$chosen_path <- "cold"
    }
  }

  chosen_idx <- match(TRUE, vapply(candidates, function(candidate) isTRUE(candidate$chosen), logical(1)))

  list(
    op = op,
    pinned_backend = pinned_backend,
    preferred = preferred,
    requested_precision = x@precision,
    chosen = candidates[[chosen_idx]]$name,
    chosen_path = candidates[[chosen_idx]]$chosen_path,
    candidates = candidates
  )
}

amatrix_backend_matrix <- function(
  x,
  ops = c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums", "solve", "chol", "qr", "svd", "eigen", "diag"),
  y_map = list()
) {
  stopifnot(is.character(ops), length(ops) >= 1L)

  rows <- lapply(ops, function(op) {
    y <- if (!is.null(y_map[[op]])) y_map[[op]] else NULL
    plan <- amatrix_backend_plan(x, op, y = y)

    data.frame(
      op = op,
      precision = plan$requested_precision,
      pinned_backend = if (is.null(plan$pinned_backend)) NA_character_ else plan$pinned_backend,
      preferred = paste(plan$preferred, collapse = " > "),
      chosen = plan$chosen,
      chosen_path = plan$chosen_path,
      resident_reuse = identical(plan$chosen_path, "resident"),
      cpu_fallback = identical(plan$chosen, "cpu") && !identical(plan$preferred[[1]], "cpu"),
      candidate_summary = paste(
        vapply(
          plan$candidates,
          function(candidate) {
            paste0(
              candidate$name,
              "[",
              if (candidate$registered) "R" else "-",
              if (candidate$available) "A" else "-",
              if (candidate$precision_compatible) "P" else "-",
              if (candidate$resident_active) "r" else "-",
              if (candidate$supported_cold) "C" else "-",
              if (candidate$supported_resident) "D" else "-",
              if (candidate$supported) "S" else "-",
              if (candidate$chosen) "C" else "-",
              "]"
            )
          },
          character(1)
        ),
        collapse = " "
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

amatrix_dispatch_op <- function(x, op, method = op, y = NULL, args = list(), fallback) {
  stopifnot(is.function(fallback))
  choice <- .amatrix_backend_for(x, op, y = y)
  backend_method <- choice$backend[[method]]

  if (!is.function(backend_method)) {
    return(fallback())
  }

  if (!is.null(.amatrix_live_resident_backend(x))) {
    .amatrix_drop_resident_binding(x)
  }

  do.call(backend_method, c(list(x = amatrix_materialize_host(x)), args))
}
