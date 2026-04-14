#' Collect full dispatch information for an aMatrix object
#'
#' Returns a snapshot of the dispatch state for an \code{aMatrix},
#' including residency, preferred backend, policy, precision, and the
#' per-operation dispatch matrix for a set of operations.
#'
#' @param x An \code{aMatrix} object.
#' @param ops Character vector of operation names to include in the
#'   dispatch matrix. Default covers the six core operations.
#' @param y_map Named list mapping operation names to right-hand-side
#'   objects used when planning binary operations such as
#'   \code{"matmul"}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{object_id}{Character. Internal object identifier.}
#'     \item{preferred_backend}{Character. Preferred backend slot value.}
#'     \item{pinned_backend}{Character or \code{NULL}. Backend to which
#'       the object is currently GPU-resident.}
#'     \item{policy}{Character. Dispatch policy slot value.}
#'     \item{precision}{Character. Precision mode (\code{"strict"} or
#'       \code{"fast"}).}
#'     \item{residency}{data.frame. Output of
#'       \code{\link{amatrix_residency_info}}.}
#'     \item{plans}{data.frame. Output of
#'       \code{\link{amatrix_backend_matrix}}.}
#'   }
#'
#' @seealso \code{\link{amatrix_backend_plan}},
#'   \code{\link{amatrix_backend_matrix}},
#'   \code{\link{amatrix_explain}}
#' @export
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
  resident_backend <- .amatrix_live_resident_backend(x)
  if (!is.null(resident_backend)) {
    backend <- tryCatch(.amatrix_get_backend(resident_backend), error = function(e) NULL)
    if (!is.null(backend) &&
        isTRUE(backend$available()) &&
        x@precision %in% unique(backend$precision_modes()) &&
        .amatrix_backend_supports_resident_op(backend, op, x = x, y = y)) {
      return(list(name = resident_backend, backend = backend))
    }
  }

  plan <- amatrix_backend_plan(x, op, y = y)
  chosen <- plan$candidates[[match(TRUE, vapply(plan$candidates, function(candidate) isTRUE(candidate$chosen), logical(1)))]]

  list(name = chosen$name, backend = .amatrix_get_backend(chosen$name))
}

#' Compute the dispatch plan for a single operation
#'
#' Evaluates each candidate backend in preference order and returns a
#' structured plan describing which backend was chosen and why each
#' candidate was accepted or rejected. The plan respects GPU residency,
#' precision compatibility, and calibration thresholds.
#'
#' @param x An \code{aMatrix} object.
#' @param op Character string naming the operation, e.g.
#'   \code{"matmul"}, \code{"crossprod"}, \code{"svd"}.
#' @param y Right-hand-side \code{aMatrix} or \code{NULL}. Used for
#'   binary operations such as \code{"matmul"} to check compatibility
#'   and calibration workload.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{op}{Character. The requested operation.}
#'     \item{pinned_backend}{Character or \code{NULL}. Backend to which
#'       \code{x} is currently GPU-resident.}
#'     \item{preferred}{Character vector. Backends evaluated in order.}
#'     \item{requested_precision}{Character. Precision mode of \code{x}.}
#'     \item{chosen}{Character. Name of the chosen backend.}
#'     \item{chosen_path}{Character. Either \code{"resident"} or
#'       \code{"cold"}.}
#'     \item{candidates}{List of per-candidate evaluation records, each
#'       a named list with logical flags for \code{registered},
#'       \code{available}, \code{precision_compatible},
#'       \code{supported_cold}, \code{supported_resident},
#'       \code{calibration_ok}, \code{supported}, and \code{chosen}.}
#'   }
#'
#' @examples
#' m <- adgeMatrix(matrix(1:6, 2, 3))
#' amatrix_backend_plan(m, "matmul")
#'
#' @seealso \code{\link{amatrix_backend_matrix}},
#'   \code{\link{amatrix_explain}},
#'   \code{\link{amatrix_execution_info}}
#' @export
amatrix_backend_plan <- function(x, op, y = NULL) {
  pinned_backend <- .amatrix_live_resident_backend(x)
  preferred <- .amatrix_backend_preference(x, op = op)
  candidates <- vector("list", length(preferred))
  found <- FALSE

  for (idx in seq_along(preferred)) {
    candidate_name <- preferred[[idx]]
    backend <- tryCatch(.amatrix_get_backend(candidate_name), error = function(e) NULL)
    entry <- list(
      name = candidate_name,
      registered = !is.null(backend),
      capabilities = character(),
      features = character(),
      precision_modes = character(),
      available = FALSE,
      precision_compatible = FALSE,
      resident_active = FALSE,
      supported_cold = FALSE,
      supported_resident = FALSE,
      calibration_ok = TRUE,
      supported = FALSE,
      chosen_path = NA_character_,
      chosen = FALSE
    )

    if (entry$registered) {
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
            .amatrix_backend_supports_resident_op(backend, op, x = x, y = y)
        )
        # Calibration gates the cold path only. Resident path is always ok
        # because the upload cost has already been paid.
        entry$calibration_ok <- (
          entry$supported_resident ||
          .amatrix_calibration_ok(x, op, candidate_name, y = y)
        )
        entry$supported <- isTRUE(
          (entry$supported_cold || entry$supported_resident) && entry$calibration_ok
        )
        if (entry$supported) {
          entry$chosen_path <- if (isTRUE(entry$supported_resident)) "resident" else "cold"
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
        calibration_ok = TRUE,
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
      candidates[[cpu_idx]]$calibration_ok <- TRUE
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

#' Tabulate dispatch plans across multiple operations
#'
#' Runs \code{\link{amatrix_backend_plan}} for each requested operation
#' and returns the results as a single data.frame, one row per
#' operation. Useful for inspecting which backend will be used across an
#' entire workload.
#'
#' @param x An \code{aMatrix} object.
#' @param ops Character vector of operation names. Defaults to the
#'   twelve standard operations.
#' @param y_map Named list mapping operation names to right-hand-side
#'   objects. Use to supply a \code{y} argument for binary operations
#'   such as \code{"matmul"}.
#'
#' @return A data.frame with one row per operation and columns:
#'   \describe{
#'     \item{op}{Character. Operation name.}
#'     \item{precision}{Character. Precision mode.}
#'     \item{pinned_backend}{Character. Backend to which \code{x} is
#'       GPU-resident, or \code{NA}.}
#'     \item{preferred}{Character. Preference order string.}
#'     \item{chosen}{Character. Selected backend.}
#'     \item{chosen_path}{Character. \code{"resident"} or
#'       \code{"cold"}.}
#'     \item{resident_reuse}{Logical. Whether the resident path is
#'       active.}
#'     \item{cpu_fallback}{Logical. Whether CPU was chosen despite not
#'       being first preference.}
#'     \item{candidate_summary}{Character. Compact flag string for all
#'       candidates.}
#'   }
#'
#' @examples
#' m <- adgeMatrix(matrix(1:6, 2, 3))
#' amatrix_backend_matrix(m, ops = c("matmul", "crossprod"))
#'
#' @seealso \code{\link{amatrix_backend_plan}},
#'   \code{\link{amatrix_execution_info}}
#' @export
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
      # cpu_fallback is TRUE iff cpu was chosen AND the user's first
      # preference was NOT cpu. Matches the authoritative semantics from
      # amatrix-15n (tests/testthat/test-regression-backend-preference-summary.R).
      cpu_fallback = identical(plan$chosen, "cpu") &&
        length(plan$preferred) >= 1L &&
        !identical(plan$preferred[[1]], "cpu"),
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
              if (isTRUE(candidate$calibration_ok)) "K" else "-",
              if (candidate$supported) "S" else "-",
              if (candidate$chosen) "X" else "-",
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

#' Low-level backend dispatch for a single operation
#'
#' Resolves the best available backend for \code{op} on \code{x},
#' attempts the GPU-resident path when applicable, and falls back to
#' the cold path (materializing \code{x} to host) if needed. If the
#' chosen backend does not implement \code{method}, the \code{fallback}
#' function is called instead.
#'
#' @param x An \code{aMatrix} object.
#' @param op Character string. Operation key used for backend selection
#'   (e.g. \code{"matmul"}, \code{"svd"}).
#' @param method Character string. Name of the backend list element to
#'   call. Defaults to \code{op}; override when the backend method name
#'   differs from the operation key.
#' @param y Right-hand-side \code{aMatrix} or \code{NULL}. Passed to
#'   the backend method and used during backend selection.
#' @param args Named list of additional arguments forwarded to the
#'   backend method on the cold path.
#' @param fallback Zero-argument function called when the chosen
#'   backend does not implement \code{method}.
#'
#' @return The result of the backend method, or the result of
#'   \code{fallback()} if the method is unavailable.
#'
#' @seealso \code{\link{amatrix_backend_plan}},
#'   \code{\link{amatrix_materialize_host}}
#' @export
amatrix_dispatch_op <- function(x, op, method = op, y = NULL, args = list(), fallback) {
  stopifnot(is.function(fallback))
  choice <- .amatrix_backend_for(x, op, y = y)
  chosen_name <- choice$name
  backend_method <- choice$backend[[method]]

  if (!is.function(backend_method)) {
    if (!identical(chosen_name, "cpu")) {
      .amatrix_log_fallback(
        op = op,
        backend = chosen_name,
        reason = sprintf("backend '%s' does not implement method '%s'",
                         chosen_name, method),
        from_backend = chosen_name
      )
    }
    return(fallback())
  }

  # Try resident path if x is GPU-resident and the backend supports this op
  # as a resident operation. This avoids dropping the binding + re-uploading.
  resident_backend_name <- .amatrix_live_resident_backend(x)
  if (!is.null(resident_backend_name)) {
    backend <- .amatrix_get_backend(resident_backend_name)
    resident_op_name <- paste0(method, "_resident")
    if (is.function(backend[[resident_op_name]]) &&
        .amatrix_backend_supports_resident_op(backend, method, x = x, y = y)) {
      lhs <- .amatrix_prepare_resident_arg(x, resident_backend_name)
      if (!is.null(lhs)) {
        out_key <- .amatrix_next_resident_key(resident_backend_name)
        result <- tryCatch(
          backend[[resident_op_name]](lhs$key, out_key),
          error = function(e) {
            .amatrix_log_fallback(
              op = op,
              backend = resident_backend_name,
              reason = sprintf("resident %s error: %s", method, conditionMessage(e)),
              from_backend = resident_backend_name
            )
            .amatrix_backend_health_mark(
              resident_backend_name, "unhealthy",
              sprintf("resident %s error: %s", method, conditionMessage(e))
            )
            NULL
          }
        )
        .amatrix_cleanup_temp_resident(list(lhs), resident_backend_name)
        if (!is.null(result)) {
          return(result)
        }
        # Clean up out_key on failure
        try(backend$resident_drop(out_key), silent = TRUE)
      }
    }

    # No resident path available or it failed — materialize but preserve
    # the binding so future ops on x don't need to re-upload.
  }

  # Wrap the cold-path backend call in tryCatch: on runtime error we
  #   (1) re-signal the ORIGINAL backend condition via signalCondition() so
  #       outer withCallingHandlers can observe the original classed error
  #       (amatrix-uu2 contract)
  #   (2) signal an amatrix_fallback condition with structured metadata so
  #       callers can observe fallback events (amatrix-hjj contract)
  #   (3) log a structured fallback event (internal telemetry)
  #   (4) mark the backend unhealthy
  #   (5) re-dispatch via the caller-supplied fallback
  # See planning_docs/quality-tracking.md §7 rule 7 (non-empty fallback log
  # is stop-ship).
  tryCatch(
    do.call(backend_method, c(list(x = amatrix_materialize_host(x)), args)),
    error = function(e) {
      if (!identical(chosen_name, "cpu")) {
        reason <- sprintf("%s runtime error: %s", method, conditionMessage(e))

        # (1) Emit a fresh observation condition carrying the ORIGINAL classes
        # of the backend error. Must not inherit from "error" because (a) we
        # have already handled the real error in this tryCatch, and (b) R's
        # signalCondition falls through to stop() on unhandled error
        # conditions and does not reach outer withCallingHandlers from inside
        # a consumed tryCatch error-handler frame. Stripping "error" gives us
        # a plain calling-handler-only signal.
        orig_classes <- setdiff(class(e), c("error", "simpleError", "condition"))
        if (length(orig_classes) > 0L) {
          observed_cond <- structure(
            class = c(orig_classes, "amatrix_observed_backend_error", "condition"),
            list(
              message = conditionMessage(e),
              call = NULL,
              op = op,
              backend = chosen_name,
              original_condition = e
            )
          )
          tryCatch(signalCondition(observed_cond), error = function(ee) NULL)
        }

        # (2) Signal a classed amatrix_fallback condition (not an error) so
        # observers can catch fallback events generically.
        fallback_cond <- structure(
          class = c("amatrix_fallback", "condition"),
          list(
            message = reason,
            call = NULL,
            op = op,
            backend = chosen_name,
            reason = reason,
            original_condition = e
          )
        )
        tryCatch(signalCondition(fallback_cond), error = function(ee) NULL)

        # (3) Internal telemetry log.
        .amatrix_log_fallback(
          op = op,
          backend = chosen_name,
          reason = reason,
          from_backend = chosen_name
        )

        # (4) Mark the backend unhealthy.
        .amatrix_backend_health_mark(chosen_name, "unhealthy", reason)

        # (5) Run the caller's fallback.
        return(fallback())
      }
      stop(e)  # CPU errors are real bugs; don't swallow them.
    }
  )
}
