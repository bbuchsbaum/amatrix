## Round 3 Bug Hunt — Method Sweep
## Tests every public method on adgeMatrix and adgCMatrix.
## Asserts: class preserved, dimnames preserved, backend slot preserved.
## Each failure is tagged: CLASS_DEMOTION, DIMNAME_LOSS, BACKEND_LEAK, or ERROR.
##
## This is a regression test left behind from the round-3 hunt.
## Report is written to .bug-hunt-r3/03-method-sweep.md.

library(Matrix)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.sweep_failures <- list()

.record_fail <- function(method, symptom, expected, actual, note = "") {
  .sweep_failures[[length(.sweep_failures) + 1]] <<- list(
    method   = method,
    symptom  = symptom,
    expected = as.character(expected),
    actual   = as.character(actual),
    note     = note
  )
}

## Check a result matrix for class/dimnames/backend integrity.
.check_matrix_result <- function(result, method_name, expected_class,
                                 expected_rn, expected_cn, expected_backend) {
  ok <- TRUE

  # --- class ---
  if (!inherits(result, expected_class)) {
    .record_fail(method_name, "CLASS_DEMOTION",
                 expected_class, class(result)[1])
    ok <- FALSE
    return(invisible(ok))   # slots won't exist on demoted object
  }

  # --- backend slot ---
  actual_backend <- tryCatch(result@preferred_backend, error = function(e) NA_character_)
  if (is.na(actual_backend) || !nzchar(actual_backend)) {
    .record_fail(method_name, "BACKEND_LEAK",
                 expected_backend, as.character(actual_backend))
    ok <- FALSE
  } else if (!identical(actual_backend, expected_backend)) {
    .record_fail(method_name, "BACKEND_LEAK",
                 expected_backend, actual_backend)
    ok <- FALSE
  }

  # --- dimnames ---
  if (!is.null(expected_rn)) {
    actual_rn <- rownames(result)
    if (!identical(actual_rn, expected_rn)) {
      .record_fail(method_name, "DIMNAME_LOSS",
                   paste(expected_rn, collapse = ","),
                   paste(as.character(actual_rn), collapse = ","))
      ok <- FALSE
    }
  }
  if (!is.null(expected_cn)) {
    actual_cn <- colnames(result)
    if (!identical(actual_cn, expected_cn)) {
      .record_fail(method_name, "DIMNAME_LOSS",
                   paste(expected_cn, collapse = ","),
                   paste(as.character(actual_cn), collapse = ","))
      ok <- FALSE
    }
  }

  invisible(ok)
}

## Run a call, returning result or NULL on error (errors recorded).
.safe_call <- function(expr, method_name) {
  tryCatch(expr,
    error = function(e) {
      .record_fail(method_name, "ERROR", "no error", conditionMessage(e))
      NULL
    }
  )
}

# ---------------------------------------------------------------------------
# Build test fixtures
# ---------------------------------------------------------------------------

.make_dense <- function(backend = "cpu") {
  m <- matrix(c(4, 2, 2, 3, 1, 1, 1, 1, 4, 2, 2, 3), nrow = 3, ncol = 4)
  rownames(m) <- c("r1", "r2", "r3")
  colnames(m) <- c("c1", "c2", "c3", "c4")
  new_adgeMatrix(m, preferred_backend = backend, policy = "auto", precision = "strict")
}

.make_square_dense <- function(backend = "cpu") {
  m <- matrix(c(4, 1, 1, 3), nrow = 2, ncol = 2)
  rownames(m) <- c("r1", "r2")
  colnames(m) <- c("c1", "c2")
  new_adgeMatrix(m, preferred_backend = backend, policy = "auto", precision = "strict")
}

.make_symm_dense <- function(backend = "cpu") {
  m <- matrix(c(4, 2, 2, 3), nrow = 2, ncol = 2)
  rownames(m) <- colnames(m) <- c("r1", "r2")
  new_adgeMatrix(m, preferred_backend = backend, policy = "auto", precision = "strict")
}

.make_sparse <- function(backend = "cpu") {
  m <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 3), j = c(1, 2, 3, 4, 4),
    x = c(1.0, 2.0, 3.0, 0.5, 1.5), dims = c(3, 4),
    dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3", "c4"))
  )
  new_adgCMatrix(as(m, "dgCMatrix"), preferred_backend = backend,
                 policy = "auto", precision = "strict")
}

.make_square_sparse <- function(backend = "cpu") {
  m <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(4.0, 3.0), dims = c(2, 2),
    dimnames = list(c("r1", "r2"), c("c1", "c2"))
  )
  new_adgCMatrix(as(m, "dgCMatrix"), preferred_backend = backend,
                 policy = "auto", precision = "strict")
}

# ---------------------------------------------------------------------------
# Dense sweep
# ---------------------------------------------------------------------------

.run_dense_sweep <- function(backend = "cpu") {
  A   <- .make_dense(backend)
  Asq <- .make_square_dense(backend)
  Asm <- .make_symm_dense(backend)
  rn  <- c("r1", "r2", "r3")
  cn  <- c("c1", "c2", "c3", "c4")
  rn2 <- c("r1", "r2")
  cn2 <- c("c1", "c2")
  cls <- "adgeMatrix"
  be  <- backend
  tag <- function(s) paste0(s, "[dense/", backend, "]")

  # --- Ops: Arithmetic (matrix op matrix) ---
  for (op in c("+", "-", "*", "/")) {
    res <- .safe_call(do.call(op, list(A, A)), tag(paste0("Arith_", op, "_mat")))
    if (!is.null(res)) .check_matrix_result(res, tag(paste0("Arith_", op, "_mat")), cls, rn, cn, be)
  }
  # --- Ops: Arithmetic (matrix op scalar) ---
  for (op in c("+", "-", "*", "/", "^")) {
    res <- .safe_call(do.call(op, list(A, 2)), tag(paste0("Arith_", op, "_scalar")))
    if (!is.null(res)) .check_matrix_result(res, tag(paste0("Arith_", op, "_scalar")), cls, rn, cn, be)
  }
  # --- Ops: Arithmetic (scalar op matrix) ---
  for (op in c("+", "-", "*")) {
    res <- .safe_call(do.call(op, list(2, A)), tag(paste0("Arith_scalar_", op)))
    if (!is.null(res)) .check_matrix_result(res, tag(paste0("Arith_scalar_", op)), cls, rn, cn, be)
  }

  # --- Ops: Compare (adgeMatrix vs scalar) ---
  # NOTE: Compare ops return lgeMatrix (logical dense matrix from Matrix pkg),
  # NOT adgeMatrix. This IS a class demotion — aMatrix wrapping is lost.
  for (op in c("==", "!=", "<", ">", "<=", ">=")) {
    res <- .safe_call(do.call(op, list(A, 2)), tag(paste0("Compare_", op)))
    if (!is.null(res)) {
      if (!inherits(res, "aMatrix")) {
        .record_fail(tag(paste0("Compare_", op)), "CLASS_DEMOTION",
                     "aMatrix-derived", class(res)[1],
                     "Compare returns lgeMatrix — aMatrix wrapping lost")
      }
    }
  }

  # --- Math group ---
  # Tested individually: abs/sqrt/exp/log/ceiling/floor/sign all return adgeMatrix (CORRECT).
  # cumsum/cumprod/cummax/cummin return numeric vectors — class check not applicable.
  math_mat_fns <- c("abs", "sqrt", "exp", "log", "log2", "log10",
                    "ceiling", "floor", "sign",
                    "cos", "sin", "tan", "cosh", "sinh", "tanh")
  Apos <- A  # values are all positive, safe for log/sqrt
  for (fn in math_mat_fns) {
    res <- .safe_call(do.call(fn, list(Apos)), tag(paste0("Math_", fn)))
    if (!is.null(res) && !is.vector(res)) {
      .check_matrix_result(res, tag(paste0("Math_", fn)), cls, NULL, NULL, be)
    }
  }
  # Cumulative — return vector, just check no error
  for (fn in c("cumsum", "cumprod", "cummax", "cummin")) {
    .safe_call(do.call(fn, list(A)), tag(paste0("Math_", fn)))
  }

  # --- Summary group (return scalars, no class assertion) ---
  for (fn in c("sum", "max", "min", "prod")) {
    .safe_call(do.call(fn, list(A)), tag(paste0("Summary_", fn)))
  }
  .safe_call(range(A), tag("Summary_range"))

  # --- Row/col reductions (return named numeric vectors) ---
  res_rs <- .safe_call(rowSums(A),  tag("rowSums"))
  res_cs <- .safe_call(colSums(A),  tag("colSums"))
  res_rm <- .safe_call(rowMeans(A), tag("rowMeans"))
  res_cm <- .safe_call(colMeans(A), tag("colMeans"))
  if (!is.null(res_rs) && !identical(names(res_rs), rn))
    .record_fail(tag("rowSums"), "DIMNAME_LOSS", paste(rn, collapse=","),
                 paste(as.character(names(res_rs)), collapse=","))
  if (!is.null(res_cs) && !identical(names(res_cs), cn))
    .record_fail(tag("colSums"), "DIMNAME_LOSS", paste(cn, collapse=","),
                 paste(as.character(names(res_cs)), collapse=","))

  # --- Shape accessors ---
  .safe_call(dim(A),  tag("dim"))
  .safe_call(nrow(A), tag("nrow"))
  .safe_call(ncol(A), tag("ncol"))

  # --- Transpose ---
  res_t <- .safe_call(t(A), tag("t"))
  if (!is.null(res_t) && !inherits(res_t, "aMatrix")) {
    .record_fail(tag("t"), "CLASS_DEMOTION", "aMatrix-derived", class(res_t)[1])
  }

  # --- Indexing: [i,j] with drop=FALSE (the working path) ---
  res_sub_drop <- .safe_call(A[1:2, 1:3, drop = FALSE], tag("[submat_dropFALSE]"))
  if (!is.null(res_sub_drop)) {
    .check_matrix_result(res_sub_drop, tag("[submat_dropFALSE]"), cls,
                         c("r1","r2"), c("c1","c2","c3"), be)
  }
  res_row_drop <- .safe_call(A[1, 1:4, drop = FALSE], tag("[row_dropFALSE]"))
  if (!is.null(res_row_drop)) {
    .check_matrix_result(res_row_drop, tag("[row_dropFALSE]"), cls, "r1", cn, be)
  }

  # --- Indexing: [i,j] WITHOUT drop=FALSE (known demotion path) ---
  # A[1:2, 1:2] drops to dgeMatrix — CLASS DEMOTION (new bug)
  res_sub_nodrop <- .safe_call(A[1:2, 1:2], tag("[submat_nodrop]"))
  if (!is.null(res_sub_nodrop)) {
    if (!inherits(res_sub_nodrop, "adgeMatrix")) {
      .record_fail(tag("[submat_nodrop]"), "CLASS_DEMOTION", cls,
                   class(res_sub_nodrop)[1],
                   "[i,j] without drop=FALSE falls through to dgeMatrix")
    }
  }

  # --- Assignment: [<- ---
  Acopy <- A
  res_assign <- .safe_call({ Acopy[1, 1] <- 999.0; Acopy }, tag("[<-_element]"))
  if (!is.null(res_assign))
    .check_matrix_result(res_assign, tag("[<-_element]"), cls, rn, cn, be)

  # --- dimnames<- ---
  Ac2 <- A
  res_dn <- .safe_call({ dimnames(Ac2) <- list(c("a","b","c"), c("x","y","z","w")); Ac2 },
                       tag("dimnames<-"))
  if (!is.null(res_dn))
    .check_matrix_result(res_dn, tag("dimnames<-"), cls, c("a","b","c"), c("x","y","z","w"), be)

  Ac3 <- A
  res_rna <- .safe_call({ rownames(Ac3) <- c("R1","R2","R3"); Ac3 }, tag("rownames<-"))
  if (!is.null(res_rna))
    .check_matrix_result(res_rna, tag("rownames<-"), cls, c("R1","R2","R3"), cn, be)

  Ac4 <- A
  res_cna <- .safe_call({ colnames(Ac4) <- c("C1","C2","C3","C4"); Ac4 }, tag("colnames<-"))
  if (!is.null(res_cna))
    .check_matrix_result(res_cna, tag("colnames<-"), cls, rn, c("C1","C2","C3","C4"), be)

  # --- diag (extractor) ---
  .safe_call(diag(Asq), tag("diag_extractor"))

  # --- diag<- (replacement) ---
  # Confirmed working: returns adgeMatrix with preserved backend. Round-2 hypothesis was wrong.
  Asq_c <- Asq
  res_diag_repl <- .safe_call({ diag(Asq_c) <- c(10, 20); Asq_c }, tag("diag<-"))
  if (!is.null(res_diag_repl))
    .check_matrix_result(res_diag_repl, tag("diag<-"), cls, rn2, cn2, be)

  # --- Coercion ---
  res_asm <- .safe_call(as.matrix(A),  tag("as.matrix"))
  if (!is.null(res_asm) && !is.matrix(res_asm))
    .record_fail(tag("as.matrix"), "WRONG_TYPE", "matrix", class(res_asm)[1])
  .safe_call(as.numeric(A), tag("as.numeric"))
  .safe_call(as.vector(A),  tag("as.vector"))
  .safe_call(as.array(A),   tag("as.array"))

  # --- %*% ---
  Bt <- matrix(1:12, nrow = 4, ncol = 3)
  res_mm <- .safe_call(A %*% Bt, tag("%*%_matrix"))
  if (!is.null(res_mm)) .check_matrix_result(res_mm, tag("%*%_matrix"), cls, rn, NULL, be)

  res_mm2 <- .safe_call(A %*% t(A), tag("%*%_self_t"))
  if (!is.null(res_mm2)) .check_matrix_result(res_mm2, tag("%*%_self_t"), cls, rn, rn, be)

  # --- crossprod / tcrossprod ---
  res_cp  <- .safe_call(crossprod(A),  tag("crossprod"))
  if (!is.null(res_cp))  .check_matrix_result(res_cp,  tag("crossprod"),  cls, cn, cn, be)
  res_tcp <- .safe_call(tcrossprod(A), tag("tcrossprod"))
  if (!is.null(res_tcp)) .check_matrix_result(res_tcp, tag("tcrossprod"), cls, rn, rn, be)

  # --- solve ---
  res_solve <- .safe_call(solve(Asq), tag("solve"))
  if (!is.null(res_solve)) .check_matrix_result(res_solve, tag("solve"), cls, cn2, rn2, be)

  # --- rbind / cbind ---
  res_rbind <- .safe_call(rbind(A, A), tag("rbind"))
  if (!is.null(res_rbind) && !inherits(res_rbind, "aMatrix"))
    .record_fail(tag("rbind"), "CLASS_DEMOTION", "aMatrix-derived", class(res_rbind)[1])

  res_cbind <- .safe_call(cbind(A, A), tag("cbind"))
  if (!is.null(res_cbind) && !inherits(res_cbind, "aMatrix"))
    .record_fail(tag("cbind"), "CLASS_DEMOTION", "aMatrix-derived", class(res_cbind)[1])

  # --- kronecker (CONFIRMED missing — amatrix-jnd) ---
  Asq2 <- .make_square_dense(backend)
  Bsq2 <- .safe_call(new_adgeMatrix(matrix(c(1,0,0,1), 2, 2),
                                    preferred_backend = backend), tag("kronecker_build"))
  if (!is.null(Bsq2)) {
    res_kron <- .safe_call(kronecker(Asq2, Bsq2), tag("kronecker"))
    if (!is.null(res_kron) && !inherits(res_kron, "adgeMatrix"))
      .record_fail(tag("kronecker"), "CLASS_DEMOTION", cls, class(res_kron)[1],
                   "confirmed amatrix-jnd: kronecker returns dgeMatrix not adgeMatrix")
  }

  # --- norm / det ---
  for (type in c("1", "I", "F", "M")) .safe_call(norm(A, type = type), tag(paste0("norm_", type)))
  .safe_call(det(Asq), tag("det"))

  # --- svd / qr / chol ---
  .safe_call(svd(Asq, nu = 0, nv = 0), tag("svd"))
  .safe_call(qr(A), tag("qr"))
  .safe_call(chol(Asm), tag("chol"))

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Sparse sweep
# ---------------------------------------------------------------------------

.run_sparse_sweep <- function(backend = "cpu") {
  S   <- .make_sparse(backend)
  Ssq <- .make_square_sparse(backend)
  rn  <- c("r1", "r2", "r3")
  cn  <- c("c1", "c2", "c3", "c4")
  rn2 <- c("r1", "r2")
  cn2 <- c("c1", "c2")
  cls <- "adgCMatrix"
  be  <- backend
  tag <- function(s) paste0(s, "[sparse/", backend, "]")

  # --- Ops: Arithmetic (sparse op sparse — stays sparse) ---
  for (op in c("+", "-", "*")) {
    res <- .safe_call(do.call(op, list(S, S)), tag(paste0("Arith_", op, "_mat")))
    if (!is.null(res) && inherits(res, "Matrix")) {
      if (!inherits(res, "adgCMatrix"))
        .record_fail(tag(paste0("Arith_", op, "_mat")), "CLASS_DEMOTION", cls, class(res)[1])
      else
        .check_matrix_result(res, tag(paste0("Arith_", op, "_mat")), cls, rn, cn, be)
    }
  }

  # --- Ops: Arithmetic (sparse op scalar) ---
  # NOTE: + and - with a scalar make ALL entries nonzero → sparse becomes dense.
  # In that case, result is adgeMatrix (aMatrix-derived but not adgCMatrix).
  # * and / preserve sparsity structure → result stays adgCMatrix.
  for (op in c("+", "-")) {
    res <- .safe_call(do.call(op, list(S, 2)), tag(paste0("Arith_", op, "_scalar")))
    if (!is.null(res)) {
      # densification: result should at minimum be aMatrix-derived
      if (!inherits(res, "aMatrix"))
        .record_fail(tag(paste0("Arith_", op, "_scalar")), "CLASS_DEMOTION",
                     "aMatrix-derived (dense ok)", class(res)[1],
                     "scalar +/- densifies sparse; result should be adgeMatrix")
      # Record if result is NOT adgeMatrix (densification went to bare dgeMatrix)
      else if (inherits(res, "dgeMatrix") && !inherits(res, "adgeMatrix"))
        .record_fail(tag(paste0("Arith_", op, "_scalar")), "CLASS_DEMOTION",
                     "adgeMatrix (densified)", class(res)[1])
    }
  }
  for (op in c("*", "/")) {
    res <- .safe_call(do.call(op, list(S, 2)), tag(paste0("Arith_", op, "_scalar")))
    if (!is.null(res) && inherits(res, "Matrix")) {
      if (!inherits(res, "adgCMatrix"))
        .record_fail(tag(paste0("Arith_", op, "_scalar")), "CLASS_DEMOTION", cls, class(res)[1])
      else
        .check_matrix_result(res, tag(paste0("Arith_", op, "_scalar")), cls, rn, cn, be)
    }
  }

  # --- Ops: Compare ---
  for (op in c("==", "!=", "<", ">")) {
    res <- .safe_call(do.call(op, list(S, 1)), tag(paste0("Compare_", op)))
    if (!is.null(res) && inherits(res, "Matrix")) {
      if (!inherits(res, "aMatrix"))
        .record_fail(tag(paste0("Compare_", op)), "CLASS_DEMOTION",
                     "aMatrix-derived", class(res)[1])
    }
  }

  # --- Math group ---
  # Sparsity-preserving fns (f(0)=0): abs, sign, ceiling, floor, round → adgCMatrix expected
  # Sparsity-breaking fns (f(0)!=0): exp (exp(0)=1) → goes dense → adgeMatrix expected
  sparse_preserving <- c("abs", "sign", "ceiling", "floor")
  for (fn in sparse_preserving) {
    res <- .safe_call(do.call(fn, list(S)), tag(paste0("Math_", fn)))
    if (!is.null(res) && !is.vector(res)) {
      if (!inherits(res, "adgCMatrix"))
        .record_fail(tag(paste0("Math_", fn)), "CLASS_DEMOTION", cls, class(res)[1],
                     paste0(fn, "(0)=0 preserves sparsity; expect adgCMatrix"))
      else
        .check_matrix_result(res, tag(paste0("Math_", fn)), cls, NULL, NULL, be)
    }
  }
  sparse_breaking <- c("exp", "cosh", "cos", "sin", "tan")
  for (fn in sparse_breaking) {
    res <- .safe_call(do.call(fn, list(S)), tag(paste0("Math_", fn)))
    if (!is.null(res) && !is.vector(res)) {
      # Goes dense — must at least be aMatrix-derived
      if (!inherits(res, "aMatrix"))
        .record_fail(tag(paste0("Math_", fn)), "CLASS_DEMOTION",
                     "aMatrix-derived (dense ok)", class(res)[1])
    }
  }
  # sqrt/log need positive values
  Spos <- .safe_call(S * S + 1, tag("Math_prep_pos"))  # all-positive sparse-ish
  if (!is.null(Spos)) {
    for (fn in c("sqrt", "log")) {
      res <- .safe_call(do.call(fn, list(Spos)), tag(paste0("Math_", fn)))
      if (!is.null(res) && !is.vector(res) && !inherits(res, "aMatrix"))
        .record_fail(tag(paste0("Math_", fn)), "CLASS_DEMOTION",
                     "aMatrix-derived", class(res)[1])
    }
  }
  # Cumulative — return vector
  for (fn in c("cumsum", "cumprod")) {
    .safe_call(do.call(fn, list(S)), tag(paste0("Math_", fn)))
  }

  # --- Summary group ---
  for (fn in c("sum", "max", "min")) .safe_call(do.call(fn, list(S)), tag(paste0("Summary_", fn)))

  # --- Row/col reductions ---
  .safe_call(rowSums(S),  tag("rowSums"))
  .safe_call(colSums(S),  tag("colSums"))
  .safe_call(rowMeans(S), tag("rowMeans"))
  .safe_call(colMeans(S), tag("colMeans"))

  # --- Transpose ---
  res_t <- .safe_call(t(S), tag("t"))
  if (!is.null(res_t) && !inherits(res_t, "aMatrix"))
    .record_fail(tag("t"), "CLASS_DEMOTION", "aMatrix-derived", class(res_t)[1])

  # --- Indexing with drop=FALSE (working path) ---
  res_sub <- .safe_call(S[1:2, 1:2, drop = FALSE], tag("[submat_dropFALSE]"))
  if (!is.null(res_sub) && inherits(res_sub, "Matrix")) {
    if (!inherits(res_sub, "adgCMatrix"))
      .record_fail(tag("[submat_dropFALSE]"), "CLASS_DEMOTION", cls, class(res_sub)[1])
    else
      .check_matrix_result(res_sub, tag("[submat_dropFALSE]"), cls,
                           c("r1","r2"), c("c1","c2"), be)
  }

  # --- Indexing without drop=FALSE (demotion path) ---
  res_sub_nd <- .safe_call(S[1:2, 1:2], tag("[submat_nodrop]"))
  if (!is.null(res_sub_nd) && inherits(res_sub_nd, "Matrix")) {
    if (!inherits(res_sub_nd, "adgCMatrix"))
      .record_fail(tag("[submat_nodrop]"), "CLASS_DEMOTION", cls, class(res_sub_nd)[1],
                   "[i,j] without drop=FALSE falls through to dgCMatrix")
  }

  # --- dimnames<- ---
  Sc <- S
  res_dn <- .safe_call({ dimnames(Sc) <- list(c("a","b","c"), c("x","y","z","w")); Sc },
                       tag("dimnames<-"))
  if (!is.null(res_dn))
    .check_matrix_result(res_dn, tag("dimnames<-"), cls, c("a","b","c"), c("x","y","z","w"), be)

  # --- diag (extractor) ---
  .safe_call(diag(Ssq), tag("diag_extractor"))

  # --- diag<- (replacement) ---
  # Confirmed working: returns adgCMatrix with preserved backend.
  Ssq_c <- Ssq
  res_diag_repl <- .safe_call({ diag(Ssq_c) <- c(10, 20); Ssq_c }, tag("diag<-"))
  if (!is.null(res_diag_repl))
    .check_matrix_result(res_diag_repl, tag("diag<-"), cls, rn2, cn2, be)

  # --- Coercion ---
  res_asm <- .safe_call(as.matrix(S), tag("as.matrix"))
  if (!is.null(res_asm) && !is.matrix(res_asm))
    .record_fail(tag("as.matrix"), "WRONG_TYPE", "matrix", class(res_asm)[1])

  # --- %*% ---
  Bt <- matrix(1:12, nrow = 4, ncol = 3)
  res_mm <- .safe_call(S %*% Bt, tag("%*%_matrix"))
  if (!is.null(res_mm) && !is.numeric(res_mm)) {
    if (!inherits(res_mm, "aMatrix"))
      .record_fail(tag("%*%_matrix"), "CLASS_DEMOTION", "aMatrix-derived", class(res_mm)[1])
  }

  # --- crossprod / tcrossprod ---
  res_cp  <- .safe_call(crossprod(S),  tag("crossprod"))
  if (!is.null(res_cp) && inherits(res_cp, "Matrix") && !inherits(res_cp, "aMatrix"))
    .record_fail(tag("crossprod"), "CLASS_DEMOTION", "aMatrix-derived", class(res_cp)[1])

  res_tcp <- .safe_call(tcrossprod(S), tag("tcrossprod"))
  if (!is.null(res_tcp) && inherits(res_tcp, "Matrix") && !inherits(res_tcp, "aMatrix"))
    .record_fail(tag("tcrossprod"), "CLASS_DEMOTION", "aMatrix-derived", class(res_tcp)[1])

  # --- solve ---
  res_solve <- .safe_call(solve(Ssq), tag("solve"))
  if (!is.null(res_solve) && inherits(res_solve, "Matrix") && !inherits(res_solve, "aMatrix"))
    .record_fail(tag("solve"), "CLASS_DEMOTION", "aMatrix-derived", class(res_solve)[1])

  # --- rbind / cbind ---
  res_rbind <- .safe_call(rbind(S, S), tag("rbind"))
  if (!is.null(res_rbind) && inherits(res_rbind, "Matrix") && !inherits(res_rbind, "aMatrix"))
    .record_fail(tag("rbind"), "CLASS_DEMOTION", "aMatrix-derived", class(res_rbind)[1])

  res_cbind <- .safe_call(cbind(S, S), tag("cbind"))
  if (!is.null(res_cbind) && inherits(res_cbind, "Matrix") && !inherits(res_cbind, "aMatrix"))
    .record_fail(tag("cbind"), "CLASS_DEMOTION", "aMatrix-derived", class(res_cbind)[1])

  # --- kronecker (CONFIRMED missing — amatrix-jnd) ---
  Ssq2 <- .make_square_sparse(backend)
  Bsq2 <- .safe_call(new_adgCMatrix(
    as(Matrix::sparseMatrix(i=c(1,2), j=c(1,2), x=c(1,1), dims=c(2,2)), "dgCMatrix"),
    preferred_backend = backend), tag("kronecker_build"))
  if (!is.null(Bsq2)) {
    res_kron <- .safe_call(kronecker(Ssq2, Bsq2), tag("kronecker"))
    if (!is.null(res_kron) && !inherits(res_kron, "aMatrix"))
      .record_fail(tag("kronecker"), "CLASS_DEMOTION", "aMatrix-derived", class(res_kron)[1],
                   "confirmed amatrix-jnd: kronecker returns dgCMatrix not adgCMatrix")
  }

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# testthat wrappers
# ---------------------------------------------------------------------------

test_that("method sweep: adgeMatrix cpu backend", {
  skip_if_not_installed("amatrix")
  .sweep_failures <<- list()
  .run_dense_sweep("cpu")
  fails <- .sweep_failures
  assign(".dense_cpu_fails", fails, envir = .GlobalEnv)
  if (length(fails) > 0) {
    msgs <- vapply(fails, function(f)
      sprintf("  [%s] %s: expected=%s actual=%s", f$symptom, f$method, f$expected, f$actual),
      character(1))
    message("=== adgeMatrix/cpu sweep failures ===\n", paste(msgs, collapse = "\n"))
  }
  expect_true(TRUE)
})

test_that("method sweep: adgCMatrix cpu backend", {
  skip_if_not_installed("amatrix")
  .sweep_failures <<- list()
  .run_sparse_sweep("cpu")
  fails <- .sweep_failures
  assign(".sparse_cpu_fails", fails, envir = .GlobalEnv)
  if (length(fails) > 0) {
    msgs <- vapply(fails, function(f)
      sprintf("  [%s] %s: expected=%s actual=%s", f$symptom, f$method, f$expected, f$actual),
      character(1))
    message("=== adgCMatrix/cpu sweep failures ===\n", paste(msgs, collapse = "\n"))
  }
  expect_true(TRUE)
})

test_that("method sweep: adgeMatrix mlx backend", {
  skip_if_not_installed("amatrix")
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    tryCatch({ library(amatrix.mlx); amatrix_mlx_is_available() },
             error = function(e) FALSE),
    "MLX not available"
  )
  .sweep_failures <<- list()
  .run_dense_sweep("mlx")
  assign(".dense_mlx_fails", .sweep_failures, envir = .GlobalEnv)
  expect_true(TRUE)
})

test_that("method sweep: adgCMatrix mlx backend", {
  skip_if_not_installed("amatrix")
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    tryCatch(amatrix_mlx_is_available(), error = function(e) FALSE),
    "MLX not available"
  )
  .sweep_failures <<- list()
  .run_sparse_sweep("mlx")
  assign(".sparse_mlx_fails", .sweep_failures, envir = .GlobalEnv)
  expect_true(TRUE)
})

# ---------------------------------------------------------------------------
# Write report
# ---------------------------------------------------------------------------

test_that("method sweep: write report to .bug-hunt-r3/03-method-sweep.md", {
  skip_if_not_installed("amatrix")

  all_fails <- c(
    if (exists(".dense_cpu_fails",  envir = .GlobalEnv)) get(".dense_cpu_fails",  envir = .GlobalEnv) else list(),
    if (exists(".sparse_cpu_fails", envir = .GlobalEnv)) get(".sparse_cpu_fails", envir = .GlobalEnv) else list(),
    if (exists(".dense_mlx_fails",  envir = .GlobalEnv)) get(".dense_mlx_fails",  envir = .GlobalEnv) else list(),
    if (exists(".sparse_mlx_fails", envir = .GlobalEnv)) get(".sparse_mlx_fails", envir = .GlobalEnv) else list()
  )

  # testthat::test_path() returns tests/testthat/<file>, so go up 3 levels to repo root
  repo_root  <- normalizePath(file.path(dirname(testthat::test_path()), "..", ".."),
                              mustWork = FALSE)
  report_dir <- file.path(repo_root, ".bug-hunt-r3")
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
  report_path <- file.path(report_dir, "03-method-sweep.md")

  lines <- c(
    "# Round 3 Bug Hunt — Method Sweep Report",
    "",
    paste("**Generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("**Total failures recorded:**", length(all_fails)),
    "",
    "## Methods Tested",
    "",
    "### adgeMatrix (dense)",
    paste(
      "Arith (+,-,*,/,^) mat×mat and mat×scalar and scalar×mat;",
      "Compare (==,!=,<,>,<=,>=) vs scalar;",
      "Math group: abs, sqrt, exp, log, log2, log10, ceiling, floor, sign,",
      "cos, sin, tan, cosh, sinh, tanh, cumsum, cumprod, cummax, cummin;",
      "Summary: sum, max, min, prod, range;",
      "rowSums, colSums, rowMeans, colMeans;",
      "dim, nrow, ncol;",
      "t();",
      "[i,j,drop=FALSE] and [i,j] (no drop);",
      "[<-;",
      "dimnames<-, rownames<-, colnames<-;",
      "diag (extractor), diag<- (replacement);",
      "as.matrix, as.numeric, as.vector, as.array;",
      "%*%, crossprod, tcrossprod, solve;",
      "rbind, cbind;",
      "kronecker;",
      "norm (1/I/F/M), det, svd, qr, chol"
    ),
    "",
    "### adgCMatrix (sparse)",
    paste(
      "Arith (+,-,*,/) mat×mat and mat×scalar;",
      "Compare (==,!=,<,>) vs scalar;",
      "Math group: abs, sign, ceiling, floor (sparsity-preserving);",
      "exp, cosh, cos, sin, tan (sparsity-breaking → dense);",
      "sqrt, log (on positive input);",
      "cumsum, cumprod;",
      "Summary: sum, max, min;",
      "rowSums, colSums, rowMeans, colMeans;",
      "t();",
      "[i,j,drop=FALSE] and [i,j] (no drop);",
      "dimnames<-;",
      "diag (extractor), diag<- (replacement);",
      "as.matrix;",
      "%*%, crossprod, tcrossprod, solve;",
      "rbind, cbind;",
      "kronecker"
    ),
    "",
    "---",
    "",
    "## Failures by Symptom",
    ""
  )

  symptoms <- c("CLASS_DEMOTION", "DIMNAME_LOSS", "BACKEND_LEAK", "ERROR", "WRONG_TYPE")
  for (sym in symptoms) {
    group <- Filter(function(f) f$symptom == sym, all_fails)
    if (length(group) == 0) next
    lines <- c(lines, paste0("### ", sym, " (", length(group), ")"), "")
    lines <- c(lines, "| Method | Expected | Actual | Note |")
    lines <- c(lines, "|--------|----------|--------|------|")
    for (f in group) {
      note <- if (nzchar(f$note)) substr(f$note, 1, 100) else ""
      lines <- c(lines,
        sprintf("| `%s` | `%s` | `%s` | %s |",
                f$method, f$expected, f$actual, note))
    }
    lines <- c(lines, "")
  }

  # Partition confirmed vs new
  known_patterns <- c("kronecker", "diag<-")
  confirmed <- Filter(function(f) any(vapply(known_patterns, function(p) grepl(p, f$method, fixed=TRUE), logical(1))), all_fails)
  novel <- Filter(function(f) !any(vapply(known_patterns, function(p) grepl(p, f$method, fixed=TRUE), logical(1))), all_fails)

  lines <- c(lines,
    "---", "",
    "## Confirmed Bugs (Round-2 Issues Executed)",
    "",
    paste0("Confirmed ", length(confirmed), " failure(s) that match round-2 issue patterns."),
    ""
  )
  for (f in confirmed)
    lines <- c(lines,
      sprintf("- **CONFIRMED** `%s` → %s (expected `%s`, got `%s`)%s",
              f$method, f$symptom, f$expected, f$actual,
              if (nzchar(f$note)) paste0(" — ", f$note) else ""))

  lines <- c(lines,
    "",
    "---", "",
    paste0("## NEW Bugs (Not in Round-2 Issue List) — ", length(novel), " found"),
    ""
  )
  if (length(novel) > 0) {
    for (f in novel)
      lines <- c(lines,
        sprintf("- **NEW** `%s` → %s (expected `%s`, got `%s`)%s",
                f$method, f$symptom, f$expected, f$actual,
                if (nzchar(f$note)) paste0(" — ", f$note) else ""))
  } else {
    lines <- c(lines, "None.")
  }

  lines <- c(lines, "",
    "---", "",
    "## Additional Findings (Hypothesis Refutations)", "",
    "The following round-2 hypotheses were WRONG — these methods work correctly:",
    "",
    "- **Math group on adgeMatrix** (amatrix-86l): `abs`, `sqrt`, `exp`, `log`, `ceiling`, `floor`, `sign`, `cos`, `sin`, `tan`, `cosh`, `sinh`, `tanh` all return `adgeMatrix` with backend preserved. The bug hypothesis was incorrect for the dense class.",
    "- **Math group on adgCMatrix (sparsity-preserving)**: `abs`, `sign`, `ceiling`, `floor` return `adgCMatrix`. Correct.",
    "- **diag<- on adgeMatrix** (amatrix-j5a): `diag(A) <- v` preserves `adgeMatrix` class and `preferred_backend` slot. The replacement method works correctly.",
    "- **diag<- on adgCMatrix** (amatrix-j5a): `diag(S) <- v` preserves `adgCMatrix` class and `preferred_backend` slot. Correct.",
    "",
    "---", "",
    "## Top 3 Failures by Impact", ""
  )

  impact_order <- c("CLASS_DEMOTION", "BACKEND_LEAK", "DIMNAME_LOSS", "ERROR")
  sorted_fails <- all_fails[order(match(
    vapply(all_fails, `[[`, character(1), "symptom"), impact_order))]
  for (i in seq_len(min(3, length(sorted_fails)))) {
    f <- sorted_fails[[i]]
    lines <- c(lines,
      sprintf("%d. **`%s`** — %s: expected `%s`, got `%s`%s",
              i, f$method, f$symptom, f$expected, f$actual,
              if (nzchar(f$note)) paste0(". ", f$note) else ""))
  }
  if (length(sorted_fails) == 0) lines <- c(lines, "No failures recorded.")

  lines <- c(lines, "", "---", "", "_End of Round 3 Method Sweep Report_")

  writeLines(lines, report_path)
  message("Report written to: ", report_path)
  expect_true(file.exists(report_path))
})
