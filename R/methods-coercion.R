setMethod("as.matrix", "adgeMatrix",      function(x, ...) as.matrix(amatrix_materialize_host(x), ...))
setMethod("as.matrix", "adgCMatrix",      function(x, ...) as.matrix(amatrix_materialize_host(x), ...))
setMethod("as.matrix", "aTransposeView",  function(x, ...) t(as.matrix(amatrix_materialize_dense(x@source), ...)))

setMethod("as.numeric", "adgeMatrix", function(x, ...) as.numeric(as.matrix(amatrix_materialize_host(x)), ...))
setMethod("as.vector",  "adgeMatrix", function(x, mode = "any") as.vector(as.matrix(amatrix_materialize_host(x)), mode))

setMethod("as.array", "adgeMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))
setMethod("as.array", "adgCMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))

setReplaceMethod("dimnames", "adgeMatrix", function(x, value) {
  am_set_dimnames(x, value)
})

setReplaceMethod("dimnames", "adgCMatrix", function(x, value) {
  am_set_dimnames(x, value)
})
