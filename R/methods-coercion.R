setMethod("as.matrix", "adgeMatrix", function(x, ...) as.matrix(amatrix_materialize_host(x), ...))
setMethod("as.matrix", "adgCMatrix", function(x, ...) as.matrix(amatrix_materialize_host(x), ...))

setMethod("as.array", "adgeMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))
setMethod("as.array", "adgCMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))

setReplaceMethod("dimnames", "adgeMatrix", function(x, value) {
  am_set_dimnames(x, value)
})

setReplaceMethod("dimnames", "adgCMatrix", function(x, value) {
  am_set_dimnames(x, value)
})
