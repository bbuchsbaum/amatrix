test_that("mlx sparse resident bridges match Matrix on small products", {
  set.seed(1)
  x_host <- Matrix::rsparsematrix(256, 128, density = 0.05)
  y_host <- matrix(rnorm(128 * 16), nrow = 128, ncol = 16)
  sp_key <- paste0("sp-", sample.int(1e6, 1))
  y_key <- paste0("y-", sample.int(1e6, 1))
  out_key <- paste0("out-", sample.int(1e6, 1))

  invisible(.Call(
    "amatrix_mlx_sparse_store_bridge",
    sp_key,
    as.double(x_host@x),
    as.integer(x_host@p),
    as.integer(x_host@i),
    as.integer(x_host@Dim),
    PACKAGE = "amatrix.mlx"
  ))
  invisible(.Call("amatrix_mlx_resident_store_bridge", y_key, y_host, PACKAGE = "amatrix.mlx"))

  ref <- as.matrix(x_host %*% y_host)
  direct <- .Call("amatrix_mlx_spmm_resident_bridge", sp_key, y_host, FALSE, PACKAGE = "amatrix.mlx")
  resident_rhs <- .Call("amatrix_mlx_spmm_resident_key_bridge", sp_key, y_key, out_key, FALSE, PACKAGE = "amatrix.mlx")

  expect_equal(direct, ref, tolerance = 1e-5)
  expect_equal(resident_rhs, ref, tolerance = 1e-5)
})
