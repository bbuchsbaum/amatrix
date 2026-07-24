# demopkg/ — not part of the amatrix build

`amatrix.demo` is an internal demonstration package showing how a package
author adopts `amatrix` to accelerate an existing statistical algorithm.
It is excluded from the `amatrix` build (`.Rbuildignore`: `^demopkg$`) and is
**not** on the R-universe registry; it exists so that the "package-author
migration" story (`planning_docs/package-author-migration.md`) has a small,
runnable, checkable proof.

What it shows:

- one classic algorithm (logistic regression via IRLS, the workhorse behind
  `glm(family = binomial)`) written **once**, with no backend-specific code;
- the hot kernels expressed as `amatrix` helpers
  (`crossprod_weighted()`, `xty_weighted()`, `solve()`);
- CPU-correct behavior on plain base matrices, with acceleration obtained by
  the *caller* wrapping the design matrix: `adgeMatrix(X, mode = "fast")`;
- tests that verify the result against `stats::glm.fit()` on every platform,
  plus a guarded cross-backend agreement test.

Try it from the repo root:

```r
# install.packages("pkgload") if needed
pkgload::load_all(".")                    # core amatrix
pkgload::load_all("demopkg/amatrix.demo") # the demo
logit_demo_benchmark()
```
