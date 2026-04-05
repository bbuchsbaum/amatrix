#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
TMPDIR="${TMPDIR:-/tmp}/af_r_repro"
mkdir -p "$TMPDIR"

echo "Building standalone C repro"
clang \
  -I/opt/homebrew/opt/arrayfire/include \
  "$ROOT/arrayfire_qr_repro.c" \
  -L/opt/homebrew/opt/arrayfire/lib \
  -Wl,-rpath,/opt/homebrew/opt/arrayfire/lib \
  -laf \
  -o "$TMPDIR/arrayfire_qr_repro"

echo "Building minimal R extension repro"
PKG_CPPFLAGS='-I/opt/homebrew/opt/arrayfire/include' \
PKG_LIBS='-L/opt/homebrew/opt/arrayfire/lib -Wl,-rpath,/opt/homebrew/opt/arrayfire/lib -laf' \
R CMD SHLIB "$ROOT/af_r_repro.c" -o "$TMPDIR/af_r_repro.so"

echo
echo "Standalone C repro"
for backend in default cpu opencl; do
  echo "=== standalone backend=$backend ==="
  "$TMPDIR/arrayfire_qr_repro" "$backend" 96
  echo "exit_code=$?"
done

echo
echo "Minimal R extension repro"
for backend in 4 1; do
  echo "=== R extension backend=$backend ==="
  set +e
  B="$backend" SO="$TMPDIR/af_r_repro.so" R -q -e '
    dyn.load(Sys.getenv("SO"))
    .Call("af_r_set_backend", as.integer(as.integer(Sys.getenv("B"))))
    print(.Call("af_r_diag"))
    set.seed(1)
    x <- matrix(rnorm(96 * 96), 96, 96)
    out <- .Call("af_r_qr", x)
    print(dim(out$q))
    print(max(abs(out$q %*% out$r - x)))
  '
  status=$?
  set -e
  echo "exit_code=$status"
done
