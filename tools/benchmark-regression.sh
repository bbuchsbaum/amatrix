#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

cd "${repo_root}"
export AMATRIX_BENCHMARK_REPO_ROOT="${repo_root}"

exec Rscript -e 'setwd(Sys.getenv("AMATRIX_BENCHMARK_REPO_ROOT")); source(file.path("tools", "benchmark-regression.R"), local = globalenv()); benchmark_regression_main(commandArgs(trailingOnly = TRUE))' --args "$@"
