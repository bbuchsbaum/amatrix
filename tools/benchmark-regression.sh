#!/usr/bin/env bash
if [[ -z "${ZSH_VERSION-}" ]]; then
  exec /bin/zsh "$0" "$@"
fi

set -eu
set -o pipefail

script_dir="$(cd "$(dirname "${(%):-%N}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

cd "${repo_root}"

benchmark_arrayfire_env="${AMATRIX_BENCHMARK_ARRAYFIRE-}"
benchmark_arrayfire_unsafe_env="${AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE-}"
unset AMATRIX_BENCHMARK_ARRAYFIRE
unset AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE

output_dir=""
args=()
for arg in "$@"; do
  if [[ "$arg" == --output-dir=* ]]; then
    output_dir="${arg#--output-dir=}"
  fi
  args+=("$arg")
done

if [[ -z "${output_dir}" ]]; then
  stamp="$(date '+%Y%m%d-%H%M%S')"
  output_dir="tools/benchmark-results/${stamp}"
  args+=("--output-dir=${output_dir}")
fi

r_args_expr="c("
for arg in "${args[@]}"; do
  escaped="${arg//\\/\\\\}"
  escaped="${escaped//\"/\\\"}"
  if [[ "${r_args_expr}" != "c(" ]]; then
    r_args_expr+=", "
  fi
  r_args_expr+="\"${escaped}\""
done
r_args_expr+=")"

if Rscript -e "setwd(\"${repo_root}\"); if (nzchar(\"${benchmark_arrayfire_env}\")) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE = \"${benchmark_arrayfire_env}\"); if (nzchar(\"${benchmark_arrayfire_unsafe_env}\")) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = \"${benchmark_arrayfire_unsafe_env}\"); options(amatrix.benchmark_regression.autorun = FALSE); source(file.path(\"tools\", \"benchmark-regression.R\"), local = globalenv()); parsed <- parse_args(${r_args_expr}); initialize_regression_benchmark_context(); if (isTRUE(parsed\$worker)) run_worker(parsed) else run_master(parsed)" ; then
  exit_code=0
else
  exit_code=$?
fi

if [[ -f "${output_dir}/benchmark-report.qmd" ]] && command -v quarto >/dev/null 2>&1; then
  (
    cd "${output_dir}"
    quarto render benchmark-report.qmd --to html --output benchmark-report.html
    quarto render benchmark-report.qmd --to pdf --output benchmark-report.pdf
  ) || echo "warning: quarto render failed for ${output_dir}/benchmark-report.qmd" >&2
fi

exit "${exit_code}"
