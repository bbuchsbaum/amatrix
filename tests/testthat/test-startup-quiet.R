# Startup-message gating (amatrix-qojr).
#
# .onAttach prints a one-line GPU-backends note. These tests exercise the
# suppression/filtering knobs for quiet CPU-only startup:
#   * options(amatrix.quiet_startup = TRUE) or AMATRIX_QUIET=1/"true"
#   * options(amatrix.optional_backends = FALSE)
#   * per-backend options(amatrix.disable_<backend> = TRUE)
#
# Contract (docs are written to this): a backend disabled via
# amatrix.disable_* is OMITTED from the note entirely -- it must not appear
# as "installed but disabled (...)". The "installed but disabled" phrasing is
# reserved for the separate probe-policy gate (e.g. an opencl backend whose
# GPU probe is not enabled). If every installed backend is disabled, no note
# is printed at all.
#
# The note only mentions *installed* backend packages, and no GPU backend
# package is guaranteed present on CI. To make the gating deterministic we run
# in a fresh child process and mock the backend spec table so it points at
# base packages ("stats"/"utils"/"tools") that are always installed. The
# "opencl" entry keeps its real name so the opencl probe-policy branch fires
# (probe not enabled -> "installed but disabled"), letting us assert the
# omission-vs-"installed but disabled" distinction.
#
# Source-tree lookup + child-process idiom mirror the other fresh-attach
# regression tests; the shared helper lives in helper-repo-source-tree.R.

.startup_quiet_repo_dir <- function() {
  .amatrix_source_tree_dir()
}

test_that("startup note is gated by quiet options, env var, and per-backend disables [amatrix-qojr]", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  repo_dir <- .startup_quiet_repo_dir()
  skip_if(is.null(repo_dir), "source tree not reachable (installed-pkg context)")

  results <- callr::r(
    function(repo_dir) {
      pkgload::load_all(repo_dir, quiet = TRUE, attach = FALSE)

      # Inject deterministic "installed" backends backed by base packages so
      # the note is non-empty regardless of which GPU backends are present.
      # "opencl" keeps its real name so the opencl probe-policy branch fires
      # (probe not enabled -> "installed but disabled").
      fake_specs <- function() {
        list(
          mlx = list(
            package = "stats",
            disable_option = "amatrix.disable_mlx",
            auto_probe = TRUE
          ),
          metal = list(
            package = "utils",
            disable_option = "amatrix.disable_metal"
          ),
          opencl = list(
            package = "tools",
            disable_option = "amatrix.disable_opencl",
            probe_env = "AMATRIX_OPENCL_PROBE_GPU"
          )
        )
      }
      fake_order <- function() c("mlx", "metal", "opencl")
      assignInNamespace(".amatrix_optional_backend_specs", fake_specs, ns = "amatrix")
      assignInNamespace(".amatrix_auto_fast_backend_order", fake_order, ns = "amatrix")

      ns <- asNamespace("amatrix")

      capture_note <- function() {
        msgs <- character()
        withCallingHandlers(
          ns$.onAttach(NULL, "amatrix"),
          packageStartupMessage = function(m) {
            msgs[[length(msgs) + 1L]] <<- conditionMessage(m)
            invokeRestart("muffleMessage")
          }
        )
        paste(msgs, collapse = "")
      }

      run_scenario <- function(opts = list(), envs = list()) {
        old_opts <- options(opts)
        on.exit(options(old_opts), add = TRUE)
        # Keep the opencl probe deterministically "not enabled" unless a
        # scenario overrides it, so opencl lands on the probe-policy branch.
        envs <- modifyList(list(AMATRIX_OPENCL_PROBE_GPU = ""), envs)
        do.call(Sys.setenv, envs)
        on.exit(Sys.unsetenv(names(envs)), add = TRUE)
        capture_note()
      }

      list(
        baseline         = run_scenario(),
        opt_quiet        = run_scenario(opts = list(amatrix.quiet_startup = TRUE)),
        env_quiet_1      = run_scenario(envs = list(AMATRIX_QUIET = "1")),
        env_quiet_true   = run_scenario(envs = list(AMATRIX_QUIET = "true")),
        env_quiet_off    = run_scenario(envs = list(AMATRIX_QUIET = "0")),
        opt_backends_off = run_scenario(opts = list(amatrix.optional_backends = FALSE)),
        disable_metal    = run_scenario(opts = list(amatrix.disable_metal = TRUE)),
        disable_opencl   = run_scenario(opts = list(amatrix.disable_opencl = TRUE)),
        disable_all      = run_scenario(opts = list(
          amatrix.disable_mlx = TRUE,
          amatrix.disable_metal = TRUE,
          amatrix.disable_opencl = TRUE
        ))
      )
    },
    args = list(repo_dir = repo_dir)
  )

  # By default the note is printed and names every installed (fake) backend.
  # opencl shows the probe-policy "installed but disabled" phrasing.
  expect_match(results$baseline, "amatrix GPU backends", fixed = TRUE)
  expect_match(results$baseline, "mlx")
  expect_match(results$baseline, "metal")
  expect_match(results$baseline, "opencl installed but disabled")

  # Quiet knobs fully suppress the note.
  expect_identical(results$opt_quiet, "")
  expect_identical(results$env_quiet_1, "")
  expect_identical(results$env_quiet_true, "")

  # An unrelated / falsey env value does not suppress it.
  expect_match(results$env_quiet_off, "amatrix GPU backends", fixed = TRUE)

  # optional_backends = FALSE suppresses the note entirely.
  expect_identical(results$opt_backends_off, "")

  # Per-backend disable OMITS that backend from the note (not "installed but
  # disabled"), leaving the others.
  expect_match(results$disable_metal, "mlx")
  expect_match(results$disable_metal, "opencl")
  expect_false(grepl("metal", results$disable_metal))

  # Disabling a probe-gated backend via its option omits it entirely rather
  # than rendering the "installed but disabled" probe-policy phrasing.
  expect_false(grepl("opencl", results$disable_opencl))
  expect_match(results$disable_opencl, "mlx")
  expect_match(results$disable_opencl, "metal")

  # Disabling every installed backend leaves nothing to announce.
  expect_identical(results$disable_all, "")
})
