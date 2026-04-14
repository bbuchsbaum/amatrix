# Agent Notes

## Quality Strategy

See `planning_docs/quality-tracking.md` for the full accuracy and performance methodology.
Summary:
- **Accuracy:** `devtools::test()` — cross-backend conformance at `1e-4` (GPU) / `1e-10` (CPU)
- **Performance:** `Rscript tools/benchmark-regression.R` — compare to `tools/baseline.csv`
- Every new exported op must be added to the coverage table in that doc.

## Bug Capture Protocol (Track 6)

Every discovered bug gets a **minimal repro first**, then the fix. This is
non-negotiable. The repro becomes a permanent regression test so the bug
cannot silently return.

**Workflow:**

1. **Write the repro before the fix.** Add a new test file under
   `tests/testthat/` named `test-regression-<short-slug>.R`. The file must
   fail on current `main` and pass after the fix.

2. **Record reproduction metadata in the test file header:**
   - Seed(s) used to construct the input
   - Dimensions and shape of the failing input
   - Backend, precision mode, and dispatch path (cold / resident)
   - R version and platform (from `sessionInfo()`)
   - Link to the beads issue or external bug report

3. **Never delete regression tests.** A closed bug can reopen. The test is
   the tombstone that says "this is fixed and we check every PR that it
   stays fixed."

4. **Empty `tests/testthat/_problems/` at release.** Extracted failing cases
   are temporary; every file there must either be fixed and reintegrated as
   `test-regression-*.R` or explicitly deleted with a commit message
   explaining the decision.

See `planning_docs/quality-tracking.md` §7 rule 5 (orphan repros are
stop-ship) and the four-tests-per-op rule in §4 Coverage Matrix.

## Issue Tracking

Use `bd` (beads) for local issue tracking in this repo.

- The project-local workspace is `.beads/`.
- The issue prefix for this repo is `amatrix`, so new issues must be `amatrix-*`.
- After changing issue state, run `bd sync` to commit beads changes to git.
- If the workspace looks inconsistent, run `bd doctor` before making further changes.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
