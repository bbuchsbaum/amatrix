# Claude Notes

## Quality Strategy

See `docs/quality-tracking.md` for the full accuracy and performance methodology.
Summary:
- **Accuracy:** `devtools::test()` — cross-backend conformance at `1e-4` (GPU) / `1e-10` (CPU)
- **Performance:** `Rscript tools/benchmark-regression.R` — compare to `tools/baseline.csv`
- Every new exported op must be added to the coverage table in that doc.

## Issue Tracking

Use `br` (`beads_rust`) for local issue tracking in this repo.

- The project-local workspace is `.beads/`.
- The issue prefix for this repo is `amatrix`, so new issues must be `amatrix-*`.
- Do not use `bd` in this repo.
- `br` is non-invasive: after changing issue state, run `br sync --flush-only`, then `git add .beads/` so the JSONL and SQLite state stay in version control.
- If the workspace looks inconsistent, run `br doctor` before making further changes.

> **N.B.** `br close`/`br update` can fail with `FOREIGN KEY constraint failed` on certain issues — this is a bug in the `fsqlite` custom runtime, not a real constraint violation. Workaround: close directly via `sqlite3`, then `br sync --flush-only --force && git add .beads/`:
> ```bash
> sqlite3 .beads/beads.db "PRAGMA foreign_keys = ON; BEGIN; \
>   UPDATE issues SET status='closed', closed_at=CURRENT_TIMESTAMP, \
>     close_reason='done', updated_at=CURRENT_TIMESTAMP, closed_by_session='manual' \
>     WHERE id='<issue-id>'; \
>   INSERT INTO events (issue_id, event_type, actor, old_value, new_value) \
>     VALUES ('<issue-id>', 'status_changed', 'bbuchsbaum', 'open', 'closed'); \
>   COMMIT;"
> br sync --flush-only --force && git add .beads/
> ```
