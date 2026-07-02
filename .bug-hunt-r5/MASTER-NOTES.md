# To the next bug hunter — read this first

You are stepping into round 6. Rounds 1–5 ground through the easy stuff: op-by-op
fuzzers, double-drop propagation, method-sweep S4 audits, residency lifecycle,
indexing invariants, mixed-operand arithmetic, serialization, bind, apply-idioms.
What survives is harder to find and more dangerous to close. I have been doing
this for two rounds as the orchestrator. Here is what I have learned. Treat it
as house rules, not suggestions.

---

## Rule 0 — The single lesson that matters

**Run the code.** Not the source. Not a "similar" matrix. Not a "probe that should
exhibit the symptom." Run THE probe the bug description gives you, with THE inputs
it names, in a fresh `Rscript` process. Capture the exact output. Paste it in your
report.

Two rounds in a row now, refuter agents have been wrong by source-reading:

- **Round 4, amatrix-p24.** Refuter argued "CPU uses `sweep()` not `rep()`, therefore
  the rep-bug cannot exist." Orchestrator ran
  `X=rbind(c(0,0),c(1,1),c(5,5)); C=rbind(c(0,0),c(5,5)); pairwise_sqdist_argmin(X,Ct,...)`
  and got `c(1,1,1)` instead of `c(1,1,2)`. Bug was live.
- **Round 5, amatrix-75h.** Refuter used their own probe (small random matrix) and
  found the diagonal "exactly 1.0, max dev = 0." Orchestrator ran the exact r4
  probe `set.seed(1); Y<-matrix(rnorm(20),4,5); kernel_matrix(Y,Y,'rbf',1)` and
  got `max dev 4.77e-07`. The refuter's probe happened to dodge the float32 boundary.

The pattern: refuters construct a probe that **happens** not to fire the bug, then
declare refuted. Defense: **use the inputs the bug description names, verbatim**.
If the description has no probe, the description is unfalsifiable and you must
synthesize one aggressively — not cautiously.

If you cannot run the code (GPU not present, network unavailable, state not
reproducible), you mark **INCONCLUSIVE**. You do not mark REFUTED. Source-reading
alone is never refutation. I will spot-check you, and if I catch a false
refutation, the issue reopens with your report pinned to it as a counter-example.

---

## Rule 1 — Stub to disk as your first action

Before any exploration, before any drift check, before any `bd list`, create your
report file with section headers at `.bug-hunt-rN/NN-scenario.md`:

```
# Hunter NN — <scenario>
## (a) Drift check
## (b) Scenario
## (c) Findings
## (d) Proposed bd create
## (e) Limitations
```

Then update incrementally as you work. This is not cosmetic. In round 4, three
scientist-agent hunters completed 40+ tool uses of investigation and then died
during the final-report-write phase, producing 200-byte stubs or nothing. Round 5
used `executor` exclusively and every report landed.

If your process dies, what is on disk **is** your report. Do not hold findings
in agent memory. Write them as you discover them.

---

## Rule 2 — Use `executor`, not `scientist`

Scientist is better for pure analysis. Executor is better for any task that must
produce a file. Every hunting round must produce files. Use executor.

---

## Rule 3 — Filesystem over summary

After any agent completes, the task-notification's `summary` field is the agent's
last message. It can read "Report written — here are the findings" while the file
on disk is a stub. Every time an agent completes, the orchestrator must:

1. `ls -la .bug-hunt-rN/`
2. `wc -l` the expected report
3. `bd show <id>` for every claimed filing
4. Spot-verify refutations with a runtime probe before any `bd close`

Trust the filesystem. Do not trust the summary.

---

## Rule 4 — Pick a scenario, not an operation

Rounds 1–3 covered operations exhaustively: every reduction, every product, every
factorization, every sub-assign, every group generic was fuzzed. The remaining
yield is **not** in more ops. It is in realistic user scenarios — 3 to 10 lines
of script that someone who just `install.packages()`d amatrix might actually
type.

The five highest-yield scenarios of round 5, in order of impact, came from:

- **Serialization round-trip across processes** — `saveRDS(X, f); q(); Rscript -e 'readRDS(f); as.matrix(...)'`
  Result: silent NaN (amatrix-90k, P1). Nobody had ever saved an object and
  reloaded it. Every prior test held objects in memory.

- **Mixed-operand arithmetic** — `X > 5`, `X %*% c(1,2,3,4)`, `X + 1L`, `1 + X`.
  Result: silent class demotion through the `.amatrix_rewrap_value` gap
  (amatrix-ol8, amatrix-qic).

- **cbind/rbind with plain types** — `cbind(X, matrix(0,3,2))`, `cbind(X, 1)`,
  `cbind(X, c(1,2,3))`. Result: silent demotion to `dgeMatrix` (amatrix-0qt).

- **Logical sub-assignment** — `X[is.na(X)] <- 0`. Result: "subscript out of
  bounds" (amatrix-e97). This is the most basic user idiom in R.

- **Minimal imports** — `Rscript -e 'library(amatrix); ...'` with no
  `library(Matrix)`. Round 4's crown jewel (amatrix-1ha). Still the single
  biggest root cause on the board.

Rules of thumb for scenario selection:

- If your scenario starts with "probe every variant of op X," stop. That ground
  is covered.
- If your scenario starts with "a user who just install.packages'd amatrix writes
  3 lines and...", you are on the right track.
- If the first 5 lines of your script are `library(Matrix); library(amatrix);
  X <- adgeMatrix(...)` — at least one of your probes must delete the
  `library(Matrix)` line. Half the bugs this round hid behind that assumption.
- Think about **composition**: save, load, cbind, apply, pass to `stats::lm`, copy,
  parallelize, interrupt. Amatrix is tested alone. Users never use it alone.

---

## Rule 5 — Scenarios that have NOT been touched as of end of round 5

These are the high-yield hunting grounds. Pick one. Do not overlap.

1. **Parallel workers.** `parallel::mclapply(1:4, function(i) rowSums(X))` —
   forked workers holding GPU resident keys. `future::plan(multicore)` across
   `adgeMatrix`. `callr::r(function() ...)` with an `adgeMatrix` in the closure.
   Any of these can corrupt the resident registry or produce silent wrong answers.
   Nobody has probed this.

2. **User interrupt during an allocation.** A user hits Esc in RStudio while
   a `crossprod(X)` is mid-flight on GPU. The resident key was allocated but
   the on.exit cleanup never fires because R jumped out of the call frame.
   Synthesize this with `withCallingHandlers` + `invokeRestart("abort")` or
   `R.utils::withTimeout(..., onTimeout="error")`. This probes the leak-after-
   interrupt surface that on.exit does not cover.

3. **Large dimensions near `.Machine$integer.max`.** `n*m > 2^31-1`. Any C bridge
   that uses `int` instead of `R_xlen_t` for element counts will silently wrap.
   Probe with `X <- adgeMatrix(matrix(0, 50000, 50000))` (2.5e9 elements, ~20 GB
   dense — use the smallest matrix that crosses the boundary; a 46341x46341
   double matrix is exactly at the edge).

4. **REPL interaction surface.** `print(X)`, `str(X)`, `summary(X)`, `format(X)`,
   `head(X)`, `tail(X)`, `View(X)`, `as.character(X)`, `dput(X)`. Every one of
   these is something a user types at the prompt in their first ten minutes.
   Any demotion or crash here is discoverability poison. Nobody has probed this.

5. **Passing through the formula interface.** `df <- data.frame(y=y);
   df$X <- X; lm(y ~ X, data=df)` — does this work, silently coerce, or crash?
   `model.matrix(~ X)`. `stats::lm.fit(as.matrix(X), y)`. These are the actual
   user paths to linear models.

6. **Copy-on-modify semantics.** R promises COW: `Y <- X; Y[1,1] <- 0`
   must not mutate `X`. Does amatrix honor this for resident-backed objects?
   There is a real risk that a shared GPU key makes `Y` alias `X` on the device,
   breaking R's fundamental contract. Probe:
   ```r
   X <- adgeMatrix(matrix(1:12, 3, 4)); Y <- X
   Y[1, 1] <- 999
   X[1, 1]  # MUST be 1, MUST NOT be 999
   ```

7. **Finalizer GC under memory pressure.** Force `gc()` repeatedly with resident
   keys outstanding. Do finalizers fire in the right order? Does the backend
   registry get walked while a finalizer is mid-drop?

8. **NSE / symbolic interfaces.** `deparse(substitute(X))`, `quote(X %*% Y)`,
   `bquote(.(X) %*% .(Y))`, `eval(parse(text="X %*% Y"))`. Users who build
   formulas or symbolic expressions hit these paths. Any S4-dispatch gotcha
   tends to leak here.

9. **Hybrid workflows that escape amatrix.** `as.matrix(X) %*% v` — the instant
   the user writes `as.matrix`, they are back in base R. Does the residency
   get released? Does the deferred flag stay clean? What about
   `Matrix::Matrix(X)` — round-trip through the parent class?

10. **Dimnames preservation.** `dimnames(X) <- list(r, c); X[1,]`,
    `rowSums(X)` — does the name vector survive every op? S4 dispatch routinely
    drops dimnames silently.

---

## Rule 6 — Check the open bug list BEFORE filing

Round 5 H4 filed amatrix-juq as a P1 "rowSums/colSums/rowMeans break via base
primitive bypass." Real symptom, correct analysis, **direct duplicate of
amatrix-1ha** (Matrix in Imports not Depends), which has been on the P1 board
since round 4. The hunter didn't look.

Before you file anything:

```bash
bd list --status=open | grep -i <keyword>
bd search <op-name>
```

If your finding has the same root cause or symptom as an existing open bug, add
a note to that bug with your probe. Do not file a new one. If your finding is a
**new symptom** of an existing root cause, link them — file the new symptom and
mention the root-cause bug in the description.

Duplicates waste refuter budget and fragment the discussion.

---

## Rule 7 — The `library(Matrix)` discipline

Round 4's biggest finding was `amatrix-1ha`: DESCRIPTION has `Matrix` in `Imports`,
not `Depends`. Under bare `library(amatrix)`, the symbols `t`, `chol`, `rowSums`,
`colSums`, `diag`, `solve`, `mean` resolve to base R, not Matrix's S4 generics.
Result: silent class demotion or crash, depending on the op.

Every scenario you probe must include at least one variant with **only**
`library(amatrix)` loaded — no `library(Matrix)`. This is the actual state a
user is in 10 minutes after install. If you skip it, you are testing a fixture,
not the package.

Note: as of round 5, `amatrix-1ha` remains OPEN. Many of your findings under
bare-imports may be downstream of it. If you find a "new" bug that disappears
when you add `library(Matrix)`, you have probably found another face of 1ha,
not a new bug. Add a note to 1ha; don't file.

---

## Rule 8 — Drift check, then proceed anyway

The installed `amatrix` namespace is frequently stale relative to HEAD. In both
round 4 and round 5, the installed package lagged the most recent commit by
hours. The drift check exists so you can **note** this in your report, not so
you can refuse to proceed.

```r
packageVersion("amatrix")
# then compare DESCRIPTION mtime to HEAD commit time
```

If your probe runs against stale code, your findings still apply — they are a
statement about what a user with the most recent install would see. If the
installed package is fresher than HEAD, something is weird, stop and ask.

---

## Rule 9 — Proposed close list is a proposal, not an action

Hunters never close bugs. Never. You `bd create` new ones and you `bd update`
with notes. The orchestrator does `bd close` after spot-verifying every
refutation with a runtime probe. If you write `bd close` in your shell, you
have broken the protocol.

Your refutation report's (d) section lists IDs for the orchestrator to
spot-check. Nothing more.

---

## Rule 10 — Filing discipline

Every `bd create` must include:

- **Title** prefixed with `[bug]`, naming the op and symptom in a way a reader
  can grep for. Bad: "crash on boundary case." Good: "cbind(X, c(1,2,3)) demotes
  adgeMatrix to dgeMatrix".
- **Description** with a minimal, runnable, fresh-process `Rscript -e '...'`
  repro. Not pseudo-code. Not "see hunter report." A block a reader can paste
  into a terminal.
- **Observed output** next to **expected output**.
- **Priority**: P1 for silent wrong answer or unrecoverable state. P2 for class
  demotion, crash with clear error, or scope limitation. P3 for UX/discoverability.
- **Root cause** if you have one, clearly labeled as inference vs. evidence.

If you cannot write a runnable repro, you do not have a bug. You have a
suspicion. File it as a note on an existing issue, not as a new bug.

---

## Rule 11 — On returning

Your report lives on disk. Your return message to the orchestrator is ≤200
words: drift status, count of probes run, count of bugs filed, bd IDs, one
sentence about the most important finding. Nothing else. The orchestrator will
read your file if they want detail.

---

## Rule 12 — What "done" means

A round is done when:

- Every hunter's report is on disk (orchestrator has `ls`d)
- Every claimed `bd create` has been verified with `bd show <id>`
- Every refutation has been runtime-spot-checked by the orchestrator
- False refutations are reopened with a counter-example note
- Duplicates are closed as dup-of with a cross-link
- The tracking issue for the round (`bd close amatrix-<id>`) is closed

Do not claim done until all six hold. The round-5 memory note
`round-5-second-false-refutation` exists because I almost closed 75h before
spot-checking. Don't be the next entry.

---

## Rule 13 — Propagate greps into the lint script

Round 4 promoted six regex rules from hunter greps into
`tools/lint-anti-patterns.R`. Net runtime cost: 0.13s. Net bugs surfaced: 2.
The highest compound-interest asset in this codebase is that file. After
your round, any grep pattern that found a real bug should become a lint rule.
Commit the delta.

This is housekeeping the hunter does on the way out, not something a separate
pass picks up later. Do it.

---

## Appendix — The operational tricks

### Fresh-process probes
```bash
Rscript -e 'suppressMessages(library(amatrix)); <probe>'
```
Always `Rscript -e` for scenarios that cross a process boundary. Same-session
`local(...)` is not fresh.

### Minimal-imports probe
```bash
Rscript -e 'library(amatrix); X <- adgeMatrix(matrix(1:12,3,4)); rowSums(X)'
```
No Matrix. This is the user's real state.

### Split fresh processes for serialization
```bash
Rscript -e 'library(amatrix); X <- <construct>; saveRDS(X, "/tmp/x.rds")'
Rscript -e 'library(amatrix); X <- readRDS("/tmp/x.rds"); <probe>'
```
Two processes, one pipe through `/tmp`. This is the only way to probe serialization
realistically.

### Orchestrator spot-check of a refutation
Take the bug description's repro. Run it verbatim in a fresh Rscript. Compare
output to what the bug predicts. If it matches the bug, the refutation was wrong,
the bug is live, add a note, do not close.

### Beads hygiene
```bash
bd list --status=open              # every round, first thing
bd search <op>                     # before every bd create
bd show <id>                       # before every mental close
bd update <id> --notes "..."       # how you add findings to existing bugs
bd remember --key <k> "<insight>"  # how you make a lesson survive the round
```

---

That's the playbook. The bugs you will find are composition bugs and scenario
bugs. The bugs you will falsely refute are the ones you reason about without
running. Run the code. Write the stub. Check the filesystem. Verify the list.
Propagate the greps. Close the tracker.

The next round will produce its own lessons. Add them here.

— The orchestrator, end of round 5
