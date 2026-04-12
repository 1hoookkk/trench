# TRENCH

## Role
- Human owns sound, taste, UX, product identity, and final musical judgment.
- Agent owns code, math, implementation, debugging, verification, and cleanup.
- Execute directly. No planning theater. No vague summaries.

## Ground truth stack
- Law (never drift): `SPEC.md`
- Live repo index (machine-readable): `REPO_TRUTH.json`
- Human jump table (readable): `REPO_TRUTH.md`
- Shipping gate + checklist: `SHIPPING.md`
- Local subtree rules (read before editing):
  - `trench-core/CLAUDE.md`
  - `pyruntime/CLAUDE.md`

## Modes
### Shape Bank mode (DEFAULT right now)
- Produce shapes (static sonic locations + explicit trajectories).
- Static targets are valid by definition.
- Do NOT score, rank, filter, prune, or “curate” outputs unless explicitly ordered.

### Trajectory mode (only when explicitly requested)
- When authoring a morph journey, judge by trajectory/motion; static plots alone are insufficient.

### Operator mode (repo hygiene / reproducibility)
- Refactors, plumbing, indexing, fixing broken paths, and making runs reproducible.
- Prefer enforcement over memory: stubs/mirrors, a single canonical truth stack, regenerated inventories.

### Shipping mode
- For release bodies, `SPEC.md` + `SHIPPING.md` win. Enforce gate/checklist.

## Working posture
- Treat the current filesystem, build wiring, imports/includes, registrations, and runtime call sites as the highest authority.
- If docs and code disagree, report the conflict explicitly and follow build-wired/runtime-wired truth.
- Treat orphaned files, experiments, graveyards, and legacy branches as non-canonical unless promoted in `REPO_TRUTH.md`.
- Prefer the shortest path to a verified result.
- Do not claim success unless you ran the relevant verification.

## Output discipline
- Return fully qualified, copy-pasteable code when changing functions or files.
- Do not use placeholders like `...existing code...` unless explicitly asked for a sketch.
- Focus on one concrete deliverable at a time unless the task clearly requires a sweep.

## Hard bans
- No RBJ cookbook for character filters.
- No compensation layers, fudge factors, or arbitrary saturation to rescue bad math.
- No parallel engine paths unless explicitly proven and approved.
- No gain baked into `c4`.
- No unapproved dependency changes.
- No shipping verbatim heritage extractions.
- No pole sanitization without a written failing test proving instability.

## Verification
- If the audible path or spectral math changes, state exactly what changed and why.
- Verify with the repo’s actual active commands from `REPO_TRUTH.md`, not stale remembered commands.
- When choosing between competing implementations, prefer the one wired into the current build or runtime.

## Repo truth maintenance
- Regenerate after moving/renaming entrypoints, paths, commands:
  - `python tools/update_repo_truth.py`
