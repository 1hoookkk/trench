# TRENCH

## Role
- I own sound, taste, UX, product identity, and final musical judgment.
- You own code, math, implementation, debugging, verification, and cleanup.
- Do the work. No planning doc. No vague summary layer between you and the task.

## Ground truth stack
- Law (never drift): `SPEC.md`
- Live repo index (machine-readable): `REPO_TRUTH.json`
- Human jump table (readable): `REPO_TRUTH.md`
- Shipping gate + checklist: `SHIPPING.md`
- Local subtree rules (read before editing):
  - `trench-core/CLAUDE.md`
  - `pyruntime/CLAUDE.md`

## What we are doing (modes)
### Shape Bank mode (DEFAULT right now)
- We are collecting **shapes**: static sonic locations + explicitly authored trajectories.
- **Static targets are valid by definition.** Do not score, rank, prune, or “optimize them away”.
- Your job is to produce correct artifacts and keep the index truthful, not to curate taste.

### Trajectory mode (only when explicitly requested)
- When authoring a morph journey, judge the **trajectory**; a single static plot is insufficient.
- Motion judgement applies **only** to trajectories, not to the whole project by default.

### Operator mode (repo hygiene / reproducibility)
- Refactors, plumbing, indexing, fixing broken paths, and making runs reproducible.
- Prefer mechanical enforcement over “remembering”: stubs/mirrors, single canonical truth files, and regenerated inventories.

### Shipping mode (when working on release bodies)
- Follow `SPEC.md` + `SHIPPING.md`. Shipping bodies must satisfy the release gate and product law.
- If Shape Bank rules conflict with Shipping rules, Shipping rules win for shipping bodies.

## Core model (always true)
- TRENCH is a cartridge-based morph filter instrument.
- MORPH = travel between authored spectral states.
- Q = aggression / emphasis behavior across that travel.
- TYPE = selects the authored body.

## Working rules
- Treat the current filesystem, build wiring, imports/includes, registrations, and runtime call sites as the highest authority.
- If `SPEC.md`, `REPO_TRUTH.json`, `REPO_TRUTH.md`, `CLAUDE.md`, and `AGENTS.md` disagree with the repo, report the conflict explicitly and follow verified wiring.
- Prefer build-wired and runtime-wired truth over orphaned files, old notes, or nearby filenames.
- Treat legacy, prototype, graveyard, and research-only paths as non-canonical unless promoted in `REPO_TRUTH.md`.
- Pick the shortest path to a verified result.
- Do not claim something works unless you ran it.

## Hard bans
- No RBJ cookbook for character filters.
- No compensation layers, fudge factors, or arbitrary saturation to rescue bad math.
- No parallel engine paths unless explicitly proven and approved.
- No gain baked into `c4`.
- No unapproved dependency changes.
- No shipping verbatim heritage extractions.
- No pole sanitization without a written failing test proving instability.

## Verification discipline
- If the audible path changes, state exactly what changed and why.
- If spectral math changes, state exactly what changed and why.
- When multiple implementations exist, determine the active one from build/runtime evidence before editing.
- When docs are stale, say so plainly and continue from verified repo truth.

## Repo truth resolution
Resolve these from the live repo before relying on them:
- active shipping plugin path
- active build/test commands
- canonical runtime path
- canonical authoring/forge path
- quarantined or research-only paths

Do not turn unresolved repo guesses into doctrine.

## Stop and ask if
- A change materially alters DSP/audible behavior (not just refactor).
- A change modifies topology, interpolation order, or cartridge format.
- A change adds/changes dependencies or introduces audio-thread allocation/locking.
- There are competing implementations and build-wired truth is unclear.
