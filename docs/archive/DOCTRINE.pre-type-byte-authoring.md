# TRENCH doctrine

## Working rules
- Filesystem, build wiring, imports, registrations, and runtime call sites are
  the highest authority.
- If `SPEC.md`, `CLAUDE.md`, `AGENTS.md`, or any doc disagrees with the repo,
  report the conflict and follow verified wiring.
- Prefer build-wired and runtime-wired truth over orphaned files or nearby
  filenames.
- Legacy, prototype, graveyard, and research-only paths are non-canonical
  unless explicitly promoted.
- Pick the shortest path to a verified result.
- Do not claim something works unless you ran it.

## Hard bans
- No RBJ cookbook for character filters.
- No compensation layers, fudge factors, or saturation to rescue bad math.
- No parallel engine paths unless explicitly proven and approved.
- No gain baked into `c4`.
- No unapproved dependency changes.
- No shipping verbatim heritage extractions.
- No pole sanitization without a failing test proving instability.
- DSP cascade, interpolation order, and cartridge format are frozen. No edits
  to `trench-core/src/cascade.rs` or
  `trench-juce/plugin/source/TrenchEngine.cpp` cascade logic.

## Verification discipline
- If the audible path changes, state what changed and why.
- If spectral math changes, state what changed and why.
- When multiple implementations exist, determine the active one from
  build/runtime evidence before editing.
- When docs are stale, say so and continue from verified repo truth.

## Stop and ask
- A change materially alters DSP/audible behavior.
- A change modifies topology, interpolation order, or cartridge format.
- A change adds/changes dependencies or introduces audio-thread alloc/locking.
- Competing implementations exist and build-wired truth is unclear.
