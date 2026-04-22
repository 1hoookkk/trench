# CLAUDE.md — Forge / Authoring Hub

Canonical guide for the authoring/generation/research half of TRENCH.
The plugin-and-filter-engine half lives at the repo-root `CLAUDE.md`;
do not touch runtime invariants from here, and do not import runtime
internals into forge tooling.

## Scope

This side owns the non-deterministic path:

    sound target  →  body spec  →  forge generation  →  candidate  →  audition  →  accept / reject

Surface:

- `authoring/` — this hub, first read target for body and design work.
- `tools/cube_authoring/` — gates, probes, triage scripts.
- `trench-juce/forge/` and `trench-juce/trench-mcp/` — the MCP server
  for measurement (`analyze_state`, `scan_vault`, `probe_frequency`,
  `scan_stability`, `find_corridors`, `compare_states`, `taste_weights`,
  `vault_seeds`).
- `vault/` — candidate bodies for triage and taste scoring.
- `docs/sonic_tables/trench_sonic_tables.ipynb` — vocabulary and target
  territory reference. Notebook view at `http://localhost:8888/` is
  convenience only; the `.ipynb` on disk is the source of truth.
- RE artifacts: `trench_re_vault/` (sibling tree), Ghidra decodes,
  Morpheus cube decodes, proven docs under
  `docs/notebooklm/packs/ground-truth/`.
- `cartridges/cubes/_morpheus/`, Morpheus pole vaults and gate verdicts.
- `BODIES.md` — the 4 launch bodies; the normative naming contract.

Out of scope (engine's job, at repo root): the frozen DSP cascade,
compile math (type 1/2/3, pole→kernel, Peak/Shelf bilinear), the
cartridge schema, plugin/runtime internals.

## First Read Order

1. `GROUND_TRUTH.md`
2. `ground_truth/README.md`
3. `sonic_tables/README.md` — when the task touches targets, vocabulary,
   or body landmarks.
4. `workbench/README.md` — when the task touches forge, foundry, MCP,
   vault, or sessions.
5. `design/README.md` — when the task touches UI, faceplate, control
   semantics, or visual references.
6. `dirty_airlock/README.md` — when the task touches RE provenance or
   cleanroom boundaries.
7. `history/README.md` — only when old solved states or historical
   precedent matter.

If the task mentions sonic tables, vowel targets, vault coverage, or
formant territory, use the local notebook at
`docs/sonic_tables/trench_sonic_tables.ipynb` and the configured
`trench-mcp` server from `.claude/settings.local.json`. In Claude Code,
prefer the server over spawning the binary manually.

## Hard Framing

- Shipping runtime truth is the live code in `trench-core/` and
  `trench-juce/plugin/`. This hub is an authoring doorway, not a
  replacement for build wiring.
- Do not treat DIRTY provenance as clean truth just because it is
  reachable from here.
- Do not treat historical runtime precedent as live truth.
- Do not use the sonic tables to justify breaking frozen DSP invariants.

## Source Hierarchy

Keep these layers separate:

1. **P2K finished bodies** — primary behavioral oracle.
2. **Atlas anchors** — evaluation baseline only.
3. **Morpheus cubes** — secondary structural priors only (pole seeds,
   not full bodies).
4. **Zero-forcing** — derived tool, not universal law.
5. **Four-law body framing** — TRENCH design doctrine, not heritage
   fact.

If you cross layers, say so explicitly as an inference.

## Hand-off to Engine

The forge emits **only** these artifacts, via the schema:

- `trench.authoring_path.pill.v1` — 4-keyframe authoring JSON.
- `trench.authoring_path.cube.v1` — 8-corner authoring JSON.
- `compiled-v1` — the wire-ready cartridge (produced by running the
  engine-side compilers on the above).

Hand-off protocol:

1. Produce an authoring-path JSON (corner kinds = `morph_designer_ref`,
   `morph_designer_inline`, `peak_shelf_ref`, `peak_shelf_inline`, or
   `native_poles`).
2. Validate with `./check`.
3. Run `tools/compile_pill.py` or `tools/compile_cube.py` to emit
   `compiled-v1`.
4. Audition via engine runtime; measure with MCP.
5. On accept, the `compiled-v1` file lands in `cartridges/engine/`.

The forge must never ship a file that bypasses the schema, never embed
runtime coefficient data inside authoring JSON, and never patch
`compiled-v1` after compile.

## Hard Rules (from authoring doctrine)

- Do not do cube-first authoring.
- Do not do stage-shuffle generation.
- Do not do direct coefficient search as the main workflow.
- Do not target Atlas terrain classes as styles to imitate.
- Do not turn one RE pattern into a universal generation rule.
- Do not make zero-forcing mandatory everywhere.
- Do not blur TRENCH product architecture with P2K heritage behavior.
- Corners are mechanism. Checkpoints are intent. Author the path, not
  the endpoints.

## Pre-Listening Gates

Reject candidates before audition if any fails:

- Low-Q identity collapse.
- Bounded gain / no trivial gain-spike cheat.
- Morph coherence through the interior, not just corners.
- Cancellation / subtractive profile when the body calls for it.
- Non-commodity terrain or behavior.
- Edge behavior that fails incoherently instead of intentionally.

Do not mix DSP correctness checks into this list. Audio-thread safety,
interpolation order, and kernel invariants are the engine's concern.

## MCP Tool Routing

- `analyze_state` — one corner or interior point vs one sonic-table
  target.
- `compare_states` — confirm a morph or Q move preserves identity.
- `scan_stability` — find cliffs before proposing a body as viable.
- `find_corridors` — locate usable morph paths; don't assume corners
  imply continuity.
- `probe_frequency` — explicit null or emphasis checks (e.g., the 12 dB
  Choke notch in Speaker Knockerz).
- `scan_vault` — rank vault bodies; territory / coverage audits.
- `vault_seeds` — top-N vault starting points for reuse.
- `taste_weights` — inspect the calibrated taste model before changing
  scoring logic.

Separate measured MCP output from notebook-derived targets. State any
cross-layer leap as an inference. Convert each result into exactly one
next move: author, reject, retune, or audit.

## Default Posture

- Prefer live shipping wiring for runtime claims (engine-side).
- Prefer this hub for authoring / design claims.
- Prefer sanitized handoffs over raw RE material.
- Prefer one canonical pointer over duplicate explanations.

## Memory

Memory lives at
`C:\Users\hooki\.claude\projects\C--Users-hooki-Trench\memory\`. Treat
it as shared across both halves; tag entries as engine-side or
forge-side in their description when the distinction matters.
