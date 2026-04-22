# CLAUDE.md — Plugin / Filter Engine

Canonical guide for the plugin-and-filter-engine half of TRENCH. The forge
half lives at `authoring/CLAUDE.md`; do not touch runtime invariants from
there, and do not import forge tooling into anything that ships.

## Scope

This side owns the deterministic path:

    heritage spec  →  compiled-v1 cartridge  →  runtime cascade  →  user audio

Surface:

- `trench-core/` — pure DSP crate (frozen cascade, cartridge loader,
  resampler, engine). See `trench-core/CLAUDE.md` for crate invariants.
- `trench-live/` — live/editor Rust.
- `trench-juce/plugin/` — JUCE host (runtime only; the `CLAUDE.md` there
  is a stub that points back here).
- `tools/compile_pill.py`, `tools/compile_cube.py`,
  `tools/bake_hedz_const.py`, `tools/pole_math.py` — offline compile
  pipeline. These live here because they are deterministic math, not
  generation.
- `cartridge.schema.json`, `schemas/` — the wire contract with forge.
- `cartridges/engine/` — shipping compiled-v1 content.

Out of scope (forge's job, under `authoring/`): body generation, taste
scoring, vault scanning, MCP measurement, notebook exploration, RE
artifacts, Morpheus decode. Consume schema-valid JSON from them; never
depend on their internals.

## Role

You own implementation, debugging, verification, code quality.
I own product, sound, UX, final judgment.

## Execution

- No preamble. No plan unless asked. Make the change.
- Use exact file paths.
- If you change code, give complete, copy-pasteable blocks.
- If a function changes, output the whole function.
- If context is missing or conflicting, say exactly what is missing.

## Audio Rules (hard)

- No allocations, locks, file I/O, or parsing on the audio thread.
- DSP topology is frozen: 12-stage DF2T cascade, 6 active + 6 passthrough,
  per-sample coefficient ramping, 32-sample control block.
- Interpolation order is frozen: Q-first, then morph. Do not reorder
  silently.
- If sound changes, say exactly what changed and why.

## Project Rules

- Follow the existing architecture. Do not invent new patterns without
  reason.
- Prefer the smallest change that preserves clarity.
- Reuse existing code before adding new code.
- Do not add dependencies without asking.
- Keep docs and code aligned. If a decision becomes durable, propose the
  doc update.

## Schema = boundary with forge

`cartridge.schema.json` is the only contract. Forge emits `compiled-v1`
or `trench.authoring_path.*.v1` JSON. Engine validates with `./check`,
then consumes. Never accept coefficient data that has not been
schema-validated. Never let forge internals (taste weights, vault
rankings, MCP outputs, notebook state) enter shipping code.

## Verification

- Run the relevant tests, not guesswork. `cargo test --lib` in
  `trench-core/`; the parity gates (`tests/hedz_cascade.rs`,
  `tests/talking_hedz_parity.rs`); schema checks
  (`tests/cartridge_schema_tests.rs`).
- Never claim something works unless you ran it.
- If bash's `/usr/bin/link` shadows MSVC `link.exe`, point rustc at MSVC
  directly via `CARGO_TARGET_X86_64_PC_WINDOWS_MSVC_LINKER`.
- If you could not run verification, say so plainly.

## Escalate Only If

- sound will materially change
- UX will materially change
- the change is hard to reverse
- the change adds architectural debt

## Filter-path reference

Which compile path turns a body's geometry into `compiled-v1` — Morph
Designer recipes, Peak/Shelf Morph two-frame sweeps, or native Morpheus
poles — is decision support, not always-on law. See the
`filter-path-selection` skill.
