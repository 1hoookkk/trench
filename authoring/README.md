# DF-II — Authoring Hub

```text
DF-II — Musical Filter.
```

The plugin is called **DF-II**. The faceplate reads "Musical Filter".
This folder is the canonical entrypoint for DF-II authoring and design work.

For filter work specifically, the repo-level contract is `FILTER_WORK.md`.
The canonical repo lives at `C:\Users\hooki\Trench` (codebase paths still
use the legacy `trench-` prefix); external trees are upstream sources, not
parallel authorities. The two rules: **DF-II is a Musical Filter** (DF-II
serial biquad cascade, body-first authoring), and **DF-II does not ship
E-mu filters**. Both at once.

Use this folder when the task is:

- authoring bodies
- shaping sonic targets
- using sonic tables
- preparing MCP-assisted authoring sessions
- consolidating design direction
- preventing E-mu/P2K clone cartridges from entering release candidates

Do **not** use this folder as the source of shipping runtime truth.
Shipping runtime remains in:

- `C:\Users\hooki\Trench\trench-juce\plugin\`
- `C:\Users\hooki\Trench\trench-core\`

## Start Here

Read in this order:

0. `SESSION_START.md` — four-project setup + launch prompts (do this once per workstation)
1. `FILTER_WORK.md`
2. `GROUND_TRUTH.md`
3. `CLAUDE.md`
4. `ground_truth/README.md`
5. `sonic_tables/README.md`
6. `workbench/README.md`
7. `design/README.md`
8. `history/README.md`

## What This Hub Solves

Before this hub, authoring truth was scattered across:

- `C:\Users\hooki\Trench`
- `C:\Users\hooki\trenchwork_clean`
- `C:\Users\hooki\trenchwork`
- `C:\Users\hooki\do-it`
- `C:\Users\hooki\trench_re_vault`

This hub gives one canonical doorway for filter work.

## Working Split

There are three practical domains:

1. `ground_truth/`
   - shipping runtime, cartridge contract, and authoring law

2. `workbench/`, `sonic_tables/`, `design/`, `history/`
   - authoring and design surfaces

3. outside references
   - useful material in other trees; never shipping truth by itself

Do not flatten references into shipped filters.

## Workstation Layout

As of 2026-05-10, all filter-work files live under this folder. The legacy
`forge/` and root-level `schemas/` directories have been folded in. Outside
this folder, the only filter-work-relevant trees are runtime
(`runtime/trench-core/`), plugin (`trench-juce/`), shipped cartridges
(`cartridges/`), research docs (`docs/`), and session state (`dev/`).

| Subdirectory | Purpose |
|---|---|
| `compilers/` | Active compilers, smoke tests, generic generators. `parity_null.py`, `compile_raw.py`, `bark.py`, `body_dsl.py`, etc. live here. |
| `tools/` | Authoring utilities that are not compilers — body authors, audition harnesses, phoneme studios, proof scripts. |
| `audition.py` | Canonical audition pipeline. `python authoring/audition.py <cartridge.json>` renders the four corner WAVs through the cascade. |
| `schemas/` | TRENCH-native JSON schemas (`trench.authoring_path.cube.v1`, `trench.compiled.cube_surface.v1`). |
| `bodies/`, `design/`, `frame_maps/`, `sonic_tables/`, `workbench/` | Original body drafts, IRs, design surfaces. |
| `ground_truth/`, `heritage/`, `imported_x3/`, `emu_designer/`, `dirty_airlock/` | Research and reference material. **Not** a shipping path. Some of these contain E-mu material pending triage. |
| `_research/` | Quarantined / archived authoring trajectories that were tried and rejected. Includes the pre-veto Peak/Shelf IRs, `HANDOFF_GPT55.md`, `SEED.md`. |
| `auditions/`, `compiled/`, `plots/` | Audition output and rendered plots. |
| `inbound/` | Staging for outside material before promotion (per FILTER_WORK.md promotion rule). |
| `history/` | Historical precedent summaries from upstream trees. |
| `q_laws/`, `diagnostics/`, `scratch/` | Engine-axis tooling and working files. |
| `atlas/`, `function_packs/`, `ui/` | Pole atlas, function-generator packs, UI work. |

## Root Launcher Pattern

User-facing pipelines may keep a thin launcher at the repo root that
delegates to the real harness in `authoring/compilers/`. Example:
`trench_eos_pipeline_compiler.py` at the root simply runpys
`authoring/compilers/trench_eos_pipeline_compiler.py`. This keeps
documented commands reachable from repo root without re-fragmenting the
workstation. Future migration sessions: do **not** fold launchers back in
or treat them as duplicates — they are the single-line UX surface for the
real harness.

## Verification

The canonical filter-work check is `./check` at the repo root. It runs the
doc-set check, P2K parity null (now reading
`authoring/compilers/parity_null.py`), Rust cascade parity, and the cartridge
schema validation.

## Current Status

Phase 2 of the consolidation is complete (all forge/ + schemas/ folded in,
E-mu clones purged). See `MIGRATION_STATE.md` for what remains.
