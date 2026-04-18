# TRENCH Authoring Hub

This is the canonical entrypoint for TRENCH authoring and design work.

Use this folder when the task is:

- authoring bodies
- shaping sonic targets
- using sonic tables
- preparing MCP-assisted authoring sessions
- consolidating design direction
- understanding clean vs DIRTY provenance

Do **not** use this folder as the source of shipping runtime truth.
Shipping runtime remains in:

- `C:\Users\hooki\Trench\trench-juce\plugin\`
- `C:\Users\hooki\Trench\trench-core\`

## Start Here

Read in this order:

0. `SESSION_START.md` — four-project setup + launch prompts (do this once per workstation)
1. `GROUND_TRUTH.md`
2. `CLAUDE.md`
3. `ground_truth/README.md`
4. `sonic_tables/README.md`
5. `workbench/README.md`
6. `design/README.md`
7. `dirty_airlock/README.md`

## What This Hub Solves

Before this hub, authoring truth was scattered across:

- `C:\Users\hooki\Trench`
- `C:\Users\hooki\trenchwork_clean`
- `C:\Users\hooki\trenchwork`
- `C:\Users\hooki\do-it`
- `C:\Users\hooki\trench_re_vault`

This hub gives one canonical doorway while keeping the internal trust split explicit.

## Internal Trust Split

There are three real domains:

1. `ground_truth/`
   - clean statements about shipping runtime, cartridge contract, and authoring law

2. `workbench/`, `sonic_tables/`, `design/`, `history/`
   - clean authoring and design surfaces

3. `dirty_airlock/`
   - DIRTY provenance and sanitized handoff boundary

Do not flatten those three into one trust level.

## Current Status

This hub is now the canonical authoring/design scaffold.
Migration is in progress; see `MIGRATION_STATE.md`.

