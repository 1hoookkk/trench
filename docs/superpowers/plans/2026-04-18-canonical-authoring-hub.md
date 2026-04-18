# Canonical Authoring Hub Migration Plan

**Date:** 2026-04-18  
**Status:** proposed  
**Goal:** collapse the current cross-repo authoring/design sprawl into one canonical hub under `C:\Users\hooki\Trench\authoring\` without losing shipping/runtime truth.

## Objective

Make one place you can point a model or an operator at and say:

- this is where authoring lives
- this is where design lives
- this is where sonic tables live
- this is where ground truth is explained
- this is where session prompts and agent guidance start

Shipping runtime truth stays where it is today:

- `C:\Users\hooki\Trench\trench-juce\plugin\`
- `C:\Users\hooki\Trench\trench-core\`

This plan creates a canonical authoring/design hub beside shipping, not inside shipping code.

## Non-negotiable rule

Single top-level hub does **not** mean single trust level.

If everything is physically collapsed into one tree with no internal boundary, the repo will immediately start confusing:

- shipping truth
- authoring truth
- DIRTY RE provenance
- historical precedent

So the hub is one root with explicit internal zones:

- `ground_truth/`
- `workbench/`
- `design/`
- `dirty_airlock/`
- `history/`

That is the minimum structure that avoids lying to yourself later.

## Proposed target structure

```text
authoring/
  README.md
  CLAUDE.md
  GROUND_TRUTH.md
  SOURCES.md
  MIGRATION_STATE.md

  ground_truth/
    shipping_runtime.md
    cartridge_contract.md
    authoring_contract.md
    provenance_map.md
    cleanroom_boundary.md

  sonic_tables/
    README.md
    tables/
    notebooks/
    exports/

  workbench/
    sessions/
    prompts/
    vault/
    mcp/
    tools/
    foundry/
    forge/

  design/
    ui_identity/
    faceplate/
    rollers/
    screenshots/
    design_system/

  history/
    do_it/
    earlier_runtime/
    solved_decisions/

  dirty_airlock/
    README.md
    manifests/
    sanitized_handoffs/
    raw_index/
```

## Canonical entrypoints

After migration, the only human/model entrypoints should be:

- `C:\Users\hooki\Trench\authoring\README.md`
- `C:\Users\hooki\Trench\authoring\CLAUDE.md`
- `C:\Users\hooki\Trench\authoring\GROUND_TRUTH.md`

Everything else is reached from there.

## What goes where

### 1. Ground truth

Purpose: teach the domain correctly in one pass.

Inputs:

- `C:\Users\hooki\Trench\TRUTH_MAP.md`
- `C:\Users\hooki\Trench\SPEC.md`
- `C:\Users\hooki\Trench\DOCTRINE.md`
- `C:\Users\hooki\Trench\CLAUDE.md`
- `C:\Users\hooki\Trench\PHONEMES.md`
- `C:\Users\hooki\Trench\BODIES.md`
- selected `docs/emu/**`
- selected `docs/architecture/**`

Output:

- `authoring\GROUND_TRUTH.md` = one-file domain brief
- `authoring\ground_truth\shipping_runtime.md` = what actually ships
- `authoring\ground_truth\authoring_contract.md` = what bodies/authoring are allowed to do
- `authoring\ground_truth\cleanroom_boundary.md` = what can never be pulled from DIRTY into clean surfaces

### 2. Sonic tables

Purpose: canonical target vocabulary and mapping surface for authoring.

Inputs:

- `C:\Users\hooki\Trench\docs\sonic_tables\`
- `C:\Users\hooki\trenchwork_clean\docs\sonic_tables\`
- relevant notebook/skill references already pointing at sonic-table MCP material

Output:

- one canonical `authoring\sonic_tables\README.md`
- one canonical table location
- one canonical notebook/export location

### 3. Workbench / forge / foundry / MCP

Purpose: put the actual authoring factory in one discoverable place.

Inputs:

- `C:\Users\hooki\trenchwork\trench-forge\`
- `C:\Users\hooki\trenchwork\trench-foundry\`
- `C:\Users\hooki\trenchwork\trench-mcp\`
- `C:\Users\hooki\trenchwork_clean\trench-mcp\`
- selected `tools/` and session docs from `C:\Users\hooki\Trench`

Output:

- canonical descriptions and links under `authoring\workbench\`
- not necessarily full code movement on day one
- day one can be an index + wrappers + migration state

### 4. Design

Purpose: stop UI/faceplate/roller work from living in random trees.

Inputs:

- `C:\Users\hooki\Trench\docs\identity\`
- `C:\Users\hooki\do-it\` roller and UI precedent assets
- `C:\Users\hooki\trenchwork_clean\assets\`
- `C:\Users\hooki\trenchwork_clean\station-editor\`

Output:

- `authoring\design\ui_identity\`
- `authoring\design\rollers\`
- `authoring\design\faceplate\`

### 5. History

Purpose: preserve prior solved runtime work without letting it masquerade as live truth.

Inputs:

- `C:\Users\hooki\do-it\`

Output:

- `authoring\history\do_it\README.md`
- selected decision summaries and validation references
- not a live dependency

### 6. DIRTY airlock

Purpose: make `trench_re_vault` reachable from the hub without pretending it is clean.

Inputs:

- `C:\Users\hooki\trench_re_vault\`

Output:

- `authoring\dirty_airlock\README.md`
- sanitized handoff index
- manifest of where raw artifacts live
- explicit “never import into shipping/runtime code” statement

This is the compromise that gives you one canonical hub without erasing the cleanroom boundary.

## Source mapping

| Source tree | Role after migration | Destination |
|---|---|---|
| `C:\Users\hooki\Trench` | shipping + root canonical host | stays put; gains `authoring/` |
| `C:\Users\hooki\trenchwork_clean` | clean authoring input | `authoring/workbench/`, `authoring/sonic_tables/`, `authoring/design/` |
| `C:\Users\hooki\trenchwork` | forge/foundry/workbench architecture | `authoring/workbench/` |
| `C:\Users\hooki\do-it` | historical runtime precedent | `authoring/history/do_it/` |
| `C:\Users\hooki\trench_re_vault` | DIRTY provenance | `authoring/dirty_airlock/` manifest + sanitized handoffs only |

## Migration phases

## Phase 1: Scaffold the hub

Create only structure and entrypoints. No content movement yet.

### Tasks

- Create `authoring/`
- Create `authoring/README.md`
- Create `authoring/CLAUDE.md`
- Create `authoring/GROUND_TRUTH.md`
- Create the six subfolders
- Create `authoring/MIGRATION_STATE.md`

### Verify

- folders exist
- all top-level entrypoint files exist
- no shipping code imports change

### Done when

- a model can be pointed at `authoring/README.md`
- a human can understand the hub layout in under 60 seconds

## Phase 2: Ground truth consolidation

Write the actual one-file domain brief and trust hierarchy.

### Tasks

- compress current repo truth into `authoring/GROUND_TRUTH.md`
- create `authoring/ground_truth/provenance_map.md`
- create `authoring/ground_truth/cleanroom_boundary.md`
- cross-link to shipping runtime files, not just docs

### Verify

- `GROUND_TRUTH.md` explicitly separates:
  - shipping runtime truth
  - authoring/factory truth
  - DIRTY provenance truth
- no statement contradicts live runtime wiring

### Done when

- you can hand one file to a model and it lands in the right domain frame

## Phase 3: Sonic tables unification

Make sonic tables one canonical surface.

### Tasks

- inventory `docs/sonic_tables/` in `Trench`
- inventory `docs/sonic_tables/` in `trenchwork_clean`
- choose one canonical storage shape
- move or mirror the winning content into `authoring/sonic_tables/`
- update references from session/skill docs to the new canonical path

### Verify

- one README
- one canonical table location
- old paths either point forward or are marked legacy

### Done when

- “where are the sonic tables?” has one answer

## Phase 4: Workbench consolidation

Make forge/foundry/MCP/session surfaces discoverable from one place.

### Tasks

- inventory active crates and tools in `trenchwork` and `trenchwork_clean`
- classify each as:
  - keep
  - wrap
  - mirror
  - archive
- create `authoring/workbench/README.md`
- add wrappers or indexes to the active workbench surfaces

### Verify

- each active workbench surface has one canonical pointer
- duplicate surfaces are called out, not silently kept

### Done when

- body-session startup no longer requires tribal memory across trees

## Phase 5: Design consolidation

Move the visual/UI reference mess into one place.

### Tasks

- collect roller precedent from `do-it`
- collect current identity/faceplate assets from `Trench`
- collect station-editor/design assets from `trenchwork_clean`
- create `authoring/design/README.md`
- classify:
  - live visual direction
  - useful precedent
  - dead experiments

### Verify

- one canonical place for screenshots, faceplate references, and control semantics

### Done when

- design work stops bouncing between `do-it`, `Trench`, and `trenchwork_clean`

## Phase 6: History and DIRTY integration

Bring the useful parts in without corrupting trust.

### Tasks

- create `authoring/history/do_it/README.md`
- capture key runtime decisions from `do-it`
- create `authoring/dirty_airlock/README.md`
- create sanitized handoff index for `trench_re_vault`
- explicitly forbid runtime imports from `dirty_airlock/`

### Verify

- DIRTY material is reachable from the hub
- DIRTY material is still unmistakably DIRTY

### Done when

- provenance is available without leaking into clean authoring/runtime claims

## Phase 7: Rewrite old entrypoints

After the hub exists, old entrypoints should stop competing with it.

### Tasks

- update root `CLAUDE.md` to point authoring/design questions at `authoring/`
- update `TRUTH_MAP.md` to mention the new hub
- update session skills and startup docs to use the canonical paths
- mark legacy paths explicitly

### Verify

- at least 80% of “where do I start?” docs point at `authoring/`

### Done when

- the hub is not just present, it is the default doorway

## What not to do

- Do not move shipping runtime code into `authoring/`
- Do not pretend `trench_re_vault` is clean just because it sits under the same top-level folder
- Do not bulk-copy entire trees on day one
- Do not keep duplicate sonic-table or session-start surfaces without naming one canonical winner
- Do not rewrite history docs before the canonical hub exists

## Recommended execution order

1. Phase 1 — scaffold
2. Phase 2 — ground truth
3. Phase 3 — sonic tables
4. Phase 5 — design
5. Phase 4 — workbench
6. Phase 6 — history + DIRTY
7. Phase 7 — rewrite entrypoints

Reason:

- ground truth and sonic tables are the highest-leverage model/operator entrypoints
- design is currently scattered and cheap to consolidate
- forge/foundry consolidation is larger and should happen after the hub language is stable

## Success criteria

- there is exactly one canonical authoring/design root: `authoring/`
- there is exactly one canonical model handoff file: `authoring/GROUND_TRUTH.md`
- there is exactly one canonical agent startup guide: `authoring/CLAUDE.md`
- sonic tables have one canonical home
- design references have one canonical home
- workbench surfaces are indexed from one place
- `do-it` is preserved as history, not mistaken for live runtime
- `trench_re_vault` is reachable but still explicitly quarantined

## First concrete slice

Do this first:

1. Create `authoring/`
2. Create `authoring/README.md`
3. Create `authoring/CLAUDE.md`
4. Create `authoring/GROUND_TRUTH.md`
5. Create `authoring/sonic_tables/README.md`
6. Create `authoring/design/README.md`
7. Create `authoring/dirty_airlock/README.md`

That gives you the canonical doorway immediately, before any large migration work.
