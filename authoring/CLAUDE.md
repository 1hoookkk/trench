# TRENCH Authoring Hub Guide

## Purpose

This file is the startup guide for models or agents entering the TRENCH authoring hub.

The target is not casual repo browsing. The target is:

- correct domain framing
- clean separation of shipping truth vs authoring truth vs DIRTY provenance
- fast orientation for body work and design work

## First Read Order

1. `GROUND_TRUTH.md`
2. `ground_truth/README.md`
3. `sonic_tables/README.md` when the task touches targets, vocabulary, or body landmarks
4. `workbench/README.md` when the task touches forge, foundry, MCP, vault, or sessions
5. `design/README.md` when the task touches UI, faceplate, control semantics, or visual references
6. `dirty_airlock/README.md` when the task touches RE provenance or cleanroom boundaries
7. `history/README.md` only when old solved states or historical precedent matter

## Hard Framing

- Shipping runtime truth is still the live code in `C:\Users\hooki\Trench\trench-juce\plugin\` and `C:\Users\hooki\Trench\trench-core\`.
- This hub is the canonical authoring/design doorway, not a replacement for build wiring.
- Do not treat DIRTY provenance as clean truth just because it is reachable from here.
- Do not treat historical runtime precedent as live truth.

## What Belongs Here

- body-session startup
- sonic-table vocabulary
- workbench and MCP entrypoints
- design direction
- source/provenance maps
- cleanroom boundary statements

## What Does Not Belong Here

- shipping DSP implementation details duplicated from runtime files
- raw DIRTY artifacts
- stale archive dumps presented as active truth

## Default Posture

- prefer live shipping wiring for runtime claims
- prefer this hub for authoring/design claims
- prefer sanitized handoffs over raw RE material
- prefer one canonical pointer over duplicate explanations

