---
name: artifact_trust_hierarchy
description: Tiered trust classification of RE artifacts — which sources are safe for CLEAN ingestion vs orientation-only
type: project
---

Trust hierarchy for DIRTY→CLEAN handoff, from NotebookLM audit (2026-03-14).

## Tier 1: High-Confidence Lookup Laws (directly ingestible)
- `semitone_table.json` — 64-entry semitone/frequency table. Cleanest numeric artifact.
- `q_radius_table.json` — 512-entry monotone radius table. Strong artifact for Q-law modeling.
- Four-corner bilinear interpolation — confirmed via FUN_1802c3d40 RE. Already implemented in trench-core.

## Tier 2: Sanitized Architectural Contracts (implementation-ready)
- **Morph Designer v2** CONTRACT.json — Tier B sanitized_behavioral. Highest trust for implementation. Template Families / Preset Wrappers / Compiled Runtime Recipes separation.
- **Morph Engine v1** CONTRACT.json — Clean design brief. Reusable filter families, immutable runtime programs.
- **Morph Engine Compiler v1** CONTRACT.json — Verified pseudocode for stage-slot compilation. Authoring→runtime boundary.

**Why:** These three contracts are the strongest assets for building the historical compiler path.
**How to apply:** Use these as the primary source when implementing any historical/compiler bridge work. Do not use raw captures directly.

## Tier 3: Patch-Specific Evidence (orientation only)
- `TalkingHedz_Complete.json` — DIRTY runtime capture. Case study, not global law. Moderate trust.
- `P2k_013_stage_plot_report.md` — Pre-reconciliation. Low trust for implementation.

## Tier 4: Unsafe / Do Not Implement From
- `morpheus_zplane_library.json` — Raw byte-level decoded payloads (f_byte, r_byte). Summary only for clean-room.
- `constants_region.json`, `wavetable1/2.json`, `uint16_table_candidates.json` — Weak/noisy, likely mis-typed memory.
- `filter_strings.json`, `filter_types.json` — DIRTY metadata with binary addresses. Orientation only.
