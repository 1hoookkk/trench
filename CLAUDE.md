# TRENCH

## Role
- I own sound, taste, UX, product identity, final musical judgment.
- You own code, math, implementation, debugging, verification, cleanup.
- I write zero code. Every task runs agent end-to-end.
- Do the work. No planning doc. No vague summary layer.

## Product (v1.0, ship 2026-07-15)
TRENCH is an authoring instrument: phoneme tokens on a time grid become the
filter trajectory. Windows JUCE 8 standalone + VST3.

## Core model
- Cartridge-based morph filter instrument.
- Cartridges are phoneme pills: morph-invariant, Q-varied spectral
  states identified by short grid labels.
- MORPH = travel between pills on the Looperator time grid, scheduled
  by the compiler. Pills themselves do not morph internally.
- Q = playback-surface tonal shaping of whichever pill the playhead
  is currently on.
- BODY = a named preset sequence of pills (Speaker Knockerz, Aluminum
  Siding, Small Talk Ah-Ee, Cul-De-Sac per BODIES.md). Compiler TBD.

## Ground truth
- Active truth + architect prompts: `TRUTH_MAP.md`
- Math + contracts: `SPEC.md`
- Cartridge wire schema: `cartridge.schema.json`
- Operating modes: `MODES.md`
- Doctrine (rules, bans, verification, escalation): `DOCTRINE.md`
- Shipping bodies + rubrics: `BODIES.md`
- Phoneme authoring contract: `PHONEMES.md`
- Shipping phoneme pills: `cartridges/engine/` baked from
  `vault/_phonemes/token_inventory_unified_v2.json` via
  `tools/bake_phoneme_pills.py`
- Hardcoded Talking Hedz ROM: `trench-core/src/hedz_rom.rs` baked from
  `vault/_phonemes/heritage_designer_sections.json` via
  `tools/bake_hedz_const.py`
- Verify workspace: `./check`
- Active shipping plugin: `trench-juce/plugin/`
- Subtree rules: `trench-core/CLAUDE.md`, `pyruntime/CLAUDE.md`
- Repo index: live filesystem. Doctrine wins; no precomputed snapshot.
