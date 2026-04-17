# TRENCH

## Role
- I own sound, taste, UX, product identity, final musical judgment.
- You own code, math, implementation, debugging, verification, cleanup.
- I write zero code. Every task runs agent end-to-end.
- Do the work. No planning doc. No vague summary layer.

## Product (v1.0, ship 2026-07-15)
TRENCH is an E-mu Z-plane filter instrument. Phoneme pills are baked
Morph Designer presets; MORPH is driven by standard DAW / synth mod
sources at playback. Windows JUCE 8 standalone + VST3.

## Core model
- Pill = one E-mu Morph Designer preset (Lo Morph / Hi Morph frames +
  6 sections with per-section type, frequency, Q/Gain).
- MORPH = interpolation between a pill's Lo and Hi frames. Driven by
  mod sources (LFO, envelope, velocity, mod wheel, DAW automation).
- Q = live resonance/gain parameter with a ±50% offset wheel, added to
  interpolated section Q/Gain.
- BODY = a named collection of pills sharing a design intent. Not a
  scheduled sequence.
- TYPE = selects the active body.

## Ground truth
- Active truth + architect prompts: `TRUTH_MAP.md`
- E-mu reference material: `docs/emu/`
- Architecture direction: `docs/architecture/zplane_truth.md`
- Math + contracts: `SPEC.md`
- Cartridge wire schema: `cartridge.schema.json`
- Operating modes: `MODES.md`
- Doctrine (rules, bans, verification, escalation): `DOCTRINE.md`
- Shipping bodies + rubrics: `BODIES.md`
- Phoneme authoring contract: `PHONEMES.md`
- Shipping phoneme pills: `cartridges/engine/` baked from
  `cartridges/engine/_source/token_inventory_unified_v2.json` +
  `cartridges/engine/_source/shapes/` via `tools/bake_phoneme_pills.py`
- Hardcoded Talking Hedz ROM: `trench-core/src/hedz_rom.rs` baked from
  `cartridges/engine/_source/heritage_designer_sections.json` via
  `tools/bake_hedz_const.py`
- Verify workspace: `./check`
- Active shipping plugin: `trench-juce/plugin/`
- Subtree rules: `trench-core/CLAUDE.md`
- Repo index: live filesystem. Doctrine wins; no precomputed snapshot.

## HARD RULE
any investigation that would read >3 files or >300 cumulative lines gets dispatched to an Explore subagent. The subagent reads the files, you only see its summary. The raw content never enters this parent context.
