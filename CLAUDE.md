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
- E-mu reference material: `docs/emu/`
- Architecture direction: `docs/architecture/zplane_truth.md`
- Math + contracts: `SPEC.md`
- Cartridge wire schema: `cartridge.schema.json`
- Operating modes: `MODES.md`
- Doctrine (rules, bans, verification, escalation): `DOCTRINE.md`
- Shipping bodies + rubrics: `BODIES.md`
- Phoneme authoring contract: `PHONEMES.md`
- Verify workspace: `./check`
- Active shipping plugin: `trench-juce/plugin/`
- Subtree rules: `trench-core/CLAUDE.md`, `pyruntime/CLAUDE.md`
- Repo index: live filesystem. Doctrine wins; no precomputed snapshot.
