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
- MORPH = travel between authored spectral states.
- Q = aggression / emphasis behavior across that travel.
- TYPE = selects the authored body.

## Ground truth
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
