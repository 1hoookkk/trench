# DF-II

**DF-II — Musical Filter.**

The plugin is called **DF-II**. The faceplate reads "Musical Filter".

A Direct Form II transposed serial biquad cascade, authored body-first and
shipped as a cartridge-based morph instrument. Not an E-mu clone. Not a
Rossum-kernel emulator. DF-II is the name. Musical Filter is what it does.

## Canonical
- Filter-work doctrine + shipping rules: `authoring/FILTER_WORK.md`
- Operating modes: `MODES.md`
- Shipping bodies + rubrics: `BODIES.md`
- Phoneme model: `PHONEMES.md`
- Cartridge wire format: `cartridge.schema.json`
- Per-session brief: `CLAUDE.md`

## Verify the workspace
    ./check

Runs doc-set sanity, cargo type-check across the workspace, shipping
pill loader smoke, vowel-formant label classification, hedz cross-
language parity, and P2K parity-null against `ref/canonical/`. The
live filesystem is the repo index — there is no precomputed snapshot.

## Shipping phoneme pills
- Flat layout: `cartridges/engine/<category>/<key>.json`
- Pill manifest: `cartridges/engine/manifest.json`
- Upstream source: `cartridges/engine/_source/`
- Rebuild from source: `python tools/bake_phoneme_pills.py`

## Entrypoints
- Rust workspace: `Cargo.toml`
- JUCE plugin repo (sibling tree): `trench-juce/plugin/`
- Authoring cockpit (sibling repo, out of shipping tree): `trench-forge/`

