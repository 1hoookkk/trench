# TRENCH

TRENCH is a cartridge-based morph filter instrument. The repo has 5 canonical
docs; everything else is implementation detail.

## Canonical
- Math + contracts (law): `SPEC.md`
- Operating modes: `MODES.md`
- Doctrine (rules, bans, verification, escalation): `DOCTRINE.md`
- Shipping bodies + rubrics: `BODIES.md`
- Per-session brief: `CLAUDE.md`

## Verify the workspace
    ./check

Runs doc-set sanity, cargo type-check across the workspace, pyruntime imports.
The live filesystem is the repo index — there is no precomputed snapshot.

## Shipping phoneme pills
- Flat layout: `cartridges/engine/<category>/<key>.json`
- Pill manifest: `cartridges/engine/manifest.json`
- Upstream source: `cartridges/engine/_source/`
- Rebuild from source: `python tools/bake_phoneme_pills.py`

## Entrypoints
- Rust workspace: `Cargo.toml`
- Python runtime API: `pyruntime/api.py`
- JUCE plugin repo (sibling tree): `trench-juce/plugin/`

