# TRENCH

TRENCH is a cartridge-based morph filter instrument. This repo is “vibe-coder friendly” by design: there is a single law doc, and a single live index you can always come back to.

## Start here (canonical)
- Law (don’t drift): `SPEC.md`
- Live repo index (paths, entrypoints, commands, inventories): `REPO_TRUTH.md`
- Shipping checklist + release gate: `SHIPPING.md`

## One habit that keeps you sane
When anything feels confusing, do this first:

1) Open `REPO_TRUTH.md`
2) If something moved, refresh it: `python tools/update_repo_truth.py`

## Shape bank (what we’re building right now)
- Shape bank docs: `vault/_shapes/README.md`
- P2K phoneme inventory docs: `vault/_phonemes/README.md`
- Regenerate the full bank: `python tools/author_sonic_bank.py --pairs --out-root vault/_shapes`

## Entrypoints
- Rust workspace: `Cargo.toml`
- Python runtime API: `pyruntime/api.py`
- JUCE plugin repo (sibling tree): `trench-juce/plugin/`

