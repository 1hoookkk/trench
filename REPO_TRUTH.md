# Repo Truth (Canonical Index)

This is the single jump-table for **what is canonical** in the repo.
If paths/entrypoints change, regenerate with: `python tools/update_repo_truth.py`.

## Start Here
- Product/DSP law: `SPEC.md`
- Shipping checklist + release gate: `SHIPPING.md`
- Agent posture + bans: `AGENTS.md`, `CLAUDE.md`
- DSP subtree rules: `trench-core/CLAUDE.md`
- Forge/compiler boundary rules: `pyruntime/CLAUDE.md`

## Entrypoints
- rust_workspace: `Cargo.toml`
- python_api: `pyruntime/api.py`
- juce_plugin_repo_dir: `trench-juce/plugin/`

## Verification (Actual Commands)
- Rust:
  - `cargo test --workspace`
  - `cargo clippy --workspace --all-targets -- -D warnings`
- Python:
  - `python -m pytest tests/test_api.py -v`
- Shape bank regen:
  - `python tools/author_sonic_bank.py --pairs --out-root vault/_shapes`

## Shape Bank (Current)
- counts_source: `vault/_shapes/manifest_unified.json`
- total_bodies: 370
- static_bodies: 50
- pair_bodies: 320
- categories: 15
- bank docs: `vault/_shapes/README.md`

## P2K Phonemes (Current)
- clusters: 32
- unmapped: 0
- heritage_fills: 7
- inventory docs: `vault/_phonemes/README.md`

## Non-canonical (Don’t Treat As Truth)
- `_GRAVEYARD/`
- `dev/tmp/`
- `target/`
- `target-codex-vizia/`

## External Workspaces (Non-canonical)
These directories are outside this repo. They may contain related work, but they are not the shipping truth for `C:\Users\hooki\Trench`.
- `trench_re_vault`: `C:\Users\hooki\trench_re_vault` (exists=True)
  - role: DIRTY airlock (reverse-engineering / capture / extraction).
  - policy: Read-only from CLEAN/shipping code. Promote only sanitized contracts into this repo.
- `trenchwork_clean`: `C:\Users\hooki\trenchwork_clean` (exists=True)
  - role: Separate CLEAN workspace snapshot (factory/tooling).
  - policy: Treat as non-canonical here; copy-in only specific files with explicit intent + provenance.

## Existence Checks (Staleness Alarm)
- `Cargo.toml`: True
- `SHIPPING.md`: True
- `SPEC.md`: True
- `pyruntime/api.py`: True
- `trench-juce/plugin`: True
- `vault/_phonemes/p2k_phoneme_inventory_v2.json`: True
- `vault/_shapes/manifest.json`: True
