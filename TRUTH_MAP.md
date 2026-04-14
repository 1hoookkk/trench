# TRENCH — Active Truth Map

One page. No nostalgia. If a path isn't listed under an active section,
it is either archive/non-canonical or queued for deletion.

## Shipping runtime

- **Plugin**: `trench-juce/plugin/` — **sibling git, not in this repo.**
  C++/JUCE 8. Windows standalone + VST3.
- **DSP source of truth**: `trench-juce/plugin/source/TrenchEngine.cpp`
  (referenced by `DOCTRINE.md`, `SPEC.md`). The 12-stage serial DF2T
  cascade and the bilinear cartridge interpolation live there.
- **Rust crates that serve the shipping runtime**: `trench-core` only.
  It exists to codify DSP math (`cascade.rs`), the cartridge schema
  (`cartridge.rs`), and heritage bakes (`hedz_rom.rs`) that the JUCE
  plugin can consume via FFI or a committed C header. If it doesn't
  serve the JUCE plugin, it doesn't belong here.

## Shipping cartridge format

- **Schema**: `cartridge.schema.json` at repo root.
- **Canonical format tag**: `"format": "compiled-v1"`.
- **Keyframe shape**: 4 keyframes (`M0_Q0`, `M100_Q0`, `M0_Q100`,
  `M100_Q100`) × 6 or 12 stages (12 = 6 active + 6 passthrough tail) ×
  `{c0, c1, c2, c3, c4}` direct-form DF2T coefficients + per-keyframe
  `boost`.
- **Rust loader**: `trench-core/src/cartridge.rs::Cartridge::from_json`
  — accepts both the keyframe format above and the array format
  (`corners.M0_Q0 = [[c0..c4], ...]`).
- **Python serializer**: `pyruntime/body.py::Body::to_compiled_json`.

## Verification command

    ./check

Runs, in order: canonical doc set check, cargo type-check the whole
workspace, `pyruntime` import check, null-target test, 41-entry
P2K parity null (hard fail at > -140 dB), Rust cascade parity vs 33
raw skins × 4 corners (hard fail at > -120 dB), and
`cartridge.schema.json` validation across `cartridges/**/*.json`.

## Authoring tools allowed to touch shipping

These write files the Rust crates or the JUCE plugin actually load.
Everything else is a research script.

| Tool | Writes what |
|---|---|
| `tools/extract_emu_filter_params.py` | `vault/_phonemes/heritage_designer_sections.json` — heritage integer grids (read-only data for forge targeting) |
| `tools/bake_hedz_const.py` | `trench-core/src/hedz_rom.rs`, `trench-core/src/hedz_golden.rs` — hardcoded Talking Hedz cartridge + golden impulse vector for the cross-language parity test |
| `pyruntime/forge_shipping.py` + `pyruntime/forge_shipping_finalists.py` | Shipping body candidates under `vault/_shipping_finalists/` |
| `tools/author_sonic_bank.py` | `vault/_shapes/` shape-bank regeneration |
| `tools/parity_null.py` | P2K parity gate consumed by `./check` |

## Archive / non-canonical paths

Kept for history or research; **not in the shipping chain**. None of
these should be imported, linked, or referenced by anything in the
active runtime or the verification command.

- `cartridges/heritage_quarantine/` — compiled biquad floats from
  trenchwork_clean. Known to diverge from the markdown-XML compile
  path (5-stage vs 6-stage runtime semantics per
  `docs/looperator/talking_hedz_x3_surfaces.md`). Do not treat as
  ground truth.
- `dev/notebooklm_sources/` — research narrative md files. Reference
  material for vocal/formant authoring, not runtime data.
- `docs/notebooklm/packs/` — imported NotebookLM bundle packs
  (ground-truth, ship-engine, frontier). Research context; the real
  ground truth is the integer grid in
  `vault/_phonemes/heritage_designer_sections.json`.
- `vault/_atlas/`, `vault/_cracked/`, `vault/_diagnostics/`,
  `vault/_plots/`, `vault/_profiles/`, `vault/_scorecards/`,
  `vault/_triage/`, `vault/_splice_candidates/`,
  `vault/_vocal_pack_2026-04-11/` — authoring scratch and diagnostics.
- `pyruntime/api.py` and everything `forge_*.py` that isn't
  `forge_shipping*` — FastAPI research cockpit. Serves the authoring
  UI, not the shipping runtime.
- `docs/archive/`, `docs/archive_symbol_inventory_2026-04-09.*` —
  snapshot-at-date inventories; historical only.

## Delete next

Claims with evidence. The purge order is safest-first.

| Path | Why it goes | Evidence |
|---|---|---|
| `trench-live/` (whole crate) | Orphan Rust plugin harness. Zero doc mentions. The shipping plugin is `trench-juce/plugin/` (JUCE, C++, sibling repo). Only references in this repo are its own files + workspace `Cargo.toml:2`. | `grep -rln "trench-live" /home/user/trench/` returns three self-references and nothing else. |
| `trench-codec/` | 119-line crate, referenced only by its own `Cargo.toml` and archive inventories. No doc owns it, no other crate depends on it. Verify before deletion with a workspace-wide import check. | Same grep pattern — only self-references plus archive dumps. |
| `xtask/` | 21 lines in `xtask/src/main.rs`. Cargo convention for workspace tooling; check whether anything in `./check` or CI invokes it before deletion. | `grep -rln xtask` shows `BODIES.md` and an archive handoff — probably historical. |
| `trench-core/src/cluster.rs`, `trench-core/src/motor.rs`, `trench-core/src/role.rs`, `trench-core/src/engine.rs` | 800+ lines of Rust that may or may not be consumed by the JUCE plugin. `engine.rs` is exported via `lib.rs:12` but only its own tests use it inside the workspace. `Cluster` / `Motor` / `Role` similarly — audit what the JUCE plugin actually FFI's into before keeping. | Workspace grep + the sibling-repo constraint. |
| `trench-live/src/editor.rs` | **Already deleted** (Phase 3 of the current surgical slice, pending commit). | `ls trench-live/src/` shows only `lib.rs`. |

## Notes

- The **sibling-repo invariant** makes everything in this repo's
  `trench-live/`, `trench-codec/`, `xtask/`, and most of `trench-core/`
  suspect until someone audits what `trench-juce/plugin/` actually
  FFI's into. That audit has to happen before the deletion list above
  can be executed safely.
- `trench-core/src/cartridge.rs` + `cascade.rs` + `agc.rs` +
  `hedz_rom.rs` + `hedz_golden.rs` are the only Rust files I can
  defend as shipping-relevant without seeing the sibling repo.
