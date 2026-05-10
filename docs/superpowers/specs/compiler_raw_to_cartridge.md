# Compiler spec - raw surface -> cartridge

**Status:** implemented at `authoring/compilers/compile_raw.py`; gated by
`runtime/trench-core/tests/compile_raw_roundtrip.rs` in `./check`.
**Authority:** `SPEC.md` section 2, `cartridge.schema.json`,
`runtime/trench-core/src/cartridge.rs::Cartridge::from_json`, `BODIES.md`.
**Scope:** in-tree Python compiler that reads an internal raw-stage-v1 authoring
surface and writes a `compiled-v1` cartridge JSON that loads cleanly into the
Rust runtime.

## 1. Purpose

`compile_raw.py` is a mechanical bridge from an internal authoring surface to
the shipping cartridge format.

It is for authored bodies, not for shipping verbatim third-party seed data.
P2K / EX3 extractions may inform the raw-stage surface, but the shipping gate
for public bodies is `BODIES.md` rubric compliance, not null parity with E-mu
references.

## 2. Input

The compiler accepts a JSON document with:

- `format: "raw-stage-v1"`
- `name`
- optional global `boost`
- exactly one of:
  - `frames`
  - `corners`
  - `keyframes`

### Native authoring surface

Preferred input is the native 2-frame surface:

- `frames.M0`
- `frames.M100`

Each frame may be written either as:

- an object with `stages` and optional per-frame `boost`
- or a direct stage-list shorthand

`compile_raw.py` duplicates `M0 -> M0_Q0/M0_Q100` and
`M100 -> M100_Q0/M100_Q100` at compile time. Q is not authored in this
surface.

### Legacy authored runtime surface

Legacy 4-corner input is still accepted for compatibility.

Required corners / keyframes for that path:

- `M0_Q0`
- `M0_Q100`
- `M100_Q0`
- `M100_Q100`

Authoring sample rate is frozen to `39062.5`. If `authoring_sample_rate_hz`
or `sampleRate` is present it must equal `39062.5`. If both are present they
must match.

Within `corners`, each corner may be written either as:

- an object with `stages` and optional per-corner `boost`
- or a direct stage-list shorthand

### Stage kinds

Each stage is authored in musical units and compiled deterministically:

- `passthrough`
- `allpole` / `resonator`
  - `pole_freq_hz`
  - `radius`
  - optional `stage_gain` (defaults to `1.0`)
- `zero_forced`
  - same as `allpole`, but with a unit-circle zero at the pole frequency
- `zero_forced_offset`
  - same as `allpole`, plus `offset_semitones`
  - the derived zero frequency is clamped into the runtime-safe range
- `explicit_zero`
  - same as `allpole`, plus `zero_freq_hz` and `zero_radius`

Maximum active stages per corner: `6`.

### Example

```json
{
  "format": "raw-stage-v1",
  "name": "aluminum_siding_seed",
  "authoring_sample_rate_hz": 39062.5,
  "boost": 4.0,
  "corners": {
    "M0_Q0": [
      {
        "kind": "explicit_zero",
        "pole_freq_hz": 11000.0,
        "radius": 0.978,
        "stage_gain": 0.78,
        "zero_freq_hz": 1000.0,
        "zero_radius": 1.0
      }
    ],
    "M0_Q100": [],
    "M100_Q0": [],
    "M100_Q100": []
  }
}
```

## 3. Output

A `compiled-v1` cartridge that passes `Cartridge::from_json`:

- `format: "compiled-v1"`
- `name`
- `sampleRate: 39062.5`
- `provenance: "compile_raw"`
- `keyframes`: 4 entries, labels `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`
- each keyframe carries `morph`, `q`, `boost`, and `stages`
- each keyframe emits `12` stages:
  - first `6` are active or passthrough-padded
  - trailing `6` are mandatory passthrough `[1,0,0,0,0]`

## 4. Contract

MUST:

- Accept only `raw-stage-v1`.
- Fail closed on missing corners, unknown stage kinds, non-finite values, or
  invalid frequency / radius ranges.
- Encode to kernel-form coefficients `[c0, c1, c2, c3, c4]`.
- Emit all 4 corners and 12 stages per keyframe.
- Preserve the frozen runtime invariants in `SPEC.md`.

MUST NOT:

- Repackage raw P2K / EX3 coefficient dumps as a shipping path by default.
- Use RBJ cookbook derivation.
- Change runtime topology, interpolation order, or cartridge format.
- Treat exact null parity with E-mu references as the shipping success metric.

## 5. Algorithm

For each corner:

1. Read up to 6 authored stages.
2. Translate stage law into raw resonator parameters:
   - pole frequency -> `a1 = -2r cos(theta)`
   - stage gain -> `val1 = stage_gain - 1`
   - zero law -> `val2`, `val3`
3. Cast the raw parameters through the f32 encode path.
4. Emit kernel coefficients:
   - `c0 = 1 + val1`
   - `c1 = a1 + val2`
   - `c2 = r^2 - val3`
   - `c3 = a1`
   - `c4 = r^2`
5. Pad active stages to 6 with passthrough.
6. Pad the tail to 12 with passthrough.
7. Emit compiled-v1 JSON.

## 6. Test strategy

Tests live in `runtime/trench-core/tests/compile_raw_roundtrip.rs` and run under
`./check`. The Rust test shells out to `python authoring/compilers/compile_raw.py`, captures
stdout, and feeds it through `Cartridge::from_json`, so any drift between the
Python emitter and the Rust loader fails the gate.

Shipped:

1. Round-trip loader acceptance for a valid raw-stage surface.
2. Round-trip loader acceptance for the direct `corners.<label> = [stage...]`
   shorthand plus `authoring_sample_rate_hz`.
3. Unknown stage-kind rejection.
4. Missing-corner rejection.

Deferred:

1. Body-rubric gates driven by `BODIES.md`.
2. "Too close to seed" tripwires against third-party references.
3. Authoring-room tooling layered on top of raw-stage-v1.

## 7. Host location

**Python `authoring/compilers/compile_raw.py`**, stdlib-only. Matches the existing in-tree
compiler pattern used by `authoring/compilers/compile_grid.py`: fast iteration in Python,
authoritative load validation in `trench-core`.

## 8. Non-goals for v1

- Legal clearance policy beyond recording the engineering boundary above.
- Automatic seed selection from the factory corpus.
- Creative optimization, scoring, or body acceptance.
- Public claims of EX3 parity.

## 9. MorphLP / LPX bridge blueprint (implementation-ready)

This subsection defines the concrete bridge from authoring controls into
`compiled-v1` for the MorphLP and MorphLPX compiler pathways.

### Inputs

- `morph_norm` in `[0,1]`.
- `morph_offset` (additive term from modulation path).
- four explicit lane bytes from authored section bytes (`+0x16,+0x17,+0x19,+0x1a`).
- two pitch bytes (`+0x18,+0x1b`) for the internal signed deltas.
- `sample_rate_slot` in `[0..3]`.
- `table_kind` in `{morphlp, morphlpx}`.

### Deterministic mapping

1. Compute row:
   `idx = clamp(floor((morph_norm + morph_offset) * 16.0), 0, 15)`.
2. Load table words:
   - MorphLP: row from `DAT_1806d73c0` (`stride=0xC`, 6 words).
   - MorphLPX: row from `DAT_1806d7480` (`stride=0x6`, 3 words) then duplicate
     to `[lane_a_0..2, lane_b_0..2]`.
3. Resolve pitch deltas exactly as decompiled:
   `pitch_a = clamp((pitch_byte_a - 0x40) * 0x80 + morph_delta_q * 0x100, -0x2000, 0x1fff)`
   and same for `pitch_b`.
4. Resolve the four pole seeds:
   `u13/u14/u11/u12 = slope[sr] * byte + base[sr]`, where
   `base=DAT_1806d73a0`, `slope=DAT_1806d73b0`.
5. Resolve the companion words:
   `r13/r14/r11/r12 = (u >> 1) + 0x6400`.
6. Build the three proved stage rows:
   - Stage 0 lane A/B from `u13/u14`, `pitch_a`, and the loaded zero words.
   - Stage 1 lane A/B from `u11/u12`, `pitch_b`, and stage-0 pole terms.
   - Stage 2 lane A/B from sentinels plus stage-1 pole terms.
7. Emit `active_stage_count = 3`, with remaining stage slots passthrough.
8. Convert stage words through the existing kernel conversion path and write
   `compiled-v1` corners.

### Shared vs corner-specific

- Shared: row index, row schema, seed tables, stage wiring, stage count.
- Lane-specific: the four byte inputs, two pitch bytes, and resulting pole/r words.
- LP only: row carries distinct M0/M100 zero triplets.
- LPX only: row triplet is morph-indexed but corner-invariant.

Blocked:
- exact binding from the two internal lanes to named `compiled-v1` corners.
  Keep these inputs explicit in the offline generator until that join is proved.

## 10. Offline preset generator: next actions

1. Add canonical constants module for:
   `DAT_1806d73a0`, `DAT_1806d73b0`, `DAT_1806d73c0`, `DAT_1806d7480`,
   and `MORPH_INDEX_SCALE=16.0`.
2. Implement `select_zero_row(table_kind, morph_norm, morph_offset)` with strict
   floor+clamp semantics and explicit LP vs LPX stride handling.
3. Implement `compile_lattice_3stage(...)` mirroring `FUN_1802c59b0` stage wiring
   and hard `active_stage_count=3`.
4. Add byte-level parity tests:
   - LPX index `0`, `7`, `15` rows produce exact expected words.
   - LP row pointer math uses `idx*0xC`, LPX uses `idx*0x6`.
   - Overflow/underflow morph values clamp to row `15`/`0`.
5. Add structure tests on emitted cartridge:
   exactly 3 active stages in this pathway and passthrough tail preserved.
6. Add offline preset emitter mode:
   given authored morph, morph offset, four lane bytes, two pitch bytes, and
   table kind, emit a `compiled-v1` cartridge snapshot for each requested
   corner projection.
7. Gate with regression fixtures built from proved rows, not inferred DSP math.

Blocked:
- exact UI ownership of the additive `morph_offset` field is unresolved; keep it
  as an explicit input parameter until that join is proved.

## 11. Unit-moving workflow (no parity gate)

For rapid product output, generate style-ready bodies directly:

```powershell
python tools/generate_emu_style_bank.py --skin 0 --intensity 1.0
```

This emits five audition-ready `compiled-v1` cartridges in:
`cartridges/engine/generated/emu_style_bank/`

Recipe set:
- `emu_soft_lp`
- `emu_talk_lp`
- `emu_vocal_push`
- `modern_hard_fold`
- `modern_air_lpx`

Use this loop to move quickly:
1. Generate bank for target skin/intensity.
2. Audition in plugin.
3. Promote winners into your shipping set.

### Bulk pack + shortlist commands

Generate a full production pack (all skins, multi-intensity, with manifest):

```powershell
python tools/generate_emu_style_bank.py --all-skins --intensities 0.85,1.0,1.2 --out-dir cartridges/engine/generated/emu_style_pack_v1 --write-manifest
```

Build a release shortlist (balanced subset, 24 files):

```powershell
python tools/build_release_shortlist.py --manifest cartridges/engine/generated/emu_style_pack_v1/manifest.json --source-dir cartridges/engine/generated/emu_style_pack_v1 --out-dir cartridges/engine/generated/emu_style_pack_v1_shortlist --target-intensity 1.0 --per-skin 6 --title "E-mu Style Pack v1 Shortlist"
```

Promote shortlist to release-candidate package folder:

```powershell
python tools/promote_release_shortlist.py --shortlist-manifest cartridges/engine/generated/emu_style_pack_v1_shortlist/manifest.shortlist.json --shortlist-dir cartridges/engine/generated/emu_style_pack_v1_shortlist --release-dir cartridges/engine/release_candidates/emu_style_pack_v1_rc1 --release-name emu_style_pack_v1_rc1
```

Single-command pipeline (generate -> shortlist -> promote):

```powershell
python tools/run_emu_pack_pipeline.py --version v2 --intensities 0.9,1.0,1.15 --target-intensity 1.0 --per-skin 6
```
