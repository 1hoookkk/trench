# Compiler spec - raw surface -> cartridge

**Status:** implemented at `tools/compile_raw.py`; gated by
`trench-core/tests/compile_raw_roundtrip.rs` in `./check`.
**Authority:** `SPEC.md` section 2, `cartridge.schema.json`,
`trench-core/src/cartridge.rs::Cartridge::from_json`, `BODIES.md`.
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
- either `corners` or `keyframes`

Required corners / keyframes:

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

Tests live in `trench-core/tests/compile_raw_roundtrip.rs` and run under
`./check`. The Rust test shells out to `python tools/compile_raw.py`, captures
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

**Python `tools/compile_raw.py`**, stdlib-only. Matches the existing in-tree
compiler pattern used by `tools/compile_grid.py`: fast iteration in Python,
authoritative load validation in `trench-core`.

## 8. Non-goals for v1

- Legal clearance policy beyond recording the engineering boundary above.
- Automatic seed selection from the factory corpus.
- Creative optimization, scoring, or body acceptance.
- Public claims of EX3 parity.
