# TRENCH spec

This file is law. It defines the math and the on-disk contract. Nothing else
lives here. Product, doctrine, modes, and shipping briefs live in their own
files.

## 1. DSP engine

- **Topology:** 12-stage serial DF2T biquad cascade (6 active stages + 6
  bypassed/passthrough SOS). Never parallel.
- **Coefficient format:** kernel-form `[c0, c1, c2, c3, c4]` interpolation
  only. Never interpolate raw biquad coefficients or frequencies.
- **Interpolation order:** 4-corner bilinear interpolation. Interpolate Q
  first, then morph.
- **Control blocks:** 32-sample control blocks with per-sample coefficient
  ramping.
- **Sample rates:** 39062.5 Hz (authoring/forge), 44100.0 Hz (plugin runtime).
- **Minifloat format:** 4-bit exponent (bias 15), 12-bit mantissa.

### Per-sample DF2T math

Per sample `n`, update each coefficient via its ramp delta:

    c_i[n] = c_i[n-1] + delta_c_i      (i in 0..4)

DF2T execution:

    y[n]  = (c0 * x[n]) + w1[n-1]
    w1[n] = (c1 * x[n]) - (c3 * y[n]) + w2[n-1]
    w2[n] = (c2 * x[n]) - (c4 * y[n])

### Minifloat decode

    decoded = ldexp((mantissa | 0x1000) * scale, exp - 15)

Sentinels: `0x0000` = 0.0; `0xDFFF` = passthrough gain; `0xFFFF` = max
constant 1.0.

### Bilinear interpolation (Q first, then morph)

    q_m0 = C_M0_Q0   + (C_M0_Q100   - C_M0_Q0)   * Q_frac
    q_m1 = C_M100_Q0 + (C_M100_Q100 - C_M100_Q0) * Q_frac
    coef = q_m0      + (q_m1        - q_m0)      * morph_frac

## 2. Cartridge format (`compiled-v1`)

The plugin consumes immutable JSON cartridges. Authoritative loader:
`trench-core/src/cartridge.rs::Cartridge::from_json`. Schema:
`cartridge.schema.json`.

```json
{
  "format": "compiled-v1",
  "name": "Body Name",
  "sampleRate": 44100,
  "keyframes": [
    {
      "label": "M0_Q0",
      "morph": 0.0, "q": 0.0,
      "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}, "... 12 stages ..."]
    },
    { "label": "M0_Q100",   "morph": 0.0, "q": 1.0, "stages": ["..."] },
    { "label": "M100_Q0",   "morph": 1.0, "q": 0.0, "stages": ["..."] },
    { "label": "M100_Q100", "morph": 1.0, "q": 1.0, "stages": ["..."] }
  ]
}
```

- `format` must equal `"compiled-v1"`.
- `keyframes` must contain all 4 corner labels: `M0_Q0`, `M0_Q100`,
  `M100_Q0`, `M100_Q100`.
- Each keyframe `stages` array contains 12 second-order sections (6 active
  + 6 passthrough). Each stage is an object with numeric `c0..c4`.
- Corners are absolute endpoints of the bilinear surface. No implicit
  midpoints, no per-stage exceptions.

### Provenance: 4-corner is a TRENCH optimization of E-mu's 2-frame Morph Designer

The 4-corner (`M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`) shape is NOT
E-mu's native Morph Designer format. E-mu stores 6 sections × 2 frames
(Lo Morph / Hi Morph) with a live Q/Gain parameter and a ±50% offset
wheel; biquad coefficients are derived at runtime from the interpolated
`(type, freq, Q)` tuple. TRENCH precomputes Q0 and Q100 snapshots at
each morph frame to avoid runtime coefficient computation — a valid
optimization on the hardware E-mu shipped this on, and free on modern
CPUs.

The long-term direction is a 2-corner runtime that matches E-mu's
native behavior exactly; see `docs/architecture/zplane_truth.md`.
Until that transition ships, the 4-corner shape defined above is
frozen.

## 3. Runtime invariants

The runtime stays:
- 6 active serial stages
- 4 corners
- kernel interpolation only

Friendlier authoring models do not grant permission to change runtime
topology, interpolation order, or cartridge format. The DSP cascade is
frozen at `trench-core/src/cascade.rs` and
`trench-juce/plugin/source/TrenchEngine.cpp`.

## 4. MorphLP / LPX Compiler Bridge (proved RE contract)

This section defines the proven authoring-to-kernel bridge for the MorphLP /
MorphLPX pathway. It is for compiler correctness, not UI design.

- **Entry functions:** `FUN_1802c5e40` (MorphLP), `FUN_1802c5f10` (MorphLPX),
  both call `FUN_1802c59b0`.
- **Morph indexing:** `index = floor((morph + offset) * 16.0)`, clamped to
  `[0,15]` using `DAT_18065a1ec`.
- **LP table:** `DAT_1806d73c0`, stride `0xC` (6 words): distinct
  `[M0_v1,M0_v2,M0_v3,M100_v1,M100_v2,M100_v3]`.
- **LPX table:** `DAT_1806d7480`, stride `0x6` (3 words): `[v1,v2,v3]` duplicated
  to both corners before the inner compiler.
- **Seed tables:** `DAT_1806d73a0` (base) and `DAT_1806d73b0` (slope), selected by
  sample-rate slot (`this + 0xc`).
- **Seed law:** `pole = slope[sr] * q_byte + base[sr]`; `r = (pole >> 1) + 0x6400`.
- **Lattice cross-coupling:** stage 0 uses table zeros; stage 1 zeros derive from
  stage 0 pole terms; stage 2 zeros derive from stage 1 pole terms.
- **Stage count:** `*(out + 0x420) = 3` hard-sets 3 active stages.

### Non-negotiable semantics for offline compilers

- Do not reinterpret LPX as a 6-word per-row table.
- Do not interpolate between table rows inside the compiler; row selection is
  discrete by the 16-step index.
- Do not emit >3 active stages for this pathway.
- Do not collapse LP and LPX into one row schema.
- Treat the byte inputs at `+0x16/+0x17/+0x19/+0x1a` as four proved lane inputs.
  Their exact binding to exported `compiled-v1` corner labels is still blocked.
