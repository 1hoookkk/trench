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

## 3. Runtime invariants

The runtime stays:
- 6 active serial stages
- 4 corners
- kernel interpolation only

Friendlier authoring models do not grant permission to change runtime
topology, interpolation order, or cartridge format. The DSP cascade is
frozen at `trench-core/src/cascade.rs` and
`trench-juce/plugin/source/TrenchEngine.cpp`.
