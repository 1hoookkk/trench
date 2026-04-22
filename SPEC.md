# TRENCH spec

This file is law. It defines the math and the on-disk contract. Nothing else
lives here. Product, doctrine, modes, and shipping briefs live in their own
files.

## 1. DSP engine

- **Topology:** 12-stage serial DF2T biquad cascade (6 active stages + 6
  bypassed/passthrough SOS). Never parallel.
- **Coefficient format:** kernel-form `[c0, c1, c2, c3, c4]` interpolation
  only. Never interpolate raw biquad coefficients or frequencies.
- **Interpolation (compiled-v1 path):** 4-corner bilinear. Interpolate Q
  first, then morph. This is a property of the compiled-v1 cartridge format,
  not a permanent engine law. See §2 (Cartridge format) and
  `docs/architecture/zplane_truth.md` §"The 4-corner → 2-corner transition".
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

### Bilinear interpolation — compiled-v1 path (Q first, then morph)

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

The following are frozen regardless of cartridge format:

- 6 active serial stages (12 total, 6 passthrough)
- kernel-form coefficient interpolation only — never raw biquad or frequency
  interpolation

The compiled-v1 path additionally uses 4-corner bilinear interpolation (§2).
This is a format-layer property, not a topology law; a future format can change
the interpolation shape without violating these invariants.

No authoring model or cartridge format change grants permission to alter cascade
topology or DF2T math. The cascade is frozen at `trench-core/src/cascade.rs` and
`trench-juce/plugin/source/TrenchEngine.cpp`.

## 4. Cube Authoring Path (`trench.authoring_path.cube.v1`)

Cube is a macro authoring surface above heritage truth. It does not create a new
canonical truth layer.

Required fields:

- `schema`
- `id`
- `name`
- `provenance`
- `exactness`
- `control_mode_default`
- `axes`
- `legacy_behavior`
- `corners`

Corner labels are fixed:

- `c000`
- `c100`
- `c010`
- `c110`
- `c001`
- `c101`
- `c011`
- `c111`

Authoring files must not store runtime coefficients, compiled corner packs,
corner caches, or opaque blobs.

## 5. Cube Compiled Surface (`trench.compiled.cube_surface.v1`)

Cube lowers into a separate compiled runtime surface while preserving backward compatibility with `compiled-v1`.

Required fields:

- `schema`
- `derived_from`
- `compiler`
- `representation`
- `exactness`
- `control_mode`
- `corner_resolution`
- `corners`

`representation` is `engine_ready_coeff_packs`.

Cube runtime interpolation is trilinear in `(x, y, z)`. `compiled-v1` bilinear
interpolation remains intact and unchanged.
