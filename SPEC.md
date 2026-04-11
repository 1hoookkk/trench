## 1. WHAT TRENCH IS

TRENCH is a cartridge-based spectral instrument for authored filter bodies. It is not an emulation, not a utility EQ, and not an open-ended modular filter lab. The product goal is not "complex filters." Complex filters are easy. The goal is filters with intent: bodies whose motion, stress points, and invariants are audible, memorable, and useful on real material.

The shipping playback surface remains minimal:
- `MORPH`
- `Q`
- `TYPE`

The native JUCE authoring surface is allowed, but it is not a generic stage editor. It exists to author and inspect body intent on an evidence surface. The primary read is the response field and the body's audible behavior. Stages remain internal runtime structure and only appear in surgical/debug contexts.

TRENCH ships a small, curated body set. Every shipping body must win on real source material by delivering:
- stable behavior across the full 4-corner surface
- a clear invariant that survives motion
- a recognizable danger zone or betrayal point
- musically legible identity on bass, vocals, drums, or full material

## 2. PRODUCT DOCTRINE

### 2.1 Intent First

Every body definition begins with four things:
- `Identity`: what the body is trying to do
- `Invariant`: what must never disappear
- `Danger Zone`: where the body becomes special or violent
- `Motion Law`: what `MORPH`, `Q`, and optional modulation mean

If a body cannot be described in those terms, it is not authored enough yet.

### 2.2 Spectrum as Evidence

The main authoring surface is a 2D evidence field:
- response trace
- poles
- zeros
- corner states

The spectrum proves the body. It is not decorative UI.

### 2.3 Stages Are Runtime, Not UX

The engine remains a serial cascade and the cartridge format stays fixed. But the user should not be forced to think in `S1..S6` just to shape tone. Stages are:
- compile targets
- serialization structure
- surgical/debug view

They are not the primary authoring abstraction.

### 2.4 No Empty Complexity

More peaks, more stages, and more weirdness do not count as progress. A stronger body has:
- clearer intent
- stronger invariant
- better motion
- better failure behavior

Do not optimize for density when the actual need is authored identity.

## 3. DSP ENGINE SPEC

*   **Topology:** 12-stage serial DF2T biquad cascade (6 active stages + 6 bypassed/passthrough SOS). Never parallel.
*   **Coefficient Format:** Kernel-form `[c0, c1, c2, c3, c4]` interpolation only. Never interpolate raw biquad coefficients or frequencies.
*   **Interpolation Order:** 4-corner bilinear interpolation. You must interpolate Q first, then morph.
*   **Control Blocks:** 32-sample control blocks with per-sample coefficient ramping.
*   **Sample Rates:** 39062.5 Hz (Authoring/Forge domain), 44100.0 Hz (Plugin runtime domain).
*   **Minifloat Format:** 4-bit exponent (bias 15), 12-bit mantissa.

**Exact DF2T Per-Sample Math:**
For each sample `n`, update the coefficients via their ramp deltas:
`c_i[n] = c_i[n-1] + delta_c_i` (for `i` in 0..4)
DF2T Execution:
`y[n] = (c0 * x[n]) + w1[n-1]`
`w1[n] = (c1 * x[n]) - (c3 * y[n]) + w2[n-1]`
`w2[n] = (c2 * x[n]) - (c4 * y[n])`

**Minifloat Decode Formula:**
`decoded_value = ldexp((mantissa | 0x1000) * scale, exp - 15)`
*Sentinels:* `0x0000` = `0.0`, `0xDFFF` = passthrough gain, `0xFFFF` = maximum constant (1.0).

**Bilinear Interpolation Formula (Q first, then morph):**
`q_m0 = C_M0_Q0 + (C_M0_Q100 - C_M0_Q0) * Q_frac`
`q_m1 = C_M100_Q0 + (C_M100_Q100 - C_M100_Q0) * Q_frac`
`interpolated_c = q_m0 + (q_m1 - q_m0) * morph_frac`

## 4. CARTRIDGE FORMAT

The plugin consumes immutable JSON cartridges in the `compiled-v1` format.
*   **4 corners:** `M0_Q0`, `M100_Q0`, `M0_Q100`, `M100_Q100`
*   **6 active stages** per corner
*   **5 coefficients** (`c0, c1, c2, c3, c4`) per stage

**JSON Schema:**
```json
{
  "version": "compiled-v1",
  "name": "Body Name",
  "corners": {
    "M0_Q0": [
      [c0, c1, c2, c3, c4],
      [c0, c1, c2, c3, c4],
      [c0, c1, c2, c3, c4],
      [c0, c1, c2, c3, c4],
      [c0, c1, c2, c3, c4],
      [c0, c1, c2, c3, c4]
    ],
    "M100_Q0": [ /* 6 stages of 5 floats */ ],
    "M0_Q100": [ /* 6 stages of 5 floats */ ],
    "M100_Q100": [ /* 6 stages of 5 floats */ ]
  }
}
```

*   `version`: Must exactly match `"compiled-v1"`.
*   `name`: The display name of the body/cartridge.
*   `corners`: The 4 absolute endpoints of the bilinear interpolation surface. Each corner array must contain exactly 6 stage arrays, and each stage array must contain exactly 5 floats.

## 5. AUTHORING SPEC

### 5.1 Primary Authoring Model

Bodies are authored against intent, then proven on the spectrum.

The default workflow is:
1. Define `Identity`, `Invariant`, `Danger Zone`, and `Motion Law`.
2. Place and edit poles/zeros directly on the response field.
3. Verify all 4 corners as one body, not as isolated snapshots.
4. Use stages only when surgical intervention is required.

### 5.2 Body Constraints

Every shipping body must have:
- one clear reason to exist
- one invariant that survives the whole morph
- one or more authored danger regions
- motion that matters between corners

If the body only sounds good at one frozen point, it is not ready.

### 5.3 Sonic Tables

The sonic tables are not decorative reference material. They provide authored targets and guide rails:
- vowel formants for intelligible cavity behavior
- nasal anti-resonances for zero placement
- instrument/body landmarks for recognizable physicality
- high-frequency and boundary landmarks for authored damage

The tables must be involved when defining:
- anchor zones
- notch territories
- vocal/formant behavior
- physical source analogies

### 5.4 Runtime Structure

The authoring layer may think in resonant objects, anchors, peaks, notches, and trajectories. The runtime remains:
- 6 active serial stages
- 4 corners
- kernel interpolation only

Do not confuse a friendlier authoring model with permission to change runtime topology.

## 6. SHIPPING BODY BRIEF

The launch set is the four named bodies in [SHIPPING.md](C:\Users\hooki\trench-juce\SHIPPING.md):
- `Speaker Knockerz`
- `Aluminum Siding`
- `Small Talk`
- `Cul-De-Sac`

These names are the public brief. Match the filter to the name, not the other way around.

For shipping generation and editing:
- keep only the four names
- generate candidates that serve those names
- reject bodies that drift from the written invariant or motion law

## 7. WHAT IS BANNED

*   No RBJ cookbook in the forge or shipping generation path.
*   No heritage branding in the shipping UI or output.
*   No verbatim P2K extractions as shipping presets. P2K and Morpheus are calibration/reference only.
*   No pole sanitization unless a failing test proves it is required.
*   No baking gain into `c4`.
*   No generic stage-first UX as the default authoring surface.
*   No complexity for its own sake.
*   Do not derive behavior from corpus statistics when a proven formula exists.
*   Do not split work across repos.
*   Do not add new dependencies without asking first.

## 8. SCOPE

**Nothing is off-limits. Measure it. Build it. Ask me if it sounds right.**

## 9. REPO STRUCTURE

*   `plugin/`: JUCE plugin, native UI, embedded cartridge assets, tests, and shipping binary targets.
*   `plugin/source/`: Editor, processor, engine, authoring surface, and runtime integration.
*   `plugin/assets/cartridges/`: Embedded compiled-v1 shipping bodies used by `TRENCH.vst3`.
*   `docs/`: Product briefs, workflow specs, architecture notes, and shipping body doctrine.
*   `build/`: Local JUCE/CMake build output.

## 10. WORKING RULES

*   **Ownership:** You own implementation, debugging, code structure, and verification. Product identity, taste, and musical judgment stay human-owned.
*   **Verification:** If the audible path changed, say exactly what changed and why. Do not claim success without running it.
*   **Measure Then Constrain:** Measure -> characterize -> validate -> only then promote to policy.
*   **No Compensation Layers:** Fix the broken domain constraint, not the symptom.
*   **Intent Beats Density:** Prefer a body with one undeniable point over a busy body with no point.
