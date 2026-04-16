---
name: Runtime Pipeline — Proven from Ghidra
description: Complete CPhantomRTFilter pipeline mapped from firmware decompilation (2026-03-14)
type: project
---

## CPhantomRTFilter — Complete Runtime Pipeline

**Class**: CPhantomRTFilter (base), CPhantomRTSSEFilter (SSE variant)

### Pipeline Order (FUN_1802c02b0)
1. Read dirty flags, dispatch vtable parameter updates
2. **FUN_1802c3d40** — Bilinear interpolation (morph × Q) on ushort banks
3. **FUN_1802c3600** — Minifloat decode (ushort → float biquad coefficients)
4. **FUN_1802c07a0** — DF2T single-stage render with per-sample ramping
5. Cascade dispatched at vtable level (not inside render function)

### Memory Layout
- `0x098`: Q parameter (smoothed float)
- `0x1e8`: Morph parameter base
- `0x200`: Morph parameter offset
- `0x2b0`: Q blend fraction (clamped, used by interpolation)
- `0x2c0`: Corner coefficient base (150 ushorts = 5 banks × 30)
- `0x420`: Stage count
- `0x428`: Current working coefficients (float, post-decode)
- `0x488`: Coefficient ramp deltas (float, per-sample)
- `0x620`: Biquad state buffer (12 floats mono, 24 stereo)

### Corner Bank Layout (at base 0x2c0)
- Offset 0x00 (0 ushorts): M0_Q0
- Offset 0x3C (30 ushorts): M0_Q100
- Offset 0x78 (60 ushorts): M100_Q0
- Offset 0xB4 (90 ushorts): M100_Q100
- Offset 0xF0 (120 ushorts): OUTPUT (5th bank, interpolation result)
- Each bank: 6 stages × 5 coefficients = 30 ushorts

### Bilinear Interpolation (FUN_1802c3d40)
- Pure linear lerp on raw ushort minifloats
- **No table lookup inside** — Q and morph are direct float fractions
- `q_m0 = C0 + (C1 - C0) * Q_frac`
- `q_m1 = C2 + (C3 - C2) * Q_frac`
- `output = q_m0 + (q_m1 - q_m0) * morph_frac`

### Minifloat Codec (FUN_1802c3600)
- 16-bit: 4-bit exponent (bias 15) + 12-bit mantissa
- Sentinels: 0x0000 → 0.0, 0xFFFF → max constant
- Decode: `ldexp((mantissa | 0x1000) * scale, exp - 15)`
- 5 values per stage, then kernel-form combining:
  - `bq[0] = decoded[0] * k + decoded[1]`
  - `bq[2] = decoded[2] * k + decoded[3]`
  - `bq[4] = state_gain * decoded[4]`

### DF2T Render (FUN_1802c07a0)
- Single-stage processor, kernel form [c0, c1, c2, c3, c4]
- DAT_18065a314 = 2.0 (kernel form identity constant)
- Per-sample: coeff[i] += delta[i] for all 5 coefficients

### Q Law — Proven Architecture
- Q is baked at **compile time** into 4 corner banks
- 512-entry damping table: `r_Q100 = 1 - T[index] * (1 - r_Q0)`
- Runtime Q knob = linear blend fraction between Q=0 and Q=100 corners
- Perceived nonlinearity from minifloat encoding, not runtime table
- Duplicated_pair policy: for historical Morph Designer filters, Q axis is degenerate
- TRENCH's 4-corner model with independent Q corners is a genuine superset

**Why:** Closes the Q law investigation. The forge can now implement Q correctly.
**How to apply:** Forge bodies with 4 unique corners already encode the correct Q behavior. Runtime just blends. No additional Q table needed in the plugin.
