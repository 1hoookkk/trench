---
name: RE Findings — Morph Designer Compiler
description: Proven RE facts about the E-mu Morpheus/P2K DSP compiler architecture
type: project
---

## Morph Designer Compiler — Proven Facts (2026-03-12)

**Source**: `FUN_1802c6590` — the Morph Designer UI-to-runtime compiler.

### Architecture

- **CASCADE, not parallel.** 6 biquad stages in series. Proven.
- **2 logical corners authored**, duplicated into **4 runtime banks**:
  - bank A (M0) at `+0x2c0`
  - bank B (M1) at `+0x2fc`
  - bank A copy at `+0x338`
  - bank B copy at `+0x374`
- The duplication is for the bilinear interpolation engine, not parallel signal paths.
- **Morph Designer has no Q-axis authoring.** Q blends between identical M0/M1 pairs. TRENCH's 4-corner model is a genuine superset.

### Stage Activation

- `type == 0` → skip (inactive)
- `type == 1, 2, 3` → compile this stage
- Loop is exactly 6 iterations

### Stage Count Rule

- Active stages 1–3: runtime count = actual count
- Active stages 4–5: runtime count **forced to 6**, remaining slots filled with bypass sentinel
- Active stages = 6: normal

### Passthrough Sentinel (`DAT_1806d7500`)

Raw `w0..w4`: `-8193, -1, -8193, -1, -8192`

After recombination (`c = w*4 + next_w`):
- `c0 = -32773, c1 = -1, c2 = -32773, c3 = -1, c4 = -8192`

Use as default stage initialization in all 6 slots.

### Per-Stage Byte Layout (6 bytes)

- byte 0: `type` (0=off, 1/2/3=topology branch)
- bytes 1, 3: feed SR base/scale tables → radius/emphasis driver (NOT frequency)
- bytes 2, 4: signed offsets `±(-64..+63)`, shifted by `sVar7` → position/spread (frequency-like)
- byte 5: padding
- `sVar7` = global shift on position bytes (NOT a gain offset)

### Type Branches

- `type == 2` and `type == 3` apply different algebraic scaling to compute 5 encoded shorts per endpoint (`local_58[0..4]` and `local_48[0..4]`)
- These branches are the unported math needed to close Talking Hedz

### Type 3 Formula (Proven)

Frequency compression for vocal/formant filters. From `FUN_1802c6590`:

```c
if (f > 0xDB && morph < 0) {
    f = ((f - 220) * (morph + 32) >> 5) + 220;
}
```

Operates in byte domain (pre-frequency-lookup). When morph is fully negative (-32), all frequencies above 0xDB compress to 220.

### AGC (Proven)

AGC verified exact against EmulatorX.dll binary (2026-02-28). All 16 table values confirmed. Algorithm confirmed from decompilation of `FUN_1802c04e0`. See `trench-core/CLAUDE.md` §AGC.

### Algebraic Projection Families

Analysis of P2K skin zero placement (2026-03-12):
- **Unit Circle**: 48% of stages — zeros on unit circle (|z|=1)
- **Near-Allpass**: 21% — zeros near poles (b≈a), minimal spectral shaping
- **Zero Inside Pole**: 22% — zeros inside pole radius
- **Stage 5**: Always carries unit-circle zeros (100%)
- **Morpheus cubes**: poles only (no numerator data). P2K skins: complete filters (val1/val2/val3).

### What's Still Unproven

- Talking Hedz specifically routes through this compiler (not proven, just assumed)
- Exact semantic names for the 4 non-type bytes
- Whether Stage 0's shelf is type 2 or type 3

### Contract Status

`contracts/cleanroom/handoffs/morph_designer/v2/CONTRACT.json`:
- `logical_corner_count = 2` ✓ proven
- `runtime_bank_count = 4` ✓ proven
- `stage_count_rule` ✓ proven
- `type-specific shaping branches` ✓ proven
- Type 2/3 math formulas: **ported** (trench-forge/src/historical.rs, all branches + tests)
