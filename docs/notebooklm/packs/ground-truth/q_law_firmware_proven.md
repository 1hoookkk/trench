---
name: q_law_firmware_proven
description: CORRECTED — CPhantomBandLimitedTable is SINC interpolation kernels for oscillator anti-aliasing, NOT Q damping. Q law source still unidentified.
type: project
---

## CPhantomBandLimitedTable — SINC Interpolation Kernels (CORRECTED 2026-03-14)

### Previous Belief (WRONG)
Believed to be 7 Q damping tables for morph filter resonance control. **Disproven** by tracing the full consumer chain through Ghidra.

### Actual Function: Polyphase SINC Resampler

CPhantomBandLimitedTable provides **band-limited interpolation kernels** for the oscillator's sample playback engine. It is consumed by FUN_1802b96c0 (CPhantomRTOscillator vtable+0x130), the high-quality interpolation method.

**Algorithm:**
1. **Band selection**: Compare phase increment (pitch ratio) against 7 frequency thresholds at table+0x50
2. **Kernel lookup**: Select the appropriate band's precomputed SINC kernel
3. **8-tap FIR**: Apply 4+4 polyphase filter coefficients to 8 neighboring samples
4. **Anti-aliasing**: Higher pitch → higher band → tighter low-pass → less aliasing

### 7 Source Tables (rdata)

| Band | Address      | T[0]     | Semitone offset | Anti-aliasing |
|------|-------------|----------|-----------------|---------------|
| 0    | 0x18065bb70 | 0.986271 | +3              | Minimal       |
| 1    | 0x18065cb70 | 0.798618 | +6              | Light         |
| 2    | 0x18065db70 | 0.608257 | +9              | Moderate      |
| 3    | 0x18065eb70 | 0.529848 | +12             | Medium        |
| 4    | 0x18065fb70 | 0.427744 | +15             | Strong        |
| 5    | 0x180660b70 | 0.331154 | +18             | Heavy         |
| 6    | 0x180661b70 | 0.283777 | +21             | Maximum       |

Source tables: 512 entries × 8 bytes (float64). Stride: 0x1000 bytes.
Pointer table base: `0x1808899d0` (7 pointers).

### Constructor Output Structure (0xa8 = 168 bytes)
- 0x00: vftable (CPhantomBandLimitedTable::vftable at 0x1806df768, destructor only)
- 0x08..0x38: 7 pointers to precomputed float32 kernel banks (16384×4 entries each)
- 0x40: 0x4000 = 16384 phases (sub-sample precision)
- 0x44: 0x1000 = 4096 stride
- 0x48: 4 (taps per phase, ×2 for forward+complement = 8-tap FIR)
- 0x50..0x68: 7 float frequency thresholds (band boundaries)

### Constructor Constants (FUN_1802e7310)
- `DAT_1806ecb18` = **3.0** — semitone step between bands
- `DAT_1806ecb58` = **511.0** — source table max index
- `DAT_1806ecb60` = **1/16383** — phase-to-source-index scale
- `DAT_1806ecb68` = **1/12** — semitone-to-octave conversion

### Runtime Chain (PROVEN)
```
phantomGlobalTableFactory (FUN_1802b21a0)
  → constructs CPhantomBandLimitedTable → DAT_1808acb88

CPhantomRTOscillator::ctor (FUN_1802b62f0)
  → stores at osc+0x18 (param_1[3])

CPhantomRTOscillator::process (FUN_1802b75e0)
  → re-injects osc[3] = DAT_1808acb88 every frame
  → dispatches to vtable methods

vtable+0x08 (FUN_1802b78e0, forward play)
  → dispatches to vtable+0x130 when crossfade enabled

vtable+0x130 (FUN_1802b96c0, SINC interpolation)
  → reads osc[3] (BandLimitedTable pointer)
  → selects band by comparing phase increment to frequency thresholds
  → 8-tap polyphase FIR convolution using band kernel data
```

### Sibling Tables (also in oscillator)
- `DAT_1808acb90` → osc+0x08 = **CPhantomCoarseTable** (128-entry semitone→frequency, FUN_1802e7820)
- `DAT_1808acb98` → osc+0x10 = **CPhantomFineTable** (401-entry sub-semitone interpolation, FUN_1802e7990)
- Both consumed by **FUN_1802bdef0** (pitch lookup, NOT damping lookup) → result at osc+0xF8 = phase increment

### Oscillator Vtable (CPhantomRTOscillator::vftable at 0x1806d62a8)
| Offset | Function     | Role |
|--------|-------------|------|
| 0x00   | FUN_1802b6520 | destructor |
| 0x08   | FUN_1802b78e0 | forward play dispatch |
| 0x10   | FUN_1802b7a80 | backward play |
| 0x18   | FUN_1802b7c20 | release |
| 0x20   | FUN_1802b7ea0 | simple forward (linear interp) |
| 0x110  | FUN_1802b9050 | stereo linear interpolation loop |
| 0x130  | FUN_1802b96c0 | **SINC interpolation loop (reads BandLimitedTable)** |
| 0x150  | FUN_1802b7800 | delay/skip handler |
| 0x160  | FUN_1802bdab0 | loop boundary guard sample copy |

### What This Means for TRENCH
- The 7 source tables in rdata are oscillator anti-aliasing data, NOT morph filter Q data
- The Q law for the morph filter has NOT been located yet
- The morph filter Q behavior must be encoded elsewhere — possibly in the morph filter runtime itself, not in the oscillator
- The cross-validation result (k matching across sample rates) is still valid but describes SINC kernel coefficients, not Q damping

### Q Law — Still Open
The actual morph filter Q law (how the Q parameter maps to pole radius/damping in the biquad cascade) remains unidentified in firmware. Next steps:
- Search for the morph filter processor (separate from CPhantomRTOscillator)
- Look for CPhantomFilter or similar classes
- The Q law may be a simple mathematical formula rather than a lookup table
