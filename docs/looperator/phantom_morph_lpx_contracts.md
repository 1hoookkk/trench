# Source
- `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_comparison_MorphLP_vs_LPX.md` (1493 bytes)
- `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_MorphLP_zero_table.md` (4237 bytes)
- `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_vtable_decompilation.md` (6979 bytes)

# What this data is
This is decompiled contract behavior for the core `CPhantomMorphLP` and `CPhantomMorphLPX` compiler classes inside the Emulator X engine. It describes how the software translates runtime morph positions (0 to 1.0) into discrete DSP zero-placement coefficients via lookup tables, frequency scaling laws, structural stages, and class hierarchies. 

# Why it's here
These contracts govern the literal math of how a "Morph" filter parameter translates into the unique "vocal" sweeps heard in the E-mu hardware. The CLEAN repo needs this table to algorithmically map a token's placement on the Looperator grid into precise DSP coefficients, accurately recreating the signature zero-dissolution behavior of hardware morph paths.

# Sanitized content

## vtable Decompilation & Class Hierarchy Structure
**Class context**: `CPhantomMorphLPX` inheriting from base `CPhantomObject`.
- **vtable Slot [0]**: `FUN_1802a38e0` (Scalar deleting destructor)
- **vtable Slot [1]**: `FUN_1802c5f10` (Compile method). This reads the 0.0-1.0 morph position, computes an index `floor((morph + offset) * 16.0)`, looks up values in a fixed zero table, and delegates compilation to the inner compiler function `FUN_1802c59b0`.
- **vtable Slot [2]**: `CharToken::getChar` (getChar method)
- The inner compile method `FUN_1802c59b0` actively controls a 3-stage active filter sequence (stages 4-6 are bypassed).
- Sample rate handling involves reading frequency biases (`DAT_1806d73a0`) and scaling factors (`DAT_1806d73b0`) based on a sample rate index [0..3].

## Comparison: CPhantomMorphLP vs. LPX
- **The Core Difference**: `CPhantomMorphLP` changes zero presence dynamically across the morph sweep, providing separate zero values for the M0 and M100 morph corners per index. `CPhantomMorphLPX` has static zeros, duplicating the same triplet across both M0 and M100.
- **MorphLP Architecture**: Uses a 16-entry × 12-byte table (`DAT_1806d73c0`) yielding 6 values per step. This allows heavy numerator character (deep LP zeros) at M0 to fade out into a near-passthrough all-pole response at M100.
- **MorphLPX Architecture**: Uses a 16-entry × 6-byte table (`DAT_1806d7480`) yielding 3 values. It provides a highly simplified, uniform LP sweep.
- **Inner Uniformity**: Both classes share the exact same morph-to-index translation logic (`16.0` scalar), index clamp range [0..15], and inner stage compiler `FUN_1802c59b0`.

## CPhantomMorphLP Zero-Placement Table
This data dictates the vocal signature. The zero pair depth decays toward a null/passthrough state at M100. Each table row (0-15) corresponds to exactly 1/16th of the user's morph range. At runtime, the DSP uses `val1` and `val2` (where `val1 == val3` representing symmetric conjugate pairs).

| Idx | M0_val1 | M0_val2 | M0_val3 | M100_val1 | M100_val2 | M100_val3 |
|---|---|---|---|---|---|---|
| 0 | 0x11ED | 0x94FC | 0x11ED | 0xDFFC | 0xFFFD | 0xDFFC |
| 1 | 0x11EE | 0x90FC | 0x11ED | 0xE2FC | 0xFDFD | 0xE2FC |
| 2 | 0x11ED | 0x8DFC | 0x11ED | 0xE5FC | 0xFC7D | 0xE5FC |
| 3 | 0x11EE | 0x8AFC | 0x11ED | 0xE8FC | 0xFAFD | 0xE8FC |
| 4 | 0x11ED | 0x87FC | 0x11ED | 0xEBFC | 0xF97D | 0xEBFC |
| 5 | 0x11ED | 0x83FC | 0x11ED | 0xEDFC | 0xF87D | 0xEDFC |
| 6 | 0x11EE | 0x80FC | 0x11ED | 0xF07D | 0xF77D | 0xF07D |
| 7 | 0x0FEE | 0x7DFC | 0x0FEF | 0xF17D | 0xF67D | 0xF17D |
| 8 | 0x0FEF | 0x7AFC | 0x0FEF | 0xF27D | 0xF47D | 0xF27D |
| 9 | 0x0FEE | 0x77FC | 0x0FEF | 0xF3FD | 0xF27D | 0xF3FD |
| 10 | 0x0FEE | 0x73FC | 0x0FEF | 0xF47D | 0xF17D | 0xF47D |
| 11 | 0x0DEE | 0x70FC | 0x0DEF | 0xF4FD | 0xF07D | 0xF4FD |
| 12 | 0x0DED | 0x6DFC | 0x0DEF | 0xF5FD | 0xEFFC | 0xF5FD |
| 13 | 0x0BEF | 0x6AFC | 0x0BEF | 0xF67D | 0xEEFC | 0xF67D |
| 14 | 0x0BEF | 0x62FC | 0x0BEF | 0xF6FD | 0xECFC | 0xF6FD |
| 15 | 0x0BEF | 0x62FC | 0x0BEF | 0xF6FD | 0xECFC | 0xF6FD |

# Integration notes
- When designing tokens leveraging the "Morpheus" hardware behavior, the compiler routine must resolve the UI grid morph parameter (0 to 1) into an explicit [0-15] index to look up this exact table.
- This represents numerator structure; it must be coupled with the denominator stage bytes pulled directly from the target token's `filter_cube_display_model.json` profile sections. 
- Apply standard minifloat translations for these u16 values (e.g. 0xFFFF saturation, 0xDFFF passthrough). Use the 16th entry (Idx=15) strictly for `morph_param >= 0.9375`.

# UNKNOWNs
- None identified. Sample rate bias and zero symmetry bounds provide enough data to directly port into the compiler routine for the Looperator.
