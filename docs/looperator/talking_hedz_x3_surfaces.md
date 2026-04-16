# Source
- `C:\Users\hooki\trench_re_vault\datasets\talking_hedz_x3_surfaces\2026-03-12\README.md` (1140 bytes)
- `C:\Users\hooki\trench_re_vault\datasets\talking_hedz_stage0\2026-03-12\README.md` (713 bytes)
- `C:\Users\hooki\trench_re_vault\datasets\talking_hedz_bracketed_residual\2026-03-12\README.md` (1778 bytes)
- `C:\Users\hooki\trench_re_vault\datasets\talking_hedz\2026-03-12\README.md` (1177 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\talking_hedz_profile.json` (7483 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\talking_hedz_null_test_2026-03-14\reports\talking_hedz_null_test_report.md` (2186 bytes)

# What this data is
This is tracking and analysis data describing the `Talking Hedz` P2K morph skin, a marquee "vocal" formant filter. Specifically, this captures efforts to decode why direct raw hardware stage structures didn't match real-world audio measurements, identifying role-mapped filter boundaries, missing stages, and providing a baseline validation profile for null-testing reconstructions. 

# Why it's here
`Talking Hedz` represents the worst-case complexity for a morphing filter body—it exhibits the foreign-pole-alignment traits (now preserved as `project_zeros_architecture` in the CLEAN repo). The Looperator compiler must understand exactly how the original runtime simplified/swapped these stages at the corners to guarantee that the user's placement of a Talking Hedz token yields the actual expected vocal character, rather than the theoretically pure but sonically incorrect 6-stage raw format.

# Sanitized content

## X3 Surface Definitions
The `talking_hedz_x3_surfaces` capture identifies three principal measurement layers built around the `[M0_Q100, M100_Q100, M100_Q0]` coordinate grid:
1. **Raw Parameter Capture** (`do-it/tools/coeff_dump.json`): Direct captured biquad coefficients straight out of the virtual engine.
2. **Impulse Derived Surface** (`do-it/tools/talking_hedz_x3_extracted.json`): A full-cascade transfer surface computed independently from impulse responses.
3. **Reference WAV Audio** (`do-it/validation/*.wav`): Shared dry-input renders to test the DSP directly.

## Sub-Dataset Variations
Rather than one monolithic model, three complementary experiments were conducted to isolate DSP discrepancies:
- **`talking_hedz` (Base Dataset Check):** Compared the raw theoretical 6-stage filter (`P2k_013`) against the actual 5-stage captured runtime model (`TalkingHedz_Complete`). Proved the real engine throws out one stage and hard-converts the final captured stage into an explicit low-pass filter (`flag = 0`). 
- **`talking_hedz_stage0`:** Investigates the purpose of Stage 0 across the morph grid. Findings isolated that Stage 0 is active and dominant at the `M100` corners, while at `M0_Q100`, Stage 0 is bypassed entirely, acting as a non-linear shaping threshold across the sweep.
- **`talking_hedz_bracketed_residual`:** Examines whether structural stages explicitly "swap roles" during the morph. Found that `stage 5` acts as a necessary structural lowpass bracket for `M0_Q0`, and source stages 1..5 genuinely cross-map roles at asymmetrical morph corners.

## Profile Validation Matrix
Reconciling against `talking_hedz_profile.json`, the definitive performance landscape surveys a 15-node grid (Morph points `[0.0, 0.5, 1.0]`, Q points `[0.2, 0.4, 0.6, 0.8, 1.0]`).
- It has a high **Movement Index** (164.4) and wide centroid sweeps mapping to a **SmoothRamp** structural terrain with a 48dB dynamic range peak, cementing the "speaking" character across varying parameter sweeps.

## Null-Test Validation
The null-test analysis (`talking_hedz_null_test_report.md` on 2026-03-14) shows progress, but evaluates to a **FAIL** rating for perfect parity:
- **`M0_Q0`** fails with high-band residual (Null vs Ref RMS: -3.10 dB).
- **`M100_Q100`** best agreement but misses the low-band slightly (Null vs Ref RMS: -9.07 dB).
- **Conclusion**: Naively plugging in the `coeff_dump.json` states does not perfectly null the reference audio. Further gain-staging and inter-stage structural limits must be integrated into the final engine implementation.

# Integration notes
- When emitting TRENCH output for Talking Hedz tokens, enforce a 5-stage architecture instead of falling back to the generic 6-stage P2K skeleton.
- Coerce the highest stage into a low-pass structural bracket implicitly dependent on the Morph assignment.
- Expect partial acoustic dissonance if directly driving generic states—the "foreign-pole alignment" zero rule must override standard DSP topology mapping to hit character accurately.

# UNKNOWNs
- **Missing Residual Mechanism**: The exact mathematical component bridging the `-3dB` to `-9dB` failure margin in the null test is undetermined. Assumed to relate to internal stage clipping, warm-up state, or undocumented gain-staging steps occurring between DSP nodes.
