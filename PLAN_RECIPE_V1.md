# Recipe-v1: kill the 4-corner bake, run E-mu recipes live

Active plan. Started 2026-04-17 after reece_stab_v1 showed +13.8 dB Hi peak
against Dillusion's +1.5 dB target. Root cause: bilinear interpolation of
pre-baked biquad coefficients does not equal E-mu's Peak/Shelf Morph math.
Fix: interpolate the raw `(freq_packed, gain_packed)` integers between Lo/Hi
frames, then run `type1/2/3_kernel` live to generate coefficients.

## Architecture calls (locked)

- **Scope:** JUCE C++ plugin only (`trench-juce/plugin/source/`). Rust
  `trench-core/` parity port deferred until after v1 ship.
- **Section count:** faithful 6 sections per E-mu. stages[0..5] driven by
  kernels, stages[6..11] passthrough.
- **Update rate:** control rate (every `kBlockSize=32` samples) with per-sample
  linear coefficient ramp. Same cadence as today.
- **Coexistence:** compiled-v1 loader stays unchanged. recipe-v1 is a parallel
  path. Speaker Knockerz and any pole-zero-authored bodies (Aluminum Siding
  via compile_raw.py) keep loading as compiled-v1.
- **Pilot body:** reece_stab_v1. Success = Hi frame peak lands at +1.5 dB in
  forge-view (Dillusion target).

## Step list

1. **Port E-mu kernels to C++.** DONE.
   - `trench-juce/plugin/source/EmuKernels.h` + `.cpp`.
   - Line-for-line port of Python reference in `tools/build_filter_reference.py`
     lines 88-153. Includes a Python-style `floor_div2` helper since C++ `/`
     truncates (diverges from `//` on odd negatives — corrupts gain_offset
     across ~50% of the gain range). `fw_words_to_kernel` combines in double
     then rounds to float to match Python oracle precision.

2. **recipe-v1 JSON schema + loader variant.** IN PROGRESS.
   - Extend `Cartridge` in `TrenchEngine.h` with a format tag enum
     (`CompiledV1` / `RecipeV1`) and a new `Recipe` member holding 6 sections
     × `(type, low_freq, low_gain, high_freq, high_gain, shift)` + lo/hi boost.
   - Extend `Cartridge::loadFromJson` in `TrenchEngine.cpp` to branch on the
     `"format"` field. compiled-v1 path unchanged. recipe-v1 path populates
     `cart.recipe` and skips `computeDeltas()`.
   - Proposed schema shape:
     ```json
     {
       "format": "recipe-v1",
       "name": "reece_stab_v1",
       "sampleRate": 39062.5,
       "lo_boost": 1.0,
       "hi_boost": 1.0,
       "sections": [
         { "index": 0, "type": 1,
           "low_freq": 30, "low_gain": 64,
           "high_freq": 90, "high_gain": 64,
           "shift": 0 }
         // up to 6 entries; type=0 means section off
       ]
     }
     ```

3. **New `Engine::updateCoefficients` branch for recipe-v1.**
   - At each control-rate tick: for each of the 6 sections, lerp
     `(low_freq, low_gain) → (high_freq, high_gain)` by current morph, apply
     Q-offset to the integer (freq shift for type 1/3, gain shift for type 2
     — see `gain_semantics` in `build_filter_reference.py`), then call
     `emu_kernels::compile_stage(type, freq_int, gain_int, shift, sr)` to get
     the 5 coefficients. Write into `targetCoeffs[0..5]`. Fill
     `targetCoeffs[6..11]` as passthrough (`c0=1, rest=0`).
   - Boost: lerp `lo_boost → hi_boost` by morph.
   - Keep the per-sample ramp logic unchanged; the cascade is frozen.

4. **Author reece_stab_v1 as recipe-v1 JSON.**
   - Map Dillusion's Peak/Shelf vocabulary to `(type, freq_packed, gain_packed)`
     per section, two frames. Save alongside compiled-v1 cartridges.
   - Register the new cartridge in `PluginProcessor.cpp` resource list.

5. **Smoke test: Hi frame peak = +1.5 dB in forge-view.**
   - Render the response with morph=1, Q=0 (or whatever corresponds to the
     Dillusion Hi). Read peak dB. Pass = ≤ ±0.5 dB of +1.5.

6. **If pass, author Aluminum Siding in recipe-v1 (if applicable).**
   - Aluminum Siding is currently pole-zero via compile_raw.py. Decision point:
     keep it on compile_raw.py or re-express as recipe-v1. Owner's call on
     sound.

## Non-goals

- No changes to the DF2T cascade (`Stage::process`, `Cascade`). Frozen.
- No changes to `trench-core/` Rust. Parity work deferred.
- No forced migration of existing compiled-v1 bodies.

## Files touched so far

- `trench-juce/plugin/source/EmuKernels.h` (new)
- `trench-juce/plugin/source/EmuKernels.cpp` (new)
