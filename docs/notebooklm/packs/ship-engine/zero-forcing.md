# Zero-Forcing — Root Cause & Fix Plan

## Root Cause (confirmed 2026-03-08)
Generated filters are spectrally flat because the entire generation path uses **RBJ cookbook biquads** (`stage_math.rs`), which lock zeros to textbook filter shapes. There is NO independent zero control — a Peak at 1kHz Q=2 always produces the same zero placement.

Original E-mu presets use **free-form resonator stages** where val1/val2/val3 can place zeros ANYWHERE relative to poles. That's what creates 30dB comb nulls, talking character, exotic interference.

## Pipeline Path
```
PresetGrammar (freq, Q, gain)
  → Node.compile() → StageRecipe (mode, freq, Q, gain)
    → stage_math::compute_stage() → RBJ cookbook biquad
      → to_kernel() → [c0..c4]
        → export::kernel_to_raw_json() → val1/val2/val3
```

## Counterpoint ZeroField (exists but unused)
- `counterpoint.rs` has `ZeroField` with `val2_delta`/`val3_delta` for independent zero placement
- Almost every callsite passes `ZeroField::ALL_POLE` — wired but never populated
- Composition struct is only used in counterpoint tests, not in the actual generation pipeline

## Zero-Forcing Formulas (from NotebookLM, verified against sources)
For a resonator stage with pole at angle θ, radius r:
- **Unit-circle zero at SAME angle** (deep notch at resonant freq):
  - `val3 = r² - 1 - val1`
  - `val2 = 2(r - 1 - val1)cos(θ)`
- **Unit-circle zero at DIFFERENT angle φ** (peak-then-notch pair):
  - `val3 = r² - 1 - val1`
  - `val2 = -2(1+val1)cos(φ) - a1`
- ARMAdillo k2 → high values naturally push zero radius → 1.0

## Fix Plan
1. Add `StageMode::Resonator` that bypasses RBJ and computes (a1, r, val1, val2, val3) directly
2. New function: `zero_forced_resonator(freq_hz, r, zero_angle, zero_radius)` → EncodedCoeffs
3. In generate.rs, for Peak/Cut archetypes: use zero-forced resonator instead of RBJ Peak/Notch
4. The `zero_angle` relative to pole angle controls spectral character:
   - Same angle = cancellation notch
   - Offset angle = peak-then-notch pair (Sculpted character)
   - Random offset = current bland behavior (avoid)
5. Alternatively: populate counterpoint.rs `ZeroField` properly and route generation through Composition

## Key Files
- `trench-forge/src/stage_math.rs` — RBJ cookbook (current path, no free zeros)
- `trench-forge/src/counterpoint.rs` — ZeroField system (unused in real generation)
- `trench-forge/src/generate.rs` — generation pipeline
- `trench-forge/src/grammar.rs` — Node.compile() → StageRecipe
- `trench-forge/src/export.rs` — kernel_to_raw_json (val1/val2/val3 output)

## NotebookLM Findings (2026-03-08)
- Original corners designed as 12-pole collective, not 6 independent stages
- Pole-zero clustering creates spike-then-null pairs
- Phase coherence between stages sharpens resonances
- Sculpted archetype [Peak, Cut, Peak, Cut, Peak, Peak] provides high spectral contrast
- Frequency overlap between adjacent stages creates complex interactions
- Multi-source vocabulary unification: all normalize to kernel-form at 39062.5Hz
- EX XML templates contain design intent (zero placement strategies)
- Morpheus 3rd axis (Depth) can map to morph species (Swap, Migrate)

## E4Ultra MorphDesigner Handoff Analysis (2026-03-08)
- Sanitized handoff at `trench_re_vault/artifacts/morph_designer_sanitized_handoff/e4ultra/morph_designer/v1/`
- **NOT the star generator**. Small, deterministic, 3-shape stage-recipe compiler (bandpass, resonant_formant, shaped_asymmetric)
- 2-corner only (low/high, duplicated to 4-corner pipeline) — inherently 1D morphs
- `global_shift = -32 + floor((freq_macro + gain_macro) * 63)` — correlated motion across all stages
- Correlated motion = coherent but boring. No antagonism, no independent actors, no crossing.
- Radius linearly derived from frequency: `radius = floor(freq * 124/256) + 118` — Q coupled to freq by construction
- 79 corpus presets × 4 sample rates = 316 exact test vectors
- **Correct use**: legacy reference lane, import lane, compatibility lane, teacher dataset, corpus validation
- **Wrong use**: main generation engine (produces competent low-entropy designer curves, not monsters)
- Replacing RBJ with MorphDesigner algebra would trade one constrained recipe for another — same ceiling
- The real fix for zero-forcing is kernel/pole-zero space with antagonistic multi-actor topologies
