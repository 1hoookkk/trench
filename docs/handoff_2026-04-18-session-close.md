# Handoff — 2026-04-18 session close

Written for the next fresh Claude session. MEMORY.md carries the permanent
facts; this doc carries the session delta and the exact next action.

## Ship doctrine (locked today)

v1 ships `trench-juce/plugin/` (Windows VST3 + standalone, 2026-07-15).
`trench-core/` Rust is DSP source of truth. JUCE C++ runtime hand-ports
each Rust stage. Render-diff harness enforces bit-accurate parity; any
drift fails a test or a CLI exit code.

`trench-live/` (Rust nih-plug) is a frozen committed scaffold (97cca17)
held for post-v1 consideration. No editor or DSP-wiring work there.

This reverses an AM-of-2026-04-18 amendment that briefly flipped ship to
Rust nih-plug. See commit 241dbd8 for the reconciled spec + plan; see the
updated `Rust Pivot 2026-04-18 (reverted)` memory for the round-trip.

## Shipped this session

- **241dbd8** — reconciled `docs/superpowers/specs/2026-04-18-space-trig-design.md` + `docs/superpowers/plans/2026-04-18-space-trig-implementation.md` + `docs/superpowers/decisions/2026-04-18-editor-visual-path.md`. Plan Task 10 (`safety_limiter` Rust) SUPERSEDED by heritage AGC in `trench-core/src/agc.rs`. Plan Task 11 now = C++ port of AGC. New Plan Task 11b = port E-mu Type 1/2/3 kernels to Rust (`trench-core/src/emu_kernels.rs`). Spec gained `## Heritage kernels (live stage compile)` section explaining why Task 11b blocks the render-diff truth claim beyond Hedz.
- **7c0666e** — render-diff harness: `trench-core/tests/render_diff_harness.rs` + `trench-core/tests/fixtures/hedz_reference_48k.wav` (1 s 48 kHz mono f32 chirp render through `FilterEngine` with baked Hedz) + `tools/render_diff.py` (bit-accurate WAV diff, handles `hound`'s `WAVE_FORMAT_EXTENSIBLE`). `cargo test --test render_diff_harness` passes; self-diff through `render_diff.py` returns max_abs_err = 0.0.

## Next action: step 2b — JUCE offline-render mode

Produce the C++ side of the render-diff baseline. Goal: bit-identical to `hedz_reference_48k.wav` for the same chirp input, same Hedz cartridge, cascade + AGC only (no SPACE, no drive, no mod_fn — matching step 2a scope).

**Landing zones:**
- `trench-juce/plugin/source/PluginProcessor.cpp` — wire an offline render entry (or separate target; see design call below).
- `trench-juce/plugin/CMakeLists.txt` — new target if option B chosen.

**Two design calls before writing any C++:**

1. **Input parity: committed input WAV vs regenerate-from-formula.**
   - Option A: commit `trench-core/tests/fixtures/chirp_input_48k.wav` (the exact buffer `generate_chirp()` builds) and load it on both sides. Zero input-side drift risk.
   - Option B: regenerate the chirp in C++ using the same formula (`CHIRP_F0_HZ`, `CHIRP_F1_HZ`, `INPUT_AMPLITUDE`, phase accumulation). Cheaper (no fixture), but `sin()` and phase arithmetic can differ across libm versions — introduces a second failure mode between the two sides.
   - **Recommended default: A.** Update `render_diff_harness.rs` to also write `chirp_input_48k.wav` on first run and read it back instead of regenerating; re-generate the golden fixture once against that input. Commit both fixtures.

2. **Entry point: `#ifdef TRENCH_OFFLINE_RENDER` in PluginProcessor vs. separate CMake target.**
   - Option A: `#ifdef` branch in `PluginProcessor.cpp::processBlock` (or a new `main()` guarded by the macro). Simple but mixes modes in one binary.
   - Option B: separate `TRENCH_OfflineRender` CMake target with its own tiny `main.cpp` that instantiates `Engine` directly (skipping JUCE host plumbing). Cleaner, doesn't pollute the plugin build, but adds a CMake target.
   - **Recommended default: B.** The render-diff path should exercise `Engine` (the code that actually ships) without being tangled with host callbacks; a dedicated target keeps it surgical.

Once both design calls are made, step 2b is ~1-2 hours of implementation: load input WAV, construct `Engine`, call `process()` equivalent on the block, write output WAV, run `python tools/render_diff.py trench-core/tests/fixtures/hedz_reference_48k.wav <cpp_out>.wav`, expect exit 0.

## After step 2b (longer-horizon)

Rough priority, ship-impact-first:

1. **Plan Task 11 — AGC C++ port** (`PostCascadeAgc.{h,cpp}`). Concrete code snippet is in the plan. Render-diff tolerance 0.0 expected.
2. **Plan Task 11b — EmuKernels Rust port** (`trench-core/src/emu_kernels.rs`). Required to extend render-diff coverage beyond Hedz. Port from `tools/bake_hedz_const.py:88-186` + cross-reference `trench-juce/plugin/source/EmuKernels.cpp`. Fixture JSON of `(type_id, freq_packed, gain_packed, shift, sr)` tuples cross-checks Rust output vs Python bake.
3. **Typed `SpatialProfile` in `cartridge.rs`** — replace the `#[serde(flatten)] raw: serde_json::Map<String, Value>` with concrete `azimuth`/`distance`/`elevation`/`itd_coeffs`/`ild_coeffs`/`band_coeffs` per spec §SPACE. Mirror in C++ `SpatialProfile`.
4. **Plan Task 2** — add `space`/`trig` host params to `PluginProcessor.cpp::createParameterLayout`.
5. **Plan Task 3** — SPACE/TRIG rollers in `PluginEditor.cpp`.
6. **Plan Tasks 4/5** — `mackie_drive` Rust + C++ (blocked on choosing Mackie 1202 model source: Cytomic circuit ref vs measured curve).
7. **Plan Tasks 6/7** — function_generator (Rust already has a scaffold at `trench-core/src/function_generator.rs`, verify against plan Task 6 spec then port).
8. **Plan Tasks 8/9** — `qsound_spatial` Rust + C++. Blocked on step 3 (typed contract).
9. **Plan Task 12** — wire everything into `processBlock` (AGC post-cascade pre-SPACE, SPACE follows AGC, no `limiter` member). Full-chain render-diff gates.
10. **Plan Task 13** — author 4 shipping bodies (Aluminum Siding, Speaker Knockerz, Small Talk Ah-Ee, Cul-De-Sac — confirm names) with authored `drive.input_gain_dB`, `mod_fn`, `spatial_profile`.

## Do NOT

- Touch `trench-live/` editor code. Scaffold stays committed; no visual/DSP work there until post-v1.
- Introduce a new `safety_limiter` or any textbook post-output peak limiter. The heritage AGC IS the tail.
- Re-surface "Rust is v1 ship path" language in any spec. The `## Ship path` block in the space-trig spec is canonical.
- Author body JSON without corresponding Rust DSP. Body content depends on Tasks 4, 6, 8, 11b shipping first.
- Build a WebView-based editor (permanent constraint per memory).

## Current build state

- `cargo test -p trench-core` — all green including `render_diff_harness`.
- `trench-juce/plugin/build/TRENCH_artefacts/Release/Standalone/TRENCH.EXE` — fresh Release build from today at 13:03. Opens; BYPASS wired; BYPASS/cascade/AGC all present in source.
- `trench-juce/plugin/` side: AGC port not yet in; `PostCascadeAgc.{h,cpp}` doesn't exist. (Plan Task 11.)

## Open decisions queued

From plan file: Mackie 1202 model source (Cytomic vs measured). Block Task 4.
From plan file: shipping-body name → file mapping (plan Task 13 Step 0).
From plan file: per-body QSound coefficient source (doc examples vs hand-author). Block Task 13 Step 3.
