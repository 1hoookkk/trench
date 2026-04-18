# SPACE + TRIG + Mackie Drive Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add SPACE (QSound post-cascade), TRIG (E-mu 8-segment function generator driving MORPH), and pre-cascade Mackie-modeled drive to TRENCH v1 without touching the frozen 12-stage DF2T cascade.

**Architecture — one source of truth:**

- **Rust `trench-core` is the canonical DSP implementation.** All new stages are written in Rust first and covered by `cargo test`.
- **C++ `TrenchEngine.cpp` ports the Rust implementations** for the JUCE runtime. The port is verified by rendering a fixed test input through both Rust and C++ implementations and requiring **sample-accurate parity** on the output WAVs.
- **Only Rust has unit tests.** The C++ plugin has no Catch2, no CPM. `PluginBasics.cpp` may live but is not re-enabled. Verification of C++ is exclusively by render-diff against Rust.
- **Cartridge schema is canonical in `trench-core/src/cartridge.rs`.** The C++ loader in `TrenchEngine.cpp` must parse the same JSON and produce equivalent state.
- The existing C++ cascade port is the pattern to follow for every new stage.

**Tech Stack:** Rust (stable), `hound` for WAV I/O in tests (already dev-dep), JUCE 8, C++17, CMake. No new external dependencies.

**Spec:** `C:/Users/hooki/Trench/docs/superpowers/specs/2026-04-18-space-trig-design.md` (v3, committed).

---

## File Structure

**New Rust files (`C:/Users/hooki/Trench/trench-core/src/`):**
- `mackie_drive.rs` — pre-cascade softclip stage
- `function_generator.rs` — 8-segment E-mu generator
- `qsound_spatial.rs` — post-cascade spatial stage
- `safety_limiter.rs` — final output peak limiter

**New Rust tests (`C:/Users/hooki/Trench/trench-core/tests/`):**
- `mackie_drive_tests.rs`
- `function_generator_tests.rs`
- `qsound_spatial_tests.rs`
- `safety_limiter_tests.rs`
- `render_diff_harness.rs` — full-chain render used by the parity tool

**New render-diff tool (`C:/Users/hooki/Trench/tools/`):**
- `render_diff.py` — runs the C++ standalone in offline mode, compares to Rust reference WAV

**Modified Rust files:**
- `trench-core/src/cartridge.rs` — add `drive`, `spatial_profile`, `mod_fn` optional blocks
- `trench-core/src/lib.rs` — register new modules

**New C++ files (`C:/Users/hooki/trench-juce/plugin/source/`):**
- `MackieDrive.h` / `MackieDrive.cpp`
- `FunctionGenerator.h` / `FunctionGenerator.cpp`
- `QSoundSpatial.h` / `QSoundSpatial.cpp`
- `SafetyLimiter.h` / `SafetyLimiter.cpp`

**Modified C++ files:**
- `PluginProcessor.{h,cpp}` — add `space`, `trig` params; wire stages in `processBlock`; add offline-render mode flag
- `PluginEditor.{h,cpp}` — add SPACE and TRIG rollers
- `TrenchEngine.{h,cpp}` — extend `Cartridge::loadFromJson` for new blocks; wire new stages into the sample loop
- `CMakeLists.txt` — add new .cpp files (no Tests re-enable)

**Modified content:**
- `plugin/assets/*.json` (4 shipping body cartridges) — author `drive`, `mod_fn`, `spatial_profile` fields

---

## Task 1: Cartridge schema bump (Rust canonical, C++ port)

**Files:**
- Modify: `trench-core/src/cartridge.rs`
- Create: `trench-core/tests/cartridge_schema_tests.rs`
- Modify: `plugin/source/TrenchEngine.h` (Cartridge struct)
- Modify: `plugin/source/TrenchEngine.cpp:135` (`Cartridge::loadFromJson`)

- [ ] **Step 1: Write failing Rust test — existing cartridge loads with defaults**

In `cartridge_schema_tests.rs`:
```rust
use trench_core::cartridge::Cartridge;

#[test]
fn existing_cartridge_without_new_blocks_loads_with_safe_defaults() {
    let json = std::fs::read_to_string("../trench-juce/plugin/assets/ah.json")
        .expect("read fixture");
    let cart = Cartridge::from_json(&json).expect("parse");
    assert_eq!(cart.drive.input_gain_db, 0.0);
    assert_eq!(cart.drive.model, "mackie_1202");
    assert!(cart.spatial_profile.is_none());
    assert!(cart.mod_fn.is_none());
}
```

- [ ] **Step 2: Write failing Rust test — cartridge with new blocks parses**

```rust
#[test]
fn cartridge_with_drive_mod_fn_spatial_parses() {
    let json = r#"{
      "format": "compiled-v1", "name": "test", "sampleRate": 48000,
      "keyframes": [],
      "drive": { "input_gain_dB": 9.0, "model": "mackie_1202" },
      "mod_fn": {
        "segments": [{ "level": 1.0, "time_ms": 50, "shape": "exp" }],
        "key-sync": 1, "tempo-sync": 1
      }
    }"#;
    let cart = Cartridge::from_json(json).expect("parse");
    assert_eq!(cart.drive.input_gain_db, 9.0);
    let mf = cart.mod_fn.as_ref().expect("mod_fn present");
    assert_eq!(mf.segments[0].level, 1.0);
    assert!(mf.key_sync);
    assert!(mf.tempo_sync);
}
```

- [ ] **Step 3: Run tests, confirm compile failure (types don't exist)**

```
cd C:/Users/hooki/Trench/trench-core && cargo test cartridge_schema
```

- [ ] **Step 4: Add types and parsing to `cartridge.rs`**

Add, inside `cartridge.rs`:
```rust
#[derive(Debug, Clone, Deserialize)]
pub struct DriveBlock {
    #[serde(rename = "input_gain_dB", default)]
    pub input_gain_db: f32,
    #[serde(default = "default_mackie_model")]
    pub model: String,
}

fn default_mackie_model() -> String { "mackie_1202".to_string() }

impl Default for DriveBlock {
    fn default() -> Self { Self { input_gain_db: 0.0, model: default_mackie_model() } }
}

#[derive(Debug, Clone, Deserialize)]
pub struct SegmentShape(pub String); // "hold" | "lin" | "exp" | "log" | "sCurve"

#[derive(Debug, Clone, Deserialize)]
pub struct FnSegment {
    pub level: f32,
    pub time_ms: f32,
    pub shape: String,
    #[serde(default)]
    pub jump: Option<i32>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ModFnBlock {
    pub segments: Vec<FnSegment>,
    #[serde(rename = "key-sync", default)]
    pub key_sync_int: i32,
    #[serde(rename = "tempo-sync", default)]
    pub tempo_sync_int: i32,
}

impl ModFnBlock {
    pub fn key_sync(&self) -> bool { self.key_sync_int != 0 }
    pub fn tempo_sync(&self) -> bool { self.tempo_sync_int != 0 }
}

#[derive(Debug, Clone, Deserialize)]
pub struct SpatialProfile {
    // Fields populated per docs/archive/qsound_spatial.md during Task 4.
    // For now, an untyped JSON blob is acceptable:
    #[serde(flatten)]
    pub raw: serde_json::Map<String, serde_json::Value>,
}

// Extend the Cartridge struct:
#[derive(Debug, Clone, Deserialize)]
pub struct Cartridge {
    // ... existing fields ...
    #[serde(default)]
    pub drive: DriveBlock,
    #[serde(default)]
    pub spatial_profile: Option<SpatialProfile>,
    #[serde(default)]
    pub mod_fn: Option<ModFnBlock>,
}
```

- [ ] **Step 5: Run tests, confirm pass**

- [ ] **Step 6: Port schema additions to C++ `Cartridge` struct**

In `TrenchEngine.h`, add:
```cpp
struct DriveBlock {
    float input_gain_dB = 0.0f;
    juce::String model = "mackie_1202";
};

struct FnSegment {
    float level = 0.0f;
    float time_ms = 0.0f;
    juce::String shape = "lin";
    int jump = -1;
};

struct ModFnBlock {
    std::vector<FnSegment> segments;
    bool keySync = true;
    bool tempoSync = false;
};

struct SpatialProfile {
    // Mirror Rust fields once populated in Task 4. Blob for now:
    juce::String rawJson;
};

struct Cartridge {
    // ... existing fields ...
    DriveBlock drive;
    std::optional<SpatialProfile> spatialProfile;
    std::optional<ModFnBlock> modFn;
};
```

- [ ] **Step 7: Parse new blocks in `TrenchEngine.cpp::Cartridge::loadFromJson`**

Read optional sub-objects; default when absent. Match the field names and rename rules from the Rust `#[serde(rename = ...)]` attributes.

- [ ] **Step 8: Verify all 8 existing embedded cartridges still load**

Build standalone, launch, cycle through all 8 body types in the selector. Expected: no crashes, same audio output as before (SPACE=0 and TRIG=0 means new blocks are inert).

- [ ] **Step 9: Commit**

```
git -C C:/Users/hooki/Trench commit -am "feat(schema): add drive/spatial/mod_fn optional blocks to cartridge (Rust canonical)"
git -C C:/Users/hooki/trench-juce commit -am "feat(schema): mirror cartridge drive/spatial/mod_fn blocks in C++"
```

---

## Task 2: Add `space` and `trig` plugin parameters with null-default behavior

**Files:**
- Modify: `plugin/source/PluginProcessor.cpp:170-190` (`createParameterLayout`)
- Modify: `plugin/source/PluginProcessor.h` (if an accessor pattern exists)

No Rust change needed — these are plugin-host parameters only. Verification is manual (parameter visible to host with default 0.0).

- [ ] **Step 1: Read `PluginProcessor.cpp` around line 170 and `PluginProcessor.h` to confirm apvts access pattern**

- [ ] **Step 2: Add parameters to `createParameterLayout`**

```cpp
params.push_back (std::make_unique<juce::AudioParameterFloat>(
    juce::ParameterID{"space", 1}, "Space",
    juce::NormalisableRange<float>(0.0f, 1.0f), 0.0f));
params.push_back (std::make_unique<juce::AudioParameterFloat>(
    juce::ParameterID{"trig", 1}, "Trig",
    juce::NormalisableRange<float>(0.0f, 1.0f), 0.0f));
```

- [ ] **Step 3: Build the standalone, verify parameters appear in host automation list**

```
cd C:/Users/hooki/trench-juce/plugin && cmake --build build --target TRENCH_Standalone
```

Launch standalone; confirm SPACE and TRIG are automatable; both default to 0.0.

- [ ] **Step 4: Commit**

```
git -C C:/Users/hooki/trench-juce commit -am "feat: add space and trig parameters, default 0.0"
```

---

## Task 3: Add SPACE and TRIG rollers to the faceplate UI

**Files:**
- Modify: `plugin/source/PluginEditor.{h,cpp}`

- [ ] **Step 1: Survey existing roller construction**

Read `PluginEditor.cpp` for the MORPH and Q roller instantiation + attachment + bounds. Record the pattern.

- [ ] **Step 2: Declare `spaceRoller`, `trigRoller`, and attachment unique_ptrs in `PluginEditor.h`**

- [ ] **Step 3: Construct and attach in the constructor using parameter IDs `"space"` and `"trig"`**

- [ ] **Step 4: Place rollers in `resized()` — commit chosen rectangles as named constants**

Example:
```cpp
static constexpr juce::Rectangle<int> kSpaceRollerBounds { x, y, w, h };
static constexpr juce::Rectangle<int> kTrigRollerBounds  { x, y, w, h };
```

Pick coords visually against the 340×510 faceplate art; document the chosen rectangles with a one-line comment explaining the layout intent.

- [ ] **Step 5: Build and visually verify**

Launch standalone, confirm four rollers are readable and movable, no overlaps.

- [ ] **Step 6: Commit**

```
git -C C:/Users/hooki/trench-juce commit -am "feat(ui): add SPACE and TRIG rollers"
```

---

## Task 4: Implement `mackie_drive` in Rust (null at 0 dB + harmonic generation)

**Files:**
- Create: `trench-core/src/mackie_drive.rs`
- Create: `trench-core/tests/mackie_drive_tests.rs`
- Modify: `trench-core/src/lib.rs` (`pub mod mackie_drive;`)

> **PLAN-PHASE BLOCKER — MACKIE MODEL SOURCE:**
> Do not start this task until one is chosen:
> (a) Cytomic-documented circuit model of the 1202 non-VLZ preamp, or
> (b) a measured curve from a real Mackie 1202.
> Do **not** substitute a tanh placeholder.
>
> Once the model source is chosen, record:
> - **Curve-match tolerance** (e.g. RMS error vs reference ≤ X) in `mackie_drive.rs` as `pub const CURVE_MATCH_RMS_TOLERANCE: f32 = ...`.
> - **Harmonic threshold** for the "nonlinearity engages" test (e.g. sum of 2f..5f bins > Y dB below fundamental at +12 dB input) as a test-file const.

- [ ] **Step 1: Write failing unit test — unity at 0 dB**

```rust
use trench_core::mackie_drive::MackieDrive;

#[test]
fn mackie_drive_is_bit_identical_at_0db() {
    let mut drive = MackieDrive::new(48000.0);
    drive.set_input_gain_db(0.0);
    let input: Vec<f32> = (0..1024).map(|i| (i as f32 * 0.01).sin()).collect();
    let output: Vec<f32> = input.iter().map(|&x| drive.process(x)).collect();
    for (a, b) in input.iter().zip(output.iter()) {
        assert_eq!(a, b); // bit-identical
    }
}
```

- [ ] **Step 2: Write failing unit test — reference curve match**

```rust
#[test]
fn mackie_drive_matches_reference_curve_within_tolerance() {
    // Reference curve is a table of (input, expected_output) pairs
    // derived from the chosen model source. Place the table at
    // trench-core/tests/fixtures/mackie_1202_reference.csv.
    let reference = load_reference_curve("tests/fixtures/mackie_1202_reference.csv");
    let mut drive = MackieDrive::new(48000.0);
    drive.set_input_gain_db(12.0);
    let rms_err = compute_rms_error(&drive, &reference);
    assert!(rms_err < MackieDrive::CURVE_MATCH_RMS_TOLERANCE);
}
```

- [ ] **Step 3: Write failing unit test — harmonics above threshold at +12 dB**

```rust
#[test]
fn mackie_drive_produces_harmonics_at_12db() {
    let mut drive = MackieDrive::new(48000.0);
    drive.set_input_gain_db(12.0);
    let sine_1khz: Vec<f32> = (0..48000)
        .map(|i| (2.0 * std::f32::consts::PI * 1000.0 * i as f32 / 48000.0).sin())
        .collect();
    let output: Vec<f32> = sine_1khz.iter().map(|&x| drive.process(x)).collect();
    let h2_h5_sum_db = measure_harmonic_sum_db(&output, 1000.0, 48000.0, 2, 5);
    assert!(h2_h5_sum_db > HARMONIC_THRESHOLD_DB_BELOW_FUNDAMENTAL);
}
```

- [ ] **Step 4: Confirm all three fail (type/function doesn't exist)**

- [ ] **Step 5: Implement `MackieDrive`**

```rust
pub struct MackieDrive {
    sample_rate: f32,
    gain_lin: f32,
}

impl MackieDrive {
    pub const CURVE_MATCH_RMS_TOLERANCE: f32 = /* set when model source chosen */;

    pub fn new(sample_rate: f32) -> Self {
        Self { sample_rate, gain_lin: 1.0 }
    }

    pub fn set_input_gain_db(&mut self, db: f32) {
        self.gain_lin = 10f32.powf(db / 20.0);
    }

    pub fn process(&mut self, x: f32) -> f32 {
        if self.gain_lin == 1.0 { return x; } // short-circuit unity
        softclip_mackie_1202(x * self.gain_lin)
    }

    pub fn reset(&mut self) { /* stateless softclip: no-op */ }
}

fn softclip_mackie_1202(x: f32) -> f32 {
    // Implementation per chosen model source.
}
```

- [ ] **Step 6: Run tests, confirm all three pass**

- [ ] **Step 7: Commit**

```
git -C C:/Users/hooki/Trench commit -am "feat(dsp): add mackie_drive Rust impl + tests"
```

---

## Task 5: Port `MackieDrive` to C++ (render-diff parity)

**Files:**
- Create: `plugin/source/MackieDrive.{h,cpp}`
- Modify: `plugin/CMakeLists.txt` (add source)

- [ ] **Step 1: Port the Rust `softclip_mackie_1202` + `MackieDrive` struct to C++**

Match the signature:
```cpp
class MackieDrive {
public:
    void prepare(double sampleRate);
    void setInputGainDb(float db);
    float processSample(float x) noexcept;
    void reset();
private:
    double sr_ = 48000.0;
    float gainLin_ = 1.0f;
};
```

- [ ] **Step 2: Build, ensure clean compile and link**

- [ ] **Step 3: Render-diff verification**

From a Rust test binary, render 1 second of pink noise through `mackie_drive::MackieDrive` at +12 dB. Save to `trench-core/tests/fixtures/mackie_ref_out.wav`.

Write a Rust standalone binary `tools/render_mackie_ref.rs` (or a one-shot test) that does this. Commit the WAV.

Then: build the C++ plugin with a temporary `#ifdef TRENCH_OFFLINE_RENDER` entry point that runs the same input through the C++ `MackieDrive` and writes `mackie_cpp_out.wav`.

Script `tools/render_diff.py` compares the two WAVs sample-for-sample and reports max absolute error. Required: max error ≤ 1 ULP (bit-identical or nearly so for a stateless softclip).

Commit the diff tool and the fixture WAV.

- [ ] **Step 4: Commit**

```
git -C C:/Users/hooki/Trench add tools/render_diff.py trench-core/tests/fixtures/mackie_ref_out.wav
git -C C:/Users/hooki/Trench commit -m "test: render-diff harness + mackie reference WAV"
git -C C:/Users/hooki/trench-juce commit -am "feat(dsp): port MackieDrive to C++ with offline-render entry"
```

---

## Task 6: Implement `function_generator` in Rust (8 segments, all shapes, jumps, REV)

**Files:**
- Create: `trench-core/src/function_generator.rs`
- Create: `trench-core/tests/function_generator_tests.rs`
- Modify: `trench-core/src/lib.rs`

- [ ] **Step 1: Write failing test — silent when no segments loaded**

```rust
#[test]
fn function_generator_outputs_0_with_no_segments() {
    let mut fg = FunctionGenerator::new(48000.0);
    for _ in 0..1000 { assert_eq!(fg.tick(), 0.0); }
}
```

- [ ] **Step 2: Write failing test — linear segment interpolation**

Test that a linear segment from 0 -> 1 over 1000ms produces ~0.5 at the midpoint.

- [ ] **Step 3: Write parametric failing tests for each shape: `hold`, `lin`, `exp`, `log`, `sCurve`**

One test per shape. For each, the midpoint output must match the expected curve value within `1e-4`.

- [ ] **Step 4: Write failing test — conditional jump loops segments**

Set up segments 0, 1, 2 with `jump: 0` on segment 2. After segment 2 end, output should reflect segment 0's start, not segment 3.

- [ ] **Step 5: Write failing test — key_sync restarts on note-on**

- [ ] **Step 6: Write failing test — REV time-reverses trajectory**

Record forward-mode output samples. Reset, enable REV, run again. Final forward sample ≈ first REV sample, first forward sample ≈ last REV sample.

- [ ] **Step 7: Write failing test — tempo_sync divides time_ms by tempo factor**

- [ ] **Step 8: Run all tests, confirm failure**

- [ ] **Step 9: Implement `FunctionGenerator`**

```rust
pub struct FunctionGenerator {
    sample_rate: f32,
    segments: Vec<FnSegment>,
    segment_idx: usize,
    sample_in_segment: u32,
    reverse: bool,
    key_sync: bool,
    tempo_sync: bool,
    tempo_bpm: f32,
    last_level: f32, // starting level for current segment's ramp
}

impl FunctionGenerator {
    pub fn new(sample_rate: f32) -> Self { /* ... */ }
    pub fn load(&mut self, block: &ModFnBlock) { /* ... */ }
    pub fn note_on(&mut self) { /* reset indices if key_sync */ }
    pub fn set_reverse(&mut self, rev: bool) { /* ... */ }
    pub fn set_tempo_bpm(&mut self, bpm: f32) { /* ... */ }
    pub fn tick(&mut self) -> f32 { /* advance one sample, honor jumps */ }
    pub fn reset(&mut self) { /* ... */ }
}

fn apply_shape(shape: &str, t01: f32, from: f32, to: f32) -> f32 { /* ... */ }
```

- [ ] **Step 10: Run all tests, confirm pass**

- [ ] **Step 11: Commit**

```
git -C C:/Users/hooki/Trench commit -am "feat(dsp): add function_generator Rust impl + tests"
```

---

## Task 7: Port `FunctionGenerator` to C++

**Files:**
- Create: `plugin/source/FunctionGenerator.{h,cpp}`
- Modify: `plugin/CMakeLists.txt`

- [ ] **Step 1: Port the Rust struct and methods to C++**

Mirror the Rust API closely:
```cpp
class FunctionGenerator {
public:
    void prepare(double sampleRate);
    void load(const ModFnBlock& block);
    void noteOn();
    void setReverse(bool rev);
    void setTempoBpm(double bpm);
    float tick() noexcept;
    void reset();
};
```

- [ ] **Step 2: Render-diff verification**

Author a small test block in JSON (fixed segments). Render 2 seconds of output through the Rust generator (save to `fg_ref_out.wav`). Render through C++ via the offline-render entry point (`fg_cpp_out.wav`). Diff: max absolute error ≤ 1e-6.

- [ ] **Step 3: Commit**

---

## Task 8: Implement `qsound_spatial` in Rust (ITD + ILD + band-law)

**Files:**
- Create: `trench-core/src/qsound_spatial.rs`
- Create: `trench-core/tests/qsound_spatial_tests.rs`
- Modify: `trench-core/src/lib.rs`
- Reference: `docs/archive/qsound_spatial.md`

- [ ] **Step 1: Read the QSound doc carefully**

Confirm: ITD law (fractional delay per ear), ILD law (gain law per ear), three-band filter coefficients, sample-rate assumptions.

- [ ] **Step 2: Populate the `SpatialProfile` Rust struct with concrete fields**

Replace the blob-of-JSON from Task 1 with typed fields:
```rust
pub struct SpatialProfile {
    pub azimuth: f32,
    pub distance: f32,
    pub elevation: f32,
    pub itd_coeffs: ItdLawCoeffs,
    pub ild_coeffs: IldLawCoeffs,
    pub band_coeffs: [BandFilterCoeffs; 3],
}
```

Update Task 1's cartridge test fixtures with real-shaped values.

- [ ] **Step 3: Write failing test — SPACE = 0 produces bit-identical bypass**

- [ ] **Step 4: Write failing test — SPACE = 1.0 widens stereo**

L/R correlation on a mono input must be lower after processing at SPACE=1.0 than at SPACE=0.0.

- [ ] **Step 5: Write failing test — ITD law applied correctly**

At known azimuth, the inter-channel delay should match the ITD coefficient table.

- [ ] **Step 6: Implement `QSoundSpatial`**

Includes:
- Fractional delay lines (L, R) with cubic interpolation
- Per-channel ILD gain
- Three biquads per channel for band-law shaping
- Explicit bypass branch when SPACE == 0.0 (no math path)
- Dry/wet crossfade linear in SPACE

- [ ] **Step 7: Run tests, confirm pass**

- [ ] **Step 8: Commit**

---

## Task 9: Port `QSoundSpatial` to C++

**Files:**
- Create: `plugin/source/QSoundSpatial.{h,cpp}`
- Modify: `plugin/CMakeLists.txt`
- Modify: `plugin/source/TrenchEngine.h` (populate `SpatialProfile` fields to match Rust)

- [ ] **Step 1: Port delay lines, gain laws, biquads**

- [ ] **Step 2: Render-diff verification**

Mono test signal → through both Rust and C++ QSound stages at SPACE=1.0 with a specific spatial profile. Diff stereo output; tolerance 1e-5 (floating-point rounding).

- [ ] **Step 3: Commit**

---

## Task 10: Implement `safety_limiter` in Rust

**Files:**
- Create: `trench-core/src/safety_limiter.rs`
- Create: `trench-core/tests/safety_limiter_tests.rs`
- Modify: `trench-core/src/lib.rs`

- [ ] **Step 1: Write failing test — transparent below threshold**

- [ ] **Step 2: Write failing test — catches peaks above ceiling**

- [ ] **Step 3: Implement peak limiter**

Simple lookahead-free. Ceiling: -0.3 dBFS. Stereo-linked gain reduction with short release (e.g. 50 ms).

- [ ] **Step 4: Run tests, confirm pass**

- [ ] **Step 5: Commit**

---

## Task 11: Port `SafetyLimiter` to C++

**Files:**
- Create: `plugin/source/SafetyLimiter.{h,cpp}`
- Modify: `plugin/CMakeLists.txt`

- [ ] **Step 1: Port limiter state machine to C++**

- [ ] **Step 2: Render-diff verification**

Render a chirp-with-peaks through both. Diff ≤ 1e-5.

- [ ] **Step 3: Commit**

---

## Task 12: Wire all stages into `processBlock` (C++ plumbing)

**Files:**
- Modify: `plugin/source/TrenchEngine.{h,cpp}`
- Modify: `plugin/source/PluginProcessor.cpp:364`

- [ ] **Step 1: Add `Engine` members: `driveL`, `driveR`, `modFn`, `spatial`, `limiter`**

- [ ] **Step 2: Apply drive in the per-sample loop before stage 1 of the cascade**

- [ ] **Step 3: Apply TRIG modulation to MORPH**

Read `trig` param per control block; tick `modFn` at the sample rate; combine:
```cpp
float fnOut = modFn.tick();
float trigAmount = *apvts.getRawParameterValue("trig");
float effectiveMorph = juce::jlimit(0.0f, 1.0f, userMorph + trigAmount * fnOut);
```

Route REV parameter to `modFn.setReverse(...)`.

- [ ] **Step 4: Apply SPACE post-cascade**

```cpp
float space = *apvts.getRawParameterValue("space");
spatial.setSpace(space);
spatial.processBlock(buffer);
```

- [ ] **Step 5: Apply limiter at final output**

- [ ] **Step 6: Wire MIDI note-on to `modFn.noteOn()`**

- [ ] **Step 7: Full-chain render-diff**

Render `tests/fixtures/noise_burst_48k.wav` through:
- Rust reference signal chain (in `trench-core/tests/render_diff_harness.rs`)
- C++ plugin offline-render mode

Require sample-accurate parity across the full chain for a known body with authored `drive`, `mod_fn`, `spatial_profile`. Diff tolerance: ≤ 1e-4 (accumulated floating-point error across many stages).

- [ ] **Step 8: Commit**

```
git -C C:/Users/hooki/trench-juce commit -am "feat: wire MackieDrive, FunctionGenerator, QSoundSpatial, SafetyLimiter into processBlock"
```

---

## Task 13: Author 4 shipping bodies (hand-written JSON)

**Files:**
- Modify: `plugin/assets/<body>.json` × 4

Bodies: Aluminum Siding, Speaker Knockerz, Small Talk Ah-Ee, Cul-De-Sac. Confirm name→file mapping with the user before starting.

For each body:

- [ ] **Step 1: Add `drive` block** with an authored `input_gain_dB` (begin 0.0 for Aluminum Siding; escalate for hotter bodies).

- [ ] **Step 2: Add `mod_fn` block** with 1–4 segments, `key-sync: 1`; loop via `jump` for repeating bodies, no jump for one-shots.

- [ ] **Step 3: Add `spatial_profile` block**

If `docs/archive/qsound_spatial.md` contains per-body or per-character example coefficient sets, copy the one matching the body's intended voice. If the doc only documents laws, hand-author azimuth/distance/elevation/band coeffs by ear: start from the doc's default/center values, audition, commit. **Do not ship identity/placeholder** — SPACE=1 should audibly differ from SPACE=0 on an authored body.

- [ ] **Step 4: Build, load body, audition**

- [ ] **Step 5: Commit per body**

---

## Task 14: Final verification

- [ ] **Step 1: `cargo test` passes in `trench-core`**

```
cd C:/Users/hooki/Trench/trench-core && cargo test
```

Expected: all tests pass.

- [ ] **Step 2: `./check` passes at repo root**

```
cd C:/Users/hooki/Trench && ./check
```

- [ ] **Step 3: Render-diff passes on all 4 shipping bodies**

Run `tools/render_diff.py` for each body with SPACE=0/TRIG=0 (pre-change parity), SPACE=1/TRIG=0, SPACE=0/TRIG=1, and SPACE=1/TRIG=1. All must meet the tolerance documented in the relevant task.

- [ ] **Step 4: Standalone audition**

For each body, in the standalone build:
- SPACE=0 + TRIG=0 matches current release character.
- SPACE>0 widens stereo audibly.
- TRIG>0 makes MORPH move on each note-on (phase reset when `key-sync=1`).
- REV reverses the TRIG trajectory.
- Hot drive bodies audibly bite harder than clean bodies.

- [ ] **Step 5: Walk spec §Verification checklist**

Any missed items become followup tasks.

- [ ] **Step 6: Final commit**

```
git -C C:/Users/hooki/Trench commit --allow-empty -m "ship: SPACE + TRIG + Mackie drive complete per spec 2026-04-18"
```

---

## Open questions held for plan execution

1. **Mackie model source** (Task 4 blocker): Cytomic circuit reference vs measured 1202 curve. Must be resolved before Task 4 starts.
2. **Curve-match tolerance + harmonic threshold** (Task 4): once model source is picked, record both in Rust source.
3. **Shipping body name→cartridge mapping** (Task 13): confirm the 4 names map to specific embedded cartridge files.
4. **QSound per-body presets availability** (Task 13 Step 3): whether `docs/archive/qsound_spatial.md` contains per-body examples or only laws.

## Advisory items

- **Task 3 roller placement**: consider invoking the `juce-lookandfeel-design` skill for disciplined pixel coords; commit chosen rectangles as named constants, not magic numbers.
- **Task 6 shape coverage**: parametric test covering `hold`, `lin`, `exp`, `log`, `sCurve` is explicitly in scope now (Step 3).
- **Task 12 generator tick rate**: per-sample `modFn.tick()` is the simplest; control-block rate is a cheap optimization if profiling ever demands it.
- **Offline-render entry point**: Task 5 adds a `#ifdef TRENCH_OFFLINE_RENDER` path to the C++ plugin. Keep it compiled-out in release; document the macro in a comment.
