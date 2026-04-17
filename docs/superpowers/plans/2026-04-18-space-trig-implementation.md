# SPACE + TRIG + Mackie Drive Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add SPACE (QSound post-cascade), TRIG (E-mu 8-segment function generator driving MORPH), and pre-cascade Mackie-modeled drive to TRENCH v1 without touching the frozen 12-stage DF2T cascade.

**Architecture:** All new DSP lives in C++ in the JUCE plugin (`C:/Users/hooki/trench-juce/plugin/source/`). The plugin runs its own C++ cascade (`TrenchEngine.cpp::processBlock`); the Rust `trench-core/` crate is reference and test impl — not runtime. New stages wrap the existing cascade without modifying its sample-loop code. Cartridge JSON schema gains three optional blocks (`drive`, `spatial_profile`, `mod_fn`) with safe defaults. Shipping bodies are hand-authored JSON.

**Tech Stack:** JUCE 8, C++17, CMake, Catch2 (currently disabled, re-enabled in Task 1). Cartridges are JSON. No new external dependencies.

**Spec:** `C:/Users/hooki/Trench/docs/superpowers/specs/2026-04-18-space-trig-design.md` (v3, committed).

---

## File Structure

**New files:**
- `C:/Users/hooki/trench-juce/plugin/source/MackieDrive.h` — pre-cascade softclip stage
- `C:/Users/hooki/trench-juce/plugin/source/MackieDrive.cpp`
- `C:/Users/hooki/trench-juce/plugin/source/FunctionGenerator.h` — 8-segment E-mu generator
- `C:/Users/hooki/trench-juce/plugin/source/FunctionGenerator.cpp`
- `C:/Users/hooki/trench-juce/plugin/source/QSoundSpatial.h` — post-cascade spatial stage
- `C:/Users/hooki/trench-juce/plugin/source/QSoundSpatial.cpp`
- `C:/Users/hooki/trench-juce/plugin/source/SafetyLimiter.h` — final output peak limiter
- `C:/Users/hooki/trench-juce/plugin/source/SafetyLimiter.cpp`
- `C:/Users/hooki/trench-juce/plugin/tests/NullBaselineTests.cpp` — verifies zero-default null
- `C:/Users/hooki/trench-juce/plugin/tests/FunctionGeneratorTests.cpp`
- `C:/Users/hooki/trench-juce/plugin/tests/MackieDriveTests.cpp`

**Modified files:**
- `C:/Users/hooki/trench-juce/plugin/source/PluginProcessor.cpp` — add `space`, `trig` params; call new stages in `processBlock`
- `C:/Users/hooki/trench-juce/plugin/source/PluginProcessor.h` — forward declarations
- `C:/Users/hooki/trench-juce/plugin/source/PluginEditor.cpp` — add two new rollers
- `C:/Users/hooki/trench-juce/plugin/source/PluginEditor.h`
- `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.h` — add fields for drive/fn-gen/spatial state; method signatures for MORPH modulation hook
- `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.cpp` — wire pre/post-cascade stages; extend `Cartridge::loadFromJson` for new blocks
- `C:/Users/hooki/trench-juce/plugin/CMakeLists.txt` — add new .cpp files; re-enable Tests.cmake inclusion
- `C:/Users/hooki/trench-juce/plugin/assets/*.json` (4 shipping body cartridges) — author `drive`, `mod_fn`, `spatial_profile` fields

---

## Task 1: Re-enable Catch2 tests and add new-file scaffolding

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/CMakeLists.txt:166`
- Modify: `C:/Users/hooki/trench-juce/plugin/tests/CMakeLists.txt` (verify exists; list new test files)

- [ ] **Step 1: Find and un-disable tests in CMakeLists.txt**

Read `C:/Users/hooki/trench-juce/plugin/CMakeLists.txt` around line 166. Comment block or `if(FALSE)` currently disables `include(tests/Tests.cmake)` or similar. Un-disable it.

- [ ] **Step 2: Verify the test target builds with the existing `PluginBasics.cpp`**

Run from `C:/Users/hooki/trench-juce/plugin/`:
```
cmake -S . -B build && cmake --build build --target TRENCH_Tests
```
Expected: build succeeds; `PluginBasics.cpp` passes.

- [ ] **Step 3: Register placeholder test files**

Add to `tests/CMakeLists.txt` target sources list:
- `tests/NullBaselineTests.cpp`
- `tests/FunctionGeneratorTests.cpp`
- `tests/MackieDriveTests.cpp`

Create each file as an empty stub:
```cpp
#include <catch2/catch_test_macros.hpp>
// Tests added in subsequent tasks.
```

- [ ] **Step 4: Rebuild and verify no regressions**

```
cmake --build build --target TRENCH_Tests && build/.../TRENCH_Tests
```
Expected: existing tests pass, empty stubs produce no failures.

- [ ] **Step 5: Commit**

```
git -C C:/Users/hooki/trench-juce/plugin add CMakeLists.txt tests/
git -C C:/Users/hooki/trench-juce/plugin commit -m "test: re-enable Catch2 and scaffold new test files"
```

---

## Task 2: Add `space` and `trig` plugin parameters with null-default behavior

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/source/PluginProcessor.cpp:170-190` (`createParameterLayout`)
- Modify: `C:/Users/hooki/trench-juce/plugin/source/PluginProcessor.h` (accessor getters if pattern exists)
- Create: `C:/Users/hooki/trench-juce/plugin/tests/NullBaselineTests.cpp`

- [ ] **Step 1: Confirm `apvts` access pattern**

Open `PluginProcessor.h` around line 38-39. Verify whether `apvts` is public, or whether the class exposes a `getApvts()` accessor. Use whichever the header provides; the test below assumes a public or friend-accessible `apvts`. If it's private with no accessor, add a test-only accessor or use `juce::AudioProcessor::getParameters()` by ID lookup.

- [ ] **Step 2: Write failing null-baseline test**

In `NullBaselineTests.cpp`:
```cpp
#include <catch2/catch_test_macros.hpp>
#include "../source/PluginProcessor.h"

TEST_CASE("space and trig parameters exist with default 0.0") {
    PluginProcessor p;
    // Use whichever access pattern Step 1 confirmed:
    auto& state = p.apvts; // or p.getApvts()
    REQUIRE(state.getRawParameterValue("space") != nullptr);
    REQUIRE(state.getRawParameterValue("trig")  != nullptr);
    REQUIRE(*state.getRawParameterValue("space") == 0.0f);
    REQUIRE(*state.getRawParameterValue("trig")  == 0.0f);
}
```

- [ ] **Step 3: Run, confirm failure**

Expected: test compiles but fails at REQUIRE (parameters don't exist yet).

- [ ] **Step 4: Add parameters to `createParameterLayout`**

In `PluginProcessor.cpp`, append inside the parameter list:
```cpp
params.push_back (std::make_unique<juce::AudioParameterFloat>(
    juce::ParameterID{"space", 1}, "Space",
    juce::NormalisableRange<float>(0.0f, 1.0f), 0.0f));
params.push_back (std::make_unique<juce::AudioParameterFloat>(
    juce::ParameterID{"trig", 1}, "Trig",
    juce::NormalisableRange<float>(0.0f, 1.0f), 0.0f));
```

- [ ] **Step 5: Run test, confirm pass**

- [ ] **Step 6: Commit**

```
git commit -am "feat: add space and trig parameters, default 0.0"
```

---

## Task 3: Add SPACE and TRIG rollers to the faceplate UI

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/source/PluginEditor.h`
- Modify: `C:/Users/hooki/trench-juce/plugin/source/PluginEditor.cpp`

- [ ] **Step 1: Survey current roller layout**

Read the existing roller construction for MORPH and Q in `PluginEditor.cpp` (grep for `RollerComponent` or `morphRoller`). Record the exact pattern (member, attachment, bounds).

- [ ] **Step 2: Declare new rollers in `PluginEditor.h`**

```cpp
RollerComponent spaceRoller;
RollerComponent trigRoller;
std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> spaceAttach;
std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> trigAttach;
```

(Match existing ownership pattern; if the plugin uses a different attachment type, mirror it.)

- [ ] **Step 3: Construct and attach in the editor constructor**

Copy the MORPH roller initialization twice, renamed. Attachment IDs = `"space"` and `"trig"`.

- [ ] **Step 4: Place rollers in `resized()`**

Lay out SPACE and TRIG below MORPH and Q, matching the faceplate bounds (340×510). Exact pixel coordinates: pick visually against the faceplate art — document the chosen rectangles in a comment.

- [ ] **Step 5: Build, run standalone, visually verify**

```
cmake --build build --target TRENCH_Standalone
```
Launch, confirm two new rollers appear, moving them changes the parameter (visible via DAW automation list or host).

- [ ] **Step 6: Commit**

```
git commit -am "feat(ui): add SPACE and TRIG rollers to faceplate"
```

---

## Task 4: Extend cartridge JSON loader with optional `drive`, `spatial_profile`, `mod_fn` blocks

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.h` (Cartridge struct)
- Modify: `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.cpp:135` (`Cartridge::loadFromJson`)
- Test: `C:/Users/hooki/trench-juce/plugin/tests/NullBaselineTests.cpp` (add cases)

- [ ] **Step 1: Write failing test — existing cartridge still loads**

```cpp
TEST_CASE("existing cartridge without new blocks loads with safe defaults") {
    auto json = loadTestFixture("ah.json"); // helper: read an embedded cart
    auto cart = trench::Cartridge::loadFromJson(json);
    REQUIRE(cart.drive.input_gain_dB == 0.0f);
    REQUIRE(cart.drive.model == "mackie_1202");
    REQUIRE(cart.spatial_profile.has_value() == false);
    REQUIRE(cart.mod_fn.has_value() == false);
}
```

Add the fixture helper if not present: read from `plugin/assets/`.

- [ ] **Step 2: Run, confirm failure (struct fields don't exist)**

- [ ] **Step 3: Add fields to `Cartridge` struct in `TrenchEngine.h`**

```cpp
struct DriveBlock {
    float input_gain_dB = 0.0f;
    juce::String model = "mackie_1202";
};

struct SpatialProfile { /* fields populated in Task 7 */ };

struct ModFunctionSegment {
    float level   = 0.0f;
    float time_ms = 0.0f;
    juce::String shape = "lin";
    int jump = -1; // -1 = no jump
};

struct ModFunctionGenerator {
    std::array<ModFunctionSegment, 8> segments{};
    bool key_sync   = true;
    bool tempo_sync = false;
};

struct Cartridge {
    // ...existing fields...
    DriveBlock drive;
    std::optional<SpatialProfile> spatial_profile;
    std::optional<ModFunctionGenerator> mod_fn;
};
```

- [ ] **Step 4: Parse in `loadFromJson`**

Read optional blocks; absent = defaults. JSON shape per spec §Body schema changes.

- [ ] **Step 5: Add test — cartridge WITH new blocks parses correctly**

```cpp
TEST_CASE("cartridge with drive/mod_fn/spatial parses") {
    auto json = R"({
      "format": "compiled-v1", "name": "test", "sampleRate": 48000,
      "keyframes": [...minimal stub...],
      "drive": { "input_gain_dB": 9.0, "model": "mackie_1202" },
      "mod_fn": { "segments": [{"level": 1.0, "time_ms": 50, "shape": "exp"}],
                  "key-sync": 1, "tempo-sync": 1 }
    })";
    auto cart = trench::Cartridge::loadFromJson(json);
    REQUIRE(cart.drive.input_gain_dB == 9.0f);
    REQUIRE(cart.mod_fn.has_value());
    REQUIRE(cart.mod_fn->segments[0].level == 1.0f);
}
```

- [ ] **Step 6: Run all tests, confirm pass, all 8 existing cartridges still load**

- [ ] **Step 7: Commit**

```
git commit -am "feat(cartridge): add optional drive/spatial/mod_fn blocks with safe defaults"
```

---

## Task 5: Implement `MackieDrive` pre-cascade stage (null at 0 dB)

**Files:**
- Create: `C:/Users/hooki/trench-juce/plugin/source/MackieDrive.h`
- Create: `C:/Users/hooki/trench-juce/plugin/source/MackieDrive.cpp`
- Test: `C:/Users/hooki/trench-juce/plugin/tests/MackieDriveTests.cpp`
- Modify: `C:/Users/hooki/trench-juce/plugin/CMakeLists.txt` (add source)

> **PLAN-PHASE DECISION — MACKIE MODEL SOURCE:**
> Pick one before writing the softclip: (a) Cytomic-documented circuit model of the 1202 non-VLZ preamp, or (b) a measured curve from a real 1202. If neither is accessible before implementation starts, escalate to the user — do not substitute a tanh placeholder. This is recorded in the spec as a non-negotiable.
>
> **PLAN-PHASE DECISION — CURVE-MATCH TOLERANCE:**
> Once the model source is chosen, pick a curve-match tolerance (e.g. RMS error vs reference ≤ X) and a harmonic-threshold for Step 2's nonlinearity test (e.g. sum of 2f..5f bins > Y dB below fundamental at +12 dB). Record both values in a comment at the top of `MackieDrive.cpp`; tests consume them as constants.

- [ ] **Step 1: Write failing test — unity at 0 dB**

```cpp
TEST_CASE("MackieDrive is bit-identical at input_gain_dB = 0") {
    MackieDrive drive;
    drive.setInputGainDb(0.0f);
    std::array<float, 512> in; // fill with a sine or noise
    std::array<float, 512> out = in;
    for (size_t i = 0; i < in.size(); ++i) out[i] = drive.processSample(in[i]);
    for (size_t i = 0; i < in.size(); ++i) REQUIRE(out[i] == in[i]);
}
```

- [ ] **Step 2: Write failing test — hot drive produces harmonics**

```cpp
TEST_CASE("MackieDrive at +12 dB introduces nonlinearity") {
    MackieDrive drive;
    drive.setInputGainDb(12.0f);
    // pump a 1 kHz sine, expect harmonic content > threshold
    // (implementation detail: FFT the output, sum 2f..5f bins)
    ...
}
```

- [ ] **Step 3: Run, confirm failures**

- [ ] **Step 4: Implement `MackieDrive`**

Skeleton:
```cpp
class MackieDrive {
public:
    void prepare(double sampleRate);
    void setInputGainDb(float db);
    float processSample(float x) noexcept; // unity when db == 0
    void reset();
private:
    float gainLin_ = 1.0f;
    // softclip curve state per selected model
};
```

Implementation uses the model chosen in the plan-phase decision above. At `db == 0`, short-circuit: return `x` unchanged.

- [ ] **Step 5: Run tests, confirm pass**

- [ ] **Step 6: Commit**

```
git commit -am "feat(dsp): add MackieDrive pre-cascade softclip stage"
```

---

## Task 6: Implement `FunctionGenerator` (E-mu 8-segment)

**Files:**
- Create: `C:/Users/hooki/trench-juce/plugin/source/FunctionGenerator.h`
- Create: `C:/Users/hooki/trench-juce/plugin/source/FunctionGenerator.cpp`
- Test: `C:/Users/hooki/trench-juce/plugin/tests/FunctionGeneratorTests.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Write failing test — silent when no segments loaded**

```cpp
TEST_CASE("FunctionGenerator outputs 0.0 with no segments") {
    FunctionGenerator fg;
    fg.prepare(48000.0);
    // no load() call
    for (int i = 0; i < 100; ++i) REQUIRE(fg.tick() == 0.0f);
}
```

- [ ] **Step 2: Write failing test — linear segment produces correct interpolation**

```cpp
TEST_CASE("FunctionGenerator linear segment interpolates 0 -> 1 over 1s") {
    FunctionGenerator fg;
    fg.prepare(48000.0);
    ModFunctionGenerator gen;
    gen.segments[0] = { .level = 0.0f, .time_ms = 0,    .shape = "hold" };
    gen.segments[1] = { .level = 1.0f, .time_ms = 1000, .shape = "lin"  };
    fg.load(gen);
    fg.noteOn(); // trigger key-sync
    // after 0.5s (24000 samples), output should be ~0.5
    for (int i = 0; i < 24000; ++i) fg.tick();
    REQUIRE(std::abs(fg.tick() - 0.5f) < 0.01f);
}
```

- [ ] **Step 3: Write failing test — conditional jump loops segments**

```cpp
TEST_CASE("FunctionGenerator with jump loops back") {
    FunctionGenerator fg;
    fg.prepare(48000.0);
    ModFunctionGenerator gen;
    // segments 0, 1, 2 with jump on seg 2 -> seg 0
    ...
    fg.load(gen);
    fg.noteOn();
    // verify output matches a looped pattern
}
```

- [ ] **Step 4: Write failing test — REV reverses time axis**

```cpp
TEST_CASE("FunctionGenerator in REV mode outputs time-reversed trajectory") {
    // run forward, record samples
    // reset, enable REV, run
    // final forward sample == first REV sample, etc.
}
```

- [ ] **Step 5: Implement `FunctionGenerator`**

Key methods:
```cpp
class FunctionGenerator {
public:
    void prepare(double sampleRate);
    void load(const ModFunctionGenerator& gen);
    void noteOn(); // key-sync reset
    void setReverse(bool rev);
    void setTempoBpm(double bpm); // for tempo-sync interpretation
    float tick() noexcept; // returns next output sample; 0 when empty
    void reset();
};
```

Internal state: current segment index, position within segment, direction (forward/reverse), sample rate, tempo bpm. No allocations after `load`.

Segment shapes: `hold`, `lin`, `exp`, `log`, `sCurve`. Compute next sample by interpolating between segment start/end levels per the shape function.

- [ ] **Step 6: Run tests, confirm pass**

- [ ] **Step 7: Commit**

```
git commit -am "feat(dsp): add FunctionGenerator E-mu 8-segment modulator"
```

---

## Task 7: Implement `QSoundSpatial` post-cascade stage

**Files:**
- Create: `C:/Users/hooki/trench-juce/plugin/source/QSoundSpatial.h`
- Create: `C:/Users/hooki/trench-juce/plugin/source/QSoundSpatial.cpp`
- Test: `C:/Users/hooki/trench-juce/plugin/tests/NullBaselineTests.cpp` (add cases)
- Modify: `TrenchEngine.h` to populate `SpatialProfile` struct with the fields the stage needs
- Reference: `C:/Users/hooki/Trench/docs/archive/qsound_spatial.md`

- [ ] **Step 1: Read the QSound reference doc**

Open `docs/archive/qsound_spatial.md`. Note the exact ITD law, ILD law, band-law coefficients, and sample-rate assumptions. Populate `SpatialProfile` struct fields accordingly.

- [ ] **Step 2: Write failing test — SPACE = 0 produces bit-identical bypass**

```cpp
TEST_CASE("QSoundSpatial at space=0 is sample-accurate null") {
    QSoundSpatial q;
    q.prepare(48000.0, 2);
    SpatialProfile profile = /* any populated profile */;
    q.setProfile(profile);
    q.setSpace(0.0f);
    // process a stereo buffer; output == input
    juce::AudioBuffer<float> buf(2, 512);
    juce::AudioBuffer<float> ref(2, 512);
    fillNoise(buf); ref.makeCopyOf(buf);
    q.processBlock(buf);
    for (int ch = 0; ch < 2; ++ch)
        for (int i = 0; i < 512; ++i)
            REQUIRE(buf.getSample(ch, i) == ref.getSample(ch, i));
}
```

- [ ] **Step 3: Write failing test — SPACE > 0 yields stereo widening**

Measure L/R correlation; at `SPACE = 1.0`, correlation should be lower than at `SPACE = 0.0` on a mono-in signal.

- [ ] **Step 4: Implement `QSoundSpatial`**

Per spec:
- ITD fractional delay lines (L and R, with interpolation)
- ILD gain law
- Three-band filter law per the doc's `band_law` table
- Bypass short-circuit when `space == 0.0f` — do **not** rely on wet math producing zero; explicit branch.
- Dry/wet linear scaling from `space`.

- [ ] **Step 5: Run tests, confirm pass**

- [ ] **Step 6: Commit**

```
git commit -am "feat(dsp): add QSoundSpatial post-cascade stage"
```

---

## Task 8: Implement `SafetyLimiter` (final output)

**Files:**
- Create: `C:/Users/hooki/trench-juce/plugin/source/SafetyLimiter.h`
- Create: `C:/Users/hooki/trench-juce/plugin/source/SafetyLimiter.cpp`
- Test: append to `NullBaselineTests.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Write failing test — transparent on normal material**

```cpp
TEST_CASE("SafetyLimiter is transparent below threshold") {
    SafetyLimiter lim;
    lim.prepare(48000.0, 2);
    // feed a buffer whose peak is below the limiter ceiling
    // expect bit-identical output
}
```

- [ ] **Step 2: Write failing test — engages on peak runaways**

```cpp
TEST_CASE("SafetyLimiter catches peaks above ceiling") {
    // feed a loud burst; expect output peak <= ceiling + small tolerance
}
```

- [ ] **Step 3: Implement peak limiter**

Simple lookahead-free hard limiter with short release. Ceiling: −0.3 dBFS (spec can refine). Stereo-linked gain reduction.

- [ ] **Step 4: Run tests, confirm pass**

- [ ] **Step 5: Commit**

```
git commit -am "feat(dsp): add SafetyLimiter at output"
```

---

## Task 9: Wire all stages into `processBlock`

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.h` — add member instances
- Modify: `C:/Users/hooki/trench-juce/plugin/source/TrenchEngine.cpp:249..477` — modify sample loop
- Modify: `C:/Users/hooki/trench-juce/plugin/source/PluginProcessor.cpp:364` — read `space`/`trig` params per block, call MIDI note-on hook

- [ ] **Step 1: Add members to `Engine`**

```cpp
MackieDrive      driveL, driveR;
FunctionGenerator modFn;
QSoundSpatial    spatial;
SafetyLimiter    limiter;
```

- [ ] **Step 2: Apply drive in the per-sample loop before the cascade**

Inside the existing sample loop in `TrenchEngine.cpp` (around line 475, where `cascadeL.stages[si].process` is called), pre-process the input through `driveL.processSample` / `driveR.processSample` before the cascade stage chain.

- [ ] **Step 3: Apply TRIG modulation to MORPH**

In the control block (where `targetMorph` is set; currently driven by `apvts` in `PluginProcessor::processBlock`), compute:

```cpp
float fnOut = modFn.tick(); // called once per sample OR once per control block, decide in step 4
float trigAmount = *apvts.getRawParameterValue("trig");
float effectiveMorph = juce::jlimit(0.0f, 1.0f, userMorph + trigAmount * fnOut);
```

Wire the REV state from the existing REV parameter into `modFn.setReverse()`.

- [ ] **Step 4: Decide generator tick rate**

Function generator ticks once per audio sample (simplest, highest resolution). Document the decision in a comment; revisit only if CPU is an issue.

- [ ] **Step 5: Apply SPACE post-cascade**

After the last cascade stage, before output write:

```cpp
float space = *apvts.getRawParameterValue("space");
spatial.setSpace(space);
spatial.processBlock(buffer); // processes in place
```

- [ ] **Step 6: Apply limiter at the very end**

```cpp
limiter.processBlock(buffer);
```

- [ ] **Step 7: Wire MIDI note-on to `modFn.noteOn()`**

In `processBlock(AudioBuffer&, MidiBuffer& midi)`, iterate MIDI messages; on note-on, call `modFn.noteOn()`.

- [ ] **Step 8: Build, run standalone, verify defaults null**

With SPACE=0, TRIG=0, all existing bodies: output should match the current release build bit-for-bit on the same input WAV. Automate this check as a test in step 9.

- [ ] **Step 9: Generate the pre-change golden WAV** *(must happen before starting Task 1's CMake changes to guarantee the fixture predates any code change — if not done yet, check out the commit BEFORE this plan's first commit, build, render, commit the WAV, then return to this task)*

From a clean checkout of the commit before Task 1 began:
```
git worktree add /tmp/pre-space-trig <pre-task-1 commit sha>
cd /tmp/pre-space-trig/trench-juce/plugin && cmake -S . -B build && cmake --build build --target TRENCH_Standalone
# Render a fixed input WAV (tests/fixtures/noise_burst_48k.wav) through the baseline body
# Save output as tests/fixtures/pre-space-trig/<body>.wav
```

Commit the fixture to the current branch under `tests/fixtures/pre-space-trig/`.

- [ ] **Step 10: Write integration test — full-chain null at all-zero**

```cpp
TEST_CASE("Full chain: SPACE=0, TRIG=0, existing body produces identical output") {
    // Render tests/fixtures/noise_burst_48k.wav through the plugin with baseline body.
    // Compare sample-by-sample against tests/fixtures/pre-space-trig/<body>.wav.
    // REQUIRE sample-accurate equality.
}
```

- [ ] **Step 11: Run all tests, commit**

```
git commit -am "feat: wire MackieDrive, FunctionGenerator, QSoundSpatial, SafetyLimiter into processBlock"
```

---

## Task 10: Author 4 shipping bodies (hand-written JSON)

**Files:**
- Modify: `C:/Users/hooki/trench-juce/plugin/assets/<body>.json` × 4

Bodies to author (names are placeholders; confirm with user):
1. Aluminum Siding
2. Speaker Knockerz
3. Small Talk Ah-Ee
4. Cul-De-Sac

For each body:

- [ ] **Step 1: Read the existing cartridge JSON**

Identify which of the 8 embedded cartridges becomes each body, or create a new file per body.

- [ ] **Step 2: Add `drive` block**

Start with `input_gain_dB` values authored by ear — begin at 0.0 for the clean body (Aluminum Siding), escalate for hotter bodies. Values are editable later; the schema is the contract.

- [ ] **Step 3: Add `mod_fn` block**

Hand-write 1–4 segments per body. Default `key-sync: 1`. For looping bodies, set a `jump` on the last active segment. For one-shot bodies, leave jumps absent.

- [ ] **Step 4: Add `spatial_profile` block**

If `docs/archive/qsound_spatial.md` contains per-body or per-character example coefficient sets, copy the one that matches the body's intended spatial voice. If the doc only documents the laws and parametric model (not per-body presets), hand-author azimuth, distance, elevation, and per-band coefficients by ear: start from the doc's default/center values and audition in the standalone build. Commit the values as authored; refinement is cheap later. Do NOT substitute identity/placeholder — SPACE=1 should audibly differ from SPACE=0 on an authored body.

- [ ] **Step 5: Build, load body, audition**

Run standalone; load the body; play; verify SPACE, TRIG, REV, ENV all behave per spec with authored content.

- [ ] **Step 6: Commit per body**

```
git commit -am "content: author <body name> with drive + mod_fn + spatial"
```

---

## Task 11: Final verification and doctrine check

- [ ] **Step 1: Run `./check` at repo root**

Run the project's verification harness per `CLAUDE.md`:

```
cd C:/Users/hooki/Trench && ./check
```

Expected: all checks pass.

- [ ] **Step 2: Run Rust reference tests**

```
cd C:/Users/hooki/Trench/trench-core && cargo test
```

Expected: all tests pass. If the Rust cascade tests regressed because of schema changes on the cartridge loader's Rust side, update `cartridge.rs` accordingly. (Per the scope decision, this plan keeps DSP changes C++-only; schema additions to the Rust loader are documentation-only unless the Rust tests load cartridges that now carry the new blocks.)

- [ ] **Step 3: Spectrum trace sanity**

Launch standalone with each authored body. Confirm:
- SPACE at 0 = null vs pre-change build.
- TRIG at 0 = null vs pre-change build.
- SPACE > 0 widens stereo.
- TRIG > 0 moves MORPH on each note-on, with phase reset when `key-sync=1`.
- REV reverses the TRIG trajectory.
- Authored hot bodies audibly bite more than clean bodies (drive verification).

- [ ] **Step 4: Audit against spec verification section**

Walk through spec §Verification checklist; confirm each item is covered by a test or an explicit manual check. Any missed items become followup tasks.

- [ ] **Step 5: Final commit**

```
git commit --allow-empty -m "ship: SPACE + TRIG + Mackie drive complete per spec 2026-04-18"
```

---

## Open questions held for plan execution

These are not code issues; they gate the plan but are not re-derivable from the spec:

1. **Mackie model source** (Task 5 blocker): Cytomic reference vs measured 1202 curve. Escalate to user if neither is available.
2. **Mackie curve-match tolerance + harmonic threshold** (Task 5 decision): once model source is picked, record both values as constants in `MackieDrive.cpp`.
3. **Golden WAV provenance** (Task 9 Step 9): the pre-change reference for null tests. Generate from the commit BEFORE Task 1 began, not from a dirty working tree. Commit under `tests/fixtures/pre-space-trig/`.
4. **Shipping body name→cartridge mapping** (Task 10): confirm the 4 names (Aluminum Siding, Speaker Knockerz, Small Talk Ah-Ee, Cul-De-Sac) map to specific embedded cartridge files before authoring.
5. **Per-body spatial coefficient availability** (Task 10 Step 4): verify whether `docs/archive/qsound_spatial.md` contains per-body presets or only parametric laws; if only laws, bodies are hand-authored by ear from center defaults.

## Advisory items (not blocking; address opportunistically during execution)

- **Task 3 roller placement**: consider invoking the `juce-lookandfeel-design` skill for disciplined pixel coords; at minimum commit chosen rectangles as named constants, not magic numbers.
- **Task 6 shape coverage**: the plan tests `lin`, `jump`, and `REV`. Add a parametric test covering `hold`, `exp`, `log`, `sCurve` shapes before shipping.
- **Task 9 generator tick rate**: spec tied to per-sample is simplest; a control-block (32-sample) tick rate is a cheap CPU optimization if profiling ever demands it. Not a v1 concern.
- **Task 11 Rust-side schema**: if any Rust test fixture JSON gains the new optional blocks, update `trench-core/src/cartridge.rs` to parse them (or skip unknown fields gracefully) rather than letting Rust tests fail.
