# Z-plane truth and TRENCH architecture direction

This document pins TRENCH's architecture to E-mu Z-plane filter ground
truth. It exists to prevent abstract inventions from drifting the
runtime away from E-mu's actual design.

## Ground truth

The canonical E-mu references live at:

- [`docs/emu/zplane_explained.md`](../emu/zplane_explained.md) — the
  Z-plane filter definition and morph interpolation behavior (Proteus
  2000-era manual).
- [`docs/emu/emulator_x3_filter_reference.md`](../emu/emulator_x3_filter_reference.md)
  — the 55-filter bank listing and Morph Designer specification
  (*Emulator X3 Reference Manual*, Chapter 6).

Any TRENCH design decision that is not either (a) provably better than
E-mu or (b) a faithful translation of E-mu truth needs justification or
removal.

## Classification framing: preserve / translate / invent

For every TRENCH concept, classify into one of three buckets:

- **(A) Better than E-mu** — modern engineering upgrade, with a stated
  reason for why it beats the original. Keep.
- **(B) E-mu truth, preserved or translated faithfully** — keep as-is,
  document provenance.
- **(C) Our invention, not yet justified** — needs a defense against
  the E-mu reference or deletion.

### Current classification

| Concept                                             | Bucket    | Notes |
|------------------------------------------------------|-----------|-------|
| Pill = Morph Designer preset                         | **B**     | Direct translation of E-mu's 6-section × 2-frame structure. |
| Phoneme vocabulary / vowel naming                    | **B**     | E-mu's own vocal-tract metaphor; see `zplane_explained.md`. |
| 30-integer heritage authoring grid                   | **B**     | The Morph Designer authoring format, preserved byte-for-byte at `cartridges/engine/_source/heritage_designer_sections.json`. |
| Type 1/2/3 firmware compile recipes                  | **B**     | E-mu's real coefficient derivation, inlined in `tools/bake_hedz_const.py`. |
| 12-stage DF2T serial cascade (6 active + 6 passthrough) | **B**  | Topology matches E-mu's 12th-order section chain. |
| Per-filter calibration skins                         | **B**     | `docs/calibration/index.json` + per-filter JSONs reverse-engineer real E-mu filter strategies (FREQUENCY_SWAP, HF_WALL_FOLD, VOWEL_FORMANT, RADIUS_ONLY_MORPH, etc.). |
| `f64` math + 32-sample control blocks + ramping      | **A**     | Modern precision and artifact control E-mu hardware couldn't afford. |
| 4-corner (`M0_Q0`/`M0_Q100`/`M100_Q0`/`M100_Q100`) bilinear shape | **A** → transition to **B** | Precomputes Q variants per frame to avoid runtime biquad computation. Justified for current runtime; transitioning to 2-corner + live Q (see below). |
| Audio-thread doctrine (no locks, no allocation)      | **A**     | Modern best practice E-mu firmware didn't need to encode. |
| `compiled-v1` cartridge JSON format                  | **A**     | Wire-safe, schema-validated, host-portable. Transitional; subsumed by 30-integer grid under the 2-corner direction. |
| Compiler (pill sequence → trajectory)                | **C dead**| E-mu has no compiler; MORPH is driven by mod sources. Concept removed from roadmap. |
| Looperator time grid as authoring surface            | **C dead**| DAW automation and synth mod sources are the surface. |
| BODIES as scheduled pill sequences                   | **C reframed** | A body is a named collection sharing a design intent (`BODIES.md`), not a scheduled trajectory. |
| "Pills do not morph internally" rule                 | **C dead**| Contradicts the definition of Z-plane filter. A pill IS the A/B frame pair. |
| `docs/hostile_authoring_workflow_spec.md`            | off-path  | Content/marketing concept (deadpan absurdist ritual framing), not engineering direction. |

## The 4-corner → 2-corner transition

### Current state

TRENCH cartridges ship with 4 bilinear corners: `M0_Q0`, `M0_Q100`,
`M100_Q0`, `M100_Q100`. Each corner stores 12 biquad stages × 5
kernel-form coefficients. At runtime, `Cartridge::interpolate` does
Q-first-then-Morph bilinear interpolation (`SPEC.md` §1) to produce the
cascade's coefficient set.

This shape is NOT E-mu's native format. E-mu's Morph Designer stores 6
sections × 2 frames × (section type, frequency, Q or Gain). At runtime,
Morph Designer interpolates section parameters between the two frames
and derives biquad coefficients live from the interpolated
`(type, freq, Q)` tuple, with a live Q Offset wheel (±50%) added on top.

TRENCH's 4-corner shape is an engineering optimization that precomputes
two Q snapshots per morph frame to avoid runtime biquad coefficient
computation. On the hardware E-mu shipped this on, runtime coefficient
derivation was expensive; on modern CPUs it is free.

### Why 2-corner is the better runtime

1. **E-mu-faithful.** Bucket (B) purity. The runtime matches the
   reference bit for bit.
2. **Shipping format equals authoring format.** The 30-integer heritage
   grid is the cartridge. No intermediate `compiled-v1` translation
   layer.
3. **Q is a live parameter** with a real offset wheel (±50%), not two
   precomputed snapshots. Full Q range, no interpolation artifacts
   between snapshot endpoints.
4. **CPU cost is negligible on modern hardware.** Six biquad coefficient
   computations per 32-sample control block ≈ 0.03% of one modern core.
5. **Authoring mental model shrinks.** Two frames per section, not four
   corners; no Q-axis confusion for morph-only filters.

### What it costs

Transition plan, in dependency order. Items marked ⚠ are parity-critical
— the `hedz_cascade` test must be green on every commit.

1. **⚠ Port E-mu type 1/2/3 firmware recipes from Python to Rust.**
   Currently inlined in `tools/bake_hedz_const.py`. New module
   proposed: `trench-core/src/emu_compile.rs`. Input:
   `(section_type, freq, q_or_gain)`. Output: 5 kernel-form biquad
   coefficients. Must be byte-for-byte against the Python reference.
2. **Define the 2-corner cartridge schema.** 30-integer grid plus
   metadata (section-type enum, Q-offset range, global boost, enable
   flags). Replaces `compiled-v1`.
3. **Rewrite `trench-core/src/cartridge.rs`** to load and interpolate
   the new schema. `Cartridge::interpolate` returns `(type, freq, Q)`
   per section; the cascade calls `emu_compile` from there.
4. **Update `trench-core/src/cascade.rs`** control-block path to call
   `emu_compile` once per section per control block. DF2T math
   untouched.
5. **⚠ Rewrite `trench-core/tests/hedz_cascade.rs`** for the new shape.
   Must stay green on every commit through the transition. Lose this
   gate = lose the byte-for-byte lock on E-mu truth.
6. **Regenerate the committed pills** from their authoring sources.
   Mechanical once steps 1–5 are green.
7. **Update authoring tools** (`tools/author_speaker_knockerz.py`,
   `tools/bake_phoneme_pills.py`, `tools/bake_hedz_const.py`) to
   target the new format.
8. **Update the JUCE plugin cartridge loader** (sibling repo). Same
   parse, same interpolate, same cascade — different struct shape.
9. **Retire the old `compiled-v1` format** once trench-core is stable
   on the new shape.

Steps 1–2 can happen on a scratch branch without disturbing anything
shipping. Steps 3–5 form a single landing. Step 6 is mechanical. Steps
7–9 fall out afterward.

The hardest risk is subtle drift in the type 1/2/3 recipe port. The
`hedz_cascade` parity gate is the enforcement mechanism.

## HTML prototype (first step before committing)

Before starting the Rust port, build a pure-JavaScript prototype of the
2-corner + live-Q architecture and run it in a browser. The E-mu
compile math is a pure function; porting to JS is a direct translation
of `tools/bake_hedz_const.py`. The DF2T cascade is similarly tight
(~20 lines of JS).

**Deliverable:** `demo/emu_zplane/{index.html, zplane.js, plot.js, heritage.json, styles.css}`.
Pure stdlib JS, no frameworks, no WASM, no libraries. Offline audio
rendering via `AudioBuffer` + `AudioBufferSourceNode`; frequency
response via Canvas 2D. Single-screen Bassbox-303 aesthetic, works on
mobile browsers.

**Why this is the right first step:**

- **Proves the architecture cheaply.** If 2-corner + live-Q sounds
  right in the browser, the Rust port has justification. If it
  doesn't, the lesson is learned for the cost of one session instead
  of a multi-week crate rewrite.
- **Becomes the reference implementation** for the Rust port's parity
  test. Same input → same output, byte-for-byte.
- **Runs on a phone.** Every pill in the heritage inventory is
  auditionable immediately.
- **Zero risk to `trench-core`** until the architecture decision is
  committed.

## What must not change without a new decision

- The DSP cascade topology (12-stage serial DF2T, 6 active + 6
  passthrough sentinel) is frozen at `trench-core/src/cascade.rs` per
  `SPEC.md` §3. The 4-corner → 2-corner transition changes the
  cartridge format and interpolation path; it does NOT change the
  cascade.
- The `hedz_cascade` parity test must stay green on every commit
  through the transition. This is the load-bearing guarantee that the
  runtime still sounds like E-mu.
- Pills authored under the current 4-corner shape are byte-preserved
  until they are re-derived from their authoring sources into the
  2-corner format. No truncation, no lossy conversion.

## Related docs

- [`docs/emu/zplane_explained.md`](../emu/zplane_explained.md)
- [`docs/emu/emulator_x3_filter_reference.md`](../emu/emulator_x3_filter_reference.md)
- [`SPEC.md`](../../SPEC.md) §1 (DSP engine), §2 (cartridge format), §3 (runtime invariants)
- [`PHONEMES.md`](../../PHONEMES.md) (authoring model)
- [`BODIES.md`](../../BODIES.md) (shipping bodies)
- [`docs/calibration/index.json`](../calibration/index.json) (E-mu filter reconstruction strategies)
