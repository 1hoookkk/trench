# Z-plane truth and TRENCH architecture direction

This document pins TRENCH's architecture to E-mu Z-plane ground truth.
It exists to stop three different things from collapsing into one lie:

- heritage authoring truth
- compiled engine representations
- runtime playback state

If these are not kept separate, the project drifts. Cache formats start
pretending to be source material, implementation shortcuts start
pretending to be authoring law, and future reverse-engineering work gets
poisoned.

## Ground truth

The canonical E-mu references live at:

- [`docs/emu/zplane_explained.md`](../emu/zplane_explained.md) — Z-plane
  definition: two complex filter frames interpolated by a Morph axis.
- [`docs/emu/emulator_x3_filter_reference.md`](../emu/emulator_x3_filter_reference.md)
  — the 55-filter bank listing, Peak/Shelf Morph description, and Morph
  Designer specification (*Emulator X3 Reference Manual*, Chapter 6).
- [`docs/emu/dillusionman_peak_shelf_morph.md`](../emu/dillusionman_peak_shelf_morph.md)
  — practical Peak/Shelf Morph vocabulary: `FREQ / SHELF / PEAK` per
  frame plus the `FilFreq` / `FilRes` control usage.
- [`docs/emu/peak_shelf_morph_reece_recipe.md`](../emu/peak_shelf_morph_reece_recipe.md)
  — worked Peak/Shelf Morph example.

Any TRENCH design decision that is not either:

1. provably better than E-mu, or
2. a faithful translation of E-mu truth

needs justification or removal.

## Provenance law

Not all references are equal.

### Canonical

These are allowed to define architecture:

- archived/manual E-mu documentation
- direct clean-room extraction artifacts
- byte-level parity references already validated in repo

### Secondary

These are useful, but cannot outrank canonical sources:

- NotebookLM exports
- NotebookLM summaries
- community writeups
- forum lore
- workflow notes
- interpretive DSP commentary

NotebookLM material is allowed to:

- corroborate a control law already present in canonical docs
- preserve phrasing or practical usage patterns
- point to likely extraction targets

NotebookLM material is **not** allowed to:

- redefine the native authoring object
- override a manual/reference statement
- silently introduce architecture law
- replace raw provenance with paraphrase

If a NotebookLM note conflicts with canonical docs, the NotebookLM note
loses.

## Classification framing: preserve / translate / invent

For every TRENCH concept, classify it into one of three buckets:

- **(A) Better than E-mu** — modern engineering upgrade with a stated
  reason for why it beats the original. Keep.
- **(B) E-mu truth, preserved or translated faithfully** — keep as-is,
  document provenance.
- **(C) Our invention, not yet justified** — needs a defense against
  the E-mu reference or deletion.

### Current classification

| Concept | Bucket | Notes |
|---|---|---|
| Pill = Morph Designer preset | **B** | Direct translation of E-mu's 6-section × 2-frame structure. |
| Phoneme vocabulary / vowel naming | **B** | E-mu's own vocal-tract metaphor; see `zplane_explained.md`. |
| 30-integer heritage authoring grid | **B** | Morph Designer authoring format, preserved byte-for-byte at `cartridges/engine/_source/heritage_designer_sections.json`. |
| Peak/Shelf authoring as left/right frame triplets | **B** | Direct translation of E-mu's `FREQ / SHELF / PEAK` per-frame object. |
| Type 1/2/3 firmware compile recipes | **B** | E-mu's real coefficient derivation, inlined in `tools/bake_hedz_const.py`. |
| 12-stage DF2T serial cascade (6 active + 6 passthrough) | **B** | Topology matches E-mu's 12th-order section chain. |
| Per-filter calibration skins | **B** | Reverse-engineered real E-mu filter strategies. |
| `f64` math + 32-sample control blocks + ramping | **A** | Modern precision and artifact control E-mu hardware could not afford. |
| 4-corner (`M0_Q0`/`M0_Q100`/`M100_Q0`/`M100_Q100`) bilinear shape | **A** constrained | Valid as a compiled/runtime representation only. Not native authoring truth. |
| Audio-thread doctrine (no locks, no allocation) | **A** | Modern best practice E-mu firmware did not need to encode. |
| `compiled-v1` cartridge JSON format | **A** transitional | Wire-safe and portable, but transitional. Must not replace heritage truth. |
| Compiler (pill sequence → trajectory) | **C dead** | E-mu has no compiler; MORPH is driven by mod sources. Remove. |
| Looperator time grid as authoring surface | **C dead** | DAW automation and synth mod sources are the surface. |
| BODIES as scheduled pill sequences | **C reframed** | A body is a named collection sharing a design intent, not a scheduled trajectory. |
| "Pills do not morph internally" rule | **C dead** | Contradicts the definition of a Z-plane filter. A pill **is** the A/B frame pair. |

## Non-negotiable architecture law

There are exactly **three layers**.

### 1. Heritage layer

This is the only canonical source of truth.

Use E-mu-native semantics only.

Examples:

- Peak/Shelf Morph as two semantic frames:
  - left `{freq_hz, shelf, peak_db}`
  - right `{freq_hz, shelf, peak_db}`
- Morph Designer as six sections, each with:
  - section shape/type
  - Lo frame `{freq_hz, gain_or_q}`
  - Hi frame `{freq_hz, gain_or_q}`

### 2. Compiled layer

This is a deterministic lowering from heritage into engine-ready data.

Allowed contents:

- 2-corner or 4-corner packs
- pole/zero packs
- biquad coefficient packs
- normalized convenience values
- stage enable flags
- interpolation metadata

This layer may optimize. It may not invent meaning.

### 3. Runtime layer

This is ephemeral playback state.

Allowed contents:

- current morph
- current resonance / gain drive
- smoothed control values
- active coefficient sets
- meter / telemetry state
- temporary buffers

Runtime state is never canonical and never hand-authored.

## Hard invariants

### Heritage layer must never contain

- `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`
- biquad coefficients
- SIMD packs
- interpolation caches
- smoothed state
- UI pixel positions
- normalized convenience values that replace source units
- anything that exists only because the engine is faster that way

### Compiled layer may contain

- corners
- precomputed coeffs
- stage packs
- interpolation tables
- normalized values
- engine-specific convenience data

But it must declare:

- what heritage object it came from
- what compiler produced it
- whether the lowering is exact or approximate

### Runtime layer may contain

- live morph state
- live resonance / peak / Q-offset state
- smoothed values
- active coeffs
- telemetry

But runtime state is disposable. It is not preset truth.

## Why four corners are not the ontology

TRENCH currently ships 4 bilinear corners:

- `M0_Q0`
- `M0_Q100`
- `M100_Q0`
- `M100_Q100`

Each corner stores 12 biquad stages × 5 kernel-form coefficients. At
runtime, `Cartridge::interpolate` does Q-first-then-Morph bilinear
interpolation to produce the cascade's coefficient set.

This shape is **not** E-mu's native authoring format.
It is a compiled/runtime representation.

### What E-mu actually authors

#### Peak/Shelf Morph

Peak/Shelf Morph is a 2-frame filter object. Each frame has:

- `FREQ`
- `SHELF`
- `PEAK`

The important law is that `FREQ` does **not** have a stable meaning by
itself. Its effect depends on the `SHELF` value.

So the canonical authoring atom is not frequency alone. It is the frame
triplet:

- `{freq_hz, shelf, peak_db}`

left and right.

#### Morph Designer

Morph Designer stores 6 sections × 2 frames ×:

- section type / shape
- frequency
- Q or Gain

Again: two frames, not four corners.

### Why four-corner thinking is dangerous

If four corners become the primary mental model, the project silently
starts assuming that:

- Q is a second authored axis
- Morph and Q form a native rectangular surface
- bilinear interpolation is the meaning of the instrument

That is implementation drift.

The correct discipline is:

- **author in frames**
- **compile to corners if needed**
- **run with whatever representation is fastest**

Never confuse the cache with the instrument.

## Peak/Shelf control law

Peak/Shelf Morph should be understood as:

- a left frame
- a right frame
- a live Morph sweep between them
- a live resonance / peak drive control acting during playback

Operationally:

- `FilFreq` is the live control that sweeps the Morph position between
  the two programmed frames.
- `FilRes` is the live control used to drive / saturate / intensify the
  filter behavior during the sweep.

This is why a Peak/Shelf preset must not be stored as a frequency list
or as four corners. The actual authored object is the two-frame posture.

## Source-of-truth schemas

### Peak/Shelf heritage schema

```json
{
  "schema": "trench.heritage.peak_shelf.v1",
  "id": "reece_stab_01",
  "name": "Reece Stab",
  "source": {
    "origin": "emu_x3",
    "family": "Peak/Shelf Morph",
    "reference": "docs/emu/dillusionman_peak_shelf_morph.md"
  },
  "left": {
    "freq_hz": 246.0,
    "shelf": -50,
    "peak_db": -24.0
  },
  "right": {
    "freq_hz": 4488.0,
    "shelf": 30,
    "peak_db": 1.5
  },
  "controls": {
    "morph_domain": [0, 255],
    "morph_driver": "FilFreq",
    "resonance_driver": "FilRes"
  },
  "notes": "Heritage object only. No derived corners."
}
```

#### Rules

- `freq_hz` is stored in Hz, not normalized.
- `shelf` preserves source semantics as a signed scalar.
- `peak_db` is stored in dB, not linear gain.
- `left` and `right` are semantic frames, not corner aliases.
- Do not add `q0`, `q100`, `corners`, or `coeffs` here.

### Morph Designer heritage schema

```json
{
  "schema": "trench.heritage.morph_designer.v1",
  "id": "talkinghedz_md_01",
  "name": "TalkingHedz",
  "source": {
    "origin": "emu_x3",
    "family": "Morph Designer",
    "reference": "docs/emu/emulator_x3_filter_reference.md"
  },
  "sections": [
    {
      "index": 0,
      "shape": "lowpass",
      "lo": { "freq_hz": 320.0, "gain_or_q": 38.0 },
      "hi": { "freq_hz": 710.0, "gain_or_q": 44.0 }
    },
    {
      "index": 1,
      "shape": "eq",
      "lo": { "freq_hz": 900.0, "gain_or_q": 12.0 },
      "hi": { "freq_hz": 1250.0, "gain_or_q": 18.0 }
    },
    {
      "index": 2,
      "shape": "eq",
      "lo": { "freq_hz": 1800.0, "gain_or_q": 10.0 },
      "hi": { "freq_hz": 2300.0, "gain_or_q": 15.0 }
    },
    {
      "index": 3,
      "shape": "off",
      "lo": { "freq_hz": 0.0, "gain_or_q": 0.0 },
      "hi": { "freq_hz": 0.0, "gain_or_q": 0.0 }
    },
    {
      "index": 4,
      "shape": "off",
      "lo": { "freq_hz": 0.0, "gain_or_q": 0.0 },
      "hi": { "freq_hz": 0.0, "gain_or_q": 0.0 }
    },
    {
      "index": 5,
      "shape": "highpass",
      "lo": { "freq_hz": 120.0, "gain_or_q": 20.0 },
      "hi": { "freq_hz": 260.0, "gain_or_q": 24.0 }
    }
  ],
  "global_control_model": {
    "morph_axis": "lo_to_hi",
    "big_wheel_role": "global_gain_or_q_offset"
  }
}
```

#### Rules

- `shape` enum is heritage-facing:
  - `off`
  - `eq`
  - `lowpass`
  - `highpass`
- `gain_or_q` preserves heritage ambiguity on purpose.
- Do not split it into engine-specific subfields in heritage.
- Store six sections even if some are `off`.

## Compiled-layer examples

The compiled layer is the only place where corners are allowed.

### Peak/Shelf compiled example

```json
{
  "schema": "trench.compiled.peak_shelf_surface.v1",
  "id": "reece_stab_01.compiled",
  "derived_from": "heritage/peak_shelf/reece_stab_01.heritage.json",
  "compiler": {
    "name": "peak_shelf_compiler",
    "version": "1.0.0"
  },
  "semantic_family": "Peak/Shelf Morph",
  "representation": "4corner_df2t",
  "exactness": "approximate",
  "derivation_notes": {
    "morph_axis": "heritage_left_to_right",
    "q_axis": "derived_runtime_axis_not_heritage_authored"
  },
  "nodes": {
    "m0_q0": { "stages": [] },
    "m1_q0": { "stages": [] },
    "m0_q1": { "stages": [] },
    "m1_q1": { "stages": [] }
  }
}
```

### Morph Designer compiled example

```json
{
  "schema": "trench.compiled.morph_designer.v1",
  "id": "talkinghedz_md_01.compiled",
  "derived_from": "heritage/morph_designer/talkinghedz_md_01.heritage.json",
  "compiler": {
    "name": "designer_compile",
    "version": "1.0.0"
  },
  "representation": "6stage_df2t_pack",
  "exactness": "approximate"
}
```

## Runtime layer contract

Runtime state is transient and disposable.

Use structs, not archival JSON, as the primary mental model.

```cpp
struct PeakShelfRuntimeState {
    float morph_norm;
    float res_norm;
    float morph_smoothed;
    float res_smoothed;
    StageCoeffs active[6];
};

struct MorphDesignerRuntimeState {
    float morph_norm;
    float global_qgain_norm;
    float morph_smoothed;
    float qgain_smoothed;
    StageCoeffs active[6];
};
```

## Transformation law

### Allowed

- `heritage -> compiled`
- `compiled -> runtime`
- `heritage -> UI labels / browser metadata`

### Forbidden

- `runtime -> heritage`
- `compiled -> heritage` as a normal authoring path
- hand-editing compiled blobs as if they were presets
- emitting heritage docs from runtime caches

If reverse recovery is ever needed, that is a forensic tool, not the
normal pipeline.

## The 4-corner -> 2-corner transition

### Current state

TRENCH cartridges currently ship with 4 bilinear corners:

- `M0_Q0`
- `M0_Q100`
- `M100_Q0`
- `M100_Q100`

Each corner stores 12 biquad stages × 5 kernel-form coefficients.

This is valid as a shipping/runtime optimization.
It is not E-mu's native source format.

### Why 2-corner + live control is the better runtime

1. **E-mu-faithful.** The runtime matches the documented authoring
   object more directly.
2. **Shipping format can match authoring format.** The heritage object
   stays legible instead of being hidden behind derived corners.
3. **Q / Peak is live control, not a fake second authored axis.**
4. **CPU cost is negligible on modern hardware.**
5. **Authoring mental model shrinks.** Two frames per section, not four
   corners.

### What must not change without a new decision

- The DSP cascade topology stays frozen unless separately justified.
- The parity gate stays green on every commit through any transition.
- Existing shipped pills remain byte-preserved until re-derived from
  their authoring sources.

## One-sentence doctrine

**Heritage defines meaning. Compiled defines machinery. Runtime defines motion.**

## Related docs

- [`docs/emu/zplane_explained.md`](../emu/zplane_explained.md)
- [`docs/emu/emulator_x3_filter_reference.md`](../emu/emulator_x3_filter_reference.md)
- [`docs/emu/dillusionman_peak_shelf_morph.md`](../emu/dillusionman_peak_shelf_morph.md)
- [`docs/emu/peak_shelf_morph_reece_recipe.md`](../emu/peak_shelf_morph_reece_recipe.md)
- [`SPEC.md`](../../SPEC.md)
- [`PHONEMES.md`](../../PHONEMES.md)
- [`BODIES.md`](../../BODIES.md)
