# TRENCH authoring model

A TRENCH pill is a baked E-mu Morph Designer preset: a single Z-plane
filter with two frames (Lo Morph / Hi Morph) interpolated via the MORPH
parameter. See [`docs/emu/`](docs/emu/) for the canonical E-mu reference
material and [`docs/architecture/zplane_truth.md`](docs/architecture/zplane_truth.md)
for TRENCH's relationship to it.

## Phoneme metaphor (E-mu origin)

The vowel / vocal-tract framing is E-mu's own language, not TRENCH's
invention. From the E-mu manual (archived at
[`docs/emu/zplane_explained.md`](docs/emu/zplane_explained.md)):

> A vowel is really a configuration of many muscles, but we consider
> it a single object. In changing from one vowel to another, we don't
> need to consider the frequencies of the resonant peaks. You remember
> the shape of your mouth for each sound and interpolate between them.

TRENCH's phoneme pill vocabulary is a faithful translation of that
framing into named tokens.

## Token vocabulary

Tokens come from the unified phoneme inventory at
`cartridges/engine/_source/token_inventory_unified_v2.json`. Each
token is a named Morph Designer preset. A token name is a stable
string; unknown token names fail validation.

The shipping flat pill layout at `cartridges/engine/<category>/<key>.json`
is a derived rebuild, produced by `tools/bake_phoneme_pills.py` from
the `_source/` subtree. Do not hand-edit the flat layout — edit
`_source/` and re-bake.

## Pill format

Pills ship as `compiled-v1` cartridges (see `SPEC.md` §2). Today the
format bakes a 4-corner shape (`M0_Q0`, `M0_Q100`, `M100_Q0`,
`M100_Q100`), which is E-mu's 2-frame Morph Designer with an added
Q-precomputation axis — a TRENCH engineering optimization over the
native 2-frame format. The long-term direction is 2-corner with live Q
computation; see
[`docs/architecture/zplane_truth.md`](docs/architecture/zplane_truth.md)
for rationale and transition plan.

## Playback model

There is no TRENCH-specific compiler. At playback, a host (the JUCE
plugin or `trench-live`) loads a pill cartridge and exposes MORPH and
Q as DAW-modulatable parameters. Users route mod sources — LFO,
envelope, velocity, mod wheel, DAW automation — to MORPH and Q the
same way they would on a hardware E-mu synth. E-mu's own manual is
explicit: *"any of the modulation sources can control the Z-plane
filter."*

A BODY is a named collection of pills sharing a design intent (see
[`BODIES.md`](BODIES.md)). Pills within a body are selected or switched
via the host's preset / pill-picker UI. There is no TRENCH-specific
trajectory compiler or time-grid scheduler.

## Authoring

Pills are authored offline via one of:

- **Direct pole-zero placement** — `tools/author_speaker_knockerz.py`
  is the reference. Hand-authored stages using Klatt-normalized
  resonators, corrected notches, and 1-pole lowpass/highpass primitives.
- **Extraction from E-mu factory filters** — `tools/extract_emu_filter_params.py`
  pulls 30-integer authoring grids from heritage bundles;
  `tools/bake_hedz_const.py` compiles them to coefficients via the
  E-mu type 1/2/3 firmware recipes.
- **Calibration-driven reconstruction** — `docs/calibration/index.json`
  plus per-filter calibration JSONs (BassBox 303, Ear Bender, Early
  Rizer, Talking Hedz, etc.) define reconstruction strategies for the
  named E-mu filters in the Emulator X3 bank.

Candidates are hand-audited against [`BODIES.md`](BODIES.md) rubrics
and the `trench-core/tests/vowel_formants.rs` classification gate.

## Constraints on authoring tools

All authoring tools MUST:

- Produce cartridges that validate against `cartridge.schema.json`
  (see `SPEC.md` §2).
- Preserve per-body invariants from [`BODIES.md`](BODIES.md) for
  shipping candidates.

All authoring tools MUST NOT:

- Expose pole/zero/frequency grids to the user as a first-class
  authoring surface.
- Accept raw coefficient edits.
- Derive stages from RBJ cookbook shapes. Direct pole-zero only (see
  [`DOCTRINE.md`](DOCTRINE.md)).
