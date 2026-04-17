# SPACE + TRIG — v1 UI and DSP Addition

**Date:** 2026-04-18
**Status:** Design, not yet implemented
**Scope:** v1.0 shipping plugin (2026-07-15)

## Intent

Add two musical gestures to TRENCH v1 without expanding the faceplate beyond
its "Musical Filter" identity:

1. **SPACE** — user-scaled QSound spatial, baked per body, off by default.
2. **TRIG** — user-scaled modulation depth from a per-body E-mu-style LFO to
   MORPH, off by default.

Both behave like character dials: zero is the null baseline, turn it up to
bring in the body's authored behavior. Nothing global, nothing synth-matrix
shaped.

## Doctrine anchor: follow E-mu

E-mu reference material (`docs/emu/`) and E-mu X3 templates define the native
modulation vocabulary. TRENCH does not invent a new modulation primitive.

- MORPH is a **modulation destination** (per `emulator_x3_filter_reference.md`
  and `zplane_explained.md`).
- Modulation **sources** follow the E-mu LFO / TwistaLoop primitive shape.
- A `PatchCord` routes source → destination with a scalar amount. TRENCH v1
  exposes one such cord on the faceplate: **TRIG** = cord amount to MORPH.

Reference templates:

- LFO template `TwistaLoop/twista.xml` — fields `speed`, `loop-select`,
  `loop-start`. Sample-based loopable oscillator used as a mod source.
- LFO template `LFO/Clocked BPM.xml` — fields `frequency`, `shape`, `delay`,
  `var`, `key-sync`, `tempo-sync`. Parametric LFO with host-tempo sync and
  per-note restart.

v1 adopts the **Clocked BPM LFO** field set as the body-level modulator.
TwistaLoop-shape authoring is out of scope for v1 and reserved for a later
milestone.

## Non-goals

- No mod matrix UI (no source pickers, no destination pickers, no per-pair
  depth grid on the face).
- No LFO rate / shape / delay knobs on the face. All LFO parameters are
  authored into the body.
- No velocity / aftertouch / mod-wheel direct parameter routing in v1 beyond
  the existing ENV→Q path.
- No tempo sync controls on the face; `tempo-sync` is a body field.
- No TRIG arm button. No SPACE enable button. Knobs at zero = off.
- No spatial azimuth / distance user controls. Spatial parameters live in the
  body.
- No changes to the frozen 12-stage DF2T cascade.
- No TwistaLoop (sample-based) modulator shape in v1.

## Final v1 faceplate surface

```
      ┌──────── BODY ────────┐
      │  ●  ●  ●  ●          │   4-body selector
      └──────────────────────┘

         MORPH          Q
          ◉             ◉

         SPACE         TRIG
          ◉             ◉

              REV   ENV
```

Two new rollers (SPACE, TRIG). Zero new buttons. REV and ENV exist; ENV keeps
its current envelope→Q routing, REV gains new semantics (see TRIG).

## SPACE

### Behavior

- Normalized plugin parameter, range `[0.0, 1.0]`, default `0.0`.
- At `0.0` the audio path short-circuits around the spatial stage; output is
  bit-identical to the current cascade output.
- At `1.0` the full QSound spatial profile for the active body is applied.
- Linear scaling of the spatial stage's wet contribution; bypass branch at
  `0.0` guarantees sample-accurate null.

### DSP placement

- Post-cascade stereo stage. The 12-stage DF2T biquad cascade stays frozen.
- Runs at audio rate inside the existing 32-sample control block.
- Uses the ITD / ILD / three-band-law coefficients already documented in
  `docs/archive/qsound_spatial.md`. No new spatial math is invented.

### Body-level data

- Each body's cartridge gains an optional `spatial_profile` block carrying the
  per-body azimuth, distance, elevation, and band coefficients.
- A body with no `spatial_profile` block defaults to identity; SPACE then
  behaves as a no-op regardless of setting.

### Why this shape

- One dial, one responsibility — matches v1 restraint doctrine.
- Off by default: user hears pure filter character first, dials space in as a
  deliberate production move.
- Spatial authoring stays inside the body; user sees a depth control.

## TRIG

### Behavior

- Normalized plugin parameter, range `[0.0, 1.0]`, default `0.0`.
- TRIG is the `PatchCord` amount from the body's authored LFO to MORPH.
  - At `0.0` the cord contributes nothing. MORPH equals the user knob.
  - At `1.0` the LFO signal is applied at full authored depth, summed with
    the user's MORPH knob and clamped to `[0.0, 1.0]`.
- LFO parameters are authored per body. User exposes depth only.

### The per-body LFO (follows `LFO/Clocked BPM.xml`)

Each body's cartridge gains an optional `mod_lfo` block with the E-mu fields:

| Field        | Type   | Meaning                                                      |
|--------------|--------|--------------------------------------------------------------|
| `frequency`  | float  | Rate. Interpretation depends on `tempo-sync`.                |
| `shape`      | int    | LFO waveform index (E-mu shape set; exact enum in plan).     |
| `delay`      | float  | Seconds of silence between key-on and LFO start.             |
| `var`        | float  | Per-voice shape variation amount (E-mu semantics).           |
| `key-sync`   | bool   | If true, LFO phase resets on note-on.                        |
| `tempo-sync` | bool   | If true, `frequency` is interpreted as a host-tempo division.|

- A body with no `mod_lfo` block defaults to a silent LFO (constant zero);
  TRIG is then inert regardless of setting.
- Default body authoring convention for v1 shipping bodies: `key-sync = 1`,
  `tempo-sync = 1`, giving the "sweep on each note" feel that motivated this
  feature. Rate, shape, delay, var are per-body creative choices.

### REV semantics

- REV inverts the polarity of the LFO contribution to MORPH (multiplies cord
  output by -1).
- REV has no effect when TRIG = 0 or when the body has no `mod_lfo`.
- REV does not affect MORPH knob behavior outside TRIG.

### Polyphony

- v1 is **monophonic** at the modulator level. The plugin runs a single LFO
  instance per active body. Key-sync restarts this single LFO on any note-on.
- Per-voice LFO instances are out of scope for v1.

### Why this shape

- Direct E-mu paradigm: LFO with key-sync + tempo-sync is exactly the "sweep
  on each note, locked to the grid" gesture trap production wants.
- All shape and personality live in the body. User scales one cord amount.
- REV reuses an existing button as a polarity switch — musical, not engineering.
- Pairs with ENV→Q for the one-two: LFO rides MORPH, ENV rides Q.

## Parameter additions

Two new plugin parameters exposed to the host and the existing parameter
layer:

| ID       | Range       | Default | Automatable |
|----------|-------------|---------|-------------|
| `space`  | 0.0 .. 1.0  | 0.0     | yes         |
| `trig`   | 0.0 .. 1.0  | 0.0     | yes         |

No other parameters change. MORPH, Q, TYPE, ENV, REV retain existing IDs.

## UI additions

- Two additional rollers matching the existing MORPH/Q roller component.
- Labels: `SPACE`, `TRIG`, in the existing faceplate label style.
- Layout fits within the current 340×510 editor bounds; exact pixel placement
  is a plan-phase decision against the active faceplate art.
- REV retains its existing button slot; its label may remain `REV`.
- No new button artwork, no new art assets beyond label text.

## Body schema changes

Additions to the body cartridge format, all optional, all default to identity
(no-op) when absent:

```
{
  // ...existing body fields...
  "spatial_profile": { /* QSound coefficients, absent = identity */ },
  "mod_lfo": {
    "frequency":  <float>,
    "shape":      <int>,
    "delay":      <float>,
    "var":        <float>,
    "key-sync":   <0|1>,
    "tempo-sync": <0|1>
  }
}
```

Existing baked cartridges in `cartridges/engine/` must load cleanly without
these blocks. Schema version bump and migration for baked pills is a
plan-phase deliverable.

Authored content for the 4 shipping bodies is a separate content task; this
spec governs the machinery, not the content.

## Verification

- Build and load the plugin standalone; both knobs appear and are automatable.
- SPACE = 0 produces sample-accurate null difference against current cascade
  output across all 4 bodies (enforced by bypass short-circuit, not by math).
- SPACE = 1 produces stereo output consistent with the body's spatial profile.
- TRIG = 0 produces identical output to current behavior across a note-on
  stream on any body.
- TRIG > 0 with an authored `mod_lfo` produces a MORPH modulation visible in
  the LCD spectrum trace, with phase reset on each note-on when
  `key-sync = 1`, and rate locked to host tempo when `tempo-sync = 1`.
- REV inverts the TRIG contribution polarity.
- Bodies without `spatial_profile` / `mod_lfo` remain functional; the
  respective control becomes inert.

## Out of scope for v1

- Multi-source mod matrix.
- TwistaLoop (sample-based) modulator shapes. Reserved for later milestone.
- Per-voice / polyphonic LFO instances.
- Velocity / aftertouch / mod-wheel direct parameter routing beyond ENV→Q.
- User-exposed LFO rate / shape / delay / var controls.
- User-exposed spatial azimuth / distance controls.
- Authoring UI for spatial profiles or LFO blocks inside the shipping plugin.
  Authoring happens in the separate authoring tool.
