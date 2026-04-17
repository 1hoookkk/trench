# SPACE + TRIG + pre-cascade Mackie drive — v1 UI and DSP Addition

**Date:** 2026-04-18
**Status:** Design, not yet implemented
**Scope:** v1.0 shipping plugin (2026-07-15)

## Intent

Three additions to TRENCH v1, two of them user-facing, one fully authored:

1. **SPACE** (user knob) — user-scaled QSound spatial, baked per body,
   off by default.
2. **TRIG** (user knob) — user-scaled modulation depth from a per-body
   E-mu-style 8-segment function generator to MORPH, off by default.
3. **Pre-cascade Mackie drive** (no UI) — per-body authored input drive into
   a Mackie-preamp-modeled soft-clip stage feeding the cascade. This is the
   voice of TRENCH's distortion character. User sees nothing; the body is
   the knob.

The faceplate stays minimal. Character lives in bodies.

## Doctrine anchor: follow E-mu

E-mu reference material (`docs/emu/`), E-mu X3 templates, and the E-mu
modulation reference define the native vocabulary. TRENCH adopts it; it does
not invent around it.

- MORPH is a **modulation destination**.
- Modulation **sources** are E-mu primitives: LFO, function generator,
  envelope, note-on controllers, realtime controllers.
- A `PatchCord` routes source → destination with a scalar amount.
- Distortion in E-mu comes from **driving into the filter** (Filter Level),
  not from post-filter shapers. TRENCH's drive stage sits pre-cascade.
- Mackie preamp saturation is the producer-historical path for the "E-mu
  sound." TRENCH bakes a Mackie-modeled soft-clip into its pre-cascade drive
  stage so that path is in the box.

v1 exposes **one cord on the faceplate**: function generator → MORPH, amount
= TRIG. All other cords, all other sources, and all other destinations
belong to later milestones.

## Non-goals

- No mod matrix UI.
- No user-exposed LFO / function-generator parameters on the face.
- No velocity / aftertouch / mod-wheel direct routing in v1 beyond existing
  ENV→Q.
- No user-facing DRIVE knob. Drive is body-authored.
- No tempo-sync controls on the face.
- No TRIG arm button. No SPACE enable button.
- No spatial azimuth / distance user controls.
- No post-cascade saturation, no waveshaper block, no tape/analog FX.
- No changes to the frozen 12-stage DF2T cascade.
- No TwistaLoop (sample-based) modulator in v1.
- No multi-voice / polyphonic modulator instances.

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

Two new rollers (SPACE, TRIG). Zero new buttons. REV and ENV exist; ENV
keeps its current envelope→Q routing, REV gains new semantics (see TRIG).

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
- Uses the ITD / ILD / three-band-law coefficients documented in
  `docs/archive/qsound_spatial.md`. No new spatial math is invented.

### Body-level data

- Each body's cartridge gains an optional `spatial_profile` block carrying
  the per-body azimuth, distance, elevation, and band coefficients.
- A body with no `spatial_profile` block defaults to identity; SPACE becomes
  a no-op regardless of setting.

## TRIG

### Behavior

- Normalized plugin parameter, range `[0.0, 1.0]`, default `0.0`.
- TRIG is the `PatchCord` amount from the body's authored **function
  generator** to MORPH.
  - At `0.0` the cord contributes nothing. MORPH equals the user knob.
  - At `1.0` the generator signal is applied at full authored depth, summed
    with the user's MORPH knob and clamped to `[0.0, 1.0]`.
- All generator parameters are authored per body. User exposes depth only.

### The per-body function generator (E-mu 8-segment)

Each body's cartridge gains an optional `mod_fn` block. The primitive is
E-mu's 8-segment function generator: levels, times, shapes, conditional
jumps. Acts as LFO, envelope, or mini-sequencer depending on authoring.

```
"mod_fn": {
  "segments": [
    { "level": <float>, "time_ms": <float>, "shape": "hold|lin|exp|log|sCurve", "jump": <int?> },
    // ... up to 8 segments total
  ],
  "key-sync":   <0|1>,
  "tempo-sync": <0|1>
}
```

| Field        | Meaning                                                          |
|--------------|------------------------------------------------------------------|
| `level`      | Target output at segment end, range `[-1.0, 1.0]`.               |
| `time_ms`    | Segment duration. Interpreted as tempo-divided if `tempo-sync=1`.|
| `shape`      | Curve from previous level to this segment's level.               |
| `jump`       | Optional. Index of segment to branch to at segment end, for loops.|
| `key-sync`   | If 1, generator restarts on note-on.                             |
| `tempo-sync` | If 1, `time_ms` values are interpreted as tempo divisions.       |

- A body with no `mod_fn` block defaults to silent output; TRIG becomes
  inert regardless of setting.
- Default shipping-body authoring convention: `key-sync = 1`. Loop via
  `jump` when a repeating modulator is desired; leave jumps absent for
  one-shot envelope behavior.
- Authoring of `mod_fn` blocks happens in the external authoring tool, not
  the shipping plugin.

### REV semantics

- REV inverts the polarity of the function generator's contribution to
  MORPH (multiplies cord output by -1).
- REV has no effect when TRIG = 0 or when the body has no `mod_fn`.
- REV does not affect MORPH knob behavior outside TRIG.

### Polyphony

- v1 is **monophonic** at the modulator level. The plugin runs a single
  function generator instance per active body. `key-sync` restarts this
  single instance on any note-on.
- Per-voice generator instances are out of scope for v1.

## Pre-cascade Mackie drive

### Behavior

- A per-body authored input-gain scalar feeds a Mackie-modeled soft-clip
  nonlinearity, whose output feeds stage 1 of the cascade.
- Hot bodies push the cascade into pole-brink territory; the Mackie curve
  carries the harmonic character on the way in.
- User sees no control. The body *is* the drive setting.

### DSP placement

- Pre-cascade, mono (or stereo if the input is stereo) stage.
- Soft-clip curve modeled on the Mackie 1202 non-VLZ preamp circuit
  (Cytomic-documented equivalent to the 24:8 / 32:8 bus desks).
- Not a tanh placeholder. The curve must match the specific Mackie
  input-stage character; plan phase selects and validates the model.
- Does not touch the frozen 12-stage cascade. Sits before it.

### Body-level data

Each body's cartridge gains a required `drive` block:

```
"drive": {
  "input_gain_dB": <float>,   // gain into the softclip, authored per body
  "model": "mackie_1202"      // reserved; v1 ships a single model
}
```

`input_gain_dB` is required for all v1 shipping bodies. Migration for
existing baked cartridges in `cartridges/engine/` defaults
`input_gain_dB = 0.0` and `model = mackie_1202`; individual shipping
bodies override to their authored hot value.

### Why baked, not a knob

- Mackie saturation is the TRENCH voice, not a producer option. Exposing it
  as a knob invites users to set it to 1 on every preset and blame the
  plugin for harshness.
- Authoring trouble, not user-bolted trouble. Hot bodies ship hot.
- Producers who want to drive harder drive into the plugin from their DAW,
  same workflow they already use with real Mackies.

## Output safety limiter

- Post-cascade, post-SPACE, hard safety limiter at the final output.
- Inaudible under normal playing. Catches peak runaways from pole-brink
  bodies pushed into regions the author did not anticipate.
- Not a tone-shaping block. Speaker protection only. Matches E-mu hardware
  headroom management on the DAC.

## Parameter additions

Two new plugin parameters exposed to the host:

| ID       | Range       | Default | Automatable |
|----------|-------------|---------|-------------|
| `space`  | 0.0 .. 1.0  | 0.0     | yes         |
| `trig`   | 0.0 .. 1.0  | 0.0     | yes         |

No other parameters change. MORPH, Q, TYPE, ENV, REV retain existing IDs.
No user-facing DRIVE parameter.

## UI additions

- Two additional rollers matching the existing MORPH/Q component.
- Labels: `SPACE`, `TRIG`.
- Layout fits within the current 340×510 editor bounds; exact pixel
  placement is a plan-phase decision against the active faceplate art.
- REV retains its existing button slot.
- No new button artwork, no new art assets beyond label text.

## Body schema changes

Additions to the body cartridge format:

```
{
  // ...existing body fields...

  // required, new in v1 schema bump
  "drive": {
    "input_gain_dB": <float>,
    "model": "mackie_1202"
  },

  // optional, absent = identity / no-op
  "spatial_profile": { /* QSound coefficients */ },
  "mod_fn":          { /* 8-segment generator */ }
}
```

Schema version bump for the cartridge format. Baked cartridges in
`cartridges/engine/` that predate this spec load with:
- `drive = { input_gain_dB: 0.0, model: "mackie_1202" }`
- `spatial_profile` absent → SPACE inert.
- `mod_fn` absent → TRIG inert.

Migration behavior is a plan-phase deliverable.

Authored content for the 4 shipping bodies (hot `input_gain_dB` values,
spatial profiles, function generator segments) is a separate content task;
this spec governs the machinery.

## Signal chain (v1 final)

```
input
  → pre-cascade: body.drive.input_gain_dB + Mackie-modeled softclip
  → 12-stage DF2T biquad cascade (frozen, with per-sample ramped coeffs)
    ├─ MORPH = user knob + (TRIG * mod_fn(t) * (REV ? -1 : 1)), clamped
    └─ Q     = user knob + envelope_follower if ENV=on
  → SPACE: bypass at 0, QSound spatial wet-scaled to SPACE at >0
  → safety limiter
  → output
```

## Verification

- Build and load plugin standalone; SPACE and TRIG appear and are
  automatable. No DRIVE parameter visible to host.
- SPACE = 0 produces sample-accurate null difference against current
  cascade output across all 4 bodies (bypass short-circuit enforced).
- SPACE = 1 produces stereo output consistent with the body's spatial
  profile.
- TRIG = 0 produces identical output to current behavior on any body.
- TRIG > 0 with an authored `mod_fn` produces MORPH modulation visible in
  the LCD spectrum trace, with phase reset on each note-on when
  `key-sync = 1`, and time axis locked to host tempo when `tempo-sync = 1`.
- REV inverts the TRIG contribution polarity.
- A body with authored hot `input_gain_dB` audibly pushes the cascade
  harder than a body with `input_gain_dB = 0.0` (measurable in harmonic
  content, not just level).
- Mackie-model soft-clip output matches the reference curve within a
  plan-phase-specified tolerance.
- Safety limiter engages only on peak excursions; transparent on normal
  material.
- Bodies without optional blocks remain functional; the respective control
  becomes inert.

## Out of scope for v1

- Multi-cord / mod matrix.
- User-facing DRIVE, LFO rate / shape, function generator authoring in
  the plugin.
- Multiple drive models. v1 ships `mackie_1202` only.
- Post-cascade saturation / waveshapers / tape FX.
- TwistaLoop (sample-based) modulator.
- Per-voice / polyphonic function generator instances.
- Velocity / aftertouch / mod-wheel direct routing beyond existing ENV→Q.
- Stereo Fuzz, TubeJam-style clipping filters.
- User-exposed spatial azimuth / distance controls.
- Authoring UI inside the shipping plugin.
