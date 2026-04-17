# SPACE + TRIG — v1 UI and DSP Addition

**Date:** 2026-04-18
**Status:** Design, not yet implemented
**Scope:** v1.0 shipping plugin (2026-07-15)

## Intent

Add two musical gestures to TRENCH v1 without expanding the faceplate beyond its
"Musical Filter" identity:

1. **SPACE** — user-scaled QSound spatial, baked per body, off by default.
2. **TRIG** — note-on MORPH sweep along a body-authored trajectory, off by default.

Both behave like character dials: zero is the null baseline, turn it up to bring
in the body's personality. Nothing global, nothing synth-matrix-shaped.

## Non-goals

- No LFO.
- No modulation matrix (no source pickers, no destination pickers, no depths per pair).
- No per-voice velocity / aftertouch / mod-wheel routing in v1.
- No tempo sync.
- No TRIG arm button. No SPACE enable button. Knobs at zero = off.
- No spatial azimuth / distance user controls. Spatial parameters live in the body.
- No changes to the frozen 12-stage DF2T cascade.

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
- At `0.0` the audio path is bit-identical to the current cascade output.
- At `1.0` the full QSound spatial profile for the active body is applied.
- Linear scaling of the spatial stage's wet contribution.

### DSP placement

- Post-cascade stereo stage. The 12-stage DF2T biquad cascade stays frozen.
- Runs at audio rate inside the existing 32-sample control block.
- Uses the ITD / ILD / three-band-law coefficients already documented in
  `docs/archive/qsound_spatial.md`. No new spatial math is invented.
- Dry/wet crossfade against the cascade output; SPACE is the wet scalar.

### Body-level data

- Each body's cartridge gains a `spatial_profile` block carrying the per-body
  azimuth, distance, elevation and band coefficients that the stage consumes.
- A body with no `spatial_profile` defaults to identity (pan = 0, distance = 0,
  flat band law); SPACE then behaves as a no-op regardless of setting.

### Why this shape

- Matches the v1 identity of restraint: one dial, one responsibility.
- "Off by default" means the user first hears the body's filter character pure,
  then dials space in as a deliberate production move — not a baseline tint.
- Keeps spatial authoring inside the body (where sound design lives), user sees
  a single depth control.

## TRIG

### Behavior

- Normalized plugin parameter, range `[0.0, 1.0]`, default `0.0`.
- On every note-on, a one-shot envelope drives MORPH from the body's authored
  start point to its authored end point, following the body's authored curve.
- TRIG is the depth / blend amount of that envelope against the user's current
  MORPH knob position.
  - At `0.0` TRIG contributes nothing. MORPH equals the user knob value.
  - At `1.0` the TRIG envelope fully overrides the MORPH knob for the duration
    of the envelope, then returns to the user's MORPH value on envelope end.
- No looping, no retrigger complexity: each note-on starts a fresh envelope.
  If a new note arrives mid-envelope, the envelope restarts from its start.

### REV semantics

- REV button flips the TRIG envelope direction: end → start instead of
  start → end. Same authored trajectory, inverse arc.
- REV has no effect when TRIG = 0 or when the body has no trajectory.
- REV does not affect MORPH knob behavior outside TRIG.

### Body-level data

Each body's cartridge gains a `trig_trajectory` block:

- `start_morph` — 0..1 baseline morph at envelope t=0.
- `end_morph` — 0..1 baseline morph at envelope t=1.
- `duration_ms` — envelope length, authored per body.
- `curve` — normalized [0,1] samples (fixed count, e.g. 64) describing the
  shape between start and end. Linear interpolation at runtime.

Exact sample count and storage format are an implementation detail for the
plan phase; the body schema must carry enough resolution to express curves,
plateaus, and mid-shot wobbles without aliasing at the authored duration.

A body with no `trig_trajectory` defaults to identity (no envelope); TRIG is
then a no-op regardless of setting.

### Why this shape

- Transient music (trap, underground) is note-on driven. Opening the filter on
  the hit is the instrument's primary expressive move.
- Pairs with ENV→Q (transient bite) for one-two physical response:
  TRIG rides MORPH across the note, ENV rides Q on the attack.
- All shape and personality live in the body. User scales one amount.
- Reuses REV as a musical button instead of an engineering toggle.

## Parameter additions

Two new plugin parameters exposed to the host and the existing parameter layer:

| ID | Range | Default | Automatable |
|----|-------|---------|-------------|
| `space`  | 0.0 .. 1.0 | 0.0 | yes |
| `trig`   | 0.0 .. 1.0 | 0.0 | yes |

No other parameters change. MORPH, Q, TYPE, ENV, REV retain existing IDs.

## UI additions

- Two additional rollers matching the existing MORPH/Q roller component.
- Labels: `SPACE`, `TRIG`, in the existing faceplate label style.
- Layout fits within the current 340×510 editor bounds; exact pixel placement
  is a plan-phase decision against the active faceplate art.
- REV retains its existing button slot; its label may remain `REV`.
- No new button artwork, no new art assets beyond label text.

## Body schema changes

Additions to the body cartridge format, both optional, both default to
identity (no-op) when absent:

```
{
  // ...existing body fields...
  "spatial_profile": { /* QSound coefficients, absent = identity */ },
  "trig_trajectory":  { /* start/end/duration/curve, absent = no envelope */ }
}
```

Authored content for the 4 shipping bodies (Speaker Knockerz, Aluminum Siding,
Small Talk Ah-Ee, Cul-De-Sac) is tracked separately as a content task; this
spec governs the machinery, not the content.

## Verification

- Build and load the plugin standalone; both knobs appear and are automatable.
- SPACE = 0 produces null difference against the current cascade output
  (sample-accurate) across all 4 bodies.
- SPACE = 1 produces stereo output consistent with the body's spatial profile.
- TRIG = 0 produces identical output to current behavior across a note-on
  stream.
- TRIG = 1 with an authored trajectory produces a MORPH sweep visible in the
  LCD spectrum trace on each note-on.
- REV flips the TRIG sweep direction.
- Disabling (setting to 0) either control restores the null baseline.
- Bodies without `spatial_profile` / `trig_trajectory` remain functional; the
  respective control becomes inert.

## Out of scope for v1

- Multi-source mod matrix.
- LFO, tempo-synced wobble, arpeggiator-style triggers.
- Velocity / aftertouch / mod-wheel direct parameter routing.
- User-exposed spatial azimuth / distance controls.
- Per-body SPACE or TRIG overrides exposed on the face.
- Authoring UI for spatial profiles or trig trajectories inside the shipping
  plugin. Authoring happens in the separate authoring tool.
