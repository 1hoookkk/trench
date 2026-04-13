# TRENCH v1 Instrument Architecture Design

Date: 2026-04-05
Revision: 2 (post adversarial review + REV decision)

## Summary

TRENCH is a cartridge-based morph filter instrument. Each body (TYPE) defines
the complete instrument: spectral architecture, waveshaper, envelope profile,
trigger gesture, and stereo spread. The producer gets five controls and a
complete sound.

## Faceplate Controls

| Control | Type | Function |
|---------|------|----------|
| MORPH | Knob (0-1) | Anchor position on the morph axis |
| Q | Knob (0-1) | Pressure applied to all slots via their pressure behaviors |
| TYPE | Selector | Loads complete instrument definition (body) |
| ENV | Toggle button | Envelope follower on/off (reactive breathing) |
| REV | Toggle button | Anatomy flip — reverses the morph axis (effective_morph = 1.0 - morph) |

### REV Behavior

REV inverts the entire body geometry. State_b becomes the starting point,
state_a becomes the destination. This is not an envelope direction flip —
it reverses the static sound, the morph sweep, and the envelope direction
simultaneously.

Implementation: `cartridge.interpolate(1.0 - morph, q)` instead of
`cartridge.interpolate(morph, q)`. One subtraction. No recompilation,
no new corners, no architecture change.

REV doubles the launch: 4 bodies become 8 distinct acoustic environments.
A cavernous sub-hole becomes a strangling high-pass. A smooth vocal opening
becomes a violent choke. REV + ENV gives the "sidechain pump" effect —
transients push morph forward, but forward now means toward the closed state.

### Hidden/DAW Parameters

| Parameter | Type | Function |
|-----------|------|----------|
| TRIG | Momentary (MIDI note-on or 0-to-1 automation) | Fires the body's one-shot morph gesture from the current MORPH anchor |
| ALIGN | Toggle (DAW-automatable) | Forces L/R stereo morph offset to zero (mono collapse) |

TRIG is a performance gesture. Producers access it via MIDI mapping or
DAW automation. Sequential retriggering creates rhythmic patterns.

ALIGN is a mixing utility for mono-compatible bass. Not on the faceplate
because bass bodies ship with near-zero spread by design. Available when
producers need to collapse any body to mono for specific bus routing.

## Signal Chain

```
Input L/R ----> [Envelope Follower] ----> morph_offset (auto-scaled)
                                              |
              effective_morph = REV ? (1.0 - MORPH) : MORPH
              morph_final = clamp(effective_morph + env_offset, 0.0, 1.0)
              azimuth = body.spatial_trajectory(morph_final)

Input L ----> [6-stage DF2T cascade @ morph_final, Q] ----> [waveshaper]
          ----> [QSound spatial L @ azimuth] ----> Out L

Input R ----> [6-stage DF2T cascade @ morph_final, Q] ----> [waveshaper]
          ----> [QSound spatial R @ azimuth] ----> Out R
```

Both channels run the cascade at the SAME morph value. Stereo differentiation
comes entirely from the QSound spatial stage (HRTF spectral coloring + ILD + ITD),
not from morph offset between channels.

### Processing Order

1. Compute effective_morph: apply REV (1.0 - MORPH if engaged, else MORPH)
2. Envelope follower reads input level, produces raw morph offset
3. Offset is auto-scaled to input level (codec-style adaptive normalization)
4. morph_final = clamp(effective_morph + env_offset * ENV, 0.0, 1.0)
5. TRIG override: when firing, TRIG gesture replaces env_offset
6. Resolve azimuth from body's spatial trajectory at morph_final
7. 4-corner bilinear interpolation (Q first, then morph) per frozen invariant
8. 6-stage serial DF2T cascade processes audio (same coefficients both channels)
9. Body-specific waveshaper applies post-cascade saturation (same curve both channels)
10. QSound spatial stage applies per-channel HRTF coloring, ILD, ITD (bypassed if ALIGN)

### TRIG + ENV + REV Interaction

| ENV | REV | TRIG | Behavior |
|-----|-----|------|----------|
| off | off | idle | morph = MORPH knob (static or DAW-automated) |
| on  | off | idle | morph = MORPH + env_offset (breathing, opens on transients) |
| off | on  | idle | morph = 1-MORPH (reversed static, body inside-out) |
| on  | on  | idle | morph = (1-MORPH) + env_offset (breathing backward, chokes on transients) |
| any | any | fire | morph = anchor + trig_gesture (one-shot overrides ENV) |

### Edge Cases

**1. Envelope follower silence spikes.**
Sparse material (drum loops with gaps) causes the adaptive range tracker to
expand toward -inf dBFS during silence. The next transient produces an
explosive morph spike. Fix: enforce a minimum dynamic range of 12 dB on the
adaptive tracker. If tracked range < 12 dB, clamp to 12 dB floor. This
prevents infinite scaling during silence while preserving responsiveness
during continuous material.

**2. TRIG to ENV handoff discontinuity.**
When TRIG gesture completes and ENV is active, the morph value could jump
from the gesture's final position to the live envelope value, producing a
click. Fix: 10ms return glide. When TRIG ends, crossfade linearly from
the gesture's final morph value to the live ENV morph value over ~441
samples (10ms at 44.1kHz). The crossfade is:
```
morph = trig_final * (1 - t) + env_morph * t, where t ramps 0 to 1 over 10ms
```

**3. ALIGN zipper noise.**
Automating ALIGN toggles ITD delay and spectral coloring instantly, causing
phase clicks. Fix: ALIGN is a 20ms crossfade between QSound-processed and
dry output, not a boolean bypass. When ALIGN transitions 0 to 1:
```
output = qsound_out * (1 - t) + dry_out * t, where t ramps over 20ms
```
The ITD delay buffer continues running during crossfade to avoid discontinuity.

**4. Fast morph jumps under high Q.**
G.726-style <1ms envelope charge with high Q and a loud transient can push
morph_final by a large delta within a single 32-sample control block. The
per-sample coefficient ramping may not be smooth enough to prevent filter
ring-up. Fix: slew limiter on morph_final. Maximum morph delta per control
block = 0.15 (prevents >15% of the morph range changing in ~0.7ms). The
slew limit applies after all morph computation (ENV, TRIG, REV) and before
interpolation. Bodies with fast envelope profiles (charge <5ms) must pass
a stability test at Q=1.0 with impulsive input.

## Body Specification (Cartridge Format)

Each body defines six authored layers:

### Layer 1: Spectral Architecture

Six SlotSpecs, one per Actor (Foundation, Mass, Throat, Bite, Air, Scar).

Per slot:
- **state_a** — morph=0 endpoint (place, focus, weight, contour, color, hue)
- **state_b** — morph=1 endpoint (same parameters)
- **pressure behavior** — what Q does to this slot (Hold, Tighten, Yield, Cross, Collapse)
- **contour / zero family** — pole-zero relationship

Compiled to 4 corners (M0_Q0, M0_Q100, M100_Q0, M100_Q100) by the macro compiler.
Runtime interpolates between corners.

### Layer 2: Morph Curve Remapping

Per-slot blend function that shapes how state_a transitions to state_b.

| Curve | Behavior |
|-------|----------|
| linear | Default. Uniform transition. |
| ease_in | Stays near state_a, accelerates toward state_b late. |
| ease_out | Moves toward state_b early, settles. |
| s_curve | Gradual start, fast middle, gradual end. |
| step_at(threshold) | Holds state_a until threshold, then snaps to state_b. |

**Implementation constraint (from adversarial review):** Morph curve remapping
is a compile-time operation, not runtime. The runtime 4-corner interpolation
uses a single morph scalar across all 6 stages and cannot accept per-slot values.

The morph curves are applied in `compile_body()`: instead of compiling corners
at morph={0, 1}, the compiler evaluates the blended state at the remapped morph
value for each slot. The per-slot curve shapes what coefficients land in each
corner, creating the "authored stations" effect. At runtime, standard linear
interpolation between the baked corners produces non-linear spectral transitions
because the corners themselves encode the curve behavior.

This means the cartridge format does not change. The curves are a forge
authoring tool, not a runtime feature. The cartridge ships pre-baked corners.

### Layer 3: Waveshaper

Body-specific nonlinear transfer curve applied after the cascade. Serial only.

Parameters per body:
- **transfer_curve** — input-to-output amplitude mapping (polynomial coefficients or lookup table)
- **curve_shape** — warm (soft-knee), hard (clipping), asymmetric, etc.

**No band_split.** The cascade's spectral emphasis naturally drives peaks harder
into the waveshaper, producing frequency-dependent saturation without a crossover.
Parallel paths are banned by doctrine.

The waveshaper is always on. It is part of the body's identity. There is no
bypass toggle. Bodies that need minimal saturation ship with a near-unity
transfer curve. Existing cartridges (Acid_Squelch_77, Velvet_Throat,
Zero_Network_Test) receive identity transfer curves (passthrough waveshaper)
to preserve their current output.

Different bodies get different curves:
- Bass body: warm soft-clip, preserves sub weight
- Glass body: sharp clip, produces brittle harmonic artifacts
- Formant body: asymmetric curve, odd-harmonic emphasis for vocal grit
- Meta body: aggressive hard-clip, destructive at high morph

### Layer 4: Envelope Profile

Codec-derived timing characteristics for the envelope follower.

Parameters per body:
- **charge_ms** — attack time constant (one-pole rise toward target)
- **discharge_ms** — release time constant (one-pole fall toward target)
- **overshoot_db** — peak detector pre-ring above the one-pole follower
- **hysteresis_db** — Schmitt trigger dead zone on the output (must exceed threshold before responding)
- **target_excursion** — normalized morph range the follower targets (0.0-1.0)

#### Envelope Follower Algorithm

```
1. Measure input level:
   - RMS over 32 samples (one control block), computed as sqrt(sum(x^2)/32)
   - Convert to dBFS: 20 * log10(rms / 1.0)

2. Adaptive normalization:
   - Track running signal range using slow peak (charge) and floor (discharge)
   - Normalize current RMS to [0, 1] within tracked range
   - This produces "same morph excursion" regardless of input level

3. Apply body timing:
   - One-pole follower: y[n] = y[n-1] + alpha * (target - y[n-1])
   - alpha_charge = 1 - exp(-1 / (charge_ms * sr / 1000))
   - alpha_discharge = 1 - exp(-1 / (discharge_ms * sr / 1000))
   - Use charge alpha when target > current, discharge when target < current

4. Apply overshoot:
   - Separate peak detector with instant attack, slow decay
   - Blend: output = follower + overshoot_weight * (peak - follower)

5. Apply hysteresis:
   - Schmitt trigger on output delta: suppress changes < hysteresis_db
   - Once threshold exceeded, pass through until signal settles

6. Scale to morph offset:
   - offset = output * target_excursion
```

Timing profiles drawn from codec behavioral research:
- CVSD-style: 39ms charge, 16ms discharge (2.4:1 asymmetry). Breathing quality.
- G.726-style: <1ms charge, 2.5ms discharge, 3.5dB overshoot. Aggressive reactive.
- Symmetric: equal charge/discharge. Smooth continuous tracking.

#### Boundary Handling

morph_final is clamped to [0.0, 1.0] after adding the envelope offset.
The adaptive normalization targets the body's `target_excursion` range,
but clamping ensures the interpolation never receives out-of-bounds values.
The Rust runtime `interpolate()` must also clamp its morph input.

### Layer 5: Trigger Gesture

One-shot morph motion law fired by TRIG (MIDI note-on or automation).

Parameters per body:
- **segments** — envelope segments (attack, hold, decay, etc.)
- **segment_curves** — shape per segment (linear, exponential, circle)
- **total_duration_ms** — gesture length
- **morph_travel** — how far from anchor the gesture travels (0.0-1.0)
- **direction** — forward (toward state_b) or backward (toward state_a)

Each body has one authored gesture:
- Bass body: fast attack bark (50ms rip, slow 200ms return)
- Glass body: instant snap, long shimmering decay (800ms)
- Formant body: moderate symmetric sweep (300ms), vowel pronunciation
- Meta body: long destabilizing tear (2000ms), slow transformation

#### Retrigger Behavior

TRIG received while gesture in flight: restart from beginning. No stacking.
The gesture always starts from the current MORPH anchor position.

#### MIDI Requirement

Plugin must enable MIDI input (change `MidiConfig::None` to
`MidiConfig::MidiCCs` or `MidiConfig::Basic`). Note-on events with
velocity > 0 fire TRIG. Implementation must handle note-on in the
`process()` method alongside audio processing.

### Layer 6: QSound Spatial Stage

HRTF-based spatial processing derived from reverse-engineered QSound
reconstruction. Not a morph offset trick — real frequency-dependent
spatial positioning using measured head-related transfer function data.

#### Source Data

Reconstruction run `qsound_spatial_v1_recon_20260305_092531`:
- ITD MAE: 7.48 us (passed, threshold 30 us)
- ILD MAE: 0.55 dB (passed, threshold 1.0 dB)
- Spectral MAE: 2.18 dB (close, threshold 1.5 dB)
- 3-band FIR model per channel: low (0.131 dB), mid (0.202 dB), high (0.178 dB) fitting error
- 9 HRTF impulse response WAVs at key spatial positions
- Training/validation metrics across 195 azimuth/elevation/distance positions

Data location: `trench_re_vault/artifacts/qsound_sanitized_handoff/qsound/spatial/v1.0.0/`

#### Architecture

Three components per channel, post-waveshaper:

```
[waveshaper output] ---> [3-band spectral coloring] ---> [ILD gain] ---> [ITD delay] ---> output
```

**1. Spectral coloring (3-band HRTF filter):**
- Three biquad stages per channel: low shelf (~300Hz), mid band (300Hz-3kHz), high shelf (>3kHz)
- Each band has a per-channel gain coefficient derived from HRTF measurements
- Gains vary with the body's authored spatial position
- This is where the "alien" quality comes from — HRTF notches at 4kHz, 8kHz, 10kHz
  interact with the body's notch architecture to create spectral collisions that
  don't exist in standard processing

**2. ILD (interaural level difference):**
- Simple per-channel gain multiplier
- 0 dB at center (0° azimuth), up to ~4.5 dB at 90°
- Monotonic with azimuth

**3. ITD (interaural time difference):**
- Short delay on one channel (max ~22 samples at 44.1kHz / ~500 us at 45° azimuth)
- Implemented as a small circular delay buffer
- Creates the phantom source positioning effect

#### Per-Body Spatial Parameters

Each body defines a **spatial trajectory**: a mapping from morph position
to QSound azimuth angle.

Parameters per body:
- **azimuth_a** — spatial position at morph=0 (degrees, -90 to +90)
- **azimuth_b** — spatial position at morph=1 (degrees, -90 to +90)
- **spatial_curve** — interpolation curve between azimuth_a and azimuth_b (linear, ease, step)
- **distance** — virtual source distance (affects spectral attenuation law: 1/(1+d^2))

At each morph position, the runtime looks up or interpolates the ITD, ILD,
and 3-band spectral gains for the current azimuth from a pre-computed table
derived from the HRTF reconstruction data.

Body-specific spatial trajectories:
- Bass body: azimuth_a=0°, azimuth_b=0°. Mono. Sub stays dead center.
- Formant body: azimuth_a=0°, azimuth_b=15°. Subtle widening as vowel opens.
- Texture body: azimuth_a=0°, azimuth_b=30°. Shimmer spreads as morph opens.
- Meta body: azimuth_a=0°, azimuth_b=45°. Full spatial tear at high morph.

REV interaction: REV flips the morph axis, so the spatial trajectory also
reverses — the body starts wide and collapses to center as morph increases.

#### ALIGN (Mono Collapse)

ALIGN (DAW-automatable parameter) bypasses the entire QSound spatial stage.
Both channels receive identical processing — no ITD, no ILD, no spectral
coloring difference. Used for mono bus routing, bass mono-compatibility
checks, or creative spatial collapse/expansion via automation.

When ALIGN is engaged, the spatial stage passes audio through unmodified.
When disengaged, the full QSound processing resumes.

#### Computational Cost

- 3 biquad stages per channel (spectral coloring) = 6 biquads total
- 1 gain multiply per channel (ILD) = 2 multiplies
- 1 delay buffer per channel (ITD) = ~22 samples max, circular buffer read/write
- Coefficient updates per control block (32 samples), not per sample
- Total: lightweight. Less than doubling the existing cascade cost.

#### Pre-computed HRTF Table

At plugin init, build a lookup table from the reconstruction data:
- Azimuth positions: 0° to 90° in 5° steps (19 entries, mirror for negative)
- Per entry: ITD_samples, ILD_dB, spectral_gain_low_L, spectral_gain_mid_L,
  spectral_gain_high_L, spectral_gain_low_R, spectral_gain_mid_R, spectral_gain_high_R
- Runtime interpolates between table entries for the current azimuth
- Table derived from `train_metrics_impulse_only.csv` and HRTF WAV analysis

#### Mono Collapse Validation

Body authors must verify: engaging ALIGN produces no audible discontinuity.
The validation gate measures level difference between spatial and collapsed
output. Bodies with azimuth_b > 30° require explicit ALIGN transition
testing during authoring.

## Frozen DSP Invariants (Unchanged)

- Topology: 12-stage serial DF2T cascade. Exactly 6 active + 6 passthrough.
- Math: kernel-form [c0, c1, c2, c3, c4] interpolation only.
- Order: 4-corner bilinear interpolation. Evaluate Q first, then morph.
- Runtime: 32-sample control blocks with per-sample coefficient ramping.
- Rates: authoring = 39062.5 Hz, plugin runtime = 44100 Hz.

The new layers (waveshaper, envelope, trigger, stereo) operate outside the
frozen cascade. They modulate its inputs (morph position) and process its
outputs (waveshaper), but the cascade itself is unchanged.

## Body Design Constraints

Every morph position must sound good. No dead zones, no level collapse,
no harsh jumps. The body is a shortcut to a good sound.

Validation criteria per body:
- Sweep test: morph 0 to 1 at Q=0 and Q=1. Every position is a usable sound.
- REV sweep test: same sweep with REV engaged. Reversed body is equally usable.
- ENV test: pink noise through body with ENV on. Breathing is musical.
- ENV+REV test: same with REV. Choking is rhythmic, not destructive.
- TRIG test: fire gesture. One-shot sounds intentional.
- Stereo test: mono collapse produces no comb artifacts or level loss.
- Midpoint audit: morph=0.5 is not a dead zone between two good endpoints.

## V1 Body Lineup

Four bodies. Three practical, one meta. All notch-heavy, all musical,
all distinct. Each body ships with all six layers authored.

REV doubles the lineup: 4 bodies x 2 orientations = 8 acoustic environments.

Body names, spectral strategies, morph narratives, waveshaper curves,
envelope profiles, trigger gestures, and stereo spreads to be designed
in a follow-up session after this architecture is approved.

## Implementation Scope

### New runtime components needed:
1. Envelope follower with adaptive scaling and body-specific timing profiles
2. Per-slot morph curve remapping in `compile_body()` (compile-time only)
3. Post-cascade serial waveshaper stage (no crossover, no parallel paths)
4. TRIG gesture generator (segment envelope, retrigger = restart)
5. QSound spatial stage: 3-band HRTF spectral coloring + ILD + ITD per channel
6. Pre-computed HRTF lookup table (azimuth to ITD/ILD/spectral gains)
7. Body spatial trajectory resolution (morph to azimuth)
8. REV morph inversion (1.0 - morph)
9. ALIGN parameter (DAW-automatable, bypasses QSound spatial stage)
10. MIDI input enable for TRIG (MidiConfig change)
11. Cartridge format extension (envelope, waveshaper, trigger, spatial trajectory, morph curves)
12. Morph + envelope offset clamping to [0.0, 1.0]

### Unchanged:
- 6-stage DF2T cascade
- 4-corner bilinear interpolation
- Kernel-form coefficient encoding
- 32-sample control blocks with per-sample ramping
- Authoring sample rate (39062.5 Hz)

### Forge extensions needed:
- Morph curve specification per SlotSpec (compile-time blend shaping)
- Waveshaper profile authoring and validation
- Envelope profile authoring (codec-style timing parameters)
- Trigger gesture authoring (segment envelope designer)
- Spatial trajectory authoring (morph to azimuth mapping per body)
- HRTF table generation from QSound reconstruction data
- ALIGN transition testing during body authoring
- Extended body validation (sweep, REV, ENV, ENV+REV, TRIG, spatial, ALIGN, midpoint)
