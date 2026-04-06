# TRENCH Modulation Exploration

Status: research / branch setup
Date: 2026-04-07

## Modulation Grammar

> **Bodies expose pressure, inversion, and gesture — not arbitrary modulation.**

All modulation in TRENCH is pre-authored per body. There is no patchcord system,
no user-facing routing matrix, no open modulation bus. A body defines what it
reacts to and how. The user's job is to play the body, not to wire it.

This rule exists to prevent TRENCH from becoming a DSP sandbox. If a modulation
feature cannot be expressed as authored body behavior, it does not ship.

## What We Have

### DSP Infrastructure (Implemented, Shipping)

| Component | Location | Status |
|-----------|----------|--------|
| 32-sample control blocks | `trench-core/src/cascade.rs` | Live |
| Per-sample coefficient ramping (linear delta) | `cascade.rs:BiquadState` | Live |
| Post-cascade boost ramping | `cascade.rs:Cascade` | Live |
| 4-corner bilinear interpolation (Q-first, morph-second) | `cartridge.rs:interpolate()` | Live, frozen |
| Sample-accurate automation flag | `trench-plugin/src/lib.rs` | Declared |
| MORPH / Q / TYPE parameters | `lib.rs:TrenchParams` | Live |

### Specified but Unimplemented (v1 Architecture Spec)

| Layer | What It Does | Timing Profile |
|-------|-------------|----------------|
| **ENV (Envelope Follower)** | Input RMS -> morph offset | CVSD: 39ms/16ms, G.726: <1ms/2.5ms, Symmetric |
| **REV (Anatomy Flip)** | `effective_morph = 1.0 - morph` | Instant (flip is arithmetic) |
| **TRIG (One-Shot Gesture)** | Multi-segment morph excursion from anchor | 50ms-2000ms per body |
| **Morph Curve Remap** | Per-slot nonlinear blend functions | Compile-time bake |
| **Waveshaper** | Post-cascade nonlinear transfer curve | Always-on, per-body |
| **QSound Spatial** | HRTF + ILD + ITD per channel | Morph-driven azimuth trajectory |
| **ALIGN** | Mono collapse (bypass QSound) | 20ms crossfade |
| **Slew Limiter** | Max morph delta = 0.15 per block | Per-block clamp |

### Heritage (E-mu Morpheus / UltraProteus / Emulator X)

From RE vault, DillusionMan documentation, and Emulator X cord architecture:

**Cord System Architecture:**
- Source -> Destination -> Amount (-128 to +127, negative = inverted)
- Multiple sources to one destination, one source to multiple destinations
- Completely flexible routing — no fixed signal path

**Note-On Sources** (measured once at key strike):
- Key Velocity, Key Number (keyboard tracking)
- Initial controller positions

**Realtime Sources** (continuously evaluated while note sustains):
- LFOs (free-running periodic)
- Auxiliary Envelopes (multi-segment)
- Function Generators (2 per voice, complex multi-segment)
- Pitch wheel, mono pressure, poly pressure
- MIDI CCs (mapped to MIDI A-D)

**Destinations:**
- Primary/secondary pitch, volume, pan
- LFO rates, envelope segment times
- Z-Plane filter: **Morph** (Transform 1), **Frequency Tracking**, **Transform 2**

**Transform 2 is significant.** In the Morpheus lineage, Transform 1 = morph axis,
Transform 2 = an orthogonal morphing axis. This is the heritage precedent for the
diagonal gesture concept — E-mu already had a second axis available as a cord
destination. TRENCH's Q parameter occupies this position but is not yet modulated
by anything except the user's knob.

**Heritage mapping to TRENCH v1:**

| E-mu Cord System | TRENCH Equivalent | Status |
|-----------------|-------------------|--------|
| Velocity -> Morph | velocity_scale on TRIG | v1.x |
| LFO -> Morph | ENV (follower, not periodic) | v1 ship |
| Pressure -> FilRes | Aftertouch -> Q | v1.x |
| FuncGen -> Morph | TRIG gesture engine | v1 ship |
| CC -> Transform 2 | Q knob (manual only) | Live |
| Negative cord amount | REV (morph inversion) | v1 ship |
| Key Number -> Morph | Not planned | — |

TRENCH v1 deliberately does NOT ship a patchcord system. Modulation is pre-authored
per body. The E-mu cord flexibility is collapsed into body-specific authored behaviors.

### Python Reference Implementations

- `pyruntime/drive.py:erode()` — One-pole envelope follower with attack/release asymmetry
- `pyruntime/macro_compile.py` — Actor-based slot compiler with pressure behaviors
- `pyruntime/analysis.py` — Morph trajectory analysis with midpoint audit

---

## What's Novel / Unexploited

### 1. Cross-Term Modulation (Exploiting the Bilinear Surface)

The interpolation is `x(m,q) = A + q*Dq + m*Dm + m*q*Dmq`.

The cross-term `Dmq` means certain spectral behaviors **only exist** when both morph AND Q move simultaneously. No current modulation source drives morph and Q in coordinated opposition or phase-locked motion.

**Idea: Diagonal Gestures.** A TRIG gesture that moves along a diagonal of the morph/Q surface instead of just the morph axis. Bodies authored with strong cross-terms would produce sounds unreachable by sequential morph-then-Q motion.

### 2. Envelope Follower Driving Q (Not Just Morph)

The spec routes ENV -> morph offset only. But Q controls pressure behavior per slot (Hold, Tighten, Yield, Cross, Collapse). An envelope follower driving Q would create **input-reactive aggression** — louder input = tighter resonance, not just different spectral position.

**Idea: Split ENV routing.** ENV drives morph at one ratio and Q at another. Per-body authored split. Creates "the harder you hit it, the more it fights back" behavior.

### 3. Per-Slot Morph Staggering

All 6 slots currently receive the same morph value. If each slot had a tiny time offset (1-4 control blocks = 0.7ms-2.8ms), the cascade would briefly exist in an impossible state where different stages are at different morph positions.

**Idea: Morph Smear.** A body-authored parameter (`slot_delay_blocks: [0, 1, 2, 1, 0, 3]`) that staggers when each slot reaches its target morph. Creates transient micro-detuning during sweeps. Zero cost at static morph positions.

### 4. REV as Modulation Source (Not Just Toggle)

REV currently flips the axis statically. If REV were a continuous 0-1 parameter instead of a toggle, it would act as a **morph inversion crossfade**. At REV=0.5, the body collapses to its midpoint from both directions simultaneously.

**Idea: REV Sweep.** Automate REV 0->1 slowly. The body literally turns inside out over time. Combined with ENV, creates "breathing backward" — transients push toward the state that was originally the resting state.

### 5. Gesture Chaining (TRIG Sequences)

Current TRIG fires one gesture and retriggers from the beginning. If a body could define 2-3 gesture variants that cycle on successive TRIG events, rapid retriggering would produce rhythmic patterns.

**Idea: TRIG Cycle.** Body defines `gestures: [bark, shimmer, choke]`. First TRIG fires bark, second fires shimmer, third fires choke, then cycles. Creates a 3-hit rhythmic morph pattern from a single MIDI note stream.

### 6. Codec-Derived Adaptive Envelope (Already Specified, Deeply Novel)

The spec already calls for CVSD-style and G.726-style envelope timing. These aren't just "fast attack / slow release" — they're modeled on **voice codec step-size adaptation algorithms**.

CVSD compressors track signal slope and increase step size when consecutive samples have the same sign. This creates an envelope follower that responds to **signal continuity**, not just level. A tone held steady gets a slow follower; a stuttering signal gets fast adaptation.

**What's unexploited:** The G.726 profile's 3.5dB overshoot creates a momentary morph spike above the steady-state position on every transient. This is a built-in accent generator. Bodies authored for dramatic mid-morph behavior (like Cul-De-Sac's "State Boundary" zone) would momentarily flash through dangerous territory on every hit.

### 7. QSound Spatial as Modulation Destination

The spatial trajectory is currently a fixed authored curve per body. If ENV or TRIG could modulate the spatial offset independently of morph, transients could cause spatial jolts without changing the spectral state.

**Idea: Spatial Punch.** On TRIG, azimuth jumps to a wide angle and decays back to the morph-derived position. The filter doesn't move, but the sound briefly shoots sideways. Physical impact sensation without spectral change.

### 8. Waveshaper Drive Modulation

The waveshaper is always-on with a fixed transfer curve. If the curve's drive amount were morph-linked (already implicit — higher morph = more spectral energy = harder waveshaper hit), an explicit drive modulation from ENV would create **compression-like behavior** where loud signals get more saturation independent of morph position.

---

## Interaction Matrix (All Sources x All Destinations)

```
              MORPH    Q       SPATIAL   WAVESHAPER-DRIVE
MORPH knob    [live]   --      [spec'd]  [implicit]
Q knob        --       [live]  --        --
ENV           [spec'd] [novel] [novel]   [novel]
TRIG          [spec'd] [novel] [novel]   --
REV           [spec'd] --      [spec'd]  --
DAW auto      [live]   [live]  --        --
MIDI vel      --       --      --        --
Aftertouch    --       --      --        --
```

`[live]` = implemented now
`[spec'd]` = in v1 architecture spec
`[novel]` = new idea from this exploration
`--` = not routed

---

## Heritage Modulation We Could Adapt

### Note-On vs Realtime Distinction
E-mu split modulation into sources evaluated once (note-on) and sources evaluated
continuously (realtime). TRENCH already mirrors this: TRIG is note-on (fires once,
plays gesture), ENV is realtime (continuous follower). The value of making this
distinction explicit in body authoring: note-on sources can scale gesture parameters
(velocity -> travel distance), while realtime sources drive continuous offsets
(envelope -> morph position). Mixing the two categories in the same authored slot
would be incoherent.

### Function Generator Shape Library
E-mu function generators supported segment shapes: linear, exponential, logarithmic,
and "circle" curves. The TRIG gesture spec already includes `segment_curves`. Heritage
shapes that create the characteristic Morpheus envelope feel are directly portable.

### Velocity-to-Morph Scaling (Note-On Source)
Heritage cord: Velocity -> FilFreq with amount scaling. TRENCH equivalent: a
`velocity_scale` float per body that scales TRIG gesture travel distance. Hard hits =
full morph excursion; soft hits = subtle movement. This is the single cheapest
heritage feature to add — one float, one multiply at note-on.

### Aftertouch-to-Q (Realtime Source)
Heritage cord: Pressure -> FilRes. The most expressive Morpheus patches used this.
TRENCH equivalent: map poly/mono aftertouch to Q offset. Gives keyboard players
continuous pressure control over resonance aggression during sustained notes.
Combined with per-slot pressure behaviors (Hold, Tighten, Yield, Cross, Collapse),
aftertouch becomes a performance axis — not just volume control.

### Negative Cord Amounts
E-mu cords scaled -128 to +127. Negative amounts inverted the modulation. TRENCH's
REV toggle is the static equivalent (morph inversion). The heritage system allowed
per-cord inversion: one cord could drive morph forward while another drove it
backward simultaneously. In TRENCH terms, this would mean ENV could push morph in
the opposite direction from TRIG gesture travel — the body breathes one way and
snaps the other. This is naturally supported if ENV offset and TRIG offset are
summed independently (already specified in v1 architecture).

---

## Priority Assessment

### Ship path (v1)

- ENV -> morph (envelope follower, one-pole, per-body timing profiles)
- REV (morph axis flip, toggle)
- TRIG (one-shot authored gesture from anchor)
- Slew limiter (max morph delta = 0.15 per block)

### v1.x extensions

- ENV split routing: morph + Q (one extra authored float per body)
- Velocity-to-morph scaling (one float per body, scales TRIG travel)
- Aftertouch-to-Q (continuous pressure control for keyboard players)
- Continuous REV — only if 2-3 bodies prove it musically; not automatic promotion

### Research path (ordered by expected yield)

1. **Diagonal gestures on morph/Q surface** — best core research bet.
   The cross-term `Dmq` is real math, not cosmetic. Diagonal motion accesses
   spectral behavior that serial morph-then-Q sweeps cannot reproduce.
   This is an instrument-level expansion, not a modulation add-on.

2. **TRIG cycling** — authored gesture sequences (bark/shimmer/choke) from
   a single MIDI note stream. Fits the grammar: bodies define the cycle,
   user just plays.

3. **Per-slot morph staggering** — research-only until null-tested and
   listen-tested. Risk: transient impossible-states may read as vague blur
   rather than musical micro-detuning. Do not promote without audio proof.

4. **Codec-derived adaptive envelope** — only after the simpler one-pole
   follower proves it earns its keep in shipped bodies.

### Explicitly deferred

- QSound as independent modulation destination (spatial punch) — secondary
  until the core spectral modulation language is undeniable.
- Waveshaper drive modulation — implicit behavior may be sufficient.
- Open patchcord routing — violates the modulation grammar. Domain 3 only.

---

## Risk Assessment

| Idea | Risk | Gate |
|------|------|------|
| Diagonal gestures | Low — exploits existing math | Prototype on one body, measure Dmq magnitude |
| ENV -> Q split | Low — one float, cheap | Prove "fights back" is audible, not just conceptual |
| Per-slot stagger | Medium — architecture churn for potentially inaudible effect | Null test + listen test before any cascade changes |
| Continuous REV | Medium — midpoint collapse can be gimmick | 2-3 bodies must prove musical midpoint behavior |
| TRIG cycling | Low — state machine, contained scope | Rhythmic test with real MIDI input |
| Codec adaptive ENV | Low technical risk, high diminishing-returns risk | Ship simple follower first, compare |
