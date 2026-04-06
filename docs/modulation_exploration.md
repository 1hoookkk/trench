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

### Heritage (E-mu Morpheus / UltraProteus)

From RE vault and DillusionMan documentation:

- **2 Function Generators**: Multi-segment envelopes (attack/hold/decay with shape control)
- **Patchcord Routing**: Source -> Destination with amount. e.g., `MidiA -> FilFreq +100%`
- **Mod Sources**: Velocity, aftertouch, MIDI CC, LFO, function generators
- **Mod Destinations**: FilFreq (morph position), FilRes (resonance), Shelf, Peak
- **Morph Control**: CC-driven sweep between Frame A (CC=0) and Frame B (CC=127)
- **Real-time**: All patchcords evaluated per MIDI event or per control block

TRENCH v1 deliberately does NOT ship a patchcord system. Modulation is pre-authored per body.

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

### Function Generator Shape Library
E-mu function generators supported segment shapes: linear, exponential, logarithmic, and "circle" curves. The TRIG gesture spec already includes `segment_curves` — but we could expand this with heritage-derived shapes that create the characteristic Morpheus envelope feel.

### Velocity-to-Morph Scaling
Morpheus patchcords allowed velocity to scale morph depth. A simple `velocity_scale` parameter per body would let MIDI velocity control how far TRIG gestures travel. Hard hits = full morph excursion; soft hits = subtle movement.

### Aftertouch-to-Q
The most expressive Morpheus patches used aftertouch -> FilRes. Mapping aftertouch to Q in TRENCH would give keyboard players continuous pressure control over resonance aggression during sustained notes.

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
