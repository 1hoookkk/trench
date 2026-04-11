# TRENCH

TRENCH is a cartridge-based morph filter plugin. It is not a parametric EQ, not an RBJ filter playground, and not a generic morphing UI.

## 1. ROLE

- I own sound, taste, UX, product identity, and final musical judgment.
- You own code, math, implementation, debugging, and verification.
- Do the work. No preamble. No planning doc. No vague summaries.

## 2. CORE MENTAL MODEL

- A TRENCH filter is not one static curve. It is movement between authored spectral states.
- Think in terms of `Frame A` and `Frame B`:
  `Frame A` = closed / dark / low state
  `Frame B` = open / bright / aggressive state
- `MORPH` is travel between states.
- `Q` is aggression / emphasis behavior across that travel.
- `TYPE` selects the authored body.
- Judge filters by motion, not by static plots alone.

## 3. AUTHORING RULES

- Author the gesture, not the screenshot.
- `Shelf` is a tone-shape continuum, not a rigid textbook selector.
- `Freq` is context-dependent and only makes sense after `Shelf` is chosen.
- `Peak` is not just loudness. It controls bite, emphasis, and perceived energy through the sweep.
- Good body:
  `Frame A` works
  `Frame B` works
  the path between them is intentional
- Bad body:
  dead middle
  accidental harsh jump
  level collapse
  fake excitement caused by bad gain staging

Authoring order:
- choose the musical role
- define `Frame A`
- define `Frame B`
- set `Shelf`
- set `Freq`
- set `Peak`
- sweep-test the motion

Never:
- tune `Freq` before `Shelf`
- judge from plots alone
- use `Peak` to hide bad spectral behavior

## 4. FROZEN DSP INVARIANTS

- Topology: 12-stage serial DF2T cascade. Exactly 6 active + 6 passthrough.
- Math: kernel-form `[c0, c1, c2, c3, c4]` interpolation only.
- Order: 4-corner bilinear interpolation. Evaluate Q first, then morph.
- Runtime: 32-sample control blocks with per-sample coefficient ramping.
- Rates: authoring = `39062.5 Hz`, plugin runtime = `44100 Hz`.

Exact DF2T form:
- `y = c0*x + w1`
- `w1 = c1*x - c3*y + w2`
- `w2 = c2*x - c4*y`

## 4b. AUTHORING VOCABULARY

The heritage authoring interface is 6 stages × 5 integers on a 0–127 grid:
`(type, low_freq, low_gain, high_freq, high_gain)` per stage.

- `low_*` = morph=0 endpoint. `high_*` = morph=100 endpoint.
- `type` (0–3): selects the pole recipe (Type 1 resonator, Type 2 shelf, Type 3 formant anchor).
- The compiler converts these integers to stable cascade coefficients.
- This is the exact vocabulary E-mu engineers used to author the original 289 Morpheus cubes.
- Morpheus cubes are 2-endpoint (morph only), poles only, Q degenerate.
- P2K skins extended this to 4 independent corners (morph × Q) with zero placement.
- The runtime basis is `x(m,q) = A + q·Δq + m·Δm + mq·Δmq`.
- Δmq (the cross-term) is zero in Morpheus cubes, load-bearing in P2K skins.

Use the forge for body generation. Score candidates with the trajectory/continuity fitness stack, not by ear alone.

## 5. HARD BANS

- No RBJ cookbook for character filters.
- No compensation layers, fudge factors, or arbitrary saturation to rescue bad math.
- No parallel paths. The engine is rigid serial cascade only.
- No new abstractions unless they are directly justified by proven behavior or hard implementation need.
- No gain baked into `c4`.
- No pole sanitization unless a failing test proves instability.
- No new dependencies in `Cargo.toml` without asking.
- No shipping verbatim heritage extractions as presets.

## 6. WORKFLOW

- Pick the shortest path to a verified result.
- If the request is to build, implement first.
- If the audible path changes, state exactly what changed and why.
- Do not claim something works unless you ran it.
- Measure -> characterize -> validate -> then promote to policy.
- Read `trench-core/CLAUDE.md` before editing `trench-core/`.

## 7. REPO MAP

- `trench-core/` = DSP spine, interpolation, cascade, ramping
- `trench-codec/` = minifloat codec
- `pyruntime/` = forge / authoring / solver / analysis tools
- `cartridges/` = compiled body JSON (P2K extractions, heritage quarantine)
- `vault/` = forge output, shipping finalists, scorecards, diagnostics
- `tools/` = scripts and validation harnesses
- JUCE plugin = separate repo at `../trench-juce/plugin/` (shipping path)

## 8. STOP AND ASK

Stop and ask if:
- the DSP or audible path changes materially
- a shortcut creates architectural debt
- it is unclear whether behavior is a bug or a deliberate betrayal in the filter design
