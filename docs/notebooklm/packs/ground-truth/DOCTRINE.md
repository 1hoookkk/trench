# TRENCH Technical Doctrine

## What TRENCH is

A filter that lies to the user in musically useful ways.

12-stage DF2T biquad cascade. Each body is a compiled program — immutable at runtime,
authored by the forge. The user has three controls: morph, Q, and body selection.

## Betrayal

The internal term for TRENCH's core behavior. Not dishonesty, not counter-intent,
not false witness. Betrayal.

Betrayal means: the control surface promises one thing, the engine delivers
something more interesting. Stable, repeatable, musically legible. Feels
intentional after use.

Four layers:
1. **Control** — knob says one thing, engine hears another.
2. **Structural** — stages do unequal jobs.
3. **Spectral** — displayed vs actual energy movement diverge.
4. **Temporal** — behavior depends on motion history.

Anti-patterns: random modulation, arbitrary saturation, unexplained gain jumps,
anything that destroys repeatability or flattens body identity.

## Zeros

Zeros are body-level, not per-stage authored. The solver owns contradiction.
Counterpoint (per-stage zero authoring) is dead. Zeros emerge from global
body-level fitting via the solver.

## Heritage data

Cubes, skins, and vocabulary extracted from heritage hardware are donor material.
Not sacred. Not end-state. They seed the forge. The forge owns what ships.

## Surfaces

- **TRENCH** — the plugin. 3 controls, loads compiled cartridges. Runtime only.
  44100 Hz. compiled-v1 format (4 corners).
- **TRENCHWORK** — the workbench. Generates, scores, curates bodies. Authoring only.
  39062.5 Hz. hyper-v1 format (32 corners).

These are separate products with separate identities. Do not blur them.

## Topology

- 12-stage serial DF2T cascade. Not parallel. Never parallel.
- Interpolation: Q first, then morph. Kernel-form [c0..c4] only.
- 32-sample control blocks, per-sample coefficient ramping.
- Minifloat: 4-bit exp (bias 15), 12-bit mantissa, ldexp decode.
- Stability by construction (r < 1). No pole sanitization unless a test proves it needed.
- Authoring compiles into immutable runtime programs. No in-place mutation.
