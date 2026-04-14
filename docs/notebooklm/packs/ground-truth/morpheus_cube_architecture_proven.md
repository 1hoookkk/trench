---
name: Morpheus Cube Architecture — Proven
description: 289 Morpheus cubes are authentically poles-only. No hidden zeros. 14 Morpheus classes are ROM-baked. P2K skins use separate Morph1/2/LP/LPX lattice compiler. Cubes are pole seeds for forge zero-forcing.
type: project
---

## Proven Facts (2026-03-16)

### Factory Architecture
- EX3 factory FUN_1802ae650 creates 46 filter types
- Types 0-13: 14 Morpheus-era classes (LP2Pole through Batman) — all ROM-baked fixed coefficients
- Types 14-45: 32 CPhantomFilterP2k instances — P2K wrapper dispatching to Morph1/2/LP/LPX/Designer/Flanger/Vocal sub-classes

### Morpheus Classes (types 0-13) — ROM-Only
- Compile methods ignore preset/cube data entirely
- Copy precomputed coefficients from ROM tables indexed by SR family
- Batman: 2 stages, Vocal1/2: 3 stages, LP/HP/BP: 2-3 stages
- These are the conventional filter section, NOT the Z-Plane morph engine

### 289 Morpheus Cubes — Poles Only By Design
- Cubes bypass the class factory. No compile step.
- Raw 4-corner × 6-stage pole data loaded directly into runtime interpolation engine
- No zeros computed at any stage. The Morpheus Z-Plane architecture is authentically poles-only.
- The cube export format isn't "missing" data — it's complete.

### P2K Lattice Compiler (FUN_1802c59b0)
- Shared inner compiler for Morph1/MorphLP/MorphLPX classes
- 3-stage cross-coupled lattice: stage N zeros = stage N-1 poles
- Numerator table (param_5) only affects stage 0; stages 1-2 always cross-coupled
- Morph1: bypass numerator [0xDFFF, 0xFFFF, 0xDFFF] = fixed zero law
- MorphLP: 16-entry LP table (per-corner zero variation) — ported to historical.rs
- MorphLPX: 16-entry LPX table (morph-invariant) — ported to historical.rs
- Morph2: separate compiler with own algebra, does NOT use the lattice

**Why:** Closes the months-long investigation into "missing" Morpheus cube zeros. They never existed. The forge's zero-forcing is additive (TRENCH feature), not restorative.

**How to apply:** Use 289 cubes as pole seeds for forge generation. Cross-breed cube poles with forge zero-forcing. The cubes are vocabulary, not incomplete data.
