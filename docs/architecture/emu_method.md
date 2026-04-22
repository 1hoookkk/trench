# The E-MU Method — Authoring Plan

> The sound designer only authored two states: **the start of the sound (M0)
> and the end of the sound (M100)**. That was their entire job.

This doc is the plan for bringing TRENCH authoring in line with how E-MU
actually worked. It is the authoring-side companion to
`zplane_truth.md` (runtime) and supersedes `cube_authoring_path.md`
(the cube premise was wrong — see §"What we got wrong"). The active
in-flight prototype for one family is `PLAN_RECIPE_V1.md` (peak/shelf
recipe-v1 in the JUCE plugin); this doc generalizes that pattern to
every family.

## The core insight

In the P2K-era decompiled Morph Designer, the authoring surface was two
morph endpoints and a filter family. The 4-corner grid that ships to
the DSP engine was a **compile artifact**, never an authoring unit. Q
was never a human-drawn axis — it was algebra baked into a named
compiler class (e.g. `CPhantomMorphLPX`) that took (M0, M100, q) and
produced coefficients at any Q value.

TRENCH inherited the *output* of that pipeline (4-corner bilinear
`compiled-v1`) and lost the *input* (2 states + named Q law). Every
authoring-side drift in the repo — Q-axis editing, cubes, MD presets
storing 4 corner dumps — traces back to editing the wrong layer.

## The proper pipeline

```
human → (M0 shape, M100 shape, family) → Q-law bake → compiled-v1 runtime
```

Three concrete things a designer touches. Nothing else.

- **M0 shape** — spectral geometry at the start of the morph sweep
- **M100 shape** — spectral geometry at the end of the sweep
- **family** — names the Q law/compiler family used to bake the runtime surface

The Q law is **code**, not data. On the current repo path it lives in
offline compiler tooling, not in `trench-core`. Each family is one
named, versioned, deterministic compiler rule:

    (m0, m100, morph ∈ [0,1], q ∈ [0,1]) → coeff_corner

For the shipping `compiled-v1` path, that family logic is baked
offline into four corners and the runtime only interpolates those
corners. Morph and Q are runtime controls; family-aware algebra is not.
Q is never peer to Morph in authoring — it is a secondary modifier the
offline family compiler responds to.

## UI consequence

Morph is the primary knob. Q sits beneath it, visually subordinate.
There is no Z. The current `trench-live` plugin has Morph + Q + Cube Z
in roughly equal visual weight — that is a symptom of the old "Q as
second axis" model and gets corrected in stage 1.

## Staged migration

The plugin must load, play, and pass tests at every commit. These
stages are ordered so nothing is broken in the middle.

### Stage 1 — Cube strip (this branch, `cube-prep-cleanup`)

**Goal:** reduce the plugin to its two real knobs.

- Strip `SurfaceFormat::CubeSurfaceV1`, `SurfaceData::Cube`, trilinear
  interp, `cube_authoring.rs`, cube integration tests, Z parameter,
  `process_block_xyz*` API.
- UI: **Morph on top, Q directly beneath it, no Z.** One column, Morph
  visually primary.
- Delete (do not commit) `CUBE_GATE.md` — under the proper method,
  cubes aren't gated, they're a category error.
- Runtime DSP, cascade, resampler, AGC, DC blocker: **untouched**.
- `cartridges/p2k/*` P2K heritage presets: **untouched**.
- `compiled-v1` format: **untouched**. Still the runtime surface.

### Stage 2 — Consolidate the offline family compilers

**Goal:** make the heritage algebra explicit and addressable.

Keep family logic out of `trench-core`. The runtime stays a pure
`compiled-v1` loader/interpolator. Family compilers live on the offline
side (`tools/` today; a future authoring subtree is acceptable only if
it preserves the same boundary).

Seed this layer from the family compilers and recipes already present
or implied by the repo:

- `morph_lp` — implemented offline at `tools/compile_morphlp.py`,
  backed by `tools/morphlp_zeros.py` and verified against
  `trench_re_vault/analysis/emulatorx_re_dumps_2026-04-21.json`
- `peak_shelf` — Peak/Shelf Morph recipe (prototype is
  `PLAN_RECIPE_V1.md`, JUCE C++ `EmuKernels` port of the Python
  reference in `tools/build_filter_reference.py`)
- `phantom_lpx` — P2K phantom family per RE findings
- `native_poles` — pole-oriented path using
  `trench-core/src/pole_block.rs` as a math helper, not as a family
  registry
- `md_generic` — heritage Morph Designer fallback

Runtime and cartridge formats are unchanged in this stage. Existing
`compiled-v1` cartridges continue to load. This stage does not add
runtime family awareness. It makes the offline compiler surface honest
and testable. `PLAN_RECIPE_V1.md`'s C++ recipe-v1 loader becomes the
`peak_shelf` family entry in the offline compiler layer.

### Stage 3 — Authoring format: 2 states + family

**Goal:** copy E-MU exactly on the authoring side.

New authoring format `*.filter.json` carrying only
`{family, m0, m100, meta}`. The `compile_filter` tool takes this + the
family's Q law → emits existing `compiled-v1` for the runtime. No
runtime change required.

This is a target authoring surface, not the current MorphLP input
contract. The current cleanroom MorphLP compiler already in-tree still
takes narrower compiler inputs (`filter_type`, `skin_idx`, `q0_byte`,
`q100_byte`) and emits `compiled-v1`. Do not describe the new source
schema as already landed until a real schema and dispatcher exist.

- Migrate `cartridges/p2k/*` to the new source format; regenerate
  `compiled-v1` outputs via the offline family compilers; verify byte-equal
  or sound-equal against current shipping.
- Rewrite `heritage_designer_sections.json` as MD-family source pairs
  under the same schema.
- Retire any tool that authors 4 corners directly.

### Stage 4 — 2-corner + live-Q runtime

**Goal:** runtime matches E-MU exactly too.

The transition described in `zplane_truth.md` §"The 4-corner → 2-corner
transition." Ship 2 corners + a Q-law identifier; the runtime calls
the Q law live per block. Byte-preserves Stage-3 authoring pairs.
4-corner cartridges remain loadable via one-shot bake path for
back-compat.

## What we got wrong (archive)

- **Cubes.** Built on the premise that humans author coefficient
  corners. E-MU designers never did. 8 corners is not "more
  expressive" than 2 states + a Q law — it is the same mistake scaled
  up. Strip, don't gate.
- **Q-axis as authoring dimension.** MD had no Q-axis authoring. Every
  tool that presented Q as a corner to edit was editing the wrong
  layer.
- **Pill as a noun.** Internally "pill" meant "4-corner compiled
  cartridge." That confused authoring with its build artifact. The
  user-facing noun is **filter**; the internal file rename is a
  mechanical follow-up.

## What stays untouched

These subsystems are already right-premise. They don't change in any
stage of this plan:

- `trench-core/src/cascade.rs` — frozen DSP topology
- `trench-core/src/resampler.rs`, `agc.rs`, `desk_drive.rs` — pure DSP
- `trench-core/src/hedz_rom.rs`, `hedz_golden.rs` — baked heritage
- `trench-core/src/pole_block.rs` — clean math helper for pole-bank
  compilation; not evidence that family registries belong in runtime
- `trench-core/src/qsound_spatial.rs`, `function_generator.rs` —
  orthogonal subsystems
- `cartridges/p2k/*` — heritage data

## References

- `docs/architecture/zplane_truth.md` — runtime side of this transition
- `PLAN_RECIPE_V1.md` — live prototype of Stage 2 for the `peak_shelf`
  family in the JUCE C++ plugin
- `tools/compile_morphlp.py` — existing cleanroom MorphLP/MorphLPX
  offline compiler
- `tools/verify_morphlp_re.py` — verification bridge from sanitized
  tables to the RE vault artifact
- `trench_re_vault/` — P2K decomp source material (filter family
  classes, Q-law algebra)
- `authoring/CLAUDE.md` — forge scope and boundary with engine
