# Canonical Filter Work

This repo is the canonical repo for TRENCH filter work:

```text
C:\Users\hooki\Trench
```

The canonical doorway is:

```text
C:\Users\hooki\Trench\authoring
```

Do not create a fourth filter repo. Do not argue about provenance taxonomy as
the main policy. The only shipping rule that matters is this:

```text
TRENCH does not ship E-mu filters.
```

References can inform taste, vocabulary, measurements, and regression checks.
The shipped cartridge must be original TRENCH work.

## Shipping Veto

A filter is blocked from shipping if any of these are true:

1. It is copied from an E-mu/P2K/Morpheus/Emulator factory preset.
2. It is a direct export, transcription, or lightly renamed version of an
   E-mu-derived cartridge.
3. It nulls like a clone against an E-mu/P2K canonical reference.
4. Its identity depends on preserving a named E-mu filter behavior instead of a
   TRENCH body goal.
5. Its source trail cannot explain how the final numbers were authored.

The goal is not to prove ancestry. The goal is to prevent shipping clones.

## What References Are For

Allowed:

- learning the shape of the design space
- understanding manual vocabulary
- debugging renderers and compilers
- measuring distance from known references
- rejecting accidental clones

Not allowed:

- shipping copied coefficients
- shipping factory-preset reconstructions
- treating a deep E-mu/P2K null as success for a product cartridge
- using E-mu names as the product identity of a shipped filter

## Authority Order

1. Shipped plugin assets and runtime wiring
2. Active authoring GUI/compiler outputs
3. Local tests that prove the candidate loads and behaves safely
4. Listening decisions against TRENCH body goals
5. Reference comparisons used as vetoes, not targets
6. Old docs, external trees, screenshots, and notes

If a reference comparison says "this is the same as E-mu," the candidate fails.
If a doc says parity is the product target, the doc is wrong.

## External Tree Roles

| Tree | Use it for | Do not use it for |
|---|---|---|
| `C:\Users\hooki\Trench` | active filter work, compiler/GUI/runtime checks, shipping candidates | dumping unreviewed material |
| `C:\Users\hooki\trenchwork_clean` | old tools, sonic tables, experiments, candidate ideas | declaring what ships |
| `C:\Users\hooki\trenchwork` | workbench architecture ideas | declaring what ships |
| `C:\Users\hooki\trench_re_vault` | research evidence and reference behavior | shipping cartridge content |
| `C:\Users\hooki\do-it` | historical precedent | product identity |

## What Belongs Here

Everything filter-work-related lives under `authoring/`. The legacy `forge/`
and root-level `schemas/` directories have been folded in.

Under `authoring/`:

- `compilers/` — active compilers, smoke tests, authoring utilities that
  generate original artifacts
- `tools/` — authoring utilities that are not compilers
  (audition harnesses, body authors, phoneme studios, proof scripts)
- `audition.py` — canonical audition pipeline
- `schemas/` — TRENCH-native JSON schemas
- `bodies/`, `design/`, `frame_maps/`, `sonic_tables/`, `workbench/` —
  original body drafts, IRs, and design surfaces
- `ground_truth/`, `heritage/`, `imported_x3/`, `_research/` — research and
  reference material that is *not* a shipping path
- `auditions/`, `compiled/`, `plots/` — audition output and plots
- `inbound/`, `dirty_airlock/` — staging for outside material
- `history/` — historical precedent summaries
- `q_laws/`, `diagnostics/`, `scratch/` — engine-axis tooling and working files

Under the runtime/plugin trees (NOT `authoring/`):

- runtime code (`runtime/trench-core/`)
- shipped plugin (`trench-juce/`)
- shipped cartridge assets (`cartridges/`)
- load/process tests (`runtime/trench-core/tests/`)
- release gates

Under `docs/` (research/reference, not authoring):

- calibration data, heritage manuals, architecture notes

## What Does Not Belong Here

- unreviewed copied preset banks
- raw research dumps pretending to be product assets
- generated coefficient blobs with no authoring story
- E-mu clone candidates in release-candidate folders
- another unindexed pile of "maybe useful" files

## Promotion Rule

Material from anywhere else becomes part of active filter work only when it has:

1. A destination in this repo.
2. A source path recorded.
3. A stated purpose: `reference`, `candidate`, `tool`, `test`, or `shipping`.
4. A verification command or reason no verification was run.
5. For `shipping`, an explicit statement that it is original TRENCH work and
   not an E-mu/P2K clone.

## Current Active Filter Path

The first verified internal slice is the EOS encoded-surface harness:

```text
trench_eos_pipeline_compiler.py                       # root launcher (delegates)
authoring/compilers/trench_eos_pipeline_compiler.py   # real harness
schema: trench.eos.encoded_surface.v1
default out: cartridges/factory/generated/live_eos_surface.json
```

The harness models the four-stage runtime shape — 6-byte stage descriptors
→ encoded uint16 quintuplets → bilinear c-domain interpolation → decoder
→ 6-stage serial biquad render kernel — and exposes explicit replacement
seams for the real lowering, decoder, and kernel code. It is **not** a
playable compiled-v1 cartridge writer. It does **not** invent E-MU
descriptor lowering, coefficient encoding, decoder behavior, or the render
equation; those must come from verified reverse engineering, not from the
harness itself. Per the doc Authority Order, this is an active
authoring/compiler output (#2), not a shipped asset (#1) or a reference
comparison (#5).

Verification commands (passing as of 2026-05-10):

```powershell
python -m py_compile trench_eos_pipeline_compiler.py authoring\compilers\trench_eos_pipeline_compiler.py
python trench_eos_pipeline_compiler.py example-fixture ...
python trench_eos_pipeline_compiler.py interp ...
python trench_eos_pipeline_compiler.py fixture ...
cargo test -p trench-core --test compile_raw_roundtrip --release -- --nocapture   # 4/4 passing
```

The previous "Peak/Shelf" slice was deleted on 2026-05-10 because that name
identifies an E-mu Morpheus filter class (DillusionMan tutorial: SHELF
blends LP↔mid-shelf↔HP, PEAK is master volume, 06-order topology) and
authoring against it is shipping-vetoed by rules #1 and #4. The cascade
plot of the disk law also inverted the authored intent at the time of
deletion. Pre-veto research IRs are archived under
`authoring/_research/peak_shelf/` for vocabulary and reject-clones use
only — see that directory's README.

The shipping rule is unchanged: until a candidate cartridge is produced
from this slice and passes the no-clone gate, nothing here is a shipping
artifact. The harness is internal scaffolding, not a release path.

## Immediate Consolidation Moves

1. Inventory outside material by usefulness, not purity.
2. Promote only files needed by active filter work.
3. Add a release gate that rejects E-mu/P2K clone candidates.
4. Make every new filter workflow start from this repo and this document.
