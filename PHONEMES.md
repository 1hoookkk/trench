# TRENCH authoring model

The only authoring flow that ships: place phoneme tokens on a time grid.
A compiler resolves the grid into a `compiled-v1` cartridge. The user never
touches stages, poles, zeros, frequencies, or Q values directly.

## Token vocabulary

Tokens come from the unified phoneme inventory at
`cartridges/engine/_source/token_inventory_unified_v2.json`. Each
token is a static sonic location whose `compiled-v1` shape file lives
under `cartridges/engine/_source/shapes/<category>/<key>.json`.

A token name is a stable string. The compiler resolves names to
spectral state via the inventory. Unknown token names fail the
compile.

The shipping flat pill layout at `cartridges/engine/<category>/<key>.json`
is a derived rebuild, produced by `tools/bake_phoneme_pills.py` from
the `_source/` subtree. Do not hand-edit the flat layout — edit
`_source/` and re-bake.

## Grid

The authoring surface is a Looperator-style grid of phoneme blocks on a
time axis. Each block is a `(token, duration, optional motion)` triple.

The grid maps onto the 4-corner morph surface of `compiled-v1`:

| grid position   | corner        |
|-----------------|---------------|
| start, low Q    | `M0_Q0`       |
| start, high Q   | `M0_Q100`     |
| end,   low Q    | `M100_Q0`     |
| end,   high Q   | `M100_Q100`   |

Between blocks the compiler schedules transitions so `MORPH` sweeps through
the sequence. `Q` and `TYPE` are playback surface controls; the cartridge
stores both Q corners.

## Compiler contract

Input: a grid (time-ordered token list + Q variants).
Output: a `compiled-v1` cartridge that validates against
`cartridge.schema.json` (see `SPEC.md` §2).

The compiler MUST:
- Pick tokens from the canonical inventory only.
- Emit 4 keyframes labeled `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`.
- Emit 12 stages per keyframe (6 active + 6 passthrough sentinel).
- Preserve per-body invariants from `BODIES.md` for shipping candidates.

The compiler MUST NOT:
- Expose pole/zero/frequency grids to the user.
- Accept raw coefficient edits.
- Derive stages from RBJ cookbook shapes. Direct pole-zero only (see
  `DOCTRINE.md`).

## What's built, what's TBD

- Token inventory: built (`cartridges/engine/_source/token_inventory_unified_v2.json`).
- Shape bank (token bodies): built (`cartridges/engine/_source/shapes/`).
- Shipping flat pill layout: built (`cartridges/engine/<category>/<key>.json`
  via `tools/bake_phoneme_pills.py`).
- Grid UI: TBD. v1 goal.
- Compiler (grid → cartridge): TBD. v1 goal.

Until the compiler ships, candidate pills come from the forge (sibling
authoring repo, out of shipping tree) and are hand-audited against
`BODIES.md` rubrics and the `trench-core/tests/vowel_formants.rs`
classification gate.
