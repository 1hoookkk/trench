# Compiler spec — grid → cartridge

**Status:** implemented at `tools/compile_grid.py`; gated by
`trench-core/tests/compile_grid_roundtrip.rs` in `./check`.
**Authority:** `PHONEMES.md` §Compiler contract, `cartridge.schema.json`,
`trench-core/src/cartridge.rs::Cartridge::from_json`.
**Scope:** in-tree Python compiler that reads an authored grid and writes a
`compiled-v1` cartridge JSON that loads cleanly into the Rust runtime.

## 1. Inputs

1. **Grid**: a JSON document with exactly 2 token entries (start + end):
   ```json
   { "name": "cul_de_sac",
     "grid": [
       { "token": "vowels.ah" },
       { "token": "vowels.ee" }
     ]
   }
   ```
   Token names use the `<category>.<key>` form that matches inventory keys.
2. **Token inventory**: `cartridges/engine/_source/token_inventory_unified_v2.json`.
   The compiler reads the `tokens` dict and resolves each token's
   `(category, key)` pair. Unknown tokens → compile error. The inventory's
   `shape_path` field is intentionally ignored (stale vault/ references from
   pre-excision cleanup); paths are computed from `category` + `key`.
3. **Shape bank**: `cartridges/engine/_source/shapes/<category>/<key>.json`.
   Each shape is itself a `compiled-v1` cartridge (morph-invariant: M0 corners
   == M100 corners; only Q varies).

## 2. Output

A `compiled-v1` cartridge that passes `Cartridge::from_json`:
- `format: "compiled-v1"`
- `name: <body name>`
- `keyframes`: 4 entries, labels `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`
- Each keyframe: `label`, optional `boost` (default 1.0), `stages` (6 or 12;
  if 12, trailing 6 must be passthrough `[1,0,0,0,0]` within 1e-10)

Validation fails closed: any schema violation or unknown token aborts the
compile with a non-zero exit and a diagnostic to stderr.

## 3. Contract (verbatim from `PHONEMES.md`)

MUST:
- Pick tokens from the canonical inventory only.
- Emit 4 keyframes labeled `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`.
- Emit 12 stages per keyframe (6 active + 6 passthrough sentinel).
- Preserve per-body invariants from `BODIES.md` for shipping candidates.

MUST NOT:
- Expose pole/zero/frequency grids to the user.
- Accept raw coefficient edits.
- Derive stages from RBJ cookbook shapes. Direct pole-zero only.

## 4. Algorithm

```
compile(grid_doc) -> compiled-v1 cartridge:
    1. require len(grid_doc.grid) == 2; else raise (see §5)
    2. resolve first_token, last_token against inventory; else raise
    3. first_shape = load cartridges/engine/_source/shapes/<cat>/<key>.json
    4. last_shape  = same, for last token
    5. emit keyframes:
         M0_Q0    ← first_shape[M0_Q0]
         M0_Q100  ← first_shape[M0_Q100]
         M100_Q0  ← last_shape[M100_Q0]
         M100_Q100← last_shape[M100_Q100]
    6. carry per-keyframe `boost` field from the source shape (default 1.0)
    7. serialize as compiled-v1 JSON
```

Because shape files are morph-invariant (each pill has M0 corners identical
to M100 corners; only Q varies), step 5 is well-defined: taking M0 corners
from the first shape and M100 corners from the last shape encodes exactly
one morph transition per output cartridge.

Reference implementation: `tools/compile_grid.py::compile_grid`.

## 5. Resolved — strict 2-token grid (Path A, no silent drops)

The cartridge format stores 4 corners only; there is **no** schedule,
sequence, or trajectory field in `cartridge.schema.json` or in the
`Cartridge` struct (`trench-core/src/cartridge.rs:66-73`).

**v1 decision:** the compiler enforces a **strict 2-token limit**. A grid
with 0, 1, 3, or more tokens is rejected with a non-zero exit code and an
error message naming the 2-token contract. Middle tokens are never silently
dropped — ambiguity about what the compiler did must be impossible.

This unblocks the 4 normative bodies in `BODIES.md` immediately without
touching the schema, the Rust loader, or any existing gate. Multi-waypoint
compilation is deferred to a future `compiled-v2` schema evolution; no
design work for (B) or (C) is in v1 scope.

Enforcement: `tools/compile_grid.py::compile_grid` (`len(grid) != 2` → raise)
and `trench-core/tests/compile_grid_roundtrip.rs::three_token_grid_is_rejected`
pin the limit end-to-end.

## 6. Test strategy

Tests live in `trench-core/tests/compile_grid_roundtrip.rs` and run under
`./check` in the `rust compile_grid round-trip` step. The Rust test shells
out to `python tools/compile_grid.py`, captures stdout, and feeds it through
`Cartridge::from_json` so any drift between the Python emitter and the Rust
loader fails the gate.

Shipped:
1. **Round-trip** (`two_token_grid_round_trips_through_rust_loader`). Compile
   `vowels.ah` → `vowels.ee`; assert loader accepts, 4 corners, 4 boosts.
2. **Strict 2-token limit** (`three_token_grid_is_rejected`). A 3-token grid
   must exit non-zero with an error naming "exactly 2 tokens".
3. **Unknown token** (`unknown_token_is_rejected`). An unresolvable token
   must exit non-zero with an error naming "unknown token".

Deferred (tracked in BODIES.md / vowel_formants.rs, not gated here yet):
4. **Per-body invariant preservation.** For each of the 4 `BODIES.md`
   bodies, compile a representative grid and assert the body-specific rubric
   (e.g., Speaker Knockerz sub anchored; Small Talk Ah-Ee two dominant
   formants) via the same classification harness that
   `trench-core/tests/vowel_formants.rs` uses.
5. **Parity against hand-authored.** When a forge-authored body cartridge
   exists, compile its grid form and byte-diff against the hand-authored
   cartridge (tolerance 1e-10 per coefficient).

## 7. Host location

**Python `tools/compile_grid.py`**, stdlib-only. Matches
`bake_phoneme_pills.py` / `bake_hedz_const.py` pattern. Round-trip parity
with the Rust loader is enforced by
`trench-core/tests/compile_grid_roundtrip.rs`, which shells out to the
Python CLI and loads the output via `Cartridge::from_json` inside a Rust
test. This preserves the single source of truth for cartridge validation
(trench-core) while keeping iteration fast on the authoring side.

## 8. Non-goals for v1

- Grid UI — separate v1 goal, not part of this spec.
- Multi-cartridge sequencing and >2-token grids — deferred to a future
  `compiled-v2` schema evolution (see §5).
- Any extension to `cartridge.schema.json` — explicitly out of scope.
- Inventory `shape_path` field cleanup — stale but harmless; the compiler
  ignores it, so fixing it is an isolated future chore.
- Pill re-authoring of `aw`/`oo` (tracked separately in `vowel_formants.rs` KNOWN_DRIFT).
