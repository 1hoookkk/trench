# MorphLP / MorphLPX Compiler Contracts (implementation surface)

This file is the authoring-surface contract for the shared lattice compiler path:
`CompileStagePackets3` (`FUN_1802c59b0`).

## Evidence anchors (proved)

- `UpdatePackets_ConstTemplateShift` (`FUN_1802c5e40`): calls
  `CompileStagePackets3(..., &DAT_1806d73c0 + index * 0xc)`.
- `UpdatePackets_LUT6` (`FUN_1802c5f10`): computes `index`, reads 3 words from
  `DAT_1806d7480 + index * 6`, duplicates them for both corners, then calls
  `CompileStagePackets3`.
- `CompileStagePackets3` (`FUN_1802c59b0`): reads `DAT_1806d73a0` and
  `DAT_1806d73b0`, writes 3 active stages, sets `*(param_3 + 0x420) = 3`.
- Xrefs to `DAT_1806d7480`: `0x1802c5fad`, `0x1802c5fb4` inside `FUN_1802c5f10`.

## Morph indexing law (proved)

From both `FUN_1802c5e40` and `FUN_1802c5f10`:

```c
index = (short)(int)((*(float *)(param_3 + 0x1e8) + *(float *)(param_3 + 0x200)) * DAT_18065a1ec);
index = clamp(index, 0, 15);
```

- `DAT_18065a1ec = 16.0f` (proven from prior dump and decompile usage).
- Mapping is 16 discrete bins with floor cast and hard clamp to `[0,15]`.
- `param_3+0x200` is a second additive morph term (strong inference: modulation
  offset source; exact UI owner is outside this function).

## Table semantics (proved)

- `DAT_1806d73c0` (MorphLP): stride `0xC` bytes, 6 words per row.
  Layout: `[M0_v1, M0_v2, M0_v3, M100_v1, M100_v2, M100_v3]`.
- `DAT_1806d7480` (MorphLPX): stride `0x6` bytes, 3 words per row.
  Layout: `[v1, v2, v3]` and duplicated to both corners before entering
  `FUN_1802c59b0`.

### LPX zero table export (`DAT_1806d7480`, 16 rows, stride `0x6`)

Each row corresponds to one morph index bin. `v1 == v3` for all rows.

| Index | v1 | v2 | v3 | M0 triplet | M100 triplet |
|---|---|---|---|---|---|
| 0 | `0x11ED` | `0x94FC` | `0x11ED` | same as row | same as row |
| 1 | `0x24F5` | `0x9BFC` | `0x24F5` | same as row | same as row |
| 2 | `0x30F9` | `0xA3FC` | `0x30F9` | same as row | same as row |
| 3 | `0x3EF9` | `0xACFC` | `0x3EF9` | same as row | same as row |
| 4 | `0x4BFB` | `0xB3FC` | `0x4BFB` | same as row | same as row |
| 5 | `0x59FC` | `0xBCFC` | `0x59FC` | same as row | same as row |
| 6 | `0x65FC` | `0xC4FC` | `0x65FC` | same as row | same as row |
| 7 | `0x73FC` | `0xCCFC` | `0x73FC` | same as row | same as row |
| 8 | `0x80FC` | `0xD4FC` | `0x80FC` | same as row | same as row |
| 9 | `0x8EFC` | `0xDDFC` | `0x8EFC` | same as row | same as row |
| 10 | `0x9AFC` | `0xE4FC` | `0x9AFC` | same as row | same as row |
| 11 | `0xA8FC` | `0xEDFC` | `0xA8FC` | same as row | same as row |
| 12 | `0xB5FC` | `0xF2FD` | `0xB5FC` | same as row | same as row |
| 13 | `0xC3FC` | `0xF6FD` | `0xC3FC` | same as row | same as row |
| 14 | `0xCFFC` | `0xFAFD` | `0xCFFC` | same as row | same as row |
| 15 | `0xDDFC` | `0xFF7D` | `0xDDFC` | same as row | same as row |

Raw byte stream (little-endian) for auditing:

```text
ed11fc94ed11f524fc9bf524f930fca3f930f93efcacf93efb4bfcb3fb4bfc59fcbcfc59fc65fcc4fc65fc73fcccfc73fc80fcd4fc80fc8efcddfc8efc9afce4fc9afca8fcedfca8fcb5fdf2fcb5fcc3fdf6fcc3fccffdfafccffcdd7dfffcdd
```

## Builder input semantics (proved unless marked)

`CompileStagePackets3(longlong this, longlong stage_bytes, longlong out, short morph_delta_q, undefined2* zero_words)`

- `this + 0xc`:
  sample-rate index used to pick seeds from `DAT_1806d73a0` and `DAT_1806d73b0`.
- `stage_bytes`:
  contains four proved byte inputs at offsets `+0x16,+0x17,+0x19,+0x1a` and
  two pitch bytes at `+0x18,+0x1b`. The exact ownership of these six bytes in
  the higher-level authoring/UI surface is not proved by this function alone.
- `morph_delta_q`:
  extra signed offset term; LP/LPX callers pass `0`.
- `zero_words`:
  6 words consumed as stage-0 lane numerators. MorphLP passes one row from
  `DAT_1806d73c0`; MorphLPX duplicates one 3-word triplet from `DAT_1806d7480`
  to fill both lane triplets.

### Seed interaction (lattice law)

First-tier and second-tier pole seed words are built as:

```c
pitch_a = clamp(((*(char *)(stage_bytes + 0x18) - 0x40) * 0x80) + morph_delta_q * 0x100, -0x2000, 0x1fff);
pitch_b = clamp(((*(char *)(stage_bytes + 0x1b) - 0x40) * 0x80) + morph_delta_q * 0x100, -0x2000, 0x1fff);

u13 = slope[sr_idx] * (short)*(char *)(stage_bytes + 0x16) + base[sr_idx];
u14 = slope[sr_idx] * (short)*(char *)(stage_bytes + 0x17) + base[sr_idx];
u11 = slope[sr_idx] * (short)*(char *)(stage_bytes + 0x19) + base[sr_idx];
u12 = slope[sr_idx] * (short)*(char *)(stage_bytes + 0x1a) + base[sr_idx];

r13 = (u13 >> 1) + 0x6400;
r14 = (u14 >> 1) + 0x6400;
r11 = (u11 >> 1) + 0x6400;
r12 = (u12 >> 1) + 0x6400;
```

Cross-coupled stage construction in exact lane form:

- Stage 0 lane A words: `[u13, pitch_a + r13, zero_words[0], zero_words[1], zero_words[2]]`
- Stage 0 lane B words: `[u14, pitch_a + r14, zero_words[3], zero_words[4], zero_words[5]]`
- Stage 1 lane A words: `[u11, pitch_b + r11, u13, r13 - pitch_a, 0xDFFF]`
- Stage 1 lane B words: `[u12, pitch_b + r12, u14, r14 - pitch_a, 0xDFFF]`
- Stage 2 lane A words: `[0xDFFF, 0xFFFF, u11, r11 - pitch_b, 0xDFFF]`
- Stage 2 lane B words: `[0xDFFF, 0xFFFF, u12, r12 - pitch_b, 0xDFFF]`
- Stage count fixed by `*(out + 0x420) = 3`.

This is the proved interaction between seed tables, zero tables, and the
4-corner kernel source words.

## Shared vs corner-specific semantics

- Shared across all four corners:
  morph index law, table row selection, stage count (=3), seed tables, lattice
  cross-coupling structure.
- Lane-specific:
  byte inputs (`+0x16/+0x17/+0x19/+0x1a`), pitch bytes (`+0x18/+0x1b`), and
  derived pole/r words.
- LP only:
  `M0` and `M100` stage-0 zero triplets differ per row (`DAT_1806d73c0`).
- LPX only:
  same triplet duplicated to both corners (`DAT_1806d7480`).

## Blocked items

- Exact semantic owner of `*(float *)(param_3 + 0x200)` (which modulator in UI
  writes this field) is blocked in this pass. The additive role in index math is
  proved; UI provenance is not.
- Exact binding from the compiler's two internal lanes plus the four byte inputs
  to named exported `compiled-v1` corners is blocked in this pass. The lane math
  is proved; the presentation-layer corner join is not.
