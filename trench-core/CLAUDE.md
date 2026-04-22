# trench-core

Pure DSP crate. No plugin dependencies.

## Layout
- `cartridge.rs` — compiled-v1 JSON loading + bilinear interpolation (Q first, then morph)
- `cascade.rs` — 12-stage serial DF2T biquad cascade with 32-sample control blocks and per-sample coefficient ramping

## Invariants (topology — frozen regardless of cartridge format)
- 6 active stages + 6 passthrough = 12 total
- DF2T math: y = c0*x + w1, w1 = c1*x - c3*y + w2, w2 = c2*x - c4*y
- Kernel-form coefficient interpolation only
- Coefficient ramping: linear over 32-sample block

## Compiled-v1 specifics (current shipping format)
- 4-corner bilinear surface: M0_Q0, M100_Q0, M0_Q100, M100_Q100
- Interpolation order: Q-axis first, then morph-axis
- See `docs/architecture/zplane_truth.md` §"The 4-corner → 2-corner transition" for the planned evolution.
