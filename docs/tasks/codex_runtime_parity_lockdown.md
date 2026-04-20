# Codex task — runtime parity lockdown

## Goal

Make `trench-core` produce bit-close output to the Python reference
implementation when the same compiled-v1 cartridge is loaded on both sides.
Today the non-resonant corners of `compiled_talking_hedz` already null at
the f32 WAV noise floor (≈ −177 dB). The resonant corners null at only
≈ −48 dB — that is the gap that stops us from shipping a plugin that
provably reproduces the reference E-mu filter math at audio fidelity.

## Definition of Done

1. `cargo test -p trench-core --test talking_hedz_parity -- --nocapture`
   reports **null ≤ −100 dB on all four corners** for
   `runtime_identity_compiled_v1_both_sides`. Raise the test's threshold
   constant from `-45.0` to `-100.0` once the fix lands — the test must
   fail on any future regression that pushes the resonant corners back up.
2. `compile_path_parity_compiled_v1_vs_calibration` still passes at its
   current −30 dB threshold. Do NOT loosen it. (It will almost certainly
   also improve once the runtime drift is gone, because the same residual
   was dominating its resonant-corner numbers.)
3. `cargo test -p trench-core` (full lib + integration) passes clean, with
   no new `#[ignore]` added. In particular `hedz_cascade`,
   `engine::tests::*`, `resampler::tests::*`, and `registry_validation`
   must all stay green.
4. `python tools/parity_null.py` still runs end to end. Rebakes of other
   manifest entries are not in scope.
5. One clean commit on the working branch with a message that explains
   what drifted, why, and how it is now prevented.

## Root cause (already diagnosed — do not re-diagnose from scratch)

`trench-core/src/cascade.rs:41-68` ramps coefficients inside
`BiquadState::process_sample` every sample:

```rust
self.coeffs[0] += self.deltas[0];
// ... through coeffs[4]
```

`trench-core/src/engine.rs:253-256`'s `snap_to_target` snaps
`output_gain` at block boundaries but does NOT snap the cascade's
coefficients:

```rust
fn snap_to_target(&mut self) {
    self.output_gain = self.target_output_gain;
    self.delta_output_gain = 0.0;
}
```

On a **static** morph/Q corner test, every block recomputes
`deltas = (target - current) / 32`. At the end of the first block
`coeffs` is supposed to equal `target`, but f64 accumulation of 32
`+= delta` operations leaves it off by a few ulp. The next block
recomputes deltas as `(target - slightly-off) / 32`, ramps again, and
the cycle repeats forever. For stages whose poles sit at radius ≈
0.999 (see M100_Q0 of `compiled_talking_hedz`), tiny ulp-scale drift in
`a1`/`a2` rotates the pole enough to de-correlate the output at
≈ −48 dB against a reference that runs fixed f64 coefficients.

The Python reference (`scipy.signal.sosfilt` inside
`tools/parity_null.py`) holds fixed f64 coefficients through the whole
buffer — no ramping, no accumulation noise. That is what the runtime
needs to match for static corners.

## Required invariants (do NOT violate)

- Ramping MUST stay for dynamic morph/Q sweeps. Removing
  `coeffs += delta` or the `set_target` delta computation will break
  any test or code path that expects coefficients to interpolate within
  a 32-sample control block during morph/Q motion.
- Block size stays at `BLOCK_SIZE = 32`.
- The passthrough cartridge, compile-v1 loader, and golden-impulse
  cross-language test (`trench-core/tests/hedz_cascade.rs`) remain
  unchanged.
- No threshold in any test may be loosened except the explicit one in
  DoD item 1 (raise `runtime_identity` threshold from −45 to −100 dB —
  tighter, not looser).

## Suggested approach (Codex may pick a different one as long as DoD holds)

Extend `engine::snap_to_target` — or introduce a sibling
`Cascade::snap_to_target` called from the same point — to:

1. For each of the 12 cascade stages, set `coeffs[i] = target_coeffs[i]`
   (you'll need to cache target coeffs in `BiquadState` at `set_target`
   time, because today only `deltas` is stored).
2. Zero out `deltas`.
3. Do the equivalent snap for `cascade.boost` (already ramped via
   `boost_delta`) — mirrors the existing `output_gain` snap.

The delta-only ramp pattern is a micro-optimisation that saves one f64
store per sample but bakes in the drift. An explicit
`target_coeffs` slot is ≈ 5 × 8 = 40 bytes per stage × 12 stages = 480
bytes of state — free.

## Verification commands (run from repo root on Windows Git Bash)

```
PATH="/c/Users/hooki/.cargo/bin:$(echo "$PATH" | tr ':' '\n' | grep -vE '^/usr/bin$|^/bin$' | paste -sd:)" \
  cargo test -p trench-core --test talking_hedz_parity -- --nocapture

PATH="/c/Users/hooki/.cargo/bin:$(echo "$PATH" | tr ':' '\n' | grep -vE '^/usr/bin$|^/bin$' | paste -sd:)" \
  cargo test -p trench-core
```

The first command must print every row `OK` with `null ≤ -100.00 dB` for
`runtime_identity` and pass the compile_path gate at its existing
threshold. The second must report `45 passed; 0 failed; 0 ignored` on
the lib and a clean bill on every integration test file.

## Boundaries

- Do NOT reshape the parity-null pipeline, rebake canonical WAVs, or
  edit `tools/rust_resampler.py`, `tools/parity_null.py`, or anything
  in `ref/canonical/`. The bug is in the Rust runtime, fix it there.
- Do NOT touch `FILTER_ARTIFACTS.json` except in one case: if the fix
  necessarily changes the `render_chain` description of `rom_hedz` or
  `compiled_*`, update the strings accordingly. No structural schema
  changes.
- Do NOT add new dependencies to `trench-core/Cargo.toml`.
- Do NOT introduce unsafe code.
- Stay on the working branch (see handoff). Do not push.

## Context pointers

- `trench-core/src/cascade.rs` — `BiquadState`, `Cascade`
- `trench-core/src/engine.rs` — `FilterEngine::snap_to_target`,
  `FilterEngine::set_parameters_xyz`, the resampled process_block path
  around line 281
- `trench-core/tests/talking_hedz_parity.rs` — the gate you are
  closing, with the runtime-identity test and the existing
  −45 dB threshold to tighten
- `trench-core/tests/hedz_cascade.rs` — cross-language golden-impulse
  parity; cannot regress
- Commit `0df3ed3` (current HEAD of `cube-prep-cleanup`) contains the
  runtime-identity gate at the current baseline
