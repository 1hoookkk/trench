# canonical_parity.rs — fixture pipeline regression

**Status:** open. Pre-existing on `cube-prep-cleanup`; not a Stage-1
regression. Logged as part of the Stage-1 cube strip commit so it
survives separately.

## Symptoms

Running:

```
cargo test -p trench-core --test canonical_parity
```

produces, among other lines:

```
! cal_Talking_Hedz          cal      M0_Q0  gain=0.7689  rel_null=    -3.9 dB
! cal_Talking_Hedz          cal    M0_Q100  gain=0.7619  rel_null=    -3.8 dB
  compiled_p2k_013: unknown source_type: compiled_v1
  compiled_talking_hedz: unknown source_type: compiled_v1

rust parity: 164 corners / 43 entries (35 raw + 6 cal), worst pass -147.4 dB, threshold -120 dB
  FAIL: P2k_013 M0_Q0 -3.9 dB
  ...
  FAIL: compiled_p2k_013 load +0.0 dB
  FAIL: compiled_talking_hedz load +0.0 dB

thread 'rust_cascade_nulls_against_python_canonical_refs' panicked:
10 corner(s) exceeded -120 dB threshold
```

## Why this is not a Stage-1 bug

1. The Stage-1 `Cartridge::interpolate(morph, q)` is byte-identical to
   the prior MorphQ path — same Q-first-then-morph lerp, variables
   renamed only (`y`→`q`, `x`→`morph`).
2. The hardware-parity gate `talking_hedz_parity` (1/1) passes after
   Stage 1. That test renders against a hardware reference and would
   regress immediately if bilinear interp shifted.
3. `engine_pills` (6/6), `cartridge_schema_tests` (7/7), `hedz_cascade`
   (2/2), and all 65 `cargo test --lib` unit tests pass after Stage 1.
4. `git diff HEAD -- trench-core/tests/canonical_parity.rs` on this
   branch shows only whitespace reformatting — no test-logic change
   relative to HEAD, so the failure is not caused by a test rewrite.

## Probable root cause

The failure log contains `unknown source_type: compiled_v1` for two
entries (`compiled_p2k_013`, `compiled_talking_hedz`). That string
originates in the Python canonical-reference generator — the fixture
pipeline that produces the golden WAVs the Rust test nulls against.
The generator doesn't recognise `compiled_v1` as a valid
`source_type`, so those entries emit placeholder/zero references and
the Rust side nulls against noise.

The -3.9 / -3.8 dB null values on the `cal_*` and raw `*` entries
suggest that a recent change to the fixture generator, WAV table, or
`AGC_TABLE` values desynced the reference set from the current Rust
cascade output. The working tree carries modifications to several
`tools/` scripts and to fixture files — any of those could be the
trigger.

## Files worth starting from

- `trench-core/tests/canonical_parity.rs` — test-side; drives Python,
  reads WAVs, computes nulls
- `tools/parity_null.py` — reference pipeline (same order as engine)
- `tools/hardware_parity_check.py`, `tools/hardware_parity_batch.py` —
  hardware-side oracle
- `ref/canonical/` — golden WAV table consumed by the Rust test
- `trench-core/tests/fixtures/` — fixture WAVs

## Remediation path

1. Reproduce: confirm the failure on a clean worktree of the current
   branch head.
2. Locate the Python source that emits `unknown source_type:
   compiled_v1` and add `compiled_v1` as a recognised source type —
   or change the test manifest to not request `compiled_v1` sources
   it cannot generate.
3. Re-bake `ref/canonical/` and re-run. If the `cal_*` nulls still
   miss -120 dB, bisect `tools/parity_null.py` changes on this branch
   against the commit that last green-gated the test
   (`65a88b4 parity: extend rust canonical_parity to 41/41`).
4. When the test passes cleanly, remove this doc.

## Not blocking

Stage 1 ships with `canonical_parity` knowingly red. All other parity
gates remain enforced. Stage 2 work on the Q-law registry is
unaffected — Q laws only change how `compiled-v1` content is *baked*,
not how the cascade runs it, and `talking_hedz_parity` continues to
guarantee the cascade is correct.
