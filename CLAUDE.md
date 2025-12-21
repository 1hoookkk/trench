# CLAUDE.md — FIELD

## Build

```bash
cmake --preset windows-release
cmake --build build --config Release
```

## Rules

1. Smallest possible diff. No drive-by refactors.
2. Keep builds green (VS2022 + CMake).
3. Audio thread: no allocations, locks, exceptions, or logging.
4. Coefficients/state in `double`; audio buffers in `float`.
5. Never interpolate biquad coefficients — interpolate poles, then rebuild.
6. Topology changes (PARALLEL↔CASCADE): crossfade two real cores.
7. Inactive parallel stages output 0 (not bypass).
8. Saturation in feedback path, not only output.
9. Rebuild coefficients at control rate (16 samples), not per-sample.

## Key Files

- `src/dsp/ZPlaneCore.h` — 14-pole engine
- `src/dsp/ZPlaneDesigner.h` — frame factory
- `docs/dsp-spec.md` — complete DSP math reference

## Parameters

- `FIELD` (0-1): filter universe position
- `TENSION` (0-1): drive + resonance intensity
- `MOTION` (0-1): morph depth
