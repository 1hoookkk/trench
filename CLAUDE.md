# CLAUDE.md — FIELD (Trench)

1) Make the smallest possible diff; no drive‑by refactors.
2) Keep builds green on Windows (VS2022 + CMake/Ninja).
3) Audio thread: no allocations, no locks, no exceptions, no logging.
4) DSP coeff/state math in **double**; audio buffers in **float**.
5) Never interpolate biquad coefficients; interpolate **poles** then rebuild (see `docs/dsp-spec.md`).
6) PARALLEL↔CASCADE changes must crossfade two real cores (no single “interpolated topology”).
7) High‑Q must stay loud: apply bandwidth/loudness compensation in PARALLEL (see `docs/dsp-spec.md`).
8) Saturation belongs in feedback/state path (Dattorro-style), not only on output.
9) Rebuild coefficients at a fixed control rate; avoid per‑sample trig where possible.
10) Every DSP change: add/extend a deterministic impulse/noise/sweep test or give exact repro steps.
11) UI: no knobs/menus/preset browser; one‑hand gestures only.
12) If something is underspecified, follow `docs/*`; otherwise ask for the latest core files.
