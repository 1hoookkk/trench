# Session 2026-04-17

## Decisions landed
- `compile_raw.py` committed. Direct pole-zero authoring via raw-stage-v1 is the authoring path for motion bodies.
- Grid-compiler path (`cartridges/engine/` pills + `compile_grid.py`) is viable only for sequence bodies. Invalidated for Speaker Knockerz, Aluminum Siding, Cul-De-Sac. Candidate for Small Talk Ah-Ee only.
- `trench-forge-dev` scaffolded as scoped doctrine exception. Firewalled, monorepo, not linked to shipping plugin.
- Existing 50 pills verified: all are static single-feature phoneme snapshots with zero morph motion. None deliver Aluminum Siding invariants. Motion bodies require authoring from scratch.

## Open bugs
- `compile_raw.py` allpole kind produces `b=(g, a1, r²)` — mirrors denominator. Wrong math. Fix: `b=(g, 0, 0)`. Firmware evidence confirms pole-only shape has flat numerator. Small edit in `encode_stage`. Blast radius zero (no authored raw files exist yet).

## Next session first task
Fix allpole. Run `./check`. Smoke-test with single-stage resonator at 440 Hz. Confirm peak lands at 440 Hz.

## Unresolved (not for tomorrow)
- Long-term authoring surface: rubric-driven search vs per-body hand authoring. Deferred until one body ships.
- `trench_re_vault` and `trenchwork_clean` are not git repos. Data preservation risk. Address next weekend.
- `trench-juce/forge/` legacy residue should move to `_research/forge_archive/`. Do before next authoring session.
