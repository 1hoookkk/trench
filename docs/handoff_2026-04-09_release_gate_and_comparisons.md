# TRENCH Handoff (2026-04-09, updated)

## Snapshot
- Release gate tooling is fail-closed at the body-gate level.
- Cul-De-Sac is currently unshippable under the current candidate set.
- Selection scripts were corrected so they no longer nominate gate-failing "least bad" candidates.
- `trench-juce` was built in Release and test suite passed after adding a runtime non-identity cartridge test.

## Current Shipping Gate Status
Source: `vault/_scorecards/shipping_release_gate.json`

- `passed: false`
- Failures:
  - `Cul-De-Sac: holdout_gate_ok=False`
  - `Cul-De-Sac: holdout_promoted=False`

This means current listening with "Cul-De-Sac C3" is known compromised against shipping doctrine.

## Corrections Applied In This Update

### 1) Hard fail-closed candidate selection (no safe fallback winners)
- `tools/promote_holdout.py`
  - Changed design selection to pick only from `gate_ok=True`.
  - If no design candidate passes gate, writes:
    - `design_candidate: null`
    - `selection_failure: "no_gate_passing_design_candidate"`
  - This prevents reporting a gate-failing candidate as a "winner".
- `tools/promote_better_than_emu.py`
  - Changed final selection to pick only from `promoted=True`.
  - If none are promoted, writes:
    - `candidate: null`
    - `selection_failure: "no_promoted_candidate"`

Verification run used:
```powershell
python tools/promote_holdout.py `
  --design-corpus datasets/wav_corpus/design_split_emotion_strat `
  --holdout-corpus datasets/wav_corpus/holdout_split_emotion_strat `
  --out vault/_scorecards/_tmp_selection_check `
  --balanced-sampling `
  --max-clips-per-class 5 `
  --allow-fallback
```

Result (`vault/_scorecards/_tmp_selection_check/better_than_emu_holdout_scorecard.json`):
- `Speaker Knockerz`: selected
- `Aluminum Siding`: selected
- `Small Talk Ah-Ee`: selected
- `Cul-De-Sac`: `design_candidate = null`, `selection_failure = no_gate_passing_design_candidate`

### 2) `trench-juce` runtime/build validation
- Added runtime test:
  - `C:/Users/hooki/trench-juce/plugin/tests/PluginBasics.cpp`
  - New case: `Shipping cartridge JSON loads into non-identity runtime cascade`
  - Asserts shipping JSON -> `Cartridge::loadFromJson` -> `Engine::loadCartridge` yields non-identity active coefficients and non-unity boost.
- Built and tested:
  - `cmake --build build --config Release`
  - `build/Release/Tests.exe --reporter compact`
  - Result: `All tests passed (17 assertions in 3 test cases)`

## Existing Pipeline Components (still valid)
- Emotional-vocal corpus build and stratified split tooling.
- Holdout protocol and release gate entrypoint:
  - `tools/run_shipping_release_gate.py`
  - `cargo xtask release-gate`
- Plot generation:
  - `tools/generate_shipping_plots.py`
  - `tools/compare_shipping_to_emu_stages.py`

## Overnight role-vocab loop (2026-04-10)

Additive generator + ranking loop for the 4 shipping labels (names-only): generates candidates from P2k-proven archetypes, validates stability, renders plots, and writes a concise report under `vault/_overnight/`.

- Script: `tools/run_overnight_role_loop.py`
- Inputs:
  - `pyruntime.forge_generator.BLUEPRINTS` (P2k-derived archetypes)
  - `datasets/role_vocab/shipping_role_targets_v1.json` (role-vocab targets)
  - `docs/sonic_tables/tables.json` (landmarks + vowel table)
- Gates / scoring:
  - Stability (fail-closed): `pyruntime.validator.validate`
  - Evidence surface: `pyruntime.forge_shipping.eval_shipping_body` (40 dB peak ceiling across 11×11 surface + talking/trajectory/continuity)
  - Role match: `pyruntime.role_vocab.body_signature` + `role_distance`
  - Small Talk vowel intent: enforce a vowel-shift gesture (start ∈ {ah, aw, uh, schwa} → end ∈ {ae, eh, ih, ee}, must cross ≥2 vowel bins). Strict literal `ah→ee` proved unreachable under current formant extraction, so this keeps the intent measurable without lying.
  - Candidate generation continues until at least `topk` candidates meet the talking/vowel gate (or attempt budget is exhausted), so winners aren’t forced to include gate-failing bodies.
- Defaults tuned to not reject the talky baseline (`Talking_Hedz` trajectory ≈ 0.398):
  - `--min-talking 0.65`
  - `--min-trajectory 0.35`
  - `--min-continuity 0.75` (musical: avoid jumpy paths)
- Outputs (additive, no deletes):
  - `vault/_overnight/<run_id>/report.md`
  - `vault/_overnight/<run_id>/manifest.json`
  - `vault/_overnight/<run_id>/candidates/<body_key>/*.json`
  - `vault/_overnight/<run_id>/plots/<body_key>/*_{freq,stage}.png`

Example:
```powershell
python tools/run_overnight_role_loop.py --seed 102 --pool 96 --topk 4
```

Latest example run: `vault/_overnight/20260410_203623_seed102/report.md`.

## Commands
### Authoritative gate
```bash
cargo xtask release-gate
```

### Faster gate (reuse corpora/splits)
```bash
cargo xtask release-gate --skip-corpus-build
```

### Selection behavior check (small, quick run)
```powershell
python tools/promote_holdout.py `
  --design-corpus datasets/wav_corpus/design_split_emotion_strat `
  --holdout-corpus datasets/wav_corpus/holdout_split_emotion_strat `
  --out vault/_scorecards/_tmp_selection_check `
  --balanced-sampling `
  --max-clips-per-class 5 `
  --allow-fallback
```

## Session 2026-04-11: P2K Corpus Analysis + Data Integrity

### Heritage cartridge cleanup
- Deleted 5 heritage compiler cartridges that had P2K extraction counterparts (BassBox_303, Early_Rizer, Fuzzi_Face, Meaty_Gizmo, Millennium). P2K extraction is authoritative.
- Quarantined 6 heritage cartridges with no extraction counterpart (types 33+: Ear_Bender, Freak_Shifta, Ooh_to_Eee, Radio_Craze, Razor_Blades, Talking_Hedz) to `cartridges/heritage_quarantine/`. These are reconstructions from `heritage_coeffs.py`, not ground truth.
- Data integrity finding: `cartridges/heritage/Talking_Hedz.json` and `cartridges/p2k/P2k_013.json` encoded the same filter (Phaser 2 = Talking Hedz) through different compilation paths. Coefficients diverged at 1e-5 in pole angles. P2K extraction matches raw skin data; heritage compiler does not.

### P2K corpus profiling
- Profile script: `tools/profile_p2k_skin.py` — extracts pole/zero frequencies, response stats at all 4 corners, motion deltas. All 33 skins profiled to `vault/_profiles/P2k_NNN_profile.json`.
- Stage-band motion matrix: `vault/_profiles/STAGE_MATRIX.md` — maps every stage's frequency band at M0_Q0 vs M100_Q100 across all 33 filters. Identifies gap motions (sub→mouth, bite→throat, etc.) that never appear in the P2K corpus.
- Codex overnight lab: `~/.codex/automations/overnight-trench-loop-tonight/automation.toml` — GPT-5.2 xhigh running every 30 min, reading profile JSONs + sonic tables + trench-mcp, writing plain-English analyses to `vault/_profiles/P2k_NNN_analysis.md`. 27/33 done as of 17:00 AEST. After analyses complete, switches to writing Forge generation briefs targeting matrix gaps.

### trench-mcp code review + fixes
- CRITICAL: `export_cartridge_json` was outputting `a1/r/val` format that the JUCE plugin can't parse. Fixed to output Direct kernel format (`c0-c4`). Hotswap was producing silence; now produces correct coefficients.
- NaN safety: `partial_cmp().unwrap()` replaced with `.unwrap_or(Ordering::Equal)`.
- Morph/Q clamping added to all tool handlers.
- `taste_weights` tool removed from code (was already commented out, stale registration in MCP).
- `find_corridors` limitation documented: horizontal morph slices only, no diagonal paths.

### Random generation — proven dead (again)
- `generate_edge_cases.py` attempted random slot generation targeting matrix gaps. 640 candidates, 2 survivors — both degenerate (ultrasonic pileup, zero talkingness, 13kHz centroid, CRITICAL reactor). Deleted the script.
- This confirms the constellation-first memory (2026-04-07): per-stage randomization cannot produce cascade interaction. Heritage splice (cross-breeding P2K corners) is the only proven generation path.

### Gate philosophy (confirmed by user)
- Hard gates (fail-closed): stability, NaN, peak ceiling, continuity, world-Q collapse prevention.
- NOT hard gates: talkingness, vowel, role-distance. These are per-target scoring only.
- Budgets, not thresholds: stability_budget, identity_budget, allowed_fracture_zones.
- Every gate failure requires a witness: exact (morph, q) cell, failing metric, dominating stage.

### Codec datasets (CVSD, G726) — status
- Raw capture data exists: `trench_re_vault/datasets/cvsd/2026-03-05/` and `trench_re_vault/datasets/g726/2026-03-05/`.
- CVSD: stimulus→capture pairs (dc_step, impulse, silence, noise) + raw probe bitstreams at 16kHz.
- G726: multiple code sizes (cs2/cs3/cs4), same stimulus set + voiced/unvoiced synth.
- `pyruntime/drive.py` has `corrode` (CVSD delta modulator) and `erode` (adaptive quantizer) as audio drive effects.
- These codecs are NOT heritage E-mu components. They are standalone codec algorithms captured for TRENCH's original use — envelope/drive behavior baked into body generation as an original design choice.
- NOT YET: codec timing extraction for envelope/body generation pipeline. The V1 spec lists "codec-derived timing" as a design goal for body layer 4 (charge/discharge/overshoot/hysteresis). User wants this baked into generation — e.g., codec timing governing per-stage or per-peak morph behavior.

### QSound spatial — status
- Capture data exists: `trench_re_vault/datasets/qsound_spatial_v1/` and `v2_live_ingress/`.
- V1 spec lists QSound as body layer 6 (HRTF: 3-band spectral + ILD + ITD).
- No implementation in trench-juce yet.

## Open Issues
- Cul-De-Sac has no gate-passing finalist right now. This is an authoring/scoring target gap, not a valid promotion case.
- Heritage cartridges for types 33-55 (Ear_Bender, Talking_Hedz, etc.) need runtime extraction from Emulator X to replace quarantined reconstructions.
- Hotswap warp function (`warp_coeffs`) assumes Rossum encoding unconditionally — needs format detection or explicit conversion layer.
