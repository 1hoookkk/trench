# TRENCH Session Handoff: 2026-04-11 (AEST)

## Critical Correction (Must Read)
- The **response stats** in `vault/_profiles/P2k_NNN_profile.json` were **wrong** due to a bug in `tools/profile_p2k_skin.py`:
  - It passed `body.boost` as the **sample rate** into `cascade_response_db(...)`.
  - It also built the frequency grid with the default plugin SR (44100) instead of authoring SR (39062.5).
- Fixed `tools/profile_p2k_skin.py` and **regenerated all 33** `P2k_NNN_profile.json` files on **2026-04-11**.
- Anything written from the old profiles (e.g. existing `P2k_NNN_analysis.md`) should be treated as **stale** until re-derived.

- Also fixed a widespread authoring bug: several places computed `freq_points()` at 44.1k but evaluated the transfer at 39.0625k (`cascade_response_db(..., sr=SR)`), which makes response plots + gates numerically wrong. Patched call sites so the **frequency grid SR matches the evaluation SR**.

## Current Truths
- All 33 P2K profiles exist in `vault/_profiles/` as `P2k_NNN_profile.json` (now regenerated with correct response stats).
- P2k_013 is labeled **“Phaser 2”** (filterType=13), but its structure reads as **vocal/formant + cavity** behavior (not a whooshy swirl).
- P2k_013 “mouth stack” at M0_Q0 (poles): ~199 Hz, ~891 Hz, ~1570 Hz, ~2348 Hz, ~4607 Hz, ~9321 Hz.
- P2k_013 motion law (from corrected profile deltas):
  - **Morph** mostly changes *pressure/stress* (peak +17–21 dB; range +20–22 dB) with **minimal centroid movement** (~0 to -63 Hz).
  - **Q** acts like a *clamp / throat shift*: range +13–16 dB and centroid +167–226 Hz (but it does **not** dominate the entire geometry the way a second morph axis would).
- Sonic tables exist at `docs/sonic_tables/tables.json` (vowel formants + singer’s formant + consonant/fricative landmarks) and should be used as the *interpretation lens*.
- Current Domain Focus: 2 (P2K anatomy → role vocabulary for “talking”)

## Identified Blockers
- Names are misleading (e.g., “Phaser 2” behaving like a throat). If we trust names, we’ll mis-author. Cost: we chase the wrong motion law.
- Any “atlas / vector DB / phoneme labeling” built on corrupted response stats will be garbage. Cost: you’ll *learn the wrong instrument*.
- “Talking” perception comes from **relative motion** (one ridge rising while another narrows; cavities swallowing energy). Copying static peak positions will produce dead, EQ-like bodies.

## Next Steps
- Re-derive any prose analyses from the **regenerated** `P2k_NNN_profile.json` files.
- Extract a compact “role signature” per body/corner: anchor band, carrier bands, void/notch territory, stress peaks, plus morph-vs-Q deltas.
- Build a corpus-grounded “phoneme atlas”: label (M0Q0, M0Q100, M100Q0, M100Q100) as `phoneme + stress + mouth-band`, then compute the smallest reusable set.

## Tonight: Vocal Pack (Ready to Record)
- Pack folder: `C:\Users\hooki\Trench\vault\_vocal_pack_2026-04-11`
  - Shipping picks (use these live first): `...\shipping\`
  - P2K references (dev/calibration, can be unruly/loud): `...\p2k_refs\`
- Hotload helper: `C:\Users\hooki\Trench\tools\hotswap_vocal_pack.ps1`

## Known Stale Artifacts (Don’t Trust Without Re-run)
- `C:\Users\hooki\Trench\vault\_shipping_finalists\manifest.json` and anything derived from it (dated 2026-04-09) predates the SR/response-stats fixes above.

## Invariants Status
- No Rust DSP/runtime changes made (topology invariants unchanged).
- Tools changed + scripts run:
  - Fixed `tools/profile_p2k_skin.py` and regenerated `vault/_profiles/P2k_NNN_profile.json`.
  - Fixed SR-consistency in authoring/gates: `pyruntime/analysis.py`, `pyruntime/api.py`, `pyruntime/forge_optimize.py`, `pyruntime/forge_shipping.py`, `pyruntime/forge_heritage_optimize.py`, `pyruntime/forge_plot.py`, `pyruntime/forge_trajectory.py`, `pyruntime/forge_roll.py`, `pyruntime/forge_session.py`.
  - No Rust builds/tests run.
