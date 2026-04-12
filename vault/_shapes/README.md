# Shape Bank (Generated)

This directory is an authored bank of **sonic targets** (static points) and **sonic trajectories**
(explicit morph journeys), written as `compiled-v1` cartridge JSONs.

## Canonical outputs
- `manifest.json` — base author run (tables-only)
- `manifest_unified.json` — base + heritage phoneme fills (P2K cluster gaps)
- `shape_bank_catalog.json` — base catalog (43 static + 320 trajectories)
- `shape_bank_catalog_unified.json` — unified catalog (50 static + 320 trajectories)
- `palette.jsonl` — base token palette (43)
- `palette_unified.jsonl` — unified token palette (50)

## Regenerate
- Base bank (363 bodies):
  - `python tools/author_sonic_bank.py --pairs --out-root vault/_shapes`
  - Static targets:
    - `vowels/`, `nasals/`, `landmarks/`, `consonants/`, `instruments/`, `moog/`, `singer/`,
      `electronic/`, `measured_bells/`
  - Trajectories:
    - `vowel_pairs/`, `vowel_nasal/`, `vowel_cons/`, `vowel_inst/`, `vowel_land/`

## Notes
- No scoring/ranking/filtering is performed here — every authored body is treated as valid because the
  target is a real sonic location.
- Catalog/palette entries include both absolute and repo-relative paths for portability.

