# P2K Phoneme Inventory (Unified)

These files map P2K profile clusters (`PH_###`) onto the unified sonic token set used by the shape
bank in `vault/_shapes`.

## Canonical
- `p2k_phoneme_inventory_v2.json` — clusters + canonical mapping (`canonical.tables_ref`,
  `canonical.shape_path`)
- `p2k_phoneme_aliases_v2.json` — one alias per cluster (tables token or heritage fill), suitable
  for UI/preset labeling
- `token_inventory_unified_v2.json` — all static tokens (tables + heritage) with attached P2K
  cluster references
- `unified_sonic_inventory_v2.json` — combined view (tables + shape bank + P2K mapping)

## Heritage gap fill
The 7 previously-unmapped clusters are represented as `heritage.*` tokens and
`vault/_shapes/heritage_phonemes/*.json` bodies:
- `heritage.ph_004`
- `heritage.ph_009`
- `heritage.ph_011`
- `heritage.ph_013`
- `heritage.ph_021`
- `heritage.ph_023`
- `heritage.ph_028`

## Legacy / raw
- `p2k_phoneme_inventory.json` — raw clustering output (no canonical mapping)
- `p2k_phoneme_aliases.json` — pre-fill aliases

