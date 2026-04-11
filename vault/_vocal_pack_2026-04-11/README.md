# Vocal Pack (2026-04-11)

Goal: a **small, reliable** set of bodies you can record vocals through *tonight* without fighting the instrument.

## Shipping Picks (start here)
Files live in `C:\Users\hooki\Trench\vault\_vocal_pack_2026-04-11\shipping\`.

1) `Small_Talk__Open_Bite.json`
   - Safe live-default. No identity cliffs in stability scan.
   - Suggested: `Q ≈ 0.15–0.35`, sweep MORPH freely.

2) `Small_Talk__Yawn_To_Bite.json`
   - More dramatic / more level. Keep output gain conservative.
   - Suggested: `Q ≈ 0.10–0.30`, MORPH mid→high for bite.

3) `Speaker_Knockerz__Rattle_Howl.json`
   - Aggressive consonant bite + throat motion.
   - Suggested: `Q ≈ 0.60–0.85`, avoid MORPH below ~0.10 (cliff zone).

4) `Speaker_Knockerz__Fracture_Hiss.json`
   - Safer/cleaner variant (lower peak), still hostile.
   - Suggested: `Q ≈ 0.60–0.90`, keep MORPH above ~0.10.

5) `Aluminum_Siding__Glass_Whistle.json`
   - “Air/foil/glass” layer: use as a **sibilance/presence sculptor**, not a full vocal replacement.
   - Suggested: `Q ≈ 0.40–0.70`, MORPH ~0.50–1.00, mix it in (or back off input).

6) `Cul_De_Sac__Comb_Dust.json`
   - Special effect shimmer / comb separation.
   - Suggested: `Q ≈ 0.60–0.90`, avoid MORPH below ~0.20 (cliffy region).

## P2K References (dev / calibration only)
Files live in `C:\Users\hooki\Trench\vault\_vocal_pack_2026-04-11\p2k_refs\`.

- These are **extraction references**, not shipping bodies.
- Expect more cliffs and less “Q obviousness”. Keep levels conservative.

## Fast Hotload (no UI friction)
If your plugin is watching `C:\Users\hooki\Trench\trench_live.json`, copy a body into it:

```powershell
Copy-Item "C:\Users\hooki\Trench\vault\_vocal_pack_2026-04-11\shipping\Small_Talk__Open_Bite.json" `
  "C:\Users\hooki\Trench\trench_live.json" -Force
```

Or run the helper picker:
- `C:\Users\hooki\Trench\tools\hotswap_vocal_pack.ps1`

