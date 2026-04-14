---
name: bat_phaser_dossier
description: Bat Phaser (BlissBatz/CPhantomBatman) — heritage filter dossier, confirmed architecture
type: project
---

## Identity

- Heritage name: Bat Phaser (P2K UI), BlissBatz (SysEx), CPhantomBatman (firmware class)
- Hardware filter type code: 66 (0x42,00h)
- Family: PHA (Phaser), 6-pole
- Origin: Emulator 4 ROM-fixed algorithm, ported to P2K/Emulator X

## Architecture (confirmed 2026-03-19)

- **NOT a Z-Plane morph surface.** Dedicated `CPhantomBatman` class with its own vftable.
- Fixed 6-pole peak/notch topology. Behavior is near-invariant across Freq/Q.
- Parameters (Freq Hz, Res %) drive specialized internal laws, not 4-corner interpolation.
- 4x larger runtime buffer than CPhantomPhaser1/2 (0x40 vs 0x10 at 44100 Hz).
- Grouped under `standard_filters` in firmware, not `morph_filters`.

## NOT in P2K skin dataset

The 33 P2K skins (P2k_000–P2k_032) are 12-pole Z-Plane morph surfaces only. Bat Phaser is a 6-pole dedicated class — not extractable as a 4-corner skin. To reproduce it, the `CPhantomBatman` compile law must be decompiled from the firmware.

## P2K skin index mapping (confirmed)

Skins are the 12-pole Z-Plane filters, starting at hardware code 131 (0x03,01h):
- P2k_000 = AceOfBass (131)
- P2k_001 = MegaSweepz (132)
- P2k_002 = EarlyRizer (133)
- P2k_003 = Millennium (134)
- P2k_004 = MeatyGizmo (135)
- P2k_005 = KlubKlassik (136)
- P2k_006 = BassBox-303 (137)
- P2k_007 = FuzziFace (138)
- P2k_008 = DeadRinger (139)
- P2k_009 = TB-OrNot-TB (140)
- P2k_010 = Ooh-To-Eee (141)
- P2k_011 = BolandBass (142)
- P2k_012 = MultiQVox (143)
- P2k_013 = TalkingHedz (144)

**Why:** Previous assumption that P2k_000=BatPhaser was wrong. The skins only contain the 12-pole morph surface filters. 6-pole standard/phaser/flanger types are firmware algorithms, not ROM data tables.

## ROM Decode (proven 2026-03-19)

- **Coefficient format: Q14 signed fixed-point** (i16 / 16384), NOT minifloat.
- True ROM base: 0x1806d6CC0 (Ghidra decompiler showed DAT_1806d6cde due to folded +0x1E offset).
- 4 loaders in vtable: slot 1 = 2-stage (0x50 stride), slots 5/9/13 = 3-stage (0x78 stride).
- Filter chain: CPhantomLP2Pole (fixed 2-pole) + Batman ROM stages = 6 or 8 pole.
- Loader 0 corner A verified: 39 dB notch at ~10 kHz. Phaser behavior confirmed.
- Previous decode failure was wrong address + wrong numeric format compounding.

**How to apply:** Use Q14 decode on the ROM tables at the corrected base addresses. Do not use minifloat decode for Batman coefficients. The 4-corner groups per stage support 2D interpolation (morph × Q).
