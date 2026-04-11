# P2K Stage-Band Motion Matrix

## Band Definitions
| Band | Range | Physical |
|---|---|---|
| sub | 20–150 Hz | floor, rumble |
| chest | 150–400 Hz | chest cavity, jaw hinge |
| throat | 400–1000 Hz | throat, F1 zone |
| mouth | 1000–2500 Hz | mouth, F2, speech |
| presence | 2500–5000 Hz | face, presence peak |
| bite | 5000–8000 Hz | teeth, sibilance |
| air | 8000+ Hz | breath, sheen |

## Matrix (M0_Q0 → M100_Q100 pole band per stage)

| Skin | Name | S1 | S2 | S3 | S4 | S5 | S6 |
|---|---|---|---|---|---|---|---|
| P2k_000 | No Filter | sub→air | pres→thro | air→pres | air | air | air→mout |
| P2k_001 | 2 Pole Lowpass | air→bite | air→mout | air→mout | air→pres | air→bite | mout→ches |
| P2k_002 | 4 Pole Lowpass | air | ches→air | air→mout | air→mout | ches→mout | air→ches |
| P2k_003 | 6 Pole Lowpass | sub→air | air→bite | air→pres | air | air→ches | air→pres |
| P2k_004 | 2 Pole Highpass | sub→bite | air→pres | air→mout | air | air | air→thro |
| P2k_005 | 4 Pole Highpass | pres→air | ches→bite | thro→pres | air | mout→pres | mouth |
| P2k_006 | 2 Pole Bandpass | sub→air | sub→ches | thro→air | air | air | air |
| P2k_007 | 4 Pole Bandpass | presence | thro→pres | thro→pres | ches→pres | ches→pres | thro→pres |
| P2k_008 | Contrary Bandpass | air→pres | mout→pres | sub→thro | air→pres | air→mout | air→thro |
| P2k_009 | Swept EQ 1 Oct | thro→air | mout→air | mout→ches | air | thro→air | ches→air |
| P2k_010 | Swept EQ 2/1 Oct | thro→mout | thro→ches | mout→pres | presence | pres→air | pres→air |
| P2k_011 | Swept EQ 3/1 Oct | sub→air | air→ches | thro→air | air | air | sub→air |
| P2k_012 | Phaser 1 | bite→air | mout→pres | presence | pres→bite | air→mout | ches→thro |
| P2k_013 | Phaser 2 | air | thro→ches | mouth | mouth | presence | ches→mout |
| P2k_014 | Bat Phaser | sub→air | ches→mout | mout→air | mout→air | mout→air | air |
| P2k_015 | Flanger Lite | air | mout→thro | pres→bite | air | bite→air | air |
| P2k_016 | Vocal Ah-Ay-Ee | air | thro→air | thro→air | air | air | mout→ches |
| P2k_017 | Vocal Oo-Ah | air | sub→thro | thro→air | air | pres→bite | air |
| P2k_018 | Dual EQ Morph | air→bite | thro→sub | mout→air | air→mout | bite→air | air→mout |
| P2k_019 | Dual EQ + LP | presence | ches→mout | mout→bite | air→bite | bite→pres | presence |
| P2k_020 | Dual EQ Morph/Expr | air→pres | ches→thro | mout→thro | air→pres | pres→mout | presence |
| P2k_021 | Peak/Shelf Morph | pres→air | throat | thro→ches | pres→mout | mout→pres | pres→mout |
| P2k_022 | Morph Designer | presence | presence | presence | mout→pres | mouth | chest |
| P2k_023 | Ace of Bass | air | throat | mouth | air | air→pres | air→bite |
| P2k_024 | MegaSweepz | presence | thro→mout | thro→mout | ches→pres | ches→sub | thro→pres |
| P2k_025 | Early Rizer | air | throat | presence | air | bite | air |
| P2k_026 | Millennium | air | air | air | air→pres | bite→mout | pres→ches |
| P2k_027 | Meaty Gizmo | pres→bite | ches→air | thro→air | air→pres | mout→pres | mout→thro |
| P2k_028 | Klub Klassik | air | sub→air | thro→air | sub→air | mouth | mout→air |
| P2k_029 | BassBox 303 | sub→air | air→mout | thro→air | air | air | sub→air |
| P2k_030 | Fuzzi Face | air→bite | thro→air | mout→air | mout→pres | presence | ches→thro |
| P2k_031 | Dead Ringer | ches→bite | pres→mout | air→mout | air→pres | air→ches | air→bite |
| P2k_032 | TB or Not TB | presence | air | air | air | bite→air | pres→bite |

## What E-mu Built Most
| Motion | Count |
|---|---|
| air→presence | 13 |
| air→mouth | 12 |
| throat→air | 11 |
| sub→air | 10 |
| air→bite | 8 |
| mouth→presence | 8 |
| mouth→air | 7 |

## Gaps (motions that never appear)
- sub→mouth (bass opening directly into speech)
- sub→presence (sub teleporting to face)
- chest→presence (chest jumping to presence)
- bite→throat (aggression collapsing into warmth)
- bite→chest (teeth folding into body)
- bite→sub (upper sizzle diving to floor)
- No filter starts 3+ stages in sub
- No filter has all 6 stages moving to different bands
