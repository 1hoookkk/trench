# Source
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_atlas.md` (13205 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_graph_atlas.md` (2398 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_graph_atlas.png` (329565 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_xml_cube_index.json` (116555 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_xml_cube_report.md` (3210 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_type_usage.md` (4501 bytes)

# What this data is
This data catalogs the usage, classification, and metadata of hardware filter models found across 896 decoded E-mu presets. It identifies which filter indices are actually actively used in presets versus which are dormant, and groups them by their structural properties (family, order, parametric controls). It also describes a visual fingerprint scatter plot of filter response measurements (`range_db` vs `peak_db`). 

# Why it's here
The Looperator token picker must allow users to select from the catalog of available filters. To do this intelligently, it needs to understand the landscape of filter types, their families, and how they relate to one another (e.g., suggesting similar filters based on proximity/adjacency). This spec provides the index of available filters and sets up the structural framework for token suggestions.

# Sanitized content

## Filter Catalog Index
Total active distinct filters catalogued in usage: **49**

Below is the definitive index of active filter types extracted from the presets, grouped by their core DSP family.

| Index | Name | Order | Param3 | Active Layers |
| :--- | :--- | :--- | :--- | :--- |
| **LPF (Low-Pass)** | | | | |
| 0 | Classic 4 LPF | 4 | q | 650 |
| 1 | Smooth 2 LPF | 2 | q | 259 |
| 2 | Steeper 6 LPF | 6 | q | 122 |
| 132 | MegaSweepz 12 LPF | 12 | q | 46 |
| 133 | EarlyRizer 12 LPF | 12 | q | 36 |
| 134 | Millennium 12 LPF | 12 | q | 27 |
| 137 | BassBox-303 12 LPF | 12 | q | 24 |
| 136 | KlubKlassik 12 LPF | 12 | q | 22 |
| **HPF (High-Pass)** | | | | |
| 9 | Deeper 4 HPF | 4 | q | 79 |
| 8 | Shallow 2 HPF | 2 | q | 52 |
| **BPF (Band-Pass)** | | | | |
| 18 | ContraBand 6 BPF | 6 | q | 47 |
| 17 | Band-pass2 4 BPF | 4 | q | 16 |
| 16 | Band-pass1 2 BPF | 2 | q | 14 |
| **PHA (Phaser)** | | | | |
| 66 | BlissBatz 6 PHA | 6 | q | 78 |
| 65 | PhazeShift2 6 PHA | 6 | q | 68 |
| 64 | PhazeShift1 6 PHA | 6 | q | 20 |
| 155 | CruzPusher 12 PHA | 12 | q | 4 |
| **EQ+ (Peak Boost)** | | | | |
| 131 | AceOfBass 12 EQ+ | 12 | q | 43 |
| 32 | Swept1oct 6 EQ+ | 6 | gain_db | 26 |
| 146 | DJAlkaline 12 EQ+ | 12 | q | 21 |
| 34 | Swept3>1oct 6 EQ+ | 6 | gain_db | 12 |
| 33 | Swept2>1oct 6 EQ+ | 6 | gain_db | 11 |
| 148 | RogueHertz 12 EQ+ | 12 | q | 8 |
| 142 | BolandBass 12 EQ+ | 12 | q | 7 |
| 147 | BassTracer 12 EQ+ | 12 | q | 6 |
| 140 | TB-OrNot-TB 12 EQ+ | 12 | q | 5 |
| **EQ- (Notch/Cut)** | | | | |
| 149 | RazorBlades 12 EQ- | 12 | q | 3 |
| **VOW (Vocal/Formant)** | | | | |
| 81 | Ooh-To-Aah 6 VOW | 6 | body_size | 17 |
| 80 | Aah-Ay-Eeh 6 VOW | 6 | body_size | 14 |
| 141 | Ooh-To-Eee 12 VOW | 12 | q | 12 |
| 153 | DeepBouche 12 VOW | 12 | q | 10 |
| 143 | MultiQVox 12 VOW | 12 | q | 8 |
| 152 | UbuOrator 12 VOW | 12 | q | 6 |
| 151 | Eeh-To-Aah 12 VOW | 12 | q | 4 |
| **REZ (Resonant/Peak)** | | | | |
| 160 | LucifersQ 12 REZ | 12 | q | 14 |
| 139 | DeadRinger 12 REZ | 12 | q | 10 |
| 159 | BassOMatic 12 REZ | 12 | q | 7 |
| 135 | MeatyGizmo 12 REZ | 12 | q | 5 |
| 145 | ZoomPeaks 12 REZ | 12 | q | 5 |
| 158 | AcidRavage 12 REZ | 12 | q | 4 |
| 161 | ToothComb 12 REZ | 12 | q | 2 |
| **FLG / DST / WAH / SFX / OFF** | | | | |
| 72 | FlangerLite 6 FLG | 6 | q | 16 |
| 138 | FuzziFace 12 DST | 12 | q | 13 |
| 127 | Off -- --- | None | off | 7 |
| 162 | EarBender 12 WAH | 12 | q | 6 |
| 156 | AngelzHairz 12 FLG | 12 | q | 4 |
| 163 | KlangKling 12 SFX | 12 | q | 4 |
| 157 | DreamWeava 12 FLG | 12 | q | 4 |
| 224 | Unknown Filter 224 | None | unknown | 1 |

## Filter Adjacency
A structural adjacency graph exists based on a "fingerprint scatter" measurement of 33 local P2K filters using `range_db` and `peak_db`. However, explicit coordinates are trapped in a PNG image artifact, so they cannot be represented as structured data here. See UNKNOWNs.

# Integration notes
- Group suggestions and lists in the picker primarily by `Family` (`LPF`, `PHA`, `EQ+`, etc.), as this represents the fundamental shape.
- Use the third parameter (`Param3`) correctly based on the index: mostly `q`, but `gain_db` for some EQ+ types and `body_size` for VOW types.
- The `Active Layers` metric indicates popularity in the original hardware; use this to bias or sort "recommended" filter tokens.

# UNKNOWNs
- **Adjacency Coordinates**: The exact numerical coordinate adjacency matrix for the 33 modeled filters cannot be extracted because it is locked within the `filter_graph_atlas.png` image capture. Without decoding the image or generating fresh measurements from the `range_db` / `peak_db` model, proximity-based suggestions ("filters near this one") are impossible to construct strictly from this dataset.
- **Coverage Gap**: Many high-index scattered filters (e.g., 131-160) are actively used in the hardware presets but have no corresponding modeled P2K response profile.
