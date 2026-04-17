# Bodies Evidence

## Scope and method

This pass uses the body contract in `C:\Users\hooki\Trench\BODIES.md`, the deterministic scaffold map in `C:\Users\hooki\Trench\docs\shipping_backbone_map.md`, and the generated shipping target file in `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json`. The workflow doc says to "`Build shipping name targets from mapped P2k bodies`" and writes `shipping_role_targets_v1.json`, so those mappings are treated as the strongest direct evidence in this report. Source: `C:\Users\hooki\Trench\docs\p2k_role_vocabulary_workflow.md:14-22`

Evidence rank used below:

- Rank A: explicit body-to-source mapping or exact runtime provenance
- Rank B: direct adjacent calibration or authoring note
- Rank C: indirect analog, gap note, or caution

## Speaker Knockerz

Body contract:

> `**Sonic identity:** highly pressurized sub-harmonic resonator that chokes`
> `and fractures low-end weight into mid-range, tearing speaker-cone cry.`
> `**Invariant:** the fundamental below 60 Hz remains anchored and phase-locked.`
> `The sub never disappears, even as upper harmonics are brutalized.`

Source: `C:\Users\hooki\Trench\BODIES.md:12-16`

### Rank A evidence

- P2K skin: `P2k_006`
  - Quote:
    > `- Speaker Knockerz -> P2k_006`
    > `Reason: closest safe BassBox lineage in current runtime.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:6-7`
  - Quote:
    > `"Speaker Knockerz": {`
    > `"source_p2k": "P2k_006",`
    > `"source_name": "2 Pole Bandpass",`
  - Source: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:5-7`

- Filter class: `CPhantomBP2Pole` (inference from the cited type-name mapping)
  - Quote:
    > `"source_p2k": "P2k_006",`
    > `"source_name": "2 Pole Bandpass",`
  - Source: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:6-7`
  - Quote:
    > `"6": "2 Pole Bandpass",`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:17`
  - Quote:
    > `- ft=6-7: BP (2/4 pole)`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:59`
  - Assessment: the documented path is `P2k_006 -> 2 Pole Bandpass -> BP family`; the strongest class-level match in the surveyed files is `CPhantomBP2Pole`. Sources: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:6-7`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:17`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:59`

- Calibration file: `C:\Users\hooki\trenchwork_clean\docs\calibration\BassBox_303.json`
  - Quote:
    > `Reason: closest safe BassBox lineage in current runtime.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:7`
  - Quote:
    > `"name": "BassBox 303",`
    > `"strategy": "FREQUENCY_SWAP",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\BassBox_303.json:2-7`
  - Quote:
    > `"role": "Anchor",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\BassBox_303.json:33`

### Rank B evidence

- RE / authoring note: bass-pressure route and Q contract
  - Quote:
    > `World: bass_pressure`
    > `Body: speaker_knockerz`
    > `Scope: 3-room vertical slice`
  - Source: `C:\Users\hooki\Trench\docs\phase1_bass_pressure_vertical_slice.md:3-6`
  - Quote:
    > `- Speaker Knockerz: squeeze bass air mass without killing the floor`
  - Source: `C:\Users\hooki\Trench\docs\body_authoring_seed.md:51`
  - Quote:
    > `- Speaker Knockerz: Q should tighten choke, rip, rattle, and cry, but preserve bass body and at least two occupied regions`
  - Source: `C:\Users\hooki\Trench\docs\world_q_behavior_gate_spec.md:241`
  - Quote:
    > `At q_max, Speaker Knockerz must satisfy all of these:`
    > `- sub anchor still audible and measured as present`
    > `- at least 2 prominent occupied regions remain`
  - Source: `C:\Users\hooki\Trench\docs\world_q_behavior_gate_spec.md:272-275`
  - Quote:
    > `- Stage 0: Sub anchor at 40-50 Hz (r=0.998, PURE contour - no zero, just pole)`
    > `- Stage 1: Low-mid resonance 150-400 Hz (r=0.990, INTERIOR_ZERO)`
    > `- Stage 2: Mid-band stress 400-1200 Hz (r=0.985, UNIT_CIRCLE zero for notch)`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:227-231`

### Rank C evidence

- Adjacent calibration gaps: bass-specific anchors still missing elsewhere in the corpus
  - Quote:
    > `"name": "Boland Bass",`
    > `"priority": "HIGH",`
    > `"reason": "No bass-synthesis calibration anchor"`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:7-9`
  - Quote:
    > `"name": "Bassomatic",`
    > `"priority": "HIGH",`
    > `"reason": "Distinct bass strategy without calibration"`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:25-27`
  - Assessment: these are relevant bass-authoring gaps, but they are weaker evidence than the direct `P2k_006` / BassBox lineage mapping. Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:7-9`, `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:25-27`

## Aluminum Siding

Body contract:

> `**Sonic identity:** brittle high-tension treble stressor that folds air and`
> `sibilance into expensive crystalline damage.`
> `**Invariant:** the midrange around 1 kHz remains scooped - a permanent`
> `acoustic void. All kinetic energy lives in the extreme highs.`

Source: `C:\Users\hooki\Trench\BODIES.md:48-52`

### Rank A evidence

- P2K skin: `P2k_026`
  - Quote:
    > `- Aluminum Siding -> P2k_026`
    > `Reason: Millennium lineage with dense HF clustering.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:8-9`
  - Quote:
    > `"Aluminum Siding": {`
    > `"source_p2k": "P2k_026",`
    > `"source_name": "Millennium",`
  - Source: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:427-429`

- Filter class: `CPhantomFilterP2k` wrapper (inference from the cited type-range mapping; underlying sub-class is not documented in the surveyed files)
  - Quote:
    > `"26": "Millennium",`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:37`
  - Quote:
    > `- Types 14-45: 32 CPhantomFilterP2k instances - P2K wrapper dispatching to Morph1/2/LP/LPX/Designer/Flanger/Vocal sub-classes`
  - Source: `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:12`
  - Assessment: `Millennium` is filter type 26, which falls inside the documented `Types 14-45` P2K-wrapper range, so `CPhantomFilterP2k` is the strongest class-level match the surveyed docs support. Sources: `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:37`, `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:12`

- Calibration file: `C:\Users\hooki\trenchwork_clean\docs\calibration\Millennium.json`
  - Quote:
    > `Reason: Millennium lineage with dense HF clustering.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:9`
  - Quote:
    > `"name": "Millennium",`
    > `"strategy": "HF_WALL_FOLD",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Millennium.json:2-7`
  - Quote:
    > `"role": "HF/Nyquist",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Millennium.json:33`

### Rank B evidence

- RE / authoring note: brittle-air contract and explicit 1 kHz void authoring
  - Quote:
    > `| 2 | Aluminum Siding | Brittle top-end damage. Sibilant fold, foil tear. |`
    > `- Aluminum Siding: weaponize brittle treble without filling the mids`
  - Source: `C:\Users\hooki\Trench\docs\body_authoring_seed.md:28`, `C:\Users\hooki\Trench\docs\body_authoring_seed.md:52`
  - Quote:
    > `- Permanent 1 kHz void (midrange scooped at all morph/Q positions)`
    > `- All energy focused in extreme highs (7-18 kHz)`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:255-256`
  - Quote:
    > `- Stage 0: High shelf/peak at 10-12 kHz (r=0.975, UNIT_CIRCLE zero at ~1 kHz for the void)`
    > `- Stage 2: Anti-resonance notch locked at 1 kHz (r=0.995, zero_forced at 1 kHz)`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:260-263`
  - Quote:
    > `Intent:`
    > `- Q should sharpen brittle shear, glass edges, foil scream, and mid-scoop tension`
  - Source: `C:\Users\hooki\Trench\docs\world_q_behavior_gate_spec.md:148-150`

### Rank C evidence

- Adjacent calibration: `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json`
  - Quote:
    > `"name": "Razor Blades",`
    > `"strategy": "NYQUIST_SHARPENING",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:2-7`
  - Quote:
    > `"role": "HF/Nyquist",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:33`
  - Assessment: this is a plausible HF-damage analog, but no surveyed file ties `Razor Blades` to `Aluminum Siding` as directly as the `Millennium` lineage note does. Sources: `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:2-7`, `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:33`

- Adjacent missing calibration: `Acid Ravage`
  - Quote:
    > `"name": "Acid Ravage",`
    > `"reason": "Aggressive acid, adjacent to Razor Blades"`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:19-21`
  - Assessment: relevant as a nearby brittle/HF family gap, but still indirect. Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\index.json:19-21`

## Small Talk Ah-Ee

Body contract:

> `**Sonic identity:** biomechanical vocal cavity that mutates from a deeply`
> `relaxed open throat into a strangled digital scream.`
> `**Invariant:** exactly two dominant formants are active at all times,`
> `mimicking human vocal cord separation, even when surrounding frequencies`
> `are mangled.`

Source: `C:\Users\hooki\Trench\BODIES.md:84-89`

### Rank A evidence

- P2K skin: `P2k_013`
  - Quote:
    > `- Small Talk Ah-Ee -> P2k_013`
    > `Reason: Talking Hedz lineage with existing vocal formant behavior.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:10-11`
  - Quote:
    > `"Small Talk Ah-Ee": {`
    > `"source_p2k": "P2k_013",`
    > `"source_name": "Phaser 2",`
  - Source: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:849-851`

- Filter class: `CPhantomPhaser2` (inference from the cited type-name mapping)
  - Quote:
    > `"13": "Phaser 2",`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:24`
  - Quote:
    > `- ft=12-13: Phaser 1/2`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:62`
  - Assessment: the direct scaffold is `P2k_013`, and the surveyed class mappings identify filter type 13 as `Phaser 2`, so `CPhantomPhaser2` is the strongest class-level match. Sources: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:850-851`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:24`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:62`

- Runtime provenance note: `Talking_Hedz` and `P2k_013` are the same filter through different paths
  - Quote:
    > `Data integrity finding: cartridges/heritage/Talking_Hedz.json and cartridges/p2k/P2k_013.json encoded the same filter (Phaser 2 = Talking Hedz) through different compilation paths.`
    > `P2K extraction matches raw skin data; heritage compiler does not.`
  - Source: `C:\Users\hooki\Trench\docs\handoff_2026-04-09_release_gate_and_comparisons.md:129`

- Calibration file: `C:\Users\hooki\trenchwork_clean\docs\calibration\Talking_Hedz.json`
  - Quote:
    > `"name": "Talking Hedz",`
    > `"strategy": "RADIUS_ONLY_MORPH",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Talking_Hedz.json:2-7`
  - Quote:
    > `"role": "HF/Nyquist",`
    > `"family": "INTERIOR_ZERO",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Talking_Hedz.json:33`, `C:\Users\hooki\trenchwork_clean\docs\calibration\Talking_Hedz.json:43`

### Rank B evidence

- Adjacent calibration file: `C:\Users\hooki\trenchwork_clean\docs\calibration\Ooh_to_Eee_(approx).json`
  - Quote:
    > `"name": "Ooh to Eee (approx)",`
    > `"strategy": "VOWEL_FORMANT",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Ooh_to_Eee_(approx).json:2-7`
  - Quote:
    > `"role": "FormantBite",`
    > `"family": "INTERIOR_ZERO",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Ooh_to_Eee_(approx).json:33`, `C:\Users\hooki\trenchwork_clean\docs\calibration\Ooh_to_Eee_(approx).json:43`
  - Quote:
    > `- cartridges/heritage/Ooh_to_Eee_(approx).json - reference vocal body (composite=0.96)`
    > `- cartridges/p2k/P2k_013.json - Talking Hedz P2K (composite=0.64)`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:18-19`

- RE / authoring note: vocal target class
  - Quote:
    > `"id": "reveal_late_vocal",`
    > `"intent": "A restrained body at rest that opens into articulate, speaking movement under slow animation.",`
    > `"internal_reference_bodies": [`
    > `"Talking Hedz",`
    > `"Ubu Orator"`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\behavior_targets\2026-03-10\launch_body_targets_v1.json:10-15`
  - Quote:
    > `- 3 active stages per corner with interior zeros`
    > `- Stage 0: F1 formant anchor`
    > `- Stage 1: F2 formant bite`
    > `- Stage 2: F3/F4 air`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:119-124`
  - Quote:
    > `- Small Talk: force an anatomical vocal read out of non-vocal material`
  - Source: `C:\Users\hooki\Trench\docs\body_authoring_seed.md:53`

### Rank C evidence

- Caution note: do not treat Talking Hedz as a computed MorphDesigner source
  - Quote:
    > `Talking Hedz is NOT a computed/compiled filter. It is a static ROM lookup.`
    > `hedz.xml is a manual Morph Designer reconstruction, NOT the runtime source.`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:36-37`
  - Assessment: this is relevant because it narrows the usable provenance path for Small Talk to `P2k_013` / extracted `Talking Hedz`, not a generic MorphDesigner reconstruction. Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:36-37`

## Cul-De-Sac

Body contract:

> `**Sonic identity:** paradoxical structure that begins as a thick blunt`
> `physical tube and physically fractures halfway through the morph into a`
> `scattered multi-band comb matrix.`
> `**Invariant:** a constant low-level resonant hum at the root note.`
> `The only tether to reality while upper harmonics undergo structural collapse.`

Source: `C:\Users\hooki\Trench\BODIES.md:121-126`

### Rank A evidence

- P2K skin: `P2k_010`
  - Quote:
    > `- Cul-De-Sac -> P2k_010`
    > `Reason: fragmented structure closest to structural collapse brief.`
  - Source: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:12-13`
  - Quote:
    > `"Cul-De-Sac": {`
    > `"source_p2k": "P2k_010",`
    > `"source_name": "Swept EQ 2/1 Octave",`
  - Source: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:1271-1273`

- Filter class: `CPhantomSweptEq2` (inference from the cited type-name mapping)
  - Quote:
    > `"10": "Swept EQ 2/1 Octave",`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:21`
  - Quote:
    > `- ft=9-11: Swept EQ (1/2:1/3:1 octave)`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:61`
  - Assessment: the direct scaffold is `P2k_010`, and the surveyed type map identifies filter type 10 as `Swept EQ 2/1 Octave`, so `CPhantomSweptEq2` is the strongest class-level match. Sources: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:1272-1273`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:21`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:61`

### Rank B evidence

- RE / authoring note: explicit fracture brief
  - Quote:
    > `| 4 | Cul-De-Sac | Structural fracture body. Tube to comb collapse. |`
    > `- Cul-De-Sac: cross a structural boundary and return as a different filter`
  - Source: `C:\Users\hooki\Trench\docs\body_authoring_seed.md:30`, `C:\Users\hooki\Trench\docs\body_authoring_seed.md:54`
  - Quote:
    > `- M0 = thick tube resonance (single dominant peak)`
    > `- M50 = State Boundary / phase cancellation null`
    > `- M100 = scattered comb filter (multiple narrow peaks)`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:282-285`
  - Quote:
    > `- The bilinear interpolation between few-stage and many-stage creates the "fracture" at midpoint`
    > `Key insight: stages that are passthrough at M0 and active at M100 naturally create the "fracture" transition`
  - Source: `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:291`, `C:\Users\hooki\Trench\docs\superpowers\plans\2026-04-06-shipping-bodies-v1.md:296`
  - Quote:
    > `Intent:`
    > `- Q may intensify discontinuity, voids, combs, and collapse`
    > `- even here, collapse must read as fracture, not accidental single-stage takeover`
  - Source: `C:\Users\hooki\Trench\docs\world_q_behavior_gate_spec.md:205-208`

- RE note: comb-null behavior is treated as a defining character source
  - Quote:
    > `Generated filters were spectrally flat because zero radii were sampled from P2K vocabulary distributions (mostly r < 0.9).`
    > `This produced shallow spectral dips instead of the 30dB comb nulls that define character in heritage presets.`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\context\zero-forcing.md:16-17`
  - Quote:
    > `Bodies authored for dramatic mid-morph behavior (like Cul-De-Sac's "State Boundary" zone) would momentarily flash through dangerous territory on every hit.`
  - Source: `C:\Users\hooki\Trench\docs\modulation_exploration.md:137`

### Rank C evidence

- Direct calibration file for `P2k_010` / `Swept EQ 2/1 Octave`
  - No evidence found in the surveyed directories.

- Adjacent fracture reference set and one matching calibration file
  - Quote:
    > `"id": "fracture_contraband",`
    > `"shipping_goal": "Original contraband body with authored instability and explicit recovery behavior.",`
    > `"internal_reference_bodies": [`
    > `"Lucifers Q",`
    > `"Meaty Gizmo",`
    > `"Razor Blades"`
  - Source: `C:\Users\hooki\trenchwork_clean\datasets\behavior_targets\2026-03-10\launch_body_targets_v1.json:178-185`
  - Quote:
    > `"name": "Razor Blades",`
    > `"strategy": "NYQUIST_SHARPENING",`
  - Source: `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:2-7`
  - Assessment: `Razor Blades` is the only fracture-adjacent reference in that target-class set that also has a surveyed calibration file; it is still weaker evidence than the direct `P2k_010` scaffold mapping. Sources: `C:\Users\hooki\trenchwork_clean\datasets\behavior_targets\2026-03-10\launch_body_targets_v1.json:178-185`, `C:\Users\hooki\trenchwork_clean\docs\calibration\Razor_Blades.json:2-7`

## Files surveyed but not cited

- `C:\Users\hooki\Trench\CLAUDE.md`
- `C:\Users\hooki\Trench\DOCTRINE.md`
- `C:\Users\hooki\trenchwork_clean\pyruntime\CLAUDE.md`
- `C:\Users\hooki\Trench\.agents\skills\trench-authoring-session\SKILL.md`
- `C:\Users\hooki\Trench\.agents\skills\trench-body-session\SKILL.md`
- `C:\Users\hooki\Trench\docs\hostile_authoring_workflow_spec.md`
- `C:\Users\hooki\Trench\docs\superpowers\specs\compiler_grid_to_cartridge.md`
- `C:\Users\hooki\Trench\docs\trench_master_prompt_library.md`
- `C:\Users\hooki\Trench\docs\world_body_candidate_schema.md`
- `C:\Users\hooki\trenchwork_clean\datasets\filter_inventory\2026-03-10\filter_asset_inventory_v1.json`
