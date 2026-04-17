# Compiler Classes Inventory

## Scope

This inventory includes `CPhantom*` classes that are filters, filter compilers, or filter runtimes. The surveyed symbol lists also contain GUI widgets and lookup tables such as `CPhantomBandLimitedTable`, `CPhantomFloatTable`, `CPhantomPanTable`, and `CPhantomGMPanTable`; those were excluded because they are not filter classes. Sources: `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:41-46`, `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:58`, `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:75`, `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:156`, `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:159`, `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:320`

Knowledge-depth labels used below:

- `parameter semantics known`: the surveyed docs expose enough math, table structure, or compile contract to describe how the class behaves
- `decompiled only`: the class/function is identified and some behavior is known, but the surveyed docs do not give a full clean-room parameter contract
- `name only`: only the class name and coarse role are present

Direct shipping-body scaffold set found in the surveyed docs:

- `Speaker Knockerz -> P2k_006 -> 2 Pole Bandpass`
- `Aluminum Siding -> P2k_026 -> Millennium`
- `Small Talk Ah-Ee -> P2k_013 -> Phaser 2`
- `Cul-De-Sac -> P2k_010 -> Swept EQ 2/1 Octave`

Sources: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:6-13`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:5-7`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:427-429`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:849-851`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:1271-1273`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:17`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:21`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:24`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:37`

## Directly implicated by the four shipping bodies

### CPhantomBP2Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:24`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:55`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:57-59`
- Shipping-body need: `Yes - Speaker Knockerz`. The direct scaffold set maps `Speaker Knockerz` to `P2k_006 -> 2 Pole Bandpass`. Sources: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:6-7`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:5-7`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:17`

### CPhantomFilterP2k

- Knowledge depth: `parameter semantics known`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:38`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:141`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:52`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:12`
- Shipping-body need: `Yes - Aluminum Siding`. `Millennium` is filter type 26, and the surveyed docs place `Types 14-45` inside the `CPhantomFilterP2k` wrapper range. Sources: `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:428-429`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:37`, `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:12`

### CPhantomPhaser2

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:31`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:324`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:62`
- Shipping-body need: `Yes - Small Talk Ah-Ee`. The direct scaffold set maps `Small Talk Ah-Ee` to `P2k_013 -> Phaser 2`. Sources: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:10-11`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:849-851`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:24`

### CPhantomSweptEq2

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:28`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:428`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:61`
- Shipping-body need: `Yes - Cul-De-Sac`. The direct scaffold set maps `Cul-De-Sac` to `P2k_010 -> Swept EQ 2/1 Octave`. Sources: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:12-13`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:1271-1273`, `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:21`

## Morph and compiler-family classes

### CPhantomMorph1

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:7`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:249`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:47`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:27-30`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomMorph1Basic

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:8`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:250`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomMorph2

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:9`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:251`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:50`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:33`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomMorphLP

- Knowledge depth: `parameter semantics known`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:10`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:253`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:48`
  - `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:37-42`
  - `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_MorphLP_zero_table.md:36-45`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomMorphLPX

- Knowledge depth: `parameter semantics known`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:11`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:254`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:49`
  - `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:12-18`
  - `C:\Users\hooki\Trench\docs\looperator\phantom_morph_lpx_contracts.md:17-26`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomMorphDesigner

- Knowledge depth: `parameter semantics known`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:12`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:252`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:51`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-findings.md:7-18`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-findings.md:31-47`
- Shipping-body need: `No direct match in the scaffold set above.`

## Vocal / ROM-compiled classes

### CPhantomVocal1

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:15`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:469`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:45`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:17`
- Shipping-body need: `No direct match in the scaffold set above.` Small Talk's direct scaffold stays `Phaser 2` / `Talking Hedz`, not `Vocal1`. Sources: `C:\Users\hooki\Trench\docs\shipping_backbone_map.md:10-11`, `C:\Users\hooki\Trench\datasets\role_vocab\shipping_role_targets_v1.json:849-851`

### CPhantomVocal2

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:16`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:470`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:46`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:17`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomBatman

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:59`
  - `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:45-46`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:6`
- Shipping-body need: `No direct match in the scaffold set above.`

## Standard named filter classes

### CPhantomLP2Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:19`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:184`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:12`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomLP4Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:20`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:185`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:13`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomLP6Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:21`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:186`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:14`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomHP2Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:22`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:168`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:15`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomHP4Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:23`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:169`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:16`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomBP4Pole

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:25`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:56`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:18`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomContraryBandpass

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:26`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:80`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:19`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomSweptEq1

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:27`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:427`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:20`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomSweptEq3

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:29`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:429`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:22`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomPhaser1

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:30`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:323`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:23`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomFlanger1

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:32`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:155`
  - `C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json:26`
  - `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:12`
- Shipping-body need: `No direct match in the scaffold set above.`

## Base and runtime filter classes

### CPhantomFilter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:37`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:133`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomRTFilter

- Knowledge depth: `decompiled only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:35`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:373`
  - `C:\Users\hooki\trench_re_vault\analysis\contradictions.md:74-75`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomRTSSEFilter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:36`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:386`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomE4Filter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\filter_types_classnames.json:39`
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:106`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomLPFilter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:187`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomHPFilter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:170`
- Shipping-body need: `No direct match in the scaffold set above.`

### CPhantomBPFilter

- Knowledge depth: `name only`
- Relevant file paths:
  - `C:\Users\hooki\trench_re_vault\analysis\symbols\classes_cphantom.txt:57`
- Shipping-body need: `No direct match in the scaffold set above.`

## Files surveyed but not cited

- `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md`
