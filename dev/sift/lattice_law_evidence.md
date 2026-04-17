# Lattice Law Evidence

## Claim under review

Claim: `P2K Lattice Compiler enforces Stage N zeros = Stage N-1 poles.`

Scope note: the surveyed files use the phrase `P2K Lattice Compiler` specifically for `FUN_1802c59b0`, the shared inner compiler for `Morph1` / `MorphLP` / `MorphLPX`. They do not use that phrase for every P2K-facing path, and several surveyed files explicitly separate ROM-only classes, `Morph2`, and `CPhantomFilterP2k` memcpy paths from that lattice family. Sources: `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:26-33`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:43-52`

## Supporting evidence

### 1. Exact statement in the notebook-ground-truth pack

- Quote:
  > `### P2K Lattice Compiler (FUN_1802c59b0)`
  > `- Shared inner compiler for Morph1/MorphLP/MorphLPX classes`
  > `- 3-stage cross-coupled lattice: stage N zeros = stage N-1 poles`
  > `- Numerator table (param_5) only affects stage 0; stages 1-2 always cross-coupled`
  > `- Morph2: separate compiler with own algebra, does NOT use the lattice`
- Source: `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:26-33`

### 2. Exact duplicate of the same statement in `trenchwork_clean`

- Quote:
  > `### P2K Lattice Compiler (FUN_1802c59b0)`
  > `- Shared inner compiler for Morph1/MorphLP/MorphLPX classes`
  > `- 3-stage cross-coupled lattice: stage N zeros = stage N-1 poles`
  > `- Numerator table (param_5) only affects stage 0; stages 1-2 always cross-coupled`
  > `- Morph2: separate compiler with own algebra, does NOT use the lattice`
- Source: `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:26-33`
- Assessment: this is the same text in a second repo location, not an independent second derivation. It strengthens availability, not independence. Sources: `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:1-4`, `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:1-4`

### 3. Adjacent contract evidence: LP / LPX definitely share one inner compiler and explicit numerator tables

- Quote:
  > `"description": "16-entry pole/gain placement tables for MorphLP and MorphLPX 3-stage LP filter classes. Extracted from DAT_1806d7480 (LPX) and DAT_1806d73c0 (LP) in EmulatorX.dll.",`
  > `"inner_compiler": "FUN_1802c59b0",`
- Source: `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:5-8`

- Quote:
  > `"note": "val1==val3 symmetry. Duplicated for both M0/M100 corners (morph-invariant zeros)."`
  > `"note": "6 distinct values: [M0_val1, M0_val2, M0_val3, M100_val1, M100_val2, M100_val3]. Per-corner zero variation."`
- Source: `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:17`, `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:42`

- Quote:
  > `vtable Slot [1]: FUN_1802c5f10 (Compile method). This reads the 0.0-1.0 morph position, computes an index floor((morph + offset) * 16.0), looks up values in a fixed zero table, and delegates compilation to the inner compiler function FUN_1802c59b0.`
  > `The inner compile method FUN_1802c59b0 actively controls a 3-stage active filter sequence (stages 4-6 are bypassed).`
- Source: `C:\Users\hooki\Trench\docs\looperator\phantom_morph_lpx_contracts.md:17-20`

- Quote:
  > `Both use the same inner compiler FUN_1802c59b0 and the same morph-to-index logic.`
  > `The difference is in the zero-placement table format:`
  > `CPhantomMorphLP changes zero presence dynamically across the morph sweep`
  > `CPhantomMorphLPX has static zeros, duplicating the same triplet across both M0 and M100.`
- Source: `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_comparison_MorphLP_vs_LPX.md:7-20`

- Quote:
  > `Unlike LPX (which duplicates the same triplet for both corners), MorphLP provides`
  > `distinct zero coefficients for M0 and M100. This means the morph axis controls`
  > `not just pole frequency but zero PRESENCE:`
  > `The cube's pole skeleton provides the denominator; this table provides the numerator.`
- Source: `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_MorphLP_zero_table.md:36-45`, `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_MorphLP_zero_table.md:89-91`

- Assessment: these contracts strongly support the existence of a shared 3-stage compiler with explicit numerator laws for `MorphLP` and `MorphLPX`, but they do not themselves restate the exact sentence `stage N zeros = stage N-1 poles`. Sources: `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:5-8`, `C:\Users\hooki\Trench\docs\looperator\phantom_morph_lpx_contracts.md:17-20`, `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_comparison_MorphLP_vs_LPX.md:7-20`, `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_MorphLP_zero_table.md:36-45`

## Contradicting or ambiguous evidence

### 1. The lattice law is not universal across all documented filter paths

- Quote:
  > `- Types 0-13: 14 Morpheus-era classes (LP2Pole through Batman) - all ROM-baked fixed coefficients`
  > `- Types 14-45: 32 CPhantomFilterP2k instances - P2K wrapper dispatching to Morph1/2/LP/LPX/Designer/Flanger/Vocal sub-classes`
- Source: `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:11-12`

- Quote:
  > `- Morph2: separate compiler with own algebra, does NOT use the lattice`
- Source: `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:33`

- Assessment: if the claim is read broadly as a law for all P2K-related filters, the surveyed docs already narrow it to a specific compiler family rather than a universal rule. Sources: `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:11-12`, `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:33`

### 2. `CPhantomFilterP2k` and some shipped filters are documented as ROM copy / static lookup, not lattice synthesis

- Quote:
  > `Talking Hedz is NOT a computed/compiled filter. It is a static ROM lookup.`
  > `There is no "Stage 0 type byte" to find. The type branching happened at E-mu's factory.`
- Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:36-36`

- Quote:
  > `| 7 | FUN_1802d3ce0 | 0x1802d3ce0 | CPhantomFilterP2k (ROM copy) | 6 | CLOSED - pure memcpy |`
- Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:52`

- Quote:
  > `- Do not conflate MorphDesigner filters (Q-degenerate, 2 corners duplicated to 4) with P2K filters (genuine 4-corner independence)`
- Source: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:168`

- Assessment: these passages contradict any broad reading that all P2K-facing filters are built by one lattice law at compile time. Sources: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:36`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:52`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:168`

### 3. Other surveyed compiler implementations shown in clean-room code compile stages independently

- Quote:
  > `freq_packed = int(round(sec.low_freq + morph * (sec.high_freq - sec.low_freq)))`
  > `gain_packed = int(round(sec.low_gain + morph * (sec.high_gain - sec.low_gain)))`
  > `return _compile_packed(sec.type, freq_packed, gain_packed, shift=shift)`
- Source: `C:\Users\hooki\trenchwork_clean\pyruntime\designer_compile.py:154-156`

- Quote:
  > `for sec in template.sections[:6]:`
  > `sp, enc = _compile_section(sec, morph, shift=shift)`
  > `stages.append(sp)`
  > `pre_encoded.append(enc)`
- Source: `C:\Users\hooki\trenchwork_clean\pyruntime\designer_compile.py:186-189`

- Quote:
  > `enc = type3_to_encoded(freq_packed, gain_packed, shift=shift)`
  > `return _encoded_to_stage_params(enc), enc`
- Source: `C:\Users\hooki\trenchwork_clean\pyruntime\heritage_coeffs.py:245-246`

- Assessment: these clean-room implementations do not show a previous-stage zero/pole coupling law. This does not disprove the `FUN_1802c59b0` claim, but it does mean the law is not surfaced across all current compiler code paths. Sources: `C:\Users\hooki\trenchwork_clean\pyruntime\designer_compile.py:154-156`, `C:\Users\hooki\trenchwork_clean\pyruntime\designer_compile.py:186-189`, `C:\Users\hooki\trenchwork_clean\pyruntime\heritage_coeffs.py:245-246`

### 4. Runtime dispatch and stage-count reports also argue against a universal one-law reading

- Quote:
  > `6 stages correct for FilterP2k. Other filter types use fewer.`
- Source: `C:\Users\hooki\trench_re_vault\analysis\contradictions.md:17`

- Quote:
  > `CPhantomRTFilter uses 8 fully-unrolled virtual methods (1/2/3/6 stages x mono/stereo).`
  > `No loop over stages at runtime - the stage count is compile-time-determined per method.`
- Source: `C:\Users\hooki\trench_re_vault\analysis\contradictions.md:74-75`

- Assessment: this does not directly contradict the narrow `FUN_1802c59b0` claim, but it reinforces that multiple stage regimes coexist and the lattice statement cannot be applied blindly across the whole filter inventory. Sources: `C:\Users\hooki\trench_re_vault\analysis\contradictions.md:17`, `C:\Users\hooki\trench_re_vault\analysis\contradictions.md:74-75`

## Assessment

The documented evidence supports the claim weakly.

Why weakly instead of strongly:

- The exact sentence `stage N zeros = stage N-1 poles` appears in two surveyed files, but they are duplicate notebook-ground-truth packs rather than independent primary derivations. Sources: `C:\Users\hooki\Trench\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:26-33`, `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:26-33`
- The adjacent clean-room contracts strongly support a shared 3-stage inner compiler with explicit numerator laws for `Morph1` / `MorphLP` / `MorphLPX`, but those contracts stop short of restating the exact cross-coupling identity in their own words. Sources: `C:\Users\hooki\trenchwork_clean\contracts\cleanroom\handoffs\lpx_zero_law\v1\CONTRACT.json:5-8`, `C:\Users\hooki\Trench\docs\looperator\phantom_morph_lpx_contracts.md:17-20`, `C:\Users\hooki\trench_re_vault\contracts\CPhantomMorphLPX\2026-03-16_comparison_MorphLP_vs_LPX.md:7-20`
- Several surveyed files explicitly show that the law is not a universal statement about all P2K-related filters, because ROM-only classes, `Morph2`, and `CPhantomFilterP2k` memcpy paths are handled separately. Sources: `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:36`, `C:\Users\hooki\trenchwork_clean\docs\context\re-proven-facts.md:52`, `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:11-12`, `C:\Users\hooki\trenchwork_clean\docs\notebooklm\packs\ground-truth\morpheus_cube_architecture_proven.md:33`

Narrow reading: `FUN_1802c59b0` for `Morph1` / `MorphLP` / `MorphLPX` probably does enforce the documented lattice cross-coupling law.

Broad reading: `P2K filters in general enforce Stage N zeros = Stage N-1 poles` is not supported by the surveyed documentation.

## Files surveyed but not cited

- `C:\Users\hooki\trench_re_vault\analysis\morph_designer_decompile.md`
- `C:\Users\hooki\trench_re_vault\analysis\morph_designer_analysis.json`
- `C:\Users\hooki\trench_re_vault\analysis\open_questions.md`
- `C:\Users\hooki\trench_re_vault\analysis\test_targets.md`
- `C:\Users\hooki\trenchwork_clean\docs\context\q-law.md`
- `C:\Users\hooki\Trench\trench-juce\forge\pyruntime\designer_compile.py`
- `C:\Users\hooki\Trench\trench-juce\forge\pyruntime\heritage_coeffs.py`
