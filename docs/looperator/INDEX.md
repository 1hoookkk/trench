# TRENCH Looperator Extraction Specifications

This directory contains sanitized hardware mapping contracts, acoustic characteristics, and display rules extracted from the Emulator X `DIRTY` vault infrastructure, rewritten for explicit integration into the TRENCH `CLEAN` repository. 

All files omit reverse-engineering process minutiae and focus solely on structural software integration requirements for building the "Looperator" module.

## Document Index

| Specification | Core Concern | Description |
|---|---|---|
| [filter_cube_display_model.md](./filter_cube_display_model.md) | **UI/UX Model** | Defines the structural deduplication rendering rule converting 6-section template clusters into singular physical geometric bounds for the grid visualizer. |
| [filter_atlas.md](./filter_atlas.md) | **Parameter Limits** | Explains the 49 baseline hardware filters from legacy presets to govern base capabilities. |
| [phantom_morph_lpx_contracts.md](./phantom_morph_lpx_contracts.md) | **DSP Constants** | Extracted explicit 16-point table driving the dynamic zero-placement rules for MorphLP vocal filters. |
| [morph_designer_behavior.md](./morph_designer_behavior.md) | **Host Boundary** | Clarifies the transition layer parsing static preset skin names against physical engine execution selection blocks. |
| [qsound_spatial.md](../archive/qsound_spatial.md) | **Spatial Matrix (archived — paper feature, no code yet)** | Derived parameters for pan ITD/ILD stereoscopic offset. |
| [talking_hedz_x3_surfaces.md](./talking_hedz_x3_surfaces.md) | **Stage Validation** | Unwinds discrepancies behind the 5-stage vs 6-stage vocal implementation, ensuring accurate filter routing over P2K hardware models. |
| [sse_filter_diff.md](./sse_filter_diff.md) | **Math Optimizations** | Confirms standard vectorization of biquad interpolations does not necessitate unique algorithmic maintenance in emulation branches. |
| [cvsd_module.md](./cvsd_module.md) | **Tonal Modeling** | Explicit timing constants ($\tau$) to simulate authentic step-based delta degradation dynamically. |

## Implementation Notes
Use these documents as immediate constraints for C++ backend compilation structures and JUCE-based UI graph evaluators in the primary TRENCH repo. When interpreting preset catalogs, `morph_designer` boundaries apply absolute decoupling semantics to protect engine execution from semantic collisions.
