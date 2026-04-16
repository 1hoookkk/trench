# Source
- `C:\Users\hooki\trench_re_vault\extract\emulatorx_morph\sse_filter_diff.json`
- `C:\Users\hooki\trench_re_vault\extract\emulatorx_morph\2026-03-12_remaining_targets_a\remaining_targets.json`

# What this data is
This file documents the binary difference between `CPhantomRTFilter` and `CPhantomRTSSEFilter` (derived via runtime type information or RTTI scanning) and records the inventory of top reverse-engineering targets remaining inside the EmulatorX codebase.

# Why it's here
In the TRENCH architecture, resolving DSP variations cleanly prevents phantom bugs. Identifying that the SSE variant has zero algorithmic divergence ensures that the core `biquad` DSP logic in the CLEAN repo doesn't need maintaining multiple parallel floating-point schemas. Furthermore, establishing the remaining target slate directs future triage efforts correctly.

# Sanitized content

## Filter Algorithmic Delta
There is **no algorithmic or behavioral difference** between the base filter and the SSE filter models. 
The differences exist strictly in the coefficient interpolation step taking advantage of basic SIMD structure optimization:
- **Base `CPhantomRTFilter`**: Computes five sequential `MOVSS`/`ADDSS`/`MOVSS` calls across the 5 coefficients of a biquad.
- **SSE `CPhantomRTSSEFilter`**: Vectorizes the first four offsets `(0x00, 0x04, 0x08, 0x0C)` into a single packed `MOVAPS`/`ADDPS`/`MOVAPS` block execution, falling back to a scalar `ADDSS` for the 5th trailing term `(0x10)`.
- **Output Validation**: Due to 100% assembly equality existing natively for the core `MULSS`/`ADDSS`/`SUBSS` operation calculating $y[n]$, the generated wave paths are bit-identical variants.

## Remaining RE Targets (Triaged List)

| Target ID | Priority | Objective | Why It Matters |
|---|---|---|---|
| `preset_name_to_numeric_selector_join` | **P1** | Map resource labels (e.g., "Talking Hedz") directly to numeric indexed DSP objects. | Closes the loop from named Skin to the actual loaded DSP logic block. |
| `runtime_bank_consumer_and_duplication_reason` | **P1** | Understand why runtime structures duplicate the two morph corners into four memory banks. | Validates if TRENCH needs a 4-bank memory footprint or only 2. |
| `stage_type_semantics_for_designer_slots` | **P1** | Demystify internal designer types IDs `(1, 2, 3)` to human readable enums. | Unlocks clean node type routing. |
| `sample_rate_family_laws_for_morph_compiler` | **P1** | Chart the mapping parameters governing frequency/gain across standard 44.1/48/96kHz variants. | Provides uniform multi-rate DSP support. |
| `macro_field_identities_in_compile_context` | P2 | Locate variables acting as the global gain/freq biases across runtime structs. | Ensures precise implementation of "Macro State" behavior models. |
| `serializer_ingress_for_designer_authoring_data` | P2 | Map the XML instantiation paths into functional compiler structures. | Permits native 6-slot serialization in standard forms. |
| `display_graph_generation_path` | P2 | Prove whether visual feedback generates from GUI cache buffer, coefficients, or raw authoring state. | Replaces guessing where to bind GUI evaluators. |

# Integration notes
- Do not maintain isolated `SSE` code paths for Emulation constraints. Unify any CPU-accelerated TRENCH math confidently using standard DSP arrays.
- The highest priority follow-up target remains bridging the UI text skin name with the precise numerical filter ID to avoid "guessing" which preset wrappers actuate which filter shapes.

# UNKNOWNs
- Why the Morph Designer outputs duplicate into four identical bank assignments. This target (`runtime_bank_consumer_and_duplication_reason`) prevents fully optimizing the local state architecture footprint.
