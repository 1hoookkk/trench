# Source
- `C:\Users\hooki\trench_re_vault\extract\morph_designer\2026-03-12_host_selection_trace_a\` (2 files)
- `C:\Users\hooki\trench_re_vault\extract\morph_designer\2026-03-12_schema_v2_a\` (4 files)
- `C:\Users\hooki\trench_re_vault\extract\morph_designer\2026-03-26_parser_boundary_b\` (6 files)
- `C:\Users\hooki\trench_re_vault\extract\morph_designer\2026-03-26_selector_boundary_a\` (3 files)

# What this data is
These are behavioral captures tracking how the original E-mu MorphDesigner tool parsed XML templates into memory, selected the correct runtime DSP classes, and mapped string resources (like "Talking Hedz") to raw mathematical shapes. They document the boundary where a serialized user token crosses into compiled runtime filter objects.

# Why it's here
The Looperator is a modern evolution of the MorphDesigner UI constraint system. Understanding the boundaries between parsed templates, named wrapper "presets", and actual compiled DSP logic allows the CLEAN repo to design a predictable selector grid that won't produce invalid object mappings or lose invariant template structures.

# Sanitized content

## `2026-03-12_host_selection_trace_a`
**Files**: `EXTRACT_MANIFEST.json`, `notes_behavioral.md`

**Description**: Maps the decoupling of the "named skin" from the runtime numeric DSP recipe.
**State-Machine Logic**: 
- The selection is decoupled into three logical layers: authoring/wrapper name surface $\to$ numeric filter-family selector $\to$ compiled runtime stage recipe.
- Numeric filter-family keys direct the actual architecture (e.g., `filterType: 13`).
- Host UI string tables (like `Talking Hedz`) act strictly as named wrappers/defaults overlaid on numeric families, preventing name collisions from disrupting DSP execution.

## `2026-03-12_schema_v2_a`
**Files**: `EXTRACT_MANIFEST.json`, `morph_designer_contract.schema.json`, `morph_designer_contract.v2.example.json`, `notes_behavioral.md`

**Description**: Redefines the contract of MorphDesigner as a family/template hierarchy, rather than a raw coefficient generator.
**Schema v2 Field Definitions**:
- `template_catalog`: Reusable structures mapping to raw stage definitions or runtime keyframes. 
- `preset_catalog`: Named wrappers containing factory macro state pointing to exact templates.
- `representation_crosswalk`: An explicit map correlating designer XML shapes, runtime banks, runtime keyframes, and visualization/display graph outputs.
- `template_equivalence_groups`: Sets of named presets that ultimately collapse down to the semantic equivalence of one template.

## `2026-03-26_parser_boundary_b`
**Files**: `EXTRACT_MANIFEST.json`, `filter_authoring_pointer_scan.json`, `filter_authoring_pointer_scan.md`, `filter_authoring_xref_map.json`, `filter_authoring_xref_map.md`, `notes_behavioral.md`

**Description**: Tracks the internal conversion of raw string-based XML tags into compiled structural memory blocks.
**Parser Entry/Exit**:
- The XML authoring tokens (e.g., `filter/type-absolute`, `filter/designer-section/1/type`) correspond directly to pointer slots inside a concrete `.data` schema aggregate spanning `0x18088a8c8` to `0x18088ab58`.
- Parsing populates object representations internally which dictates vtable assignments, dispatching the layout cleanly through `CPhantomMorphDesigner` or `CPhantomMorph2` objects to memory structures like `MorphDesignerCompile36B`.

## `2026-03-26_selector_boundary_a`
**Files**: `EXTRACT_MANIFEST.json`, `notes_behavioral.md`, `selector_boundary_candidates.json`

**Description**: Evaluates the explicit class selection triggers dictating boundary validity for incoming parameter fields.
**Selector Boundary Conditions**:
- A valid selection relies on the `<type-absolute>` XML field.
- **`type-absolute = 108`**: Maps identically to host registry key `0x6c` (Dec 108), dynamically spawning the `CPhantomMorphDesigner` class and hitting the `MorphDesignerCompile36B` compiler routines.
- **`type-absolute = 144`**: Matches registry key `0x90` (Dec 144) resulting in an implicit `CPhantomFilterP2k` class object constructed under subtype `13`.
- The system rejects or ignores raw configuration blobs missing compatible registry/class parity, validating based on class constructor slot lookup natively.

# Integration notes
- Token definitions must completely isolate descriptive labels (e.g. presets) from functional `type-absolute` keys in code models.
- Apply the Schema v2 logic inside Looperator logic architecture: parse raw structures $\to$ assign normalized Template $\to$ assign optional preset alias $\to$ render Display surface visually.
- Never feed variable templates directly into fixed index parsers; guarantee the correct underlying compiler maps to the right DSP family struct shape (108 mapped exclusively onto Morph behavior).

# UNKNOWNs
- **Direct String Entry Function**: A complete unbroken trace from the raw strings being ingested via a physical parse function into the memory struct was missing/blocked in reverse engineering tools. 
- **Filter Name Join Mechanism**: How the exact explicit string table keys link their references to numeric ID selection blocks logic is missing the load/transport vector trace.
