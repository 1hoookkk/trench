# TRENCH Ground Truth

This is the one-file domain brief for authoring and design work.

## The System In One Pass

TRENCH is a cartridge-based Z-plane filter instrument with a frozen shipping runtime and a separate authoring/design world around it.

The mistake to avoid is flattening everything into one truth layer.

There are three:

1. **Shipping runtime truth**
2. **Authoring/design truth**
3. **DIRTY provenance truth**

## 1. Shipping Runtime Truth

This is what actually makes sound in the current shipping tree.

- Runtime plugin path: `C:\Users\hooki\Trench\trench-juce\plugin\`
- Runtime DSP path: `C:\Users\hooki\Trench\trench-juce\plugin\source\PluginProcessor.cpp` -> `TrenchEngine.cpp`
- Runtime cartridge contract: `compiled-v1`
- Runtime math/cartridge law: `C:\Users\hooki\Trench\SPEC.md`
- Runtime doctrine: `C:\Users\hooki\Trench\DOCTRINE.md`

When a statement about runtime conflicts with repo prose, build wiring wins.

## 2. Authoring / Design Truth

This is the clean factory around the shipping runtime.

It includes:

- body-session method
- sonic tables
- vault and MCP-assisted analysis
- forge/foundry/workbench surfaces
- design references for faceplate, controls, and rollers

Primary sources:

- `C:\Users\hooki\Trench\PHONEMES.md`
- `C:\Users\hooki\Trench\docs\body_authoring_seed.md`
- `C:\Users\hooki\Trench\docs\sonic_tables\`
- `C:\Users\hooki\trenchwork_clean`
- `C:\Users\hooki\trenchwork`

This layer informs what should be built and shipped. It does not override live runtime wiring.

## 3. DIRTY Provenance Truth

This is reverse-engineering and capture truth.

Primary source:

- `C:\Users\hooki\trench_re_vault`

DIRTY truth is used to learn target behavior and produce sanitized contracts.
It is not a runtime dependency and not a clean authoring surface.

## Current External Tree Roles

- `C:\Users\hooki\trenchwork_clean`
  - strongest clean authoring/factory tree
  - sonic tables, MCP, runtime and plugin experiments, vault, tools

- `C:\Users\hooki\trenchwork`
  - broader workbench architecture
  - forge, foundry, MCP, plugin/runtime shell

- `C:\Users\hooki\do-it`
  - historical runtime implementation and solved-state precedent
  - useful for old validation logic and prior runtime decisions

- `C:\Users\hooki\trench_re_vault`
  - DIRTY airlock
  - capture, extraction, analysis, sanitization source

## Canonical Authoring Rule

Use this hub for:

- "how do I author?"
- "where are the sonic tables?"
- "what is the workbench?"
- "where do design references live?"
- "where does DIRTY end and clean begin?"

Use live runtime files for:

- "what ships?"
- "what code runs?"
- "what is the actual cartridge/runtime contract?"

## Immediate Canonical Entry Files

- `README.md`
- `CLAUDE.md`
- `ground_truth/README.md`
- `sonic_tables/README.md`
- `workbench/README.md`
- `design/README.md`
- `dirty_airlock/README.md`

