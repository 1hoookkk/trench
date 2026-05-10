# TRENCH Ground Truth

This is the one-file domain brief for authoring and design work.

## The System In One Pass

TRENCH is a cartridge-based Z-plane filter instrument with a frozen shipping runtime and a separate authoring/design world around it.

The mistake to avoid is treating reference material as shippable product.

There are three working layers:

1. **Shipping runtime truth**
2. **Authoring/design truth**
3. **Reference material**

## 1. Shipping Runtime Truth

This is what actually makes sound in the current shipping tree.

- Runtime plugin path: `C:\Users\hooki\Trench\trench-juce\plugin\`
- Runtime DSP path: `C:\Users\hooki\Trench\trench-juce\plugin\source\PluginProcessor.cpp` -> `TrenchEngine.cpp`
- Runtime cartridge contract: `compiled-v1`
- Runtime math/cartridge law: `C:\Users\hooki\Trench\SPEC.md`
- Runtime doctrine: `C:\Users\hooki\Trench\DOCTRINE.md`

When a statement about runtime conflicts with repo prose, build wiring wins.

## 2. Authoring / Design Truth

This is the factory around the shipping runtime.

It includes:

- body-session method
- sonic tables
- vault and MCP-assisted analysis
- forge/foundry/workbench surfaces
- design references for faceplate, controls, and rollers

Primary sources:

- `C:\Users\hooki\Trench\authoring\FILTER_WORK.md`
- `C:\Users\hooki\Trench\authoring\`
- `C:\Users\hooki\Trench\authoring\compilers\`
- `C:\Users\hooki\Trench\authoring\tools\`
- `C:\Users\hooki\Trench\PHONEMES.md`
- `C:\Users\hooki\Trench\docs\body_authoring_seed.md`
- `C:\Users\hooki\Trench\docs\sonic_tables\`

Reference sources:

- `C:\Users\hooki\trenchwork_clean`
- `C:\Users\hooki\trenchwork`

This layer informs what should be built and shipped. It does not override live runtime wiring.

## 3. Reference Material

This is material that can inform taste, measurements, debugging, or rejection
of clones.

Examples:

- `C:\Users\hooki\trench_re_vault`
- `C:\Users\hooki\trenchwork_clean`
- `C:\Users\hooki\trenchwork`
- `C:\Users\hooki\do-it`

Reference material is not a runtime dependency and does not become a shipping
filter by being useful.

The hard shipping rule:

```text
TRENCH does not ship E-mu filters.
```

E-mu/P2K references can be used to reject accidental clones. They are not the
target for product cartridges.

## Current External Tree Roles

- `C:\Users\hooki\trenchwork_clean`
  - older authoring/factory reference tree
  - sonic tables, MCP, runtime and plugin experiments, vault, tools
  - not canonical unless an artifact is imported, indexed, or cited here

- `C:\Users\hooki\trenchwork`
  - broader workbench architecture
  - forge, foundry, MCP, plugin/runtime shell
  - reference/import source only

- `C:\Users\hooki\do-it`
  - historical runtime implementation and solved-state precedent
  - useful for old validation logic and prior runtime decisions

- `C:\Users\hooki\trench_re_vault`
  - research evidence and reference behavior
  - not a shipping cartridge source

## Canonical Authoring Rule

For filter work, `C:\Users\hooki\Trench` is the canonical repo and
`authoring/` is the canonical doorway.

Use this hub for:

- "how do I author?"
- "where are the sonic tables?"
- "what is the workbench?"
- "where do design references live?"
- "how do I prevent E-mu/P2K clone cartridges from shipping?"

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
- `history/README.md`
