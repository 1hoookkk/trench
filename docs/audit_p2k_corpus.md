# P2K Corpus Audit — Session Brief

## Objective

Analyze 10+ P2K corpus bodies to extract the structural rules that make them work. Rewrite CLAUDE.md authoring rules based on findings.

## Method

Use `trench-mcp` tools:
- `scan_vault` — find all available P2K bodies
- `analyze_state` — examine each body at M0/M50/M100 at Q0/Q50/Q100
- `compare_states` — measure morph and Q deltas
- `probe_frequency` — check specific frequency bands
- `find_corridors` — identify morph sweep patterns

For each body, extract:
1. Per-stage frequency placement at all 4 corners
2. Radius distribution (min, max, spread)
3. Zero placement strategy (allpole vs interior vs unit-circle)
4. Shared b0 values per corner (from c0 field)
5. Morph strategy: swap, fan, cluster-scatter, radius-only
6. Q behavior: compression target, radius modulation
7. Whether any stage is passthrough at any corner
8. DC gain at each corner
9. Cascade peak frequency migration across morph sweep
10. Number of visible formant peaks at each morph position

## Target Bodies

Start with these known-good bodies from the P2K corpus:
- BassBox 303 (Sub Physicalizer) — frequency swap, sub-focused
- Talking Hedz (Alien Formant) — radius-only morph, vocal
- Ear Bender (Pad Breather) — HF wall fold
- Meaty Gizmo (Circuit Breaker) — one-vs-five suppression
- Fuzzi Face — resonance phalanx
- Any 5 additional P2K skins with high taste scores

## Questions to Answer

1. Do ANY heritage bodies use passthrough stages at any corner? Or are all 6 always active?
2. What is the typical b0 range? Does it ever go below 0.01?
3. What radius range is typical? Do bodies regularly use r > 0.995?
4. How do bodies handle sub-bass without DC explosion?
5. What zero placement patterns appear most often?
6. How wide is the frequency spread at M0 vs M100?
7. Do stages cross frequencies during morph (frequency swap)?
8. What makes the morph sweep sound like movement vs a static shift?

## Output

Write findings to `docs/p2k_audit_findings.md` with:
- Raw data tables per body
- Extracted rules (numbered, falsifiable)
- Recommended CLAUDE.md rewrites
- List of assumptions from current CLAUDE.md that are wrong

## Constraints

- Read-only audit. Do not modify any code.
- Do not run the optimizer or forge.
- Do not deploy anything to the plugin.
- Trust the heritage data. If a body does something unexpected, the body is right.
