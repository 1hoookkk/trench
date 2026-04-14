---
name: learning_loop_status
description: Status of the 3-layer self-learning loop — what's wired, what's broken, how to use MCP + RLS for generation feedback
type: project
---

Three-layer learning loop for forge generation quality.

## Layer 1 — RLS verdict learning (WIRED, 2026-03-16)
- Every Y/N verdict in TRENCHWORK updates `rls_state.json` via online RLS.
- Pipeline reads RLS weights (not stale OLS `taste_model.json`) and blends 30% into trust_score.
- Feature vector: [movement_norm, emergence_norm, fitness, dead_zones, novelty].

### Current RLS weight interpretation (66 updates as of 2026-03-16)
- movement: ~0 (neutral signal)
- emergence: -1.61 (user kills high-emergence — raw spectral energy != quality)
- fitness: +0.35 (mild positive — composite fitness is a weak keep signal)
- dead_zones: +1.04 (user keeps bodies with dead zones — these are character features, not flaws)
- novelty: -1.02 (user kills novel-for-novelty's-sake — prefers coherent familiar patterns)

**How to apply:** When tuning generation parameters or scorer weights, respect these biases. Don't optimize for emergence or novelty as quality signals — the user's verdict history says they're anti-correlated with keeps.

## Layer 2 — CLAUDE.md / project memory (PARTIALLY WIRED)
- This memory file captures RLS interpretations for future sessions.
- Future sessions should: check `rls_state.json` weights, run `scan_vault` via MCP to see current vault quality, and use both to inform forge guidance.
- Not automated — requires Claude to read this file and act on it.

## Layer 3 — MCP → forge feedback (NOT WIRED)
- `scan_vault` MCP tool ranks vault bodies by composite (talk×0.4 + ridge×0.3 + coherence×0.3).
- These rankings don't feed back into the forge pipeline. Vault survivors don't seed next generation.
- Planned work: "vault as vocabulary" — vault survivors feed back as seeds.
- Gap: No mechanism to translate MCP scan results into `FitnessWeights` adjustments or vocabulary bias.

**Why:** Without Layer 3, generation is memoryless across sessions. RLS biases ranking of what's already generated, but doesn't change WHAT gets generated. The forge keeps producing the same distribution of bodies and then RLS re-ranks them. True learning would shift the generation distribution toward body types that survive triage.

**How to apply:** When implementing vault-as-vocabulary, prioritize bodies with high scan_vault composite scores as seeds. When adjusting FitnessWeights, use RLS weight signs as directional guidance (e.g., reduce emergence weight since user kills high-emergence).
