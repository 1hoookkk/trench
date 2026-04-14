---
name: Forge Direction — Body-Level Generation
description: User's stated priorities for forge improvement. Body-level refactor is #1. Low-Q gate is fastest win. Historical contract is separate lane.
type: project
---

## User's forge direction (2026-03-14)

**Most impact**: Body-level generation refactor — replace `generate_stage() x 6` with body-level generation. Generate one body skeleton with internal anatomy, then instantiate stages inside role bounds.

**Fastest win**: Low-Q identity gate — hard rejection for bodies that only come alive at extreme Q. Now implemented in `expressiveness.rs::low_q_identity()`.

**Important but separate**: Historical handoff contract — finish sanitized vector contract for historical stage-code seam.

## Key principles (user-stated)

- Soft anatomy priors, hard rejection
- Structural spectral spread (forced, not emergent)
- First-class zero/subtractive authoring (carve behavior, not only pushes)
- Low-Q identity as a modern forge criterion
- Good role vocabulary: anchor, mass, carve, pressure, detail, boundary, cleanup
- Do NOT hardcode rigid stage ontology ("Stage 0 always X")
- Do NOT resurrect Bank thesis as main engine architecture

## Working hypothesis

The forge is weak because it generates six cousin stages instead of one coherent body with differentiated internal roles.

**Why:** User has repeatedly confirmed this diagnosis across multiple messages. The skeleton anatomy tightening (Phase 3) was a step, but the root cause requires architectural change in `generate.rs`.

**How to apply:** When working on forge generation, always think body-first, not stage-first. Any new generation code should start from a body plan and derive stages from it.
