---
name: dishonesty-catalogue
description: Full ranked catalogue of musical dishonesty acts for TRENCH runtime — 12 categories, implementation tiers, anti-patterns
type: project
---

**Why:** "Don't lie randomly. Lie in ways that preserve body identity while increasing stress, motion, and character."

## Ranked Top 10 Experiments
1. Morph warping by body region (Cat 2)
2. Split-Q law — q_body / q_edge / q_boundary (Cat 1, 4)
3. Anchor preservation compensation (Cat 5, 7, 11)
4. Carve hysteresis — zeros lag/overshoot relative to poles (Cat 6, 9)
5. Interior-vs-boundary unequal scaling (Cat 7)
6. Stress-zone split-code generalization beyond Type 3 (Cat 3)
7. Frequency crowding / repulsion laws (Cat 3)
8. Input-dependent carve depth (Cat 10)
9. Direction-dependent interpolation (Cat 8)
10. Boundary stage soft topology shift (Cat 12)

## 12 Categories
1. Split one control into two internal truths (f_intent vs f_realized)
2. Non-uniform morph travel (sigmoid, sticky attractors, hysteresis)
3. Frequency-space dishonesty (region warping, identity-preserving, stress zones)
4. Q dishonesty (contextual, asymmetric, persistence law, overstated display)
5. Gain / posture preservation (preserve weight not level)
6. Zero dishonesty (zeros lag poles, overshoot, become sticky near carve bands)
7. Stage-role dishonesty (6 political actors not 6 democratic stages)
8. Interpolation dishonesty (role-space interp, region-sensitive weights)
9. Time-domain dishonesty (directional memory, pressure accumulation, resonance fatigue)
10. Input-dependent dishonesty (subtle internal self-adjustment)
11. Spectral-center dishonesty (preserve perceived center of gravity)
12. Topology dishonesty (stage law changes across morph invisibly)

## 4-Layer Framework
1. Control — knob says one thing, engine hears another
2. Structural — stages do unequal jobs
3. Spectral — displayed vs actual energy diverge
4. Temporal — behavior depends on motion history

## Implementation: DishonestyProfile struct with toggles
- `morph_warp_mode`, `split_q_mode`, `anchor_preserve_amount`
- `carve_hysteresis_ms`, `stress_profile_mode`, `boundary_soft_shift`

## Anti-patterns
- No random modulation, arbitrary saturation, unexplained gain jumps
- No lies that destroy repeatability or flatten body identity
- Best dishonesty: stable, repeatable, musically legible, feels intentional

## How to apply
Build as post-processing toggles, not hardcoded. Audition in foundry REPL.
Must be deterministic — same coordinate eventually settles to same state.
